/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "assert.h"
#include "physics.h"
#include "macros.h"
#include "vec.h"
#include "force.h"
#include "invblock.h"
#include "confio.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "copyrite.h"
#include "qmmm.h"
#include <stdio.h>
#include <string.h>
#include "gmx_fatal.h"
#include "typedefs.h"
#include <stdlib.h>

#include <mkl_types.h>
#include <mkl_lapack.h>


#define PS2AU 41341.3733365586

/* diagonalize hermitian matrix using MKL routine zheev */
static void 	diag(int n, double *w, double **V, double **M)
{
	char jobz, uplo;
	int i, j, lwork, info;
	double *rwork, wkopt;
	double *M_mkl, *work;

	jobz = 'V';
	uplo = 'U';
    	snew(M_mkl, n*n);
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			M_mkl[j + i*n] = M[i][j];
		}
	}

        lwork = -1;
        dsyev(&jobz, &uplo, &n, M_mkl, &n, w, &wkopt, &lwork, &info );
        lwork = (int)wkopt;
	snew(work, lwork);
	dsyev(&jobz, &uplo, &n, M_mkl, &n, w, work, &lwork, &info);
	if (info != 0)
	{
		gmx_fatal(FARGS, "Lapack returned error code: %d in dsyev", info);
	}

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			V[i][j] = M_mkl[j + i*n];
		}
	}

	sfree(M_mkl);
	sfree(work);
}


//Switching function for MM charges
double swf(double r, double cutoff)
{
	double start = 0.8 * cutoff;

	if (r < start)
		return 1.0;

	if (r > cutoff)
		return 0.0;

	double a = 1.0 / (cutoff * cutoff - start * start);
	double a3 = a * a * a;
	double a4 = a3 * a;
	double a5 = a4 * a;

	double b = r * r - start * start;
	double b3 = b * b * b;
	double b4 = b3 * b;
	double b5 = b4 * b;

	return 1.0;//- 10.0 * a3 * b3 + 15.0 * a4 * b4 - 6.0 * a5 * b5;
}

real inproduct(real *a, real *b, int n)
{
	int
		i;
	real
		dot=0.0;

	/* computes the inner product between two vectors (a.b), both of
	 * which have length n.
	 */  
	for(i=0;i<n;i++){
		dot+=a[i]*b[i];
	}
	return(dot);
}

int hop(int step, int state_i, int state_j, t_QMrec *qm)
{
	int
		swap = 0;
	real
		d11=0.0,d12=0.0,d21=0.0,d22=0.0;

	/* calculates the inproduct between the current Ci vector and the
	 * previous CI vector. A diabatic hop will be made if d12 and d21
	 * are much bigger than d11 and d22. In that case hop returns true,
	 * otherwise it returns false.
	 */  
	if(step){ /* only go on if more than one step has been done */
		d11 = inproduct(qm->CIvec[state_i],qm->CIvecold[state_i],qm->CIdim);
		d12 = inproduct(qm->CIvec[state_i],qm->CIvecold[state_j],qm->CIdim);
		d21 = inproduct(qm->CIvec[state_j],qm->CIvecold[state_i],qm->CIdim);
		d22 = inproduct(qm->CIvec[state_j],qm->CIvecold[state_j],qm->CIdim);
	}
	fprintf(stderr,"-------------------\n");
	fprintf(stderr,"Overlap between states:\t%i\t%i\n", state_i+1, state_j+1);
	fprintf(stderr,"d11 = %13.8f\n",d11);
	fprintf(stderr,"d12 = %13.8f\n",d12);
	fprintf(stderr,"d21 = %13.8f\n",d21);
	fprintf(stderr,"d22 = %13.8f\n",d22);
	fprintf(stderr,"-------------------\n");

	if((fabs(d12)>sqrt(0.7))&&(fabs(d21)>sqrt(0.7)))
		swap = 1;

	return(swap);
}

/* Terachem interface routines */

void init_msdsh_tc(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
{
	FILE    
		*rffile=NULL,*out=NULL;
	char
		*buf;
	int
		i;

	if(!qm->nQMcpus){

		/* init gradually switching on of the SA */
		qm->SAstep = 0;

		/* we read the number of cpus and environment from the environment
		 * if set.  
		 */
		snew(buf,1000);
		buf = getenv("SLURM_NPROCS");
		if (buf)
		{
			sscanf(buf,"%d",&qm->nQMcpus);
		}
		else
		{
			qm->nQMcpus=1;
		}
		fprintf(stderr,"number of CPUs for Terachem = %d\n",qm->nQMcpus);
		
		/* Firefly input file */
		snew(buf,1000);    
		buf = getenv("TERACHEM_INP");
		if (buf)
		{
			sscanf(buf,"%s",qm->QMinp);
		}
		else
		{
			gmx_fatal(FARGS,"no $TERACHEM_INP input file defined\n");
		}
		fprintf(stderr,"Terachem input file = %s\n",qm->QMinp);

		/* Firefly batch script*/
		snew(buf,1000);    
		buf = getenv("TERACHEM_EXE");
		if (buf)
		{
			snew(qm->gauss_exe,1000);
			sscanf(buf,"%s",qm->gauss_exe);
			fprintf(stderr,"Terachem Batch sript = %s\n",qm->gauss_exe);
		}
		else
		{
			gmx_fatal(FARGS,"no $TERACHEM_EXE batch script defined\n");
		}

		qm->step = 0;
		qm->swap = FALSE;
	}


	fprintf(stderr,"Terachem initialised...\n");
}  
                                                 
void write_msdsh_input_tc(int step,gmx_bool swap,
		t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, const double dw)
{
	int
		i, j, closed, occ;
	t_QMMMrec
		*QMMMrec;
	char
		*buf, *buf1, *pos;
	char
		periodic_system[37][3]={"XX","H ","He","Li","Be","B ","C ","N ",
			"O ","F ","Ne","Na","Mg","Al","Si","P ",
			"S ","Cl","Ar","K ","Ca","Sc","Ti","V ",
			"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga",
			"Ge","As","Se","Br","Kr"};
	FILE
		*out, *lat, *inp;
	char datex;

	double swcut, rm, dist;
	double *mmr, *mmc;
	
	QMMMrec = fr->qr;
	out = fopen("tc.inp","w");
	if (out==NULL) 
	{
		gmx_fatal(FARGS, "Error while trying to write tc.inp file");
	}
	inp = fopen(qm->QMinp,"r");
	if (inp==NULL) 
	{
		gmx_fatal(FARGS, "qm_msdsh_ff.c", 148, "Error while trying to read %s file", qm->QMinp);
	}

      	/* Now put qm.xyz and mm.xyz into tc.inp */
        fprintf(out,"coordinates             qm.xyz\n");
        if (mm->nrMMatoms > 0) 
		fprintf(out,"pointcharges            mm.xyz\n");


	/* Now reading input and copy it to tc.inp, everything except of coordinates and pointcharges*/
	snew(buf,1000);
	snew(buf1,1000);
	while (!feof(inp))
	{
		fgets(buf, 1000, inp);
		if (strncmp(buf,"coordinates",11) != 0 && strncmp(buf,"pointcharges",12) != 0)
		{
			fputs(buf, out);
		}
	}


	fclose(out);
	fclose(inp);

	/* Now write QM System into qm.xyz */
	fprintf(stderr,"Number of QM Atoms= %d\n",qm->nrQMatoms);

	out = fopen("qm.xyz","w");
	if (out==NULL) 
	{
		gmx_fatal(FARGS, "Error while trying to write qm.xyz file");
	}

       	fprintf(out,"%d\n",qm->nrQMatoms);
       	fprintf(out,"\n");

	for (i=0;i<qm->nrQMatoms;i++)
	{
		fprintf(out,"%s %.10lf %.10lf %.10lf\n",
				periodic_system[qm->atomicnumberQM[i]],
				qm->xQM[i][XX]*10,
				qm->xQM[i][YY]*10,
				qm->xQM[i][ZZ]*10);
	}

	fclose(out);

	/* Now write QM System into mm.xyz */
	if (mm->nrMMatoms > 0)
	{
		/* Now write MM point charges into mm.xyz */
		fprintf(stderr,"Number of MM Atoms= %d\n",mm->nrMMatoms);
              	out = fopen("mm.xyz","w");
		if (out==NULL) 
		{
			gmx_fatal(FARGS, "Error while trying to write qm.xyz file");
		}

       		fprintf(out,"%d\n",mm->nrMMatoms);
       		fprintf(out,"\n");

		for (i=0;i<mm->nrMMatoms;i++)
		{
			fprintf(out,"%.10lf %.10lf %.10lf %.10lf\n",
				mm->MMcharges[i],
				mm->xMM[i][XX]*10,
				mm->xMM[i][YY]*10,
				mm->xMM[i][ZZ]*10);
		}

		fclose(out);
	}

	sfree(buf);
	sfree(buf1);
	fflush(stderr);

} /* write_msdsh_SH_input */

void read_msdsh_output_tc(rvec QMgrad[],rvec MMgrad[],int step,
		gmx_bool swapped,t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, real *QMener, const gmx_bool readCI)
{
	int
		i,j,atnum,k,tots;
	char
		buf[300], buf1[14];
	char gradstring[] = "Gradient units are Hartree/Bohr";
	char gradmmstring[] = "------- MM / Point charge part -------";
        char enerstring[] = "terachem gradient in au";
	char cistring[] = " CI vector";
	char ncsfstring[] = " Number of CSFs:";
	char convfailstring[] = " ** WVFN ****  MAXIMUM NUMBER OF ITERATIONS REACHED";
	char tmpstr[300];
	char *pch;
	FILE
		*in, *chgrad, *inb;
	int *CSS;
	double fbuf;
	double a[81], b[9], sln[9];
	double xa,ya,za, ra;
	double tf[3], tq[3];
	double r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,rmx,rmy,rmz;
	double r21x,r21y,r21z,r31x,r31y,r31z,r32x,r32y,r32z;
	double dx, dy, dz;
	double xx,xy,xz,yx,yy,yz,zx,zy,zz;
	double tx,ty,tz,l;

	double mmener;
	rvec *mmgrd;

	int n = 7;
	int ndim = 9;

	in=fopen("grad.xyz","r");
	if (in == NULL) 
	{
		gmx_fatal(FARGS, "Error while reading Terchem output grad.xyz file, please copy it from terachem scratch directory");
	}

	while (NULL != fgets(buf,300,in))
	{
		if (!strncmp(enerstring, buf, strlen(enerstring)))
		{
			fprintf(stderr, "Reading QM Energy Information... ");
			/* Line contains information about energy in case of single reference methods eg. RHF, DFT */
			pch = strtok(buf, " ,");
			pch = strtok(NULL, " ,");
			pch = strtok(NULL, " ,");
			pch = strtok(NULL, " ,");
			pch = strtok(NULL, " ,");
			pch = strtok(NULL, " ,");
		 	QMener[0] = atof(pch);

			fprintf(stderr, "DONE\n");
			fflush(stderr);

			continue;
		}
	}

	fclose(in);

	in=fopen("tc-qmmm.out","r");
	if (in == NULL) 
	{
		gmx_fatal(FARGS, "Error while reading Terchem output tc-qmmm.out");
	}


	/* Parsing qmmm.out */
	while (NULL != fgets(buf,300,in))
	{
		if (!strncmp(gradstring, buf, strlen(gradstring)))
		{
	        	fprintf(stderr, "Reading QM Gradient Information... ");

			//skip next two lines
		        fgets(buf,300,in);
		        fgets(buf,300,in);

			//Read QM gradients from the next lines
			for(i=0;i<qm->nrQMatoms;i++)
			{
			        fgets(buf,300,in);
               			pch = strtok(buf, " ");
			 	QMgrad[i][XX] = atof(pch);
 				pch = strtok(NULL, " ");
			 	QMgrad[i][YY] = atof(pch);
				pch = strtok(NULL, " ");
			 	QMgrad[i][ZZ] = atof(pch);

//              			fprintf(stderr,"%d %.9lf %.9lf %.9lf\n", i, QMgrad[i][XX], QMgrad[i][YY], QMgrad[i][ZZ]);

			}

			fprintf(stderr, "DONE\n");
			fflush(stderr);

			continue;
		}

		if (!strncmp(gradmmstring, buf, strlen(gradmmstring)))
		{
	        	fprintf(stderr, "Reading MM Gradient Information... ");

			//Read MM gradients from the next lines
			for(i=0;i<mm->nrMMatoms;i++)
			{
			        fgets(buf,300,in);
               			pch = strtok(buf, " ");
			 	MMgrad[i][XX] = atof(pch);
 				pch = strtok(NULL, " ");
			 	MMgrad[i][YY] = atof(pch);
				pch = strtok(NULL, " ");
			 	MMgrad[i][ZZ] = atof(pch);
			}

			fprintf(stderr, "DONE\n");
			fflush(stderr);

			continue;
		}
	}

//	fprintf(stderr,"%.9lf %.9lf %.9lf\n", MMgrad[0][XX], MMgrad[0][YY], MMgrad[0][ZZ]);



	/* Strangely Terachem also put MM-MM point charges interactions into energy and gradient,
           so wee need to cleanup it */
        snew(mmgrd, mm->nrMMatoms);

	for(i = 0; i < mm->nrMMatoms; i++)
	{
		mmgrd[i][XX] = 0.0;
		mmgrd[i][YY] = 0.0;
		mmgrd[i][ZZ] = 0.0;	
	}

        mmener = 0.0;

        //Cycle over all pair of MM point charges
        for(i = 0; i < mm->nrMMatoms-1; i++)
	        for(j = i+1; j < mm->nrMMatoms; j++)
		{
		  //Build vector connecting two charges in Bohr
		  xa = (mm->xMM[j][XX] - mm->xMM[i][XX]) / BOHR2NM;
		  ya = (mm->xMM[j][YY] - mm->xMM[i][YY]) / BOHR2NM;
		  za = (mm->xMM[j][ZZ] - mm->xMM[i][ZZ]) / BOHR2NM;
		  ra = sqrt(xa*xa + ya*ya + za*za);

		  //Contribution to MM-MM energy in Hartree
		  mmener += mm->MMcharges[i] * mm->MMcharges[j] / ra;

		  //Contribution to MM-MM gradient in Hartree/Bohr 
                  //Gradient should be in the i->j direction for i
                  //               and in the j->i direction for j
                  //               with Qi*Qj/rij^2 length in both cases
                 
                  mmgrd[i][XX] += (xa/ra) * mm->MMcharges[i] * mm->MMcharges[j] / (ra*ra);
                  mmgrd[i][YY] += (ya/ra) * mm->MMcharges[i] * mm->MMcharges[j] / (ra*ra);
                  mmgrd[i][ZZ] += (za/ra) * mm->MMcharges[i] * mm->MMcharges[j] / (ra*ra);

                  mmgrd[j][XX] -= (xa/ra) * mm->MMcharges[i] * mm->MMcharges[j] / (ra*ra);
                  mmgrd[j][YY] -= (ya/ra) * mm->MMcharges[i] * mm->MMcharges[j] / (ra*ra);
                  mmgrd[j][ZZ] -= (za/ra) * mm->MMcharges[i] * mm->MMcharges[j] / (ra*ra);
        	}

	//Now substract MM-MM interactions from actuall energy and gradient
	fprintf(stderr,"QMener= %.10lf MMener= %.10lf\n", QMener[0], mmener);
 
	QMener[0] -= mmener;

	for(i = 0; i < mm->nrMMatoms; i++)
	{
		MMgrad[i][XX] -= mmgrd[i][XX];
		MMgrad[i][YY] -= mmgrd[i][YY];
		MMgrad[i][ZZ] -= mmgrd[i][ZZ];
	}

//	fprintf(stderr,"%.9lf %.9lf %.9lf\n", MMgrad[0][XX], MMgrad[0][YY], MMgrad[0][ZZ]);


        sfree(mmgrd);

        fclose(in);
}

int do_terachem(int step,char *exe)
{
	char *buf;

	/* make the call to the molpro binary through system()
	 * The location of the binary will be picked up from the 
	 * environment using getenv().
	* */
	snew(buf,100);
	sprintf(buf,"%s %s",exe,"tc.inp");
	fprintf(stderr,"Calling '%s'\n",buf);
#ifdef GMX_NO_SYSTEM
	printf("Warning-- No calls to system(3) supported on this platform.");
	gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#else
	if ( system(buf) != 0 )	return 1;
#endif
	//sprintf(buf,"copy firefly.out step-%d.out",step);
	//system(buf);
	sfree(buf);
        return 0;
}

real call_msdsh_tc(t_commrec *cr, t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, 
		rvec f[], rvec fshift[])
{ 
	/* a molpor call routine intended for doing diabatic surface
	 * "sliding". See the manual for the theoretical background of this
	 * TSH method.  
	 */
	int
		state,i,j;
	real
		*QMener, *DeltaE, res;
	gmx_bool
		swap=FALSE; /* the actual swap */
	rvec
		*QMgrad,*MMgrad, *f0, *f1, *f2;
	char
		*exe=NULL, buf[1000];

	snew(QMener, qm->nstates);
	snew(QMgrad, qm->nrQMatoms);
	snew(MMgrad, mm->nrMMatoms);

        fprintf(stderr, "--------------------\n");
        fprintf(stderr, "BEGIN OF STEP %d\n", qm->step);
        fprintf(stderr, "--------------------\n");

	/* Random failures prtoection. This tries to launch Terachem second time, if first fails */
	write_msdsh_input_tc(qm->step,qm->swap,fr,qm,mm,0.0);
	if (do_terachem(qm->step,qm->gauss_exe) != 0)
	{
                fprintf(stderr, "First call to Terachem failed, trying once more...\n");
                if (do_terachem(qm->step,qm->gauss_exe) != 0) gmx_fatal(FARGS,"Call to Terachem failed\n");
	}
	read_msdsh_output_tc(QMgrad,MMgrad,qm->step,qm->swap,fr,qm,mm,QMener,TRUE);

	/*output energies */
	fprintf(stderr,"The energes are ");
	for (i=0; i<qm->nstates; i++)
	{
		fprintf(stderr, "%f ", QMener[i]);
	}
	fprintf(stderr, "\n");
       	fprintf(fr->qr->QMoutput, "%i\t%i\t", qm->step, qm->currstate);

	for (i=0; i<qm->nstates; i++)
	{
		fprintf(fr->qr->QMoutput, "%lf\t", QMener[i]);
	}

	fprintf(fr->qr->QMoutput, "\n");
	fflush(fr->qr->QMoutput);

	/* add the QMMM forces to the gmx force array and fshift */
	for(i=0;i<qm->nrQMatoms;i++){
		for(j=0;j<DIM;j++){
			f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
			fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
		}
	}
	for(i=0;i<mm->nrMMatoms;i++){
		for(j=0;j<DIM;j++){
			f[i+qm->nrQMatoms][j]      = HARTREE_BOHR2MD*MMgrad[i][j];
			fshift[i+qm->nrQMatoms][j] = HARTREE_BOHR2MD*MMgrad[i][j];
		}
	}
	for (i=0; i<qm->nstates; i++)
	{
		QMener[i] = QMener[i]*HARTREE2KJ*AVOGADRO;
	}
        fprintf(stderr, "END OF STEP %d\n", qm->step);

	/* increment step */
	qm->step++;

	sfree(MMgrad);
	sfree(QMgrad);
	res = QMener[qm->currstate-1];
	sfree(QMener);

	return res;
} /* call_msdsh */

/* end of molpro sub routines */

