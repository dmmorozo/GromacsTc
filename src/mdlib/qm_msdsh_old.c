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

#define PS2AU 41341.3733365586

/* TODO: this should be made thread-safe */

/* Firefly interface routines */


void init_msdsh_ff(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
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
		fprintf(stderr,"number of CPUs for Firefly = %d\n",qm->nQMcpus);
		
		/* Firefly input file */
		snew(buf,1000);    
		buf = getenv("FIREFLY_INP");
		if (buf)
		{
			sscanf(buf,"%s",qm->QMinp);
		}
		else
		{
			sscanf("pck_watbox.inp","%s",qm->QMinp);
		}
		fprintf(stderr,"Firefly input file = %s\n",qm->QMinp);

		/* Firefly batch script*/
		snew(buf,1000);    
		buf = getenv("FIREFLY_EXE");
		if (buf)
		{
			snew(qm->gauss_exe,1000);
			sscanf(buf,"%s",qm->gauss_exe);
			fprintf(stderr,"Firefly Batch sript = %s\n",qm->gauss_exe);
		}
		else
		{
			gmx_fatal(FARGS,"no $FIREFLY_EXE batch script defined\n");
		}

		qm->step = 0;
		qm->swap = FALSE;
	}


	fprintf(stderr,"Firefly initialised...\n");
}  

void write_msdsh_input_ff(int step,gmx_bool swap,
		t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
	int
		i, closed, occ;
	t_QMMMrec
		*QMMMrec;
	char
		*buf, *buf1;
	char
		periodic_system[37][3]={"XX","H ","He","Li","Be","B ","C ","N ",
			"O ","F ","Ne","Na","Mg","Al","Si","P ",
			"S ","Cl","Ar","K ","Ca","Sc","Ti","V ",
			"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga",
			"Ge","As","Se","Br","Kr"};
	FILE
		*out, *lat, *inp;
	char datex;
	
	//checking existence of previous run firefly.dat file
	out = fopen("firefly.dat","r");
	if (out) 
	{
		datex = TRUE;
		fclose(out);
	} else
	{
		datex = FALSE;
	}

	QMMMrec = fr->qr;
	out = fopen("ff.inp","w");
	if (out==NULL) 
	{
		gmx_fatal(FARGS, "Error while trying to write firefly.inp file");
	}
	inp = fopen(qm->QMinp,"r");
	if (inp==NULL) 
	{
		gmx_fatal(FARGS, "qm_msdsh_ff.c", 148, "Error while trying to read %s file", qm->QMinp);
	}

	/* Print out Step Number */
	fprintf(out, "!step %i\n", step);

	/* Now reading input and copy it to firefly.inp, while $data group not found */
	snew(buf,1000);
	snew(buf1,1000);
	while (strncmp(buf," $DATA",6) != 0 && strncmp(buf," $data",6) != 0 && !feof(inp) )
	{
		fgets(buf, 1000, inp);
		fputs(buf, out);
	}

	/* Job Title and Symmetry String */
	fgets(buf, 1000, inp);
	fputs(buf, out);
	fgets(buf, 1000, inp);
	fputs(buf, out);
	fprintf(stderr,"Number of QM Atoms= %d\n",qm->nrQMatoms);
	/* Now write QM System */
	for (i=0;i<qm->nrQMatoms;i++)
	{
#ifdef GMX_DOUBLE
		fprintf(out,"%s   %d.0   %10.7lf   %10.7lf   %10.7lf\n",
				periodic_system[qm->atomicnumberQM[i]],
				qm->atomicnumberQM[i],
				qm->xQM[i][XX]*10,
				qm->xQM[i][YY]*10,
				qm->xQM[i][ZZ]*10);
#else
		fprintf(out,"%s,%10.7f,%10.7f,%10.7f\n",
				periodic_system[qm->atomicnumberQM[i]],
				qm->xQM[i][XX]*10,
				qm->xQM[i][YY]*10,
				qm->xQM[i][ZZ]*10);
#endif
	}
	fprintf(out, " $END\n");

	if (mm->nrMMatoms > 0)
	{
	/* Now write EFP data */
		fprintf(out, " $EFRAG\n");
		fprintf(out, " coord=cart\n");
		fprintf(stderr,"Number of MM Atoms= %d\n",mm->nrMMatoms);
	/*write EFRAG */
		fprintf(out, " FRAGNAME=PROT\n");
		for (i=0; i < 3; i++)
		{
#ifdef GMX_DOUBLE
			fprintf(out," O%d          %15.10lf %15.10lf %15.10lf\n",
                                        i,
					mm->xMM[i][XX]*10,
					mm->xMM[i][YY]*10,
					mm->xMM[i][ZZ]*10);
#else
			fprintf(out," O%d          %15.10f %15.10f %15.10f\n",
                                        i
					mm->xMM[i][XX]*10,
					mm->xMM[i][YY]*10,
					mm->xMM[i][ZZ]*10);
#endif          
		}
		fprintf(out," $END\n");
	/*Now write PROT fragment with charges*/
		fprintf(out," $PROT\n");
		fprintf(out,"TITLE\n");
		fprintf(out," COORDINATES\n");
        	for (i=0; i < mm->nrMMatoms; i++)
		{
#ifdef GMX_DOUBLE
			fprintf(out," O%d          %15.10lf %15.10lf %15.10lf\n",
                                        i,
					mm->xMM[i][XX]/BOHR2NM,
					mm->xMM[i][YY]/BOHR2NM,
					mm->xMM[i][ZZ]/BOHR2NM);
#else
			fprintf(out," O%d          %15.10f %15.10f %15.10f\n",
                                        i
					mm->xMM[i][XX]/BOHR2NM,
					mm->xMM[i][YY]/BOHR2NM,
					mm->xMM[i][ZZ]/BOHR2NM);
#endif          
		}
		fprintf(out," STOP\n");
		fprintf(out," MONOPOLES\n");
        	for (i=0; i < mm->nrMMatoms; i++)
		{
#ifdef GMX_DOUBLE
			fprintf(out," O%d          %15.10lf\n",
                                        i,
					mm->MMcharges[i]);
#else
			fprintf(out," O%d          %15.10f\n",
                                        i
					mm->MMcharges[i]);
		
#endif          
		}

       		fprintf(out," STOP\n");
		fprintf(out," $END\n");
	}
	/*End with coordinates */
	fclose(inp);
	fclose(out);

	/*Do some dirty script job*/
	/* Changing iroot parameter in $det to actual qm->currstate */
	//Windows call
	sprintf(buf,"sed -e 's/iroot=[0-9]*/iroot=%d/' %s > %s", qm->currstate, "ff.inp", "firefly.inp");
	fprintf(stderr,"Calling '%s'\n",buf);
#ifdef GMX_NO_SYSTEM
	printf("Warning-- No calls to system(3) supported on this platform.");
	gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#else
	if ( system(buf) != 0 )
		gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#endif
	//Windows call
	system("rm -f ff.inp");
	/*Copy $vec group*/
	if (qm->step == 0 || !datex) 
	{
		//Windows call
		sprintf(buf,"sed -n '/$VEC/,/$END/p' %s >> %s", qm->QMinp, "firefly.inp");
		fprintf(stderr,"Calling '%s'\n",buf);
#ifdef GMX_NO_SYSTEM
		printf("Warning-- No calls to system(3) supported on this platform.");
		gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#else
		if ( system(buf) != 0 )
		gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#endif
	}else
	{
		sprintf(buf,"sed -n '/--- OPTIMIZED MCSCF MO-S ---/,/$END/p' %s >> %s", "firefly.dat", "firefly.inp");
		fprintf(stderr,"Calling '%s'\n",buf);
#ifdef GMX_NO_SYSTEM
		printf("Warning-- No calls to system(3) supported on this platform.");
		gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#else
		if ( system(buf) != 0 )
		gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#endif
	}
	sfree(buf);
	sfree(buf1);
	fflush(stderr);
} /* write_msdsh_SH_input */

void read_msdsh_output_ff(rvec QMgrad[],rvec MMgrad[],int step,
		gmx_bool swapped,t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, real *QMener)
{
	int
		i,j,atnum,k,tots;
	char
		buf[300], buf1[14];
	char gradstring[] = " $GRAD";
	char cistring[] = " CI vector";
	char ncsfstring[] = " Number of CSFs:";
	char convfailstring[] = " ** WVFN ****  MAXIMUM NUMBER OF ITERATIONS REACHED";
	char tmpstr[300];
	char *pch;
	FILE
		*in, *chgrad;
	int *CSS;
	double fbuf;
	double a[81], b[9], sln[9];
	double xa,ya,za;
	double tf[3], tq[3];
	double r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,rmx,rmy,rmz;
	double r21x,r21y,r21z,r31x,r31y,r31z,r32x,r32y,r32z;
	double dx, dy, dz;
	double xx,xy,xz,yx,yy,yz,zx,zy,zz;
	double tx,ty,tz,l;
	int n = 7;
	int ndim = 9;
	in=fopen("firefly.dat","r");
	if (in == NULL) 
	{
		gmx_fatal(FARGS, "Error while reading Firefly punch file firefly.dat");
	}
	fprintf(stderr, "Reading QM Gradient Information... ");
	while (NULL != fgets(buf,300,in))
	{
		if (!strncmp(gradstring, buf, 6))
		{
			fgets(buf,300,in);
			/* next line contains information about energy in case of single reference methods eg. RHF, DFT
			CASSCF Energies will be red later from output file*/
			if (qm->QMmethod != eQMmethodCASSCF) 
			{
				pch = strtok(buf, " ");
				pch = strtok(NULL, " ");
				QMener[0] = atof(pch);
			}
			/* next lines contain the gradients of the QM atoms */
			for(i=0;i<qm->nrQMatoms;i++)
			{
				fgets(buf, 300, in);
#ifdef GMX_DOUBLE
				sscanf(buf,"%15c %lf %lf %lf",
					    tmpstr,
						&QMgrad[i][XX],
						&QMgrad[i][YY],
						&QMgrad[i][ZZ]);
				
#else
				sscanf(buf,"%15c %f %f %f",
					    tmpstr,
						&QMgrad[i][XX],
						&QMgrad[i][YY],
						&QMgrad[i][ZZ]);
#endif     
			}
		}	
	}
	for(i=0;i<qm->nrQMatoms;i++){
		fprintf(stderr, "%.8lf %.8lf %.8lf\n", QMgrad[i][XX], QMgrad[i][YY], QMgrad[i][ZZ]);
	}
	fprintf(stderr, "DONE\n");
	fclose(in);
	/*Reading gradients on MM atoms/point charges */
	if (mm->nrMMatoms > 0) 
	{
		in=fopen("firefly.out","r");
		if (in == NULL) 
		{
			gmx_fatal(FARGS, "Error while reading Firefly output file firefly.out");
		}
		fprintf(stderr, "Reading Point Charges Gradient Information... ");
                while (!feof(in))
		{
			fgets(buf,300,in);
			while ((strncmp(buf," READY TO DO 1E- PART OF GRAD WRT FRAGMENT", 42) != 0) && (!feof(in)))
			{
				fgets(buf,300,in);
			}
			fgets(buf,300,in);
			for (i=0; i<mm->nrMMatoms; i++)
			{
#ifdef GMX_DOUBLE
				sscanf(buf,"%9c %lf %lf %lf",
					    tmpstr,
						&MMgrad[i][XX],
						&MMgrad[i][YY],
						&MMgrad[i][ZZ]);
				
#else
				sscanf(buf,"%9c %f %f %f",
					    tmpstr,
						&MMgrad[i][XX],
						&MMgrad[i][YY],
						&MMgrad[i][ZZ]);
#endif     
				fgets(buf,300,in);
			}
		}
		fclose(in);

		fprintf(stderr, "DONE\n");
	}

	/* In Case of CASSCF Energies should be read from firefly.out and ci vectors from CIVECTR */
	if (qm->QMmethod == eQMmethodCASSCF) 
	{
		snew(CSS, qm->nstates);
		in=fopen("firefly.out","r");
		if (in == NULL) 
		{
			gmx_fatal(FARGS, "Error while reading Firefly output file firefly.out");
		}
		fprintf(stderr, "Reading energy of CASSCF states... ");
		while (!feof(in))
		{
			fgets(buf,300,in);
			while ((strncmp(buf," THE DENSITIES ARE STATE AVERAGED OVER", 38) != 0) && (!feof(in)))
			{
				fgets(buf,300,in);
			}
                        fgets(buf,300,in);
			for (i=0; i<qm->nstates; i++)
			{
				if (!feof(in))
				{
                                        sscanf(buf,"%7c %d",tmpstr,&CSS[i]);
					sscanf(buf,"%21c %lf",tmpstr,&QMener[i]);
					fgets(buf,300,in);
				}
			}
		}
		if (QMener[0] == 0.0)
		{
			gmx_fatal(FARGS, "CASSCF did not converge");
		}
		fprintf(stderr, "DONE\n");
		fclose(in);
		/*CI vectors should be red from binary file CIVECTR*/
		in=fopen("CIVECTR","rb");
		if (in == NULL) 
		{
			gmx_fatal(FARGS, "Error while reading CIVECTR file. This file should be copied from firefly TMPDIR to WRKDIR.");
		}
		fprintf(stderr, "Reading CI vectors... \n");
       		fprintf(stderr, "Requested pure-spin roots numbers: ");
		for (i=0; i<qm->nstates; i++)
		{
			fprintf(stderr, "%d ", CSS[i]);	
		}
		fprintf(stderr, "\n");
		/*Header*/
		fread(tmpstr, 1, 4, in); 
		/*Full Number of states*/
		fread(&tots,4,1,in);
		/*Number of determinants*/
		fread(tmpstr, 1, 4, in); 
		fread(&qm->CIdim,4,1,in);
		fread(tmpstr, 1, 8, in); 
		fprintf(stderr, "Number of states= %d; Number of determinants in CI= %d\n",tots,qm->CIdim);
		if (!step)
		{
			snew(qm->CIvec, qm->nstates);
			snew(qm->CIvecold, qm->nstates);
			for (i=0; i<qm->nstates; i++) {
				snew(qm->CIvec[i], qm->CIdim);
				snew(qm->CIvecold[i], qm->CIdim);
			}
		}
		if (step) {
			for (j=0; j<qm->nstates; j++) {
				for (i=0; i<qm->CIdim; i++) {
					qm->CIvecold[j][i] = qm->CIvec[j][i];
				}
			}
		}
		i=0;
		for (j=0; j<tots; j++) 
		{
			if (i<qm->nstates) 
			{
				fread(tmpstr, 1, 4, in);	
				fread(qm->CIvec[i], 8	, qm->CIdim, in);
				fread(tmpstr, 1, 4, in);
        	                if (j==CSS[i]-1)
				{
				 i++;
				}
			}
		}
		fprintf(stderr, "DONE\n");
		fclose(in);
                sfree(CSS);
	}

	fprintf(fr->qr->QMoutput, "%i\t%i\t", step, qm->currstate);
#ifdef GMX_DOUBLE
	for (i=0; i<qm->nstates; i++)
	{
		fprintf(fr->qr->QMoutput, "%lf\t", QMener[i]);
	}
#else
	for (i=0; i<qm->nstates; i++)
	{
		fprintf(fr->qr->QMoutput, "%f\t", QMener[i]);
	}
#endif
	fprintf(fr->qr->QMoutput, "\n");
	fflush(fr->qr->QMoutput);
	fflush(stderr);
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

void do_firefly(int step,char *exe)
{
	char *buf;

	/* make the call to the molpro binary through system()
	 * The location of the binary will be picked up from the 
	 * environment using getenv().
	* */
	snew(buf,100);
	sprintf(buf,"%s %s",exe,"firefly");
	fprintf(stderr,"Calling '%s'\n",buf);
#ifdef GMX_NO_SYSTEM
	printf("Warning-- No calls to system(3) supported on this platform.");
	gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#else
	if ( system(buf) != 0 )	gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#endif
	//sprintf(buf,"copy firefly.out step-%d.out",step);
	//system(buf);
	sfree(buf);
}

real call_msdsh_ff(t_commrec *cr, t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, 
		rvec f[], rvec fshift[])
{ 
	/* a molpor call routine intended for doing diabatic surface
	 * "sliding". See the manual for the theoretical background of this
	 * TSH method.  
	 */
	int
		state,i,j;
	real
		*QMener, DeltaE, res;
	gmx_bool
		swap=FALSE; /* the actual swap */
	rvec
		*QMgrad,*MMgrad;
	char
		*exe=NULL;

// File for Landau-Zener probabilities output
	FILE                
		*lzfile=NULL;
	real 
		lzl,lzh,d12=0.0,d21=0.0;


	snew(QMener, qm->nstates);
	snew(exe,1000);
	sprintf(exe,"%s/%s",qm->gauss_dir,qm->gauss_exe);

	/* copy the QMMMrec pointer */
	snew(QMgrad,qm->nrQMatoms);
	snew(MMgrad,mm->nrMMatoms);
	write_msdsh_input_ff(qm->step,qm->swap,fr,qm,mm);
	do_firefly(qm->step,qm->gauss_exe);
	read_msdsh_output_ff(QMgrad,MMgrad,qm->step,qm->swap,fr,qm,mm, QMener);
	fprintf(stderr,"The energes are ");
	for (i=0; i<qm->nstates; i++)
	{
		fprintf(stderr, "%f ", QMener[i]);
	}
	fprintf(stderr, "\n");

	/* Calculation of Landau-Zener probabilities for neigbour states*/
	if (qm->step) {
	//First for higher state 
        	if (qm->currstate < qm->nstates)
		{
			d12 = inproduct(qm->CIvec[qm->currstate],qm->CIvecold[qm->currstate - 1],qm->CIdim);
			lzh = exp(-0.25*M_PI*(QMener[qm->currstate] - QMener[qm->currstate - 1])*qm->delta_t*PS2AU/fabs(d12));
		} else lzh = 0.0;
	//Second for lower state 
        	if (qm->currstate > 1)
		{
			d21 = inproduct(qm->CIvec[qm->currstate - 2],qm->CIvecold[qm->currstate - 1],qm->CIdim);
			lzl = exp(-0.25*M_PI*(QMener[qm->currstate - 1] - QMener[qm->currstate - 2])*qm->delta_t*PS2AU/fabs(d21));
		} else lzl = 0.0;
	//output
		lzfile = fopen("lz.log","a");
		fprintf(lzfile, "%i\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n", qm->step, fabs(d12), fabs(d21), lzh, lzl);
		fclose(lzfile);

	}
	/* hops to the first surface meeting the energy criterion
	 * what if three surfaces come close? => ignored atm
	 */
	for (i=0; i<qm->nstates; i++)
   	{
		if (i != qm->currstate - 1) 
		{
			DeltaE = fabs(QMener[i] - QMener[qm->currstate - 1]);
			fprintf(stderr, "DeltaE= %.5lf\n", DeltaE);
			if (DeltaE < 0.01)
			{
				swap = (qm->step && hop(qm->step, i, qm->currstate - 1, qm));
				if (swap) 
				{
					fprintf(stderr, "HOP: state %i -> state %i\n", qm->currstate, i+1);
					qm->currstate = i + 1;
					break;
				}
			}
		}
	}
	if (swap){/* change surface, so do another call */
		write_msdsh_input_ff(qm->step,qm->swap,fr,qm,mm);
		do_firefly(qm->step,qm->gauss_exe);
		read_msdsh_output_ff(QMgrad,MMgrad,qm->step,qm->swap,fr,qm,mm,QMener);
	}
           
	/* add the QMMM forces to the gmx force array and fshift
	*/
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
	fprintf(stderr,"step %5d, SA = %5d, state = %5d, energy=%.5lf\n",
			qm->step,(qm->SAstep>0),qm->currstate,QMener[qm->currstate-1]);
	qm->step++;
	sfree(MMgrad);
	sfree(QMgrad);
	sfree(exe);
	res = QMener[qm->currstate-1];
	sfree(QMener);
	return(res);
	
} /* call_msdsh */

/* end of molpro sub routines */

