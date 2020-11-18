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


/* TODO: this should be made thread-safe */

/* Molpro interface routines */

void init_msdsh(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
{
	FILE    
		*rffile=NULL,*out=NULL;
	char
		*buf;
	int
		i;

	/* per layer we make a new subdir for integral file, checkpoint
	 * files and such. These dirs are stored in the QMrec for
	 * convenience 
	 */


	if(!qm->nQMcpus){
		for(i=0;i<DIM;i++)
		{
			//qm->SHbasis[i]=basissets[qm->QMbasis][i];
		} 

		/* init gradually switching on of the SA */
		qm->SAstep = 0;

		/* we read the number of cpus and environment from the environment
		 * if set.  
		 */
		snew(buf,1000);
		buf = getenv("NCPUS");
		if (buf)
		{
			sscanf(buf,"%d",&qm->nQMcpus);
		}
		else
		{
			qm->nQMcpus=1;
		}
		fprintf(stderr,"number of CPUs for molpro = %d\n",qm->nQMcpus);

		snew(buf,1000);
		buf = getenv("MEM");
		if (buf)
		{
			sscanf(buf,"%d",&qm->QMmem);
		}
		else
		{
			qm->QMmem=50000000;
		}
		fprintf(stderr,"memory for molpro = %d\n",qm->QMmem);

		snew(buf,1000);
		buf = getenv("ACC");
		if (buf)
		{
			sscanf(buf,"%d",&qm->accuracy);
		}
		else
		{
			qm->accuracy=8;
		}  
		fprintf(stderr,"accuracy in MCSCF = %d\n",qm->accuracy); 

		snew(buf,1000);
		/* molpro always solves CPHF/CPMCSCF */
		qm->cpmcscf = TRUE;
		//fprintf(stderr, "NOT using CP-MCSCF\n");
		/*buf = getenv("CPMCSCF");
		  if (buf)
		  {
		  sscanf(buf,"%d",&i);
		  qm->cpmcscf = (i!=0);
		  }
		  else
		  qm->cpmcscf=FALSE;
		  if (qm->cpmcscf)
		  fprintf(stderr,"using cp-mcscf in l1003\n");
		  else
		  fprintf(stderr,"NOT using cp-mcscf in l1003\n");*/

		snew(buf,1000);
		buf = getenv("SASTEP");
		if (buf)
		{
			sscanf(buf,"%d",&qm->SAstep);
		}
		else
		{
			/* init gradually switching on of the SA */
			qm->SAstep = 0;
		}
		/* we read the number of cpus and environment from the environment
		 * if set.  
		 */
		fprintf(stderr,"Level of SA at start = %d\n",qm->SAstep);

		if (qm->bTS || qm->bOPT)
		{
			gmx_fatal(FARGS, "Error optimization of QM subsystem or transition state not suported with molpro");
		}

		/* molpro settings on the system */
		snew(buf,1000);
		/* recycle the variables for molpro here */
		buf = getenv("MOLPRO_DIR");

		if (buf)
		{
			snew(qm->gauss_dir,1000);
			sscanf(buf,"%s",qm->gauss_dir);
		}
		else
		{
			gmx_fatal(FARGS,"no $MOLPRO_DIR, check molpro manual\n");
		}

		snew(buf,1000);    
		buf = getenv("MOLPRO_EXE");
		if (buf)
		{
			snew(qm->gauss_exe,1000);
			sscanf(buf,"%s",qm->gauss_exe);
		}
		else
		{
			gmx_fatal(FARGS,"no $MOLPRO_EXE, check molpro manual\n");
		}

		qm->step = 0;
		qm->swap = FALSE;
	}

	/* allows the more obvious input for the user:
	 * e.g. to run on the third excited state
	 * startstate = 3
	 * store the current state in c array convention
	 */
	qm->currstate -= 1;

	fprintf(stderr,"molpro initialised...\n");
}  

void write_msdsh_input(int step,gmx_bool swap,
		t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
	int
		i, closed, occ;
	t_QMMMrec
		*QMMMrec;
	char
		periodic_system[37][3]={"XX","H ","He","Li","Be","B ","C ","N ",
			"O ","F ","Ne","Na","Mg","Al","Si","P ",
			"S ","Cl","Ar","K ","Ca","Sc","Ti","V ",
			"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga",
			"Ge","As","Se","Br","Kr"};
	FILE
		*out, *lat;

	QMMMrec = fr->qr;
	out = fopen("molpro.inp","w");
	if (out==NULL) 
	{
		gmx_fatal(FARGS, "Error while opening molpro.inp");
	}
	lat = fopen("molpro.lat","w");
	if (lat==NULL) 
	{
		gmx_fatal(FARGS, "Error while opening molpro.lat");
	}

	fprintf(out, "***, step %i\n", step);
	fprintf(out, "memory,%d\n", qm->QMmem);
	fprintf(out, "file,2,molpro.wfu,new\n");
	fprintf(out, "gthresh,printci=0.0\n");
	fprintf(out, "gprint,civector\n");
	fprintf(out, "geomtyp=xyz\n");
	fprintf(out, "set,zsymel=nosym\n");
	fprintf(out, "set,charge=%2i\n",qm->QMcharge);
	fprintf(out, "geometry={\n");
	fprintf(out, "%d\n\n", qm->nrQMatoms);
	for (i=0;i<qm->nrQMatoms;i++)
	{
#ifdef GMX_DOUBLE
		fprintf(out,"%s,%10.7lf,%10.7lf,%10.7lf\n",
				periodic_system[qm->atomicnumberQM[i]],
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
	fprintf(out, "}\n");
	fprintf(out, "gdirect\n");

	if (qm->QMbasis != eQMbasisINPUT)
	{
		fprintf(out, "basis,%s\n", eQMbasis_names[qm->QMbasis]);
	}
	else
	{
		fprintf(out, "include,molpro.basis\n");
	}

	/* MM point charge data */
	if(QMMMrec->QMMMscheme!=eQMMMschemeoniom && mm->nrMMatoms)
	{
		fprintf(stderr,"nr mm atoms in molpro.c = %d\n",mm->nrMMatoms);
		fprintf(out, "lattice,molpro.lat,molpro.chgrad\n");
		fprintf(lat, "MM point charges from gromacs\n");
		fprintf(lat, "%d\n", mm->nrMMatoms);
		for (i=0;i<mm->nrMMatoms;i++) 
		{
#ifdef GMX_DOUBLE
			fprintf(lat,"%10.7lf,%10.7lf,%10.7lf,%8.4lf,1\n",
					mm->xMM[i][XX]*10,
					mm->xMM[i][YY]*10,
					mm->xMM[i][ZZ]*10,
					mm->MMcharges[i]);
#else
			fprintf(lat,"%10.7f,%10.7f,%10.7f,%8.4f,1\n",
					mm->xMM[i][XX]*10,
					mm->xMM[i][YY]*10,
					mm->xMM[i][ZZ]*10,
					mm->MMcharges[i]);
#endif
		}
	}

	closed = qm->nelectrons/2 - qm->CASelectrons/2;
	occ = closed + qm->CASorbitals;
	fprintf(out, "{matrop\nread,Cstart,orbitals,natural,file=molpro.orb\nsave,Cstart,3000.2\n}\n");
	fprintf(out, "multi;\n");
	fprintf(out, "maxiter,40;\n");
	fprintf(out, "occ,%i;\n", occ);
	fprintf(out, "closed,%i;\n", closed);
	fprintf(out, "config,csf;\n");
	fprintf(out, "wf,%i,1,%i;\n", qm->nelectrons, qm->multiplicity-1);
	fprintf(out, "state,%i;\n", qm->nstates);
	fprintf(out, "start,3000.2,natural;\nnatorb,3001.2;\n");
	fprintf(out, "cpmcscf,grad,%i.1,save=5105.1;\n", qm->currstate+1);
	fprintf(out, "show,energy\n");
	fprintf(out, "{matrop\nload,Cend,orbitals,3001.2,natural\nwrite,Cend,molpro.orb,new\n}\n");

	fprintf(out, "forces\n");
	fprintf(out, "---\n");

	fclose(lat);
	fclose(out);
} /* write_msdsh_SH_input */

void read_msdsh_output(rvec QMgrad[],rvec MMgrad[],int step,
		gmx_bool swapped,t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, real *QMener)
{
	int
		i,j,atnum;
	char
		buf[300];
	char gradstring[] = " Atom          dE/dx               dE/dy               dE/dz";
	char enerstring[19]; 
	char cistring[] = " CI vector";
	char ncsfstring[] = " Number of CSFs:";
	char convfailstring[] = " ** WVFN ****  MAXIMUM NUMBER OF ITERATIONS REACHED";
	char tmpstr[300];
	char *pch;
	FILE
		*in, *chgrad;
	sprintf(enerstring, " ENERGY(1:%i)      =", qm->nstates);

	in=fopen("molpro.out","r");
	if (in == NULL) 
	{
		gmx_fatal(FARGS, "Error while reading molpro ouput");
	}

	while (NULL != fgets(buf,300,in))
	{
		if (!strncmp(convfailstring, buf, 51))
		{
			gmx_fatal(FARGS, "CASSCF did not converge");
		}
		if (!strncmp(gradstring, buf, 60))
		{
			fgets(buf,300,in);
			/* next lines contain the gradients of the QM atoms */
			for(i=0;i<qm->nrQMatoms;i++)
			{
				fgets(buf, 300, in);
#ifdef GMX_DOUBLE
				sscanf(buf,"%d %lf %lf %lf\n",
						&atnum,
						&QMgrad[i][XX],
						&QMgrad[i][YY],
						&QMgrad[i][ZZ]);
#else
				sscanf(buf,"%d %f %f %f\n",
						&atnum, 
						&QMgrad[i][XX],
						&QMgrad[i][YY],
						&QMgrad[i][ZZ]);
#endif     
			} 
		} 
		if (!strncmp(enerstring, buf, 19))
		{
			pch = strtok(buf, " =[]");
			for (i=0; i<qm->nstates; i++) 
			{
				pch = strtok(NULL, " =[]");
				QMener[i] = atof(pch);
			}
		}
		if (!step)
		{
			if (!strncmp(ncsfstring, buf, 16))
			{
				snew(qm->CIvec, qm->nstates);
				snew(qm->CIvecold, qm->nstates);
				sscanf(buf,"%s %s %s %d",tmpstr,tmpstr,tmpstr,&qm->CIdim);
				for (i=0; i<qm->nstates; i++) {
					snew(qm->CIvec[i], qm->CIdim);
					snew(qm->CIvecold[i], qm->CIdim);
				}
			}
		}
		if (!strncmp(cistring, buf, 10))
		{
			fgets(buf,300,in);
			fgets(buf,300,in);
			if (step) {
				for (j=0; j<qm->nstates; j++) {
					for (i=0; i<qm->CIdim; i++) {
						qm->CIvecold[j][i] = qm->CIvec[j][i];
					}
				}
			}
			for (i=0; i<qm->CIdim; i++) {
				fgets(buf,300,in);
				/* skip the configuration info */
				pch = strtok(buf, " ");
				j = 0;
				while ( (pch = strtok(NULL, " ")) != NULL) {
					qm->CIvec[j][i] = atof(pch);
					j++;
				}
			}
		}
	}

	if (mm->nrMMatoms > 0) 
	{
		chgrad = fopen("molpro.chgrad", "r");
		if (chgrad == NULL)
		{
			gmx_fatal(FARGS, "Error while reading molpro gradient ouput");
		}
		fgets(buf,300,chgrad);
		fgets(buf,300,chgrad);
		/* the next lines are the gradients of the MM atoms */
		for(i=0;i<mm->nrMMatoms;i++)
		{
			fgets(buf,300,chgrad);
#ifdef GMX_DOUBLE
			sscanf(buf,"%lf %lf %lf\n",
					&MMgrad[i][XX],
					&MMgrad[i][YY],
					&MMgrad[i][ZZ]);
#else
			sscanf(buf,"%f %f %f\n",
					&MMgrad[i][XX],
					&MMgrad[i][YY],
					&MMgrad[i][ZZ]);
#endif	
		}
		fclose(chgrad);
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
	fclose(in);
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

	if((fabs(d12)>0.5)&&(fabs(d21)>0.5))
		swap = 1;

	return(swap);
}

void do_molpro(int step,char *exe)
{
	char
		buf[100];

	/* make the call to the molpro binary through system()
	 * The location of the binary will be picked up from the 
	 * environment using getenv().
	 */
	if(step) /* hack to prevent long inputfiles */
		/* need to think about another way to handle
		 * the location of the scratch files...
		 */
		sprintf(buf,"%s %s",
				exe,
				"molpro.inp");
	else
		sprintf(buf,"%s %s",
				exe,
				"molpro.inp");
	fprintf(stderr,"Calling '%s'\n",buf);
#ifdef GMX_NO_SYSTEM
	printf("Warning-- No calls to system(3) supported on this platform.");
	gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#else
	if ( system(buf) != 0 )
		gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#endif
}

real call_msdsh(t_commrec *cr, t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, 
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
		*buf=NULL;
	char
		*exe=NULL;

	snew(QMener, qm->nstates);
	snew(exe,1000);
	sprintf(exe,"%s/%s",qm->gauss_dir,qm->gauss_exe);
	/* hack to do ground state simulations */
	if(!qm->step){
		snew(buf,1000);
		buf = getenv("STATE");
		if (buf)
			sscanf(buf,"%d",&state);
		else
			state=2;
		if(state==1)
			qm->swap=TRUE;
	}
	/* end of hack */


	/* copy the QMMMrec pointer */
	snew(QMgrad,qm->nrQMatoms);
	snew(MMgrad,mm->nrMMatoms);
	/* at step 0 there should be no SA */
	/*  if(!step)
	 * qr->bSA=FALSE;*/
	/* temporray set to step + 1, since there is a chk start */
	write_msdsh_input(qm->step,qm->swap,fr,qm,mm);

	do_molpro(qm->step,exe);
	read_msdsh_output(QMgrad,MMgrad,qm->step,qm->swap,fr,qm,mm, QMener);
	fprintf(stderr,"The energes are ");
	for (i=0; i<qm->nstates; i++)
	{
		fprintf(stderr, "%f ", QMener[i]);
	}
	fprintf(stderr, "\n");
	/* hops to the first surface meeting the energy criterion
	 * what if three surfaces come close? => ignored atm
	 */
	for (i=0; i<qm->nstates; i++)
   	{
		if (i != qm->currstate) 
		{
			DeltaE = fabs(QMener[i] - QMener[qm->currstate]);
			if (DeltaE < 0.01)
			{
				swap = (qm->step && hop(qm->step, i, qm->currstate, qm));
				if (swap) 
				{
					fprintf(stderr, "HOP: state %i -> state %i\n", qm->currstate+1, i+1);
					qm->currstate = i;
					break;
				}
			}
		}
	}
	if (swap){/* change surface, so do another call */
		write_msdsh_input(qm->step,qm->swap,fr,qm,mm);
		do_molpro(qm->step,exe);
		read_msdsh_output(QMgrad,MMgrad,qm->step,qm->swap,fr,qm,mm,&DeltaE);
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
	fprintf(stderr,"step %5d, SA = %5d, state = %5d\n",
			qm->step,(qm->SAstep>0),qm->currstate+1);
	qm->step++;
	sfree(MMgrad);
	sfree(QMgrad);
	sfree(exe);
	sfree(buf);
	res = QMener[qm->currstate];
	sfree(QMener);
	return(res);
} /* call_msdsh */

/* end of molpro sub routines */

