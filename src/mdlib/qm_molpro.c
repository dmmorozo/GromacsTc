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

#ifdef GMX_QMMM_MOLPRO

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
#include "typedefs.h"


/* TODO: this should be made thread-safe */

/* Molpro interface routines */

void init_molpro(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
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
    fprintf(stderr,"molpro initialised...\n");
}  



void write_molpro_SH_input(int step,gmx_bool swap,
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

    fprintf(out, "***, input-file generated by gromacs\n");
    fprintf(out, "memory,%d\n", qm->QMmem);
    fprintf(out, "file,2,molpro.wfu,new\n");
    // XXX 
    // fprintf(out, "gthresh,energy=%d\n", qm->accuracy);
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

    if (qm->QMmethod != eQMmethodCASSCF)
    {
        gmx_fatal(FARGS, "TSH only possible with CASSCF");
    }
    else
    {
        closed = qm->nelectrons/2 - qm->CASelectrons/2;
        occ = closed + qm->CASorbitals;
        fprintf(out, "{matrop\nread,Cstart,orbitals,natural,file=molpro.orb\nsave,Cstart,3000.2\n}\n");
        fprintf(out, "multi;\n");
        fprintf(out, "occ,%i;\n", occ);
        fprintf(out, "closed,%i;\n", closed);
        fprintf(out, "config,csf;\n");
        fprintf(out, "wf,%i,1,%i;\n", qm->nelectrons, qm->multiplicity-1);
        fprintf(out, "state,2;\nweight,0.5,0.5;\nstart,3000.2,natural;\nnatorb,3001.2;\n");
        /* swap == true: back in the ground state */
        if (swap)
        {
            fprintf(out, "cpmcscf,grad,1.1,save=5105.1;\n");
        }
        else
        {
            fprintf(out, "cpmcscf,grad,2.1,save=5105.1;\n");
        }
        fprintf(out, "show,energy\n");
        fprintf(out, "{matrop\nload,Cend,orbitals,3001.2,natural\nwrite,Cend,molpro.orb,new\n}\n");

    }

    fprintf(out, "forces\n");
    fprintf(out, "---\n");

    fclose(lat);
    fclose(out);
} /* write_molpro_SH_input */

void write_molpro_input(int step ,t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
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

    fprintf(out, "***, input-file generated by gromacs\n");
    fprintf(out, "memory,%d\n", qm->QMmem);
    fprintf(out, "file,2,molpro.wfu,new\n");
    // XXX
    //fprintf(out, "gthresh,energy=%d\n", qm->accuracy);
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

	if (step > 0)
	{
		fprintf(out, "{matrop\nread,Cstart,orbitals,file=molpro.orb\nsave,Cstart,3000.2\n}\n");
	}
    switch(qm->QMmethod)
    {
        case eQMmethodAM1:
            gmx_fatal(FARGS, "Error AM1 not available with molpro");
            break;
        case eQMmethodPM3:
            gmx_fatal(FARGS, "Error PM3 not available with molpro");
            break;
        case eQMmethodRHF:
            fprintf(out, "{rhf;wf,%d,1,%d;", qm->nelectrons, qm->multiplicity-1);
			if (step > 0)
			{
				fprintf(out, "start,3000.2;");
			}
			fprintf(out, "save,3001.2}\n");
            break;
        case eQMmethodUHF:
            fprintf(out, "{uhf;wf,%d,1,%d;", qm->nelectrons, qm->multiplicity-1);
			if (step > 0)
			{
				fprintf(out, "start,3000.2;");
			}
			fprintf(out, "save,3001.2}\n");
            break;
        case eQMmethodDFT:
            gmx_fatal(FARGS, "Error only B3LYP available with molpro");
            break;
        case eQMmethodB3LYP:
            fprintf(out, "df=[b3lyp]\n");
            if (qm->multiplicity == 1) 
            {
                fprintf(out, "{rks;wf,%d,1,%d;", qm->nelectrons, qm->multiplicity-1);
            }
            else
            {
                fprintf(out, "{uks;wf,%d,1,%d;", qm->nelectrons, qm->multiplicity-1);
            }
			if (step > 0)
			{
				fprintf(out, "start,3000.2;");
			}
			fprintf(out, "save,3001.2}\n");
            break; 
        case eQMmethodCASSCF:
            if (qm->nelectrons%2 == 1 || qm->CASelectrons%2 == 1) 
            {
                gmx_fatal(FARGS, "Error only even numbers of electrons and CASelectrons are supported. Use QMmethod INPUT instead");
            }
            closed = qm->nelectrons/2 - qm->CASelectrons/2;
            occ = closed + qm->CASorbitals;
            fprintf(out, "{casscf;frozen,0;occ,%i;closed,%i;wf,%i,1,%i", occ, closed, qm->nelectrons, qm->multiplicity-1);
			if (step > 0)
			{
				fprintf(out, "start,3000.2;");
			}
			fprintf(out, "natorb,3001.2}\n");
            /* add cpmscf gradients here, otherwise molpro crashes */
            break;
        case eQMmethodINPUT:
            fprintf(out, "include,molpro.method\n");
            break;
        default:
            gmx_fatal(FARGS, "Error illegal QMmethod name\n");
            break;
    }
	fprintf(out, "{matrop\nload,Cend,orbitals,3001.2\nwrite,Cend,molpro.orb,new\n}\n");
    fprintf(out, "forces\n");
    fprintf(out, "show,energy\n");
    fprintf(out, "---\n");

    fclose(lat);
    fclose(out);
}  /* write_molpro_input */

real read_molpro_output(rvec QMgrad[],rvec MMgrad[],int step,
        t_QMrec *qm, t_MMrec *mm)
{
    int
        i,j,atnum;
    char
        buf[300];
    char gradstring[] = " Atom          dE/dx               dE/dy               dE/dz";
    char enerstring[] = " ENERGY           =";
    char tmpstr[300];
    real
        QMener;
    FILE
        *in, *chgrad;

    in=fopen("molpro.out","r");
    if (in == NULL) 
    {
        gmx_fatal(FARGS, "Error while reading molpro ouput");
    }

    while (NULL != fgets(buf,300,in))
    {
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
#ifdef GMX_DOUBLE
            sscanf(buf,"%s %s %lf %s\n",tmpstr, tmpstr, &QMener, tmpstr);
#else
            sscanf(buf,"%s %s %f %s\n",tmpstr, tmpstr, &QMener, tmpstr);
#endif
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
    fclose(in);
    return(QMener);  
}

real read_molpro_SH_output(rvec QMgrad[],rvec MMgrad[],int step,
        gmx_bool swapped,t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, real *DE)
{
    int
        i,j,atnum;
    char
        buf[300];
    char gradstring[] = " Atom          dE/dx               dE/dy               dE/dz";
    char enerstring[] = " ENERGY(1:2)      =";
    char cistring[] = " CI vector";
    char ncsfstring[] = " Number of CSFs:";
    char tmpstr[300];
    real
        QMener1, QMener2, DeltaE;
    FILE
        *in, *chgrad;

    in=fopen("molpro.out","r");
    if (in == NULL) 
    {
        gmx_fatal(FARGS, "Error while reading molpro ouput");
    }

    while (NULL != fgets(buf,300,in))
    {
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
            /* the multi programm writes HARTREE instead of AU
             * should work anyway...
             * */
#ifdef GMX_DOUBLE
            sscanf(buf,"%s %s %s %lf %lf%s %s\n", tmpstr, tmpstr, tmpstr, &QMener1, &QMener2, tmpstr, tmpstr);
#else
            sscanf(buf,"%s %s %s %f %f%s %s\n", tmpstr, tmpstr, tmpstr, &QMener1, &QMener2, tmpstr, tmpstr);
#endif
        }
        if (!step)
        {
            if (!strncmp(ncsfstring, buf, 16))
            {
                sscanf(buf,"%s %s %s %d",tmpstr,tmpstr,tmpstr,&qm->CIdim);
                snew(qm->CIvec1,qm->CIdim);
                snew(qm->CIvec1old,qm->CIdim);
                snew(qm->CIvec2,qm->CIdim);
                snew(qm->CIvec2old,qm->CIdim);
            }
        }
        if (!strncmp(cistring, buf, 10))
        {
            fgets(buf,300,in);
            fgets(buf,300,in);
            if (step) {
                for (i=0; i<qm->CIdim; i++) {
                    qm->CIvec1old[i] = qm->CIvec1[i];
                    qm->CIvec2old[i] = qm->CIvec2[i];
                }
            }
            for (i=0; i<qm->CIdim; i++) {
                fgets(buf,300,in);
#ifdef GMX_DOUBLE
                sscanf(buf, "%s %lf %lf\n", tmpstr, &qm->CIvec1[i], &qm->CIvec2[i]);
#else
                sscanf(buf, "%s %f %f\n", tmpstr, &qm->CIvec1[i], &qm->CIvec2[i]);
#endif
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

#ifdef GMX_DOUBLE
    fprintf(fr->qr->QMoutput, "%i\t%i\t%lf\n", step, swapped, fabs(QMener1-QMener2));
#else
    fprintf(fr->qr->QMoutput, "%i\t%i\t%f\n", step, swapped, fabsf(QMener1-QMener2));
#endif
    fflush(fr->qr->QMoutput);
    fclose(in);
    *DE = QMener2-QMener1;
    if (swapped)
    {
        return QMener1;
    }
    return(QMener2);  
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

int hop(int step, t_QMrec *qm)
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
        d11 = inproduct(qm->CIvec1,qm->CIvec1old,qm->CIdim);
        d12 = inproduct(qm->CIvec1,qm->CIvec2old,qm->CIdim);
        d21 = inproduct(qm->CIvec2,qm->CIvec1old,qm->CIdim);
        d22 = inproduct(qm->CIvec2,qm->CIvec2old,qm->CIdim);
    }
    fprintf(stderr,"-------------------\n");
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

real call_molpro(t_commrec *cr,  t_forcerec *fr, 
        t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
    /* normal molpro jobs */
    int
        i,j;
    real
        QMener=0.0;
    rvec
        *QMgrad,*MMgrad;
    char
        *exe;

    snew(exe,1000);
    sprintf(exe,"%s/%s",qm->gauss_dir,qm->gauss_exe);
    snew(QMgrad,qm->nrQMatoms);
    snew(MMgrad,mm->nrMMatoms);

    write_molpro_input(qm->step,fr,qm,mm);
    do_molpro(qm->step,exe);
    QMener = read_molpro_output(QMgrad,MMgrad,qm->step,qm,mm);
    /* put the QMMM forces in the force array and to the fshift
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
    QMener = QMener*HARTREE2KJ*AVOGADRO;
    qm->step++;
    sfree(QMgrad);
    sfree(MMgrad);
    sfree(exe);
    return(QMener);
} /* call_molpro */

real call_molpro_SH(t_commrec *cr, t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, 
        rvec f[], rvec fshift[])
{ 
    /* a molpor call routine intended for doing diabatic surface
     * "sliding". See the manual for the theoretical background of this
     * TSH method.  
     */
    int
        state,i,j;
    real
        QMener=0.0, DeltaE=0.0;
    gmx_bool
        swap=FALSE; /* the actual swap */
    rvec
        *QMgrad,*MMgrad;
    char
        *buf=NULL;
    char
        *exe=NULL;

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
    write_molpro_SH_input(qm->step,qm->swap,fr,qm,mm);

    do_molpro(qm->step,exe);
    QMener = read_molpro_SH_output(QMgrad,MMgrad,qm->step,qm->swap,fr,qm,mm, &DeltaE);
    fprintf(stderr,"The energy difference is %f\n",DeltaE);
    if  (DeltaE < 0.01 ){
        if(!qm->swap){
            swap    = (qm->step && hop(qm->step,qm));
            qm->swap = swap;
        } 
        else { /* already on the other surface, so check if we go back */
            swap    = (qm->step && hop(qm->step,qm));
            qm->swap =!swap; /* so qm->swap shoud be false again */
        }
        if (swap){/* change surface, so do another call */
            write_molpro_SH_input(qm->step,qm->swap,fr,qm,mm);
            do_molpro(qm->step,exe);
            QMener = read_molpro_SH_output(QMgrad,MMgrad,qm->step,qm->swap,fr,qm,mm,&DeltaE);
        }
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
    QMener = QMener*HARTREE2KJ*AVOGADRO;
    fprintf(stderr,"step %5d, SA = %5d, swap = %5d\n",
            qm->step,(qm->SAstep>0),qm->swap);
    qm->step++;
    sfree(MMgrad);
    sfree(QMgrad);
    sfree(exe);
    sfree(buf);
    return(QMener);
} /* call_molpro_SH */

/* end of molpro sub routines */

#else
int
gmx_qmmm_molpro_empty;
#endif

