#ifdef WRAD

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include <mpi.h>
#include "atomic_data/Atomic.h"


#ifdef STARS
#include "stars.h"
#endif

// ============================================================================================
int putsource(struct CELL *cell,struct RUNPARAMS *param,int level,REAL aexp, struct OCT *curoct){
  REAL X0=1./pow(2,param->lcoarse);
  int flag;

#ifdef WRADTEST
  // ========================== FOR TESTS ============================
  // =================================================================

  int igrp;
  REAL  dxcur=pow(0.5,curoct->level);
  int icell=cell->idx;
  REAL xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
  REAL yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
  REAL zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;
	      
  
#ifndef TESTCLUMP
  // ============= STROMGREN SPHERE CASE =======================
  if((fabs(xc-0.)<=X0)*(fabs(yc-0.)<=X0)*(fabs(zc-0.)<=X0)){
    if((xc>0.)*(yc>0.)*(zc>0.)){
      cell->rfield.src=param->srcint/pow(X0,3)*param->unit.unit_t/param->unit.unit_n*pow(aexp,2)/8.;///pow(1./16.,3);
      cell->rfieldnew.src=param->srcint/pow(X0,3)*param->unit.unit_t/param->unit.unit_n*pow(aexp,2)/8.;///pow(1./16.,3);
      flag=1;
    }
    else{
      flag=0;
    }
  }
  else{
    cell->rfield.src=0.;
    cell->rfieldnew.src=0.;
    flag=0;
  }
#else
  REAL factgrp[NGRP];
  FACTGRP; //defined in Atomic.h

  cell->rfield.src=0.;
  cell->rfieldnew.src=0.;
  if(xc<=X0){
    for(igrp=0;igrp<NGRP;igrp++){
      curoct->cell[icell].rfield.e[igrp]=factgrp[igrp]*1e10*param->unit.unit_t*pow(param->unit.unit_l,2)/param->clight; 
      curoct->cell[icell].rfield.fx[igrp]=factgrp[igrp]*1e10*param->unit.unit_t*pow(param->unit.unit_l,2);
      curoct->cell[icell].rfield.fy[igrp]=0.; 
      curoct->cell[icell].rfield.fz[igrp]=0.; 
      
      curoct->cell[icell].rfieldnew.e[igrp] =curoct->cell[icell].rfield.e[igrp];
      curoct->cell[icell].rfieldnew.fx[igrp]=curoct->cell[icell].rfield.fx[igrp]; 
      curoct->cell[icell].rfieldnew.fy[igrp]=curoct->cell[icell].rfield.fy[igrp]; 
      curoct->cell[icell].rfieldnew.fz[igrp]=curoct->cell[icell].rfield.fz[igrp]; 
    }
  }
#endif

#else
  // ========================== FOR COSMOLOGY CASES ============================
  // ===========================================================================



#ifdef STARS
	struct PART *nexp;
	struct PART *curp;

	cell->rfield.src   =0.;
	cell->rfieldnew.src=0.;
	flag=0;

	nexp=cell->phead;
	if(nexp==NULL) return 0;
 	do{ 	curp=nexp;
		nexp=curp->next;

		if (curp->isStar){
			if ( (param->cosmo->tphy - curp->age) < param->stars->tlife  ) {
				cell->rfield.src +=  param->srcint/pow(X0,3)*param->unit.unit_t/param->unit.unit_n*pow(aexp,2); // switch to code units 
				flag++;
			}
		}
	}while(nexp!=NULL);

	cell->rfieldnew.src=cell->rfield.src;
#else
#ifdef WHYDRO2
  if((cell->field.d>param->denthresh)&&(cell->rfield.temp<param->tmpthresh)){
    cell->rfield.src=param->srcint*cell->field.d/pow(X0,3)*param->unit.unit_t/param->unit.unit_n*pow(aexp,2); // switch to code units 
    cell->rfieldnew.src=cell->rfield.src;
    flag=1;
  }
  else{
    cell->rfield.src=0.;
    cell->rfieldnew.src=0.;
    flag=0;
  }

#endif
#endif




#endif


  return flag;
}

// ============================================================================================

int FillRad(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, int init, REAL aexp){
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dx;
  int icomp,icell;
  int nread;
  REAL d;
  int nc=0;
  int igrp;
  //Anouk Inputs: aaazzzzzzzzzeeeeeeeeeeeeeeeeeeeeeeeerrrrrtyyuuuuuuuuuklmsqqqqqqqqqqqqqqqwqsqsdfghjklmwxccvb  nn&é
  int flag;

  curoct=firstoct[level-1];
  if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
    nextoct=curoct;
    do{
      curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
      for(icell=0;icell<8;icell++){	

#ifndef TESTCLUMP
	if(curoct->cell[icell].child!=NULL){
	  curoct->cell[icell].rfield.src=0.;
	  continue; // src are built on the finest level
	}
#endif

        flag=putsource(&(curoct->cell[icell]),param,level,aexp,curoct); // creating sources if required
	
	// filling the temperature, nh, and xion
#ifdef WCHEM

#ifdef WRADHYD
	d=curoct->cell[icell].field.d; // baryonic density [unit_mass/unit_lenght^3]
	curoct->cell[icell].rfield.nh=d/(PROTON_MASS/param->unit.unit_mass); // switch to atom/unit_length^3

	
	if (curoct->cell[icell].rfield.nh <= 0) {
		printf("densité negative FILLRAD");
		printf("%e\t%e\t%e\t%e\t%e\t%e\t", curoct->cell[icell].field.d,curoct->cell[icell].field.u,curoct->cell[icell].field.v,curoct->cell[icell].field.w,curoct->cell[icell].field.p,curoct->cell[icell].field.E);		
		abort();	
	}

	curoct->cell[icell].rfieldnew.nh=curoct->cell[icell].rfield.nh;

	curoct->cell[icell].rfield.eint=curoct->cell[icell].field.p/(GAMMA-1.); // 10000 K for a start
	curoct->cell[icell].rfieldnew.eint=curoct->cell[icell].field.p/(GAMMA-1.); 
	curoct->cell[icell].rfieldnew.xion=curoct->cell[icell].field.X;
	curoct->cell[icell].rfield.xion=curoct->cell[icell].field.X;
	/* if((curoct->cell[icell].rfieldnew.eint==0.)) { */
	/*   printf("zero eint\n"); */
	/*   abort(); */
	/* } */
#endif

#endif
	if(init==1){
#ifdef WCHEM

#ifndef COOLING
	  REAL temperature=1e4;
	  REAL xion=1.2e-3;
#else
	  REAL temperature=1e2;
	  REAL xion=1e-5;
#endif
	  REAL eint;

	  curoct->cell[icell].rfield.xion=xion; 
	  curoct->cell[icell].rfieldnew.xion=xion; 
	  
#ifndef WRADHYD
	  // note below the a^5 dependance is modified to a^2 because a^3 is already included in the density
	  eint=(1.5*curoct->cell[icell].rfield.nh*KBOLTZ*(1.+xion)*temperature)*pow(aexp,2)/pow(param->unit.unit_v,2)/param->unit.unit_mass;
	  curoct->cell[icell].rfield.eint=eint; // 10000 K for a start
	  curoct->cell[icell].rfieldnew.eint=eint; // 10000 K for a start
	  E2T(&curoct->cell[icell].rfieldnew,aexp,param);
	  E2T(&curoct->cell[icell].rfield,aexp,param);
#endif

#endif

#ifndef TESTCLUMP
	  for(igrp=0;igrp<NGRP;igrp++){
	    curoct->cell[icell].rfield.e[igrp]=0.+EMIN; 
	    curoct->cell[icell].rfield.fx[igrp]=0.;
	    curoct->cell[icell].rfield.fy[igrp]=0.; 
	    curoct->cell[icell].rfield.fz[igrp]=0.; 

	    curoct->cell[icell].rfieldnew.e[igrp] =curoct->cell[icell].rfield.e[igrp];
	    curoct->cell[icell].rfieldnew.fx[igrp]=curoct->cell[icell].rfield.fx[igrp]; 
	    curoct->cell[icell].rfieldnew.fy[igrp]=curoct->cell[icell].rfield.fy[igrp]; 
	    curoct->cell[icell].rfieldnew.fz[igrp]=curoct->cell[icell].rfield.fz[igrp]; 
	  }
#endif
	}

#ifdef WRADHYD
	E2T(&curoct->cell[icell].rfieldnew,aexp,param);
	E2T(&curoct->cell[icell].rfield,aexp,param);
#endif

	nc+=flag;
      }
    }while(nextoct!=NULL);
  }

#ifdef WMPI
  int nctot;
  MPI_Allreduce(&nc,&nctot,1,MPI_INT,MPI_SUM,cpu->comm);
  nc=nctot;
#endif



  if(cpu->rank==0) printf("== SRC STAT === > Found %d sources \n",nc);
    
  return nc;
}


#endif
