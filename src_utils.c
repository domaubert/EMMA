#ifdef WRAD

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include <mpi.h>
#include "atomic_data/Atomic.h"


#ifdef WHYDRO2
// ============================================================================================
int putsource(struct CELL *cell,struct RUNPARAMS *param,int level,REAL aexp){
  REAL X0=1./pow(2,param->lcoarse);
  int flag;
  if(cell->field.d>param->srcthresh){
    cell->rfield.src=param->srcint/pow(X0,3)*param->unit.unit_t/param->unit.unit_n*pow(aexp,2); // switch to code units
    cell->rfieldnew.src=cell->rfield.src;
    flag=1;
  }
  else{
    cell->rfield.src=0.;
    flag=0;
  }


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
  int flag;

  curoct=firstoct[level-1];
  if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
    nextoct=curoct;
    do{
      curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
      for(icell=0;icell<8;icell++){	

	if(curoct->cell[icell].child!=NULL){
	  curoct->cell[icell].rfield.src=0.;
	  continue; // src are built on the finest level
	}
	
        flag=putsource(&(curoct->cell[icell]),param,level,aexp); // creating sources if required
	
	// filling the temperature, nh, and xion
#ifdef WCHEM
	d=curoct->cell[icell].field.d; // baryonic density [unit_mass/unit_lenght^3]
	curoct->cell[icell].rfield.nh=d/(PROTON_MASS*MOLECULAR_MU/param->unit.unit_mass); // switch to atom/unit_length^3
	curoct->cell[icell].rfieldnew.nh=d/(PROTON_MASS*MOLECULAR_MU/param->unit.unit_mass); // switch to atom/unit_length^3

#ifdef WRADHYD
	curoct->cell[icell].rfield.eint=curoct->cell[icell].field.p/(GAMMA-1.); // 10000 K for a start
	curoct->cell[icell].rfieldnew.eint=curoct->cell[icell].field.p/(GAMMA-1.); 
#endif

#endif
	if(init==1){
#ifdef WCHEM
	  REAL temperature=1e4;
	  REAL xion=0.2e-3;
	  REAL eint;

	  curoct->cell[icell].rfield.xion=xion; 
	  curoct->cell[icell].rfieldnew.xion=xion; 
	  
#ifndef WRADHYD
	  // note below the a^5 dependance is modified to a^2 because a^3 is already included in the density
	  eint=(1.5*curoct->cell[icell].rfield.nh*KBOLTZ*(1.+xion)*temperature)*pow(aexp,2)/pow(param->unit.unit_v,2)/param->unit.unit_mass;
	  curoct->cell[icell].rfield.eint=eint; // 10000 K for a start
	  curoct->cell[icell].rfieldnew.eint=eint; // 10000 K for a start
#endif

#endif

	  for(igrp=0;igrp<NGRP;igrp++){
	    curoct->cell[icell].rfield.e[igrp]=0.+EMIN; 
	    curoct->cell[icell].rfield.fx[igrp]=0.; 
	    curoct->cell[icell].rfield.fy[igrp]=0.; 
	    curoct->cell[icell].rfield.fz[igrp]=0.; 

	    curoct->cell[icell].rfieldnew.e[igrp]=0.+EMIN; 
	    curoct->cell[icell].rfieldnew.fx[igrp]=0.; 
	    curoct->cell[icell].rfieldnew.fy[igrp]=0.; 
	    curoct->cell[icell].rfieldnew.fz[igrp]=0.; 
	  }
	}

#ifdef WRADHYD
	E2T(&curoct->cell[icell].rfieldnew,aexp,param);
	//printf("temp=%e\n",curoct->cell[icell].rfieldnew.temp);
#endif

	nc+=flag;
      }
    }while(nextoct!=NULL);
  }

  printf("== SRC STAT === > Found %d sources \n",nc);

  return nc;
}


#endif
#endif
