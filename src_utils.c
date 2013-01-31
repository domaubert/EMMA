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
    //cell->rfield.src=param->srcint/pow(X0,3)*param->unit.unit_t/param->unit.unit_n; // switch to code units
    
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
	curoct->cell[icell].rfield.nh=d*(1.-HELIUM_MASSFRACTION)/(PROTON_MASS/param->unit.unit_mass); // switch to atom/unit_length^3
	curoct->cell[icell].rfieldnew.nh=d*(1.-HELIUM_MASSFRACTION)/(PROTON_MASS/param->unit.unit_mass); // switch to atom/unit_length^3
#endif
	if(init==1){
#ifdef WCHEM
	  curoct->cell[icell].rfield.temp=1e4; // 10000 K for a start
	  curoct->cell[icell].rfield.xion=1.2e-3; 
	  curoct->cell[icell].rfieldnew.temp=1e4; // 10000 K for a start
	  curoct->cell[icell].rfieldnew.xion=1.2e-3; 
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


	nc+=flag;
      }
    }while(nextoct!=NULL);
  }

  printf("== SRC STAT === > Found %d sources \n",nc);

  return nc;
}


#endif
#endif
