#ifdef WRAD

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include <mpi.h>
#include "atomic_data/Atomic.h"
#include "stars.h"


#ifdef COARSERAD
#ifdef STARS

void collectstars(struct CELL *cellcoarse, struct CELL *cell,struct RUNPARAMS *param, REAL dxcoarse, REAL aexp, REAL tcur, int *ns){

  REAL srcint = param->srcint;
  struct PART *nexp;
  struct PART *curp;
  
  // found a root
  if(cell->child==NULL){
    nexp=cell->phead;
    if(nexp!=NULL){
      do{ 	
	curp=nexp;
	nexp=curp->next;
	
	if ((curp->isStar)){
	  (*ns)++;
	  //printf("coucou %e %e %e\n",tcur,curp->age,param->stars->tlife);
	  if(tcur>= curp->age){ // for inter-level communications
	    //printf("hehey\n");
	    if ( (tcur - curp->age) < param->stars->tlife  ) {
	      //printf("hihi %e\n",dxcoarse);
	      cellcoarse->rfield.src +=  (curp->mass*param->unit.unit_mass)*srcint/pow(dxcoarse,3)*param->unit.unit_t/param->unit.unit_n*pow(aexp,2); // switch to code units 
	    }
	  }
	}
      }while(nexp!=NULL);
	
    }
  }
  else{
    // recursive call
    int icell;
    for(icell=0;icell<8;icell++){
      collectstars(cellcoarse,&(cell->child->cell[icell]),param,dxcoarse,aexp,tcur,ns);
    }
  }
}


void collectE(struct CELL *cellcoarse, struct CELL *cell,struct RUNPARAMS *param,int level, REAL aexp, REAL tcur, int *ns){

  REAL fact=pow(2.,3*(param->lcoarse-level));
  // found a root
  if(cell->child==NULL){
    int igrp;
    for(igrp=0;igrp<NGRP;igrp++){
      cellcoarse->rfield.e[igrp]   +=cell->rfield.e[igrp] *fact;
      cellcoarse->rfield.fx[igrp]  +=cell->rfield.fx[igrp]*fact;
      cellcoarse->rfield.fy[igrp]  +=cell->rfield.fy[igrp]*fact;
      cellcoarse->rfield.fz[igrp]  +=cell->rfield.fz[igrp]*fact;
    }
  }
  else{
    // recursive call
    int icell;
    for(icell=0;icell<8;icell++){
      collectE(cellcoarse,&(cell->child->cell[icell]),param,level+1,aexp,tcur,ns);
      (*ns)++;
    }
  }
}


// --------------------------------------------------------------------

int putsource2coarse(struct CELL *cell,struct RUNPARAMS *param,int level,REAL aexp, REAL tcur){
  
  
  cell->rfield.src   =0.;
  cell->rfieldnew.src=0.;
  

  REAL dxcoarse=pow(0.5,level);
  int ns=0;
  collectstars(cell,cell,param,dxcoarse,aexp,tcur,&ns);

  cell->rfieldnew.src=cell->rfield.src;

  return ns;
}

// --------------------------------------------------------------------

int putE2coarse(struct CELL *cell,struct RUNPARAMS *param,int level,REAL aexp, REAL tcur){
  
  
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    cell->rfield.e[igrp]   =0.;
    cell->rfield.fx[igrp]   =0.;
    cell->rfield.fy[igrp]   =0.;
    cell->rfield.fz[igrp]   =0.;
  }
  
  int ns=0;
  collectE(cell,cell,param,level,aexp,tcur,&ns);

  for(igrp=0;igrp<NGRP;igrp++){
    cell->rfieldnew.e[igrp]    =cell->rfield.e[igrp]  ;
    cell->rfieldnew.fx[igrp]   =cell->rfield.fx[igrp] ;
    cell->rfieldnew.fy[igrp]   =cell->rfield.fy[igrp] ;
    cell->rfieldnew.fz[igrp]   =cell->rfield.fz[igrp] ;

  }

  return ns;
}


#endif
#endif

// ============================================================================================
int putsource(struct CELL *cell,struct RUNPARAMS *param,int level,REAL aexp, REAL tcur, struct OCT *curoct,  struct CPUINFO *cpu){
  REAL X0=1./pow(2,param->lcoarse);
  REAL  dxcur=pow(0.5,curoct->level);
  int flag;

#ifdef WRADTEST
  // ========================== FOR TESTS ============================
  // =================================================================

  int igrp;
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

  REAL srcint = param->srcint;
  
  struct PART *nexp;
  struct PART *curp;
  
  cell->rfield.src   =0.;
  cell->rfieldnew.src=0.;
  flag=0;
  
  nexp=cell->phead;
  if(nexp==NULL) return 0;
  int nss=0;
  do{ 	curp=nexp;
    nexp=curp->next;
    
    if ((curp->isStar)){
      nss++;
      //printf("star found ! t=%e age=%e agelim=%e idx=%d\n",tcur,curp->age,param->stars->tlife,curp->idx);
      if(tcur>= curp->age){ // for inter-level communications
	if ( (tcur - curp->age) < param->stars->tlife  ) {
	  cell->rfield.src +=  (curp->mass*param->unit.unit_mass)*srcint/pow(dxcur,3)*param->unit.unit_t/param->unit.unit_n*pow(aexp,2); // switch to code units 
	  flag=1;
	  
	}
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
  REAL tcur=a2t(param,aexp);

  //if(cpu->rank==RANK_DISP) printf("Building Source field\n");

  curoct=firstoct[level-1];
  if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
    nextoct=curoct;
    do{
      curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
      for(icell=0;icell<8;icell++){	


#ifdef COARSERAD
	// if RT is restricted to coarse grid
	// gather E and F from the mother cell
	/* if(level>param->lcoarse){ */

	/*   putcoarse2E(&(curoct->cell[icell]),param,level,aexp,tcur); */
	/*   for(igrp=0;igrp<NGRP;igrp++){ */
	/*     curoct->cell[icell].rfield.e[igrp]=curoct->parent->rfield.e[igrp]; */
	/*     curoct->cell[icell].rfield.fx[igrp]=curoct->parent->rfield.fx[igrp]; */
	/*     curoct->cell[icell].rfield.fy[igrp]=curoct->parent->rfield.fy[igrp]; */
	/*     curoct->cell[icell].rfield.fz[igrp]=curoct->parent->rfield.fz[igrp]; */
	/*     curoct->cell[icell].rfieldnew.e[igrp]=curoct->parent->rfield.e[igrp]; */
	/*     curoct->cell[icell].rfieldnew.fx[igrp]=curoct->parent->rfield.fx[igrp]; */
	/*     curoct->cell[icell].rfieldnew.fy[igrp]=curoct->parent->rfield.fy[igrp]; */
	/*     curoct->cell[icell].rfieldnew.fz[igrp]=curoct->parent->rfield.fz[igrp]; */
	/*   } */
	/* } */
	
	if(level==param->lcoarse){
	  // WE pump up the high resolution conserved RT quantities (if the cell is split)
	  if(curoct->cell[icell].child!=NULL){
	    int nE;
	    nE=putE2coarse(&(curoct->cell[icell]),param,level,aexp,tcur);
	    /* if(isnan(curoct->cell[icell].rfieldnew.e[0])) */
	    /*   { */
	    /* 	printf("WTF\n"); */
	    /* 	abort(); */
	    /*   } */
	    //printf("PUMPING UP FROM %d CELLs\n",nE);
	  }
	}
#endif


	// filling the temperature, nh, and xion
#ifdef WCHEM
#ifdef WRADHYD
	d=curoct->cell[icell].field.d; // baryonic density [unit_mass/unit_lenght^3]
	curoct->cell[icell].rfield.nh=d/(PROTON_MASS/param->unit.unit_mass); // switch to atom/unit_length^3
	curoct->cell[icell].rfieldnew.nh=curoct->cell[icell].rfield.nh;
	curoct->cell[icell].rfield.eint=curoct->cell[icell].field.p/(GAMMA-1.); // 10000 K for a start
	curoct->cell[icell].rfieldnew.eint=curoct->cell[icell].field.p/(GAMMA-1.); 
	curoct->cell[icell].rfieldnew.nhplus=curoct->cell[icell].field.dX/(PROTON_MASS/param->unit.unit_mass);
	curoct->cell[icell].rfield.nhplus=curoct->cell[icell].field.dX/(PROTON_MASS/param->unit.unit_mass);
	/* if(curoct->cell[icell].rfield.nhplus>curoct->cell[icell].rfield.nh){ */
	/*   printf("MEGA WTF \n"); */
	/* } */
#endif
#endif

#ifndef TESTCLUMP
	if(curoct->cell[icell].child!=NULL){
	  curoct->cell[icell].rfield.src=0.;
	  continue; // src are built on the finest level
	}
#endif

        flag=putsource(&(curoct->cell[icell]),param,level,aexp,tcur,curoct,cpu); // creating sources if required
	
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

	  curoct->cell[icell].rfield.nhplus=xion*curoct->cell[icell].rfield.nh; 
	  curoct->cell[icell].rfieldnew.nhplus=xion*curoct->cell[icell].rfieldnew.nh; 
	  
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



  //if(cpu->rank==RANK_DISP) printf("== SRC STAT === > Found %d sources \n",nc);
  return nc;
}


#endif
