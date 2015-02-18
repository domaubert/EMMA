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
#include "chem_utils.h"

#ifdef COARSERAD
#ifdef STARS
#ifndef SNTEST
void collectstars(struct CELL *cellcoarse, struct CELL *cell,struct RUNPARAMS *param, REAL dxcoarse, REAL aexp, REAL tcur, int *ns){
#ifdef PIC
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
	      cellcoarse->rfield.src +=  (curp->mass*param->unit.unit_mass)*srcint/POW(dxcoarse*param->unit.unit_l,3)*param->unit.unit_t/param->unit.unit_N*POW(aexp,2); // switch to code units
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
#endif // PIC
}




// --------------------------------------------------------------------

int putsource2coarse(struct CELL *cell,struct RUNPARAMS *param,int level,REAL aexp, REAL tcur){


  cell->rfield.src   =0.;
  cell->rfieldnew.src=0.;


  REAL dxcoarse=POW(0.5,level);
  int ns=0;
  collectstars(cell,cell,param,dxcoarse,aexp,tcur,&ns);

  cell->rfieldnew.src=cell->rfield.src;

  return ns;
}

#endif //SNTEST
#endif //STARS

// --------------------------------------------------------------------

void collectE(struct CELL *cellcoarse, struct CELL *cell,struct RUNPARAMS *param,int level, REAL aexp, REAL tcur, int *ns){

  REAL fact=POW(2.,3*(param->lcoarse-level));
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


#endif //COARSERAD

// ============================================================================================
int putsource(struct CELL *cell,struct RUNPARAMS *param,int level,REAL aexp, REAL tcur, struct OCT *curoct,  struct CPUINFO *cpu){
  REAL X0=1./POW(2,param->lcoarse);
  REAL  dxcur=POW(0.5,curoct->level);
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

int lifetime_test = 1;
#ifdef SNTEST
REAL tcur_in_yrs = tcur*param->unit.unit_t/MYR *1e6;
if ( tcur_in_yrs >= LIFETIME_OF_STARS_IN_TEST) lifetime_test = 0;
#endif // SNTEST

  //if((FABS(xc-0.5)<=X0)*(FABS(yc-0.5)<=X0)*(FABS(zc-0.5)<=X0) && lifetime_test){
  if(curoct->x==0.5 && curoct->y==0.5 && curoct->z==0.5 && icell==0 && lifetime_test){
    if((xc>0.)*(yc>0.)*(zc>0.)){
      //cell->rfield.src=param->srcint/POW(X0,3)*param->unit.unit_t/param->unit.unit_n*POW(aexp,2)/8.;///8.;///8.;///POW(1./16.,3);
      cell->rfield.src=param->srcint/POW(X0*param->unit.unit_l,3)*param->unit.unit_t/param->unit.unit_N*POW(aexp,2);///8.;///8.;///POW(1./16.,3);
      cell->rfieldnew.src=cell->rfield.src;
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
  // ============= END STROMGREN SPHERE CASE =======================

#else //ifdef TESTCLUMP
  REAL factgrp[NGRP];
  FACTGRP; //defined in Atomic.h

  cell->rfield.src=0.;
  cell->rfieldnew.src=0.;
  flag=0;
  if(fabs(xc)<=X0){
    for(igrp=0;igrp<NGRP;igrp++){
      curoct->cell[icell].rfield.e[igrp]=factgrp[igrp]*1e10*param->unit.unit_t/param->unit.unit_l/param->unit.unit_N/param->clight;
      curoct->cell[icell].rfield.fx[igrp]=factgrp[igrp]*1e10*param->unit.unit_t/param->unit.unit_l/param->unit.unit_N;
      curoct->cell[icell].rfield.fy[igrp]=0.;
      curoct->cell[icell].rfield.fz[igrp]=0.;
      curoct->cell[icell].rfieldnew.e[igrp] =curoct->cell[icell].rfield.e[igrp];
      curoct->cell[icell].rfieldnew.fx[igrp]=curoct->cell[icell].rfield.fx[igrp];
      curoct->cell[icell].rfieldnew.fy[igrp]=curoct->cell[icell].rfield.fy[igrp];
      curoct->cell[icell].rfieldnew.fz[igrp]=curoct->cell[icell].rfield.fz[igrp];
    }
  }
#endif //TESTCLUMP

#else //ifndef WRADTEST
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

    if ((curp->isStar==1)){
      nss++;
      //printf("star found ! t=%e age=%e agelim=%e idx=%d COUNT=%d\n",tcur,curp->age,param->stars->tlife,curp->idx,(( (tcur - curp->age) < param->stars->tlife  )&&(tcur>= curp->age)));
      cell->rfield.src +=  (curp->mass*param->unit.unit_d)*srcint/POW(dxcur,3)*param->unit.unit_t/param->unit.unit_N*POW(aexp,2); // switch to code units
      //printf("SRC= %e\n",cell->rfield.src);
      flag=1;
    }

#ifdef DECREASE_EMMISIVITY_AFTER_TLIFE
    if (curp->isStar==3 || curp->isStar==4){ //Supernovae + decreasing luminosity OR decreasing luminosity
      nss++;
      REAL slope = -4.; //slope derived from http://www.stsci.edu/science/starburst99/figs/fig77.html
      //printf("star found ! t=%e age=%e agelim=%e idx=%d COUNT=%d\n",tcur,curp->age,param->stars->tlife,curp->idx,(( (tcur - curp->age) < param->stars->tlife  )&&(tcur>= curp->age)));
      REAL src = (curp->mass*param->unit.unit_d)*srcint/POW(dxcur,3)*param->unit.unit_t/param->unit.unit_N*POW(aexp,2);
      REAL t = (param->cosmo->tphy - curp->age) / param->stars->tlife;
      cell->rfield.src +=  src*POW(t,slope);;
      //printf("SRC= %e\n",cell->rfield.src);
      flag=1;
    }
#endif

  }while(nexp!=NULL);

  cell->rfieldnew.src=cell->rfield.src;

#else //ifndef STARS
#ifdef WHYDRO2
  if((cell->field.d>param->denthresh)&&(cell->rfield.temp<param->tmpthresh)){
    cell->rfield.src=param->srcint*cell->field.d/POW(X0*param->unit.unit_l,3)*param->unit.unit_t/param->unit.unit_N*POW(aexp,2); // switch to code units
    cell->rfieldnew.src=cell->rfield.src;
    flag=1;
  }
  else{
    cell->rfield.src=0.;
    cell->rfieldnew.src=0.;
    flag=0;
  }

#endif //WHYDRO
#endif //STARS
#endif

  return flag;
}

// ==================================================================================================================
#ifdef HOMOSOURCE
void homosource(struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, int levext){
  struct OCT *nextoct;
  struct OCT *curoct;
  int level;
  int icell;
  REAL bkg=0.;
  int N=0;

  for(level=param->lcoarse;level<=param->lmax;level++){
    REAL  dx3=POW(0.5,3*level);
    curoct=firstoct[level-1];
    if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
      nextoct=curoct;
      do{
	curoct=nextoct;
	nextoct=curoct->next;
	if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
	for(icell=0;icell<8;icell++){
	  if(curoct->cell[icell].child==NULL){
	    bkg+=curoct->cell[icell].rfield.src*dx3;
	    if(curoct->cell[icell].rfield.src>0) {
	      N++;
	      //printf("lev=%d src=%e N=%d cpu=%d\n",level,curoct->cell[icell].rfield.src,N,cpu->rank);
	    }
	  }
	}
      }while(nextoct!=NULL);
    }
  }

#ifdef WMPI
  REAL bkgtot;
  MPI_Allreduce(&bkg,&bkgtot,1,MPI_REEL,MPI_SUM,cpu->comm);
  bkg=bkgtot;
  int Ntot;
  MPI_Allreduce(&N,&Ntot,1,MPI_INT,MPI_SUM,cpu->comm);
  N=Ntot;
#endif

  //if(cpu->rank==RANK_DISP) printf("Call lev=%d Homogeneous field found = %e with %d src\n",levext,bkg,N);
  param->bkg=bkg; // saving the value

}
#endif


// ==================================================================================================================
void cleansource(struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu){
  struct OCT *nextoct;
  struct OCT *curoct;
  int level;
  int icell;
  REAL bkg=0.;
  int N=0;

  for(level=param->lcoarse;level<param->lmax;level++){
    curoct=firstoct[level-1];
    if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
      nextoct=curoct;
      do{
	curoct=nextoct;
	nextoct=curoct->next;
	if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
	for(icell=0;icell<8;icell++){
	  curoct->cell[icell].rfield.src=0.;
	  curoct->cell[icell].rfieldnew.src=0.;
	}
      }while(nextoct!=NULL);
    }
  }

}


// ============================================================================================

int FillRad(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, int init, REAL aexp, REAL tloc){
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

#ifdef TESTCOSMO
  REAL tcur=a2t(param,aexp);
#else
  REAL tcur=tloc; // PROBABLY NEEDS TO BE FIXED FOR NON COSMO RUN WITH STARS (LATER....)
#endif
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
	if(level==param->lcoarse){
	  // WE pump up the high resolution conserved RT quantities (if the cell is split)
	  if(curoct->cell[icell].child!=NULL){
	    int nE;
	    nE=putE2coarse(&(curoct->cell[icell]),param,level,aexp,tcur);
	  }
	}
#endif


	// filling the temperature, nh, and xion
#ifdef WCHEM
#ifdef WRADHYD
	d=curoct->cell[icell].field.d; // baryonic density [unit_mass/unit_lenght^3]
	//curoct->cell[icell].rfield.nh=d/(PROTON_MASS*MOLECULAR_MU)*param->unit.unit_d; // switch to atom/m^3
	curoct->cell[icell].rfield.nh=d; // [unit_N] note d in unit_d and nh in unit_N are identical
	curoct->cell[icell].rfieldnew.nh=curoct->cell[icell].rfield.nh;
	curoct->cell[icell].rfield.eint=curoct->cell[icell].field.p/(GAMMA-1.); // 10000 K for a start
	curoct->cell[icell].rfieldnew.eint=curoct->cell[icell].field.p/(GAMMA-1.);
 	/* curoct->cell[icell].rfieldnew.nhplus=curoct->cell[icell].field.dX/(PROTON_MASS*MOLECULAR_MU)*param->unit.unit_d; */
	/* curoct->cell[icell].rfield.nhplus=curoct->cell[icell].field.dX/(PROTON_MASS*MOLECULAR_MU)*param->unit.unit_d; */
	curoct->cell[icell].rfield.nhplus=curoct->cell[icell].field.dX; // [unit_N] note d in unit_d and nh in unit_N are identical
	curoct->cell[icell].rfieldnew.nhplus=curoct->cell[icell].rfield.nhplus; // [unit_N] note d in unit_d and nh in unit_N are identical

#endif
#endif

#ifndef TESTCLUMP
	if(curoct->cell[icell].child!=NULL){
	  curoct->cell[icell].rfield.src=0.;
	  curoct->cell[icell].rfieldnew.src=0.;
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
	  eint=(1.5*curoct->cell[icell].rfield.nh*KBOLTZ*(1.+xion)*temperature)*POW(aexp,2)/POW(param->unit.unit_v,2)/param->unit.unit_mass;
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
