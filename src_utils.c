/**
  * \file src_utils.c
  */

#ifdef WRAD

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "prototypes.h"
#include "oct.h"
#include "atomic_data/Atomic.h"
#include "chem_utils.h"
#include "tools.h"
#include "convert.h"

#ifdef WMPI
#include <mpi.h>
#endif // WMPI

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
    // below is essentially for outputs
    cellcoarse->rfield.nhplus +=cell->rfield.nhplus*fact;
     cellcoarse->rfield.eint +=cell->rfield.eint*fact;

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
  REAL  dvcur=POW(dxcur,3);
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
      cell->rfield.src=param->srcint/POW(X0*param->unit.unit_l,3)*param->unit.unit_t/param->unit.unit_N*POW(aexp,2);//8.;///8.;///POW(1./16.,3);
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



#ifdef UVBKG
 /*
  * set uvbkg according to the value in uvbkg.dat
  */
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    //(comoving photons/s/m3)

    cell->rfield.src += param->uv.value[igrp] * param->unit.unit_t*POW(aexp,5);
    cell->rfieldnew.src=cell->rfield.src;
  }
#else
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
   /* --------------------------------------------------------------------------
    * decreasing luminosity state, at the and of their life star still radiate
    * with a luminosity function of a decreasing power law of time.
    * Slope derived from http://www.stsci.edu/science/starburst99/figs/fig77.html
    * --------------------------------------------------------------------------
    */

      nss++;
      REAL slope = -4.;
      //printf("star found ! t=%e age=%e agelim=%e idx=%d COUNT=%d\n",tcur,curp->age,param->stars->tlife,curp->idx,(( (tcur - curp->age) < param->stars->tlife  )&&(tcur>= curp->age)));
      REAL src = (curp->mass*param->unit.unit_d)*srcint/POW(dxcur,3)*param->unit.unit_t/param->unit.unit_N*POW(aexp,2);
      REAL t = (param->cosmo->tphy - curp->age) / param->stars->tlife;
      cell->rfield.src +=  src*POW(t,slope);
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
#endif //UVBKG
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
// ==================================================================================================================
// ===============UV background ====================================================================================
// ==================================================================================================================

#ifdef UVBKG
// ----------------------------------------------------------//
/**
  * \brief Read a UV background from uvbkg.dat and store it.
  *
  * input values must be in phot/s/m^3
  *
  * uvbkg.dat format:
  * 1 int N : number of samples
  * N lines of 2 floats: (redshift) and (comoving photons/s/m3)
  */
// ----------------------------------------------------------//
void setUVBKG(struct RUNPARAMS *param, char *fname){

  //TODO consider NGRP>1
  //printf("Reading UVBKG from file :%s\n",fname);

  FILE *buf;
  buf=fopen(fname,"r");
  if(buf==NULL){
    printf("ERROR : cannot open the file %s, please check\n",fname);
    abort();
  }else{

    size_t rstat=fscanf(buf,"%d",&param->uv.N);

    param->uv.redshift=(REAL*)calloc(param->uv.N,sizeof(REAL));
    param->uv.Nphot=(REAL*)calloc(param->uv.N,sizeof(REAL));
    param->uv.value=(REAL*)calloc(NGRP,sizeof(REAL));

    int i;
    for(i=0; i<param->uv.N; i++){
      rstat=fscanf(buf,"%lf %lf",&param->uv.redshift[i],&param->uv.Nphot[i]);
    }
  }
}

// ----------------------------------------------------------//
/**
  * \brief Linear fit of UV background
  *
  * Made a linear interpolation of the data from uvbkg.dat
  * to get a value at current aexp.
  */
// ----------------------------------------------------------//
void setUVvalue(struct RUNPARAMS *param, REAL aexp){

  //TODO consider NGRP>1
  if(NGRP>1) printf("WARNING BAD BEHAVIOR FOR BKG with NGRP>1 !\n");

  REAL z = 1./aexp - 1.;
  if (z > param->uv.redshift[0]){
    int igrp;
    for(igrp=0;igrp<NGRP;igrp++){
      int iredshift;
      for(iredshift=1; iredshift<param->uv.N; iredshift++){
        if(z <= param->uv.redshift[iredshift]){

          REAL y1 = param->uv.Nphot[iredshift];
          REAL y2 = param->uv.Nphot[iredshift-1];
          REAL x1 = param->uv.redshift[iredshift];
          REAL x2 = param->uv.redshift[iredshift-1];

          param->uv.value[igrp] = (z-x1) * (y2-y1)/(x2-x1) + y1;
          break;
        }
      }
    }
  }else{
    int igrp;
    for(igrp=0;igrp<NGRP;igrp++) param->uv.value[igrp]=0.;
  }
}
#endif


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

#ifdef UVBKG
  setUVvalue(param, aexp);
#endif

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
