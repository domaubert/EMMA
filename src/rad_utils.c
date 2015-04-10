
#ifdef WRAD
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>

#include <mpi.h>
#include "hydro_utils.h"

#ifdef WCHEM
#include "chem_utils.h"
#endif

#ifdef GPUAXL
#include "rad_utils_gpu.h"
#endif

#include <omp.h>

//================================================================================
void distribE(struct CELL *cellcoarse, struct CELL *cell,struct RUNPARAMS *param,int level){

  // found a root
  if(cell->child==NULL){
    int igrp;
    for(igrp=0;igrp<NGRP;igrp++){
      cell->rfield.e[igrp] =cellcoarse->rfield.e[igrp];
      cell->rfield.fx[igrp]=cellcoarse->rfield.fx[igrp];
      cell->rfield.fy[igrp]=cellcoarse->rfield.fy[igrp];
      cell->rfield.fz[igrp]=cellcoarse->rfield.fz[igrp];
    }
  }
  else{
    // recursive call
    int icell;
    for(icell=0;icell<8;icell++){
      distribE(cellcoarse,&(cell->child->cell[icell]),param,level+1);
    }
  }
}


//================================================================================
void diffR(struct Rtype *W2, struct Rtype *W1, struct Rtype *WR){
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    WR->e[igrp]=W2->e[igrp]- W1->e[igrp];
    WR->fx[igrp]=W2->fx[igrp]- W1->fx[igrp];
    WR->fy[igrp]=W2->fy[igrp]- W1->fy[igrp];
    WR->fz[igrp]=W2->fz[igrp]- W1->fz[igrp];
  }
    WR->src=W2->src- W1->src;
#ifdef STARS
    WR->snfb=W2->snfb- W1->snfb;
#endif
#ifdef WCHEM
    WR->nhplus=W2->nhplus-W1->nhplus;
    WR->eint=W2->eint-W1->eint;
    WR->nh=W2->nh-W1->nh;
#endif
}

//================================================================================
void minmod_R(struct Rtype *Wm, struct Rtype *Wp, struct Rtype *Wr){

  REAL beta=1.; // 1. for MINBEE 2. for SUPERBEE
  // FLUX LIMITER
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){

    if(Wp->e[igrp]>0){
      Wr->e[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->e[igrp],Wp->e[igrp])),FMIN(Wm->e[igrp],beta*Wp->e[igrp]));
    }
    else{
      Wr->e[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->e[igrp],Wp->e[igrp])),FMAX(Wm->e[igrp],beta*Wp->e[igrp]));
    }


    if(Wp->fx[igrp]>0){
      Wr->fx[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->fx[igrp],Wp->fx[igrp])),FMIN(Wm->fx[igrp],beta*Wp->fx[igrp]));
    }
    else{
      Wr->fx[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->fx[igrp],Wp->fx[igrp])),FMAX(Wm->fx[igrp],beta*Wp->fx[igrp]));
    }


    if(Wp->fy[igrp]>0){
      Wr->fy[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->fy[igrp],Wp->fy[igrp])),FMIN(Wm->fy[igrp],beta*Wp->fy[igrp]));
    }
    else{
      Wr->fy[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->fy[igrp],Wp->fy[igrp])),FMAX(Wm->fy[igrp],beta*Wp->fy[igrp]));
    }


    if(Wp->fz[igrp]>0){
      Wr->fz[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->fz[igrp],Wp->fz[igrp])),FMIN(Wm->fz[igrp],beta*Wp->fz[igrp]));
    }
    else{
      Wr->fz[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->fz[igrp],Wp->fz[igrp])),FMAX(Wm->fz[igrp],beta*Wp->fz[igrp]));
    }

  }

  if(Wp->src>0){
      Wr->src=FMAX(FMAX(0.,FMIN(beta*Wm->src,Wp->src)),FMIN(Wm->src,beta*Wp->src));
    }
    else{
      Wr->src=FMIN(FMIN(0.,FMAX(beta*Wm->src,Wp->src)),FMAX(Wm->src,beta*Wp->src));
    }


#ifdef STARS
  if(Wp->snfb>0){
      Wr->snfb=FMAX(FMAX(0.,FMIN(beta*Wm->snfb,Wp->snfb)),FMIN(Wm->snfb,beta*Wp->snfb));
    }
    else{
      Wr->snfb=FMIN(FMIN(0.,FMAX(beta*Wm->snfb,Wp->snfb)),FMAX(Wm->snfb,beta*Wp->snfb));
    }
#endif

#ifdef WCHEM
  if(Wp->nhplus>0){
    Wr->nhplus=FMAX(FMAX(0.,FMIN(beta*Wm->nhplus,Wp->nhplus)),FMIN(Wm->nhplus,beta*Wp->nhplus));
  }
  else{
    Wr->nhplus=FMIN(FMIN(0.,FMAX(beta*Wm->nhplus,Wp->nhplus)),FMAX(Wm->nhplus,beta*Wp->nhplus));
  }

  if(Wp->eint>0){
    Wr->eint=FMAX(FMAX(0.,FMIN(beta*Wm->eint,Wp->eint)),FMIN(Wm->eint,beta*Wp->eint));
  }
  else{
    Wr->eint=FMIN(FMIN(0.,FMAX(beta*Wm->eint,Wp->eint)),FMAX(Wm->eint,beta*Wp->eint));
  }

  if(Wp->nh>0){
    Wr->nh=FMAX(FMAX(0.,FMIN(beta*Wm->nh,Wp->nh)),FMIN(Wm->nh,beta*Wp->nh));
  }
  else{
    Wr->nh=FMIN(FMIN(0.,FMAX(beta*Wm->nh,Wp->nh)),FMAX(Wm->nh,beta*Wp->nh));
  }
#endif

}

//================================================================================
void interpminmod_R(struct Rtype *W0, struct Rtype *Wp, struct Rtype *Dx, struct Rtype *Dy, struct Rtype *Dz,REAL dx,REAL dy,REAL dz){
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    Wp->e[igrp] =W0->e[igrp] +dx*Dx->e[igrp] +dy*Dy->e[igrp] +dz*Dz->e[igrp];
    Wp->fx[igrp] =W0->fx[igrp] +dx*Dx->fx[igrp] +dy*Dy->fx[igrp] +dz*Dz->fx[igrp];
    Wp->fy[igrp] =W0->fy[igrp] +dx*Dx->fy[igrp] +dy*Dy->fy[igrp] +dz*Dz->fy[igrp];
    Wp->fz[igrp] =W0->fz[igrp] +dx*Dx->fz[igrp] +dy*Dy->fz[igrp] +dz*Dz->fz[igrp];
  }
    Wp->src =W0->src +dx*Dx->src +dy*Dy->src +dz*Dz->src;

#ifdef WCHEM
    Wp->nhplus =W0->nhplus + dx*Dx->nhplus + dy*Dy->nhplus + dz*Dz->nhplus;
    Wp->eint =W0->eint +dx*Dx->eint +dy*Dy->eint +dz*Dz->eint;
    Wp->nh =W0->nh +dx*Dx->nh +dy*Dy->nh +dz*Dz->nh;
#endif
#ifdef STARS
    Wp->snfb =W0->snfb +dx*Dx->snfb +dy*Dy->snfb +dz*Dz->snfb;
#endif
}

//================================================================================

void coarse2fine_rad2(struct CELL *cell, struct Rtype *Wi, REAL cloc){

  struct OCT * oct;
  struct Rtype *W0;
	  struct Rtype *Wp;
	  struct Rtype *Wm;
	  struct Rtype Wint;
	  struct Rtype Dp,Dm;
	  struct Rtype D[3];
	  struct Rtype *W;
	  int inei2;
	  int vcell[6],vnei[6];
	  int dir;
	  REAL dxcur;

	  oct=cell2oct(cell);
	  getcellnei(cell->idx, vnei, vcell); // we get the neighbors
	  dxcur=POW(0.5,oct->level);

	  W0=&(cell->rfield);
	  // Limited Slopes
	  for(dir=0;dir<3;dir++){

	    inei2=2*dir;
	    if(vnei[inei2]==6){
	      Wm=&(oct->cell[vcell[inei2]].rfield);
	    }
	    else{
	      Wm=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].rfield);
	    }

#ifdef TRANSXM
	      if((oct->x==0.)&&(inei2==0)){
		Wm=&(cell->rfield);
	      }
#endif

#ifdef TRANSYM
	      if((oct->y==0.)&&(inei2==2)){
		Wm=&(cell->rfield);
	      }
#endif

#ifdef TRANSZM
	      if((oct->z==0.)&&(inei2==4)){
		Wm=&(cell->rfield);
	      }
#endif



	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      Wp=&(oct->cell[vcell[inei2]].rfield);
	    }
	    else{
	      Wp=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].rfield);
	    }
#ifdef TRANSXP
	      if(((oct->x+2.*dxcur)==1.)&&(inei2==1)){
		Wp=&(cell->rfield);
	      }
#endif

#ifdef TRANSYP
	      if(((oct->y+2.*dxcur)==1.)&&(inei2==3)){
		Wp=&(cell->rfield);
	      }
#endif

#ifdef TRANSZP
	      if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
		Wp=&(cell->rfield);
	      }
#endif


	    diffR(Wp,W0,&Dp);
	    diffR(W0,Wm,&Dm);


	    minmod_R(&Dm,&Dp,D+dir);


	  }

	  // Interpolation
	  int ix,iy,iz;
	  int icell;

	  // =================================================
	  // =================================================
	  // =================================================
	  // =================================================

	  int tag=0;
	  for(iz=0;iz<2;iz++){
	    for(iy=0;iy<2;iy++){
	      for(ix=0;ix<2;ix++){
		icell=ix+iy*2+iz*4;
		interpminmod_R(W0,&Wint,D,D+1,D+2,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Wp contains the interpolation

		REAL F=SQRT(Wint.fx[0]*Wint.fx[0]+Wint.fy[0]*Wint.fy[0]+Wint.fz[0]*Wint.fz[0]);
		REAL E=Wint.e[0];
		REAL f=F/cloc/E;
		REAL F0=SQRT(W0->fx[0]*W0->fx[0]+W0->fy[0]*W0->fy[0]+W0->fz[0]*W0->fz[0]);
		REAL E0=W0->e[0];
		REAL f0=F0/cloc/E0;
		if(f>1.){
		  //printf("flux error f=%e F=%e E=%e c=%e/ f0=%e F0=%e E0=%e CORRECTION\n",f,F,E,cloc,f0,F0,E0);
		  //tag=1;
		  /* Wint.fx[0]/=1.02*f; */
		  /* Wint.fy[0]/=1.02*f; */
		  /* Wint.fz[0]/=1.02*f; */

		  memcpy(Wi+icell,W0,sizeof(struct Rtype));
		}
		else{
		  memcpy(Wi+icell,&Wint,sizeof(struct Rtype));
		}
	      }
	    }
	  }

	  /* if(tag){ */
	  /*   for(iz=0;iz<8;iz++){ */
	  /*     REAL F0=SQRT(Wi[iz].fx[0]*Wi[iz].fx[0]+Wi[iz].fy[0]*Wi[iz].fy[0]+Wi[iz].fz[0]*Wi[iz].fz[0]); */
	  /*     REAL E0=Wi[iz].e[0]; */
	  /*     REAL f0=F0/cloc/E0; */
	  /*     printf("i=%d  f=%e F=%e E=%e\n",iz,f0,F0,E0); */
	  /*   } */
	  /*   abort(); */
	  /* } */


}


// =====================================================================
// =====================================================================

void coarse2fine_radlin(struct CELL *cell, struct Rtype *Wi){


	  struct OCT * oct;

	  struct Rtype *W0;
	  struct Rtype *Wp;
	  struct Rtype *Wm;
	  struct Rtype Wint;
	  struct Rtype *W;

	  struct Rtype D[6];

	  int inei2;
	  int vcell[6],vnei[6];
	  int dir;
	  REAL dxcur;

	  oct=cell2oct(cell);
	  getcellnei(cell->idx, vnei, vcell); // we get the neighbors
	  dxcur=POW(0.5,oct->level);

	  W0=&(cell->rfield);
	  // Limited Slopes
	  for(dir=0;dir<3;dir++){

	    inei2=2*dir;
	    if(vnei[inei2]==6){
	      Wm=&(oct->cell[vcell[inei2]].rfield);
	    }
	    else{
	      Wm=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].rfield);
	    }
#ifdef TRANSXM
	      if((oct->x==0.)&&(inei2==0)){
		Wm=&(cell->rfield);
	      }
#endif

#ifdef TRANSYM
	      if((oct->y==0.)&&(inei2==2)){
		Wm=&(cell->rfield);
	      }
#endif

#ifdef TRANSZM
	      if((oct->z==0.)&&(inei2==4)){
		Wm=&(cell->rfield);
	      }
#endif


	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      Wp=&(oct->cell[vcell[inei2]].rfield);
	    }
	    else{
	      Wp=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].rfield);

	    }
#ifdef TRANSXP
	      if(((oct->x+2.*dxcur)==1.)&&(inei2==1)){
		Wp=&(cell->rfield);
	      }
#endif

#ifdef TRANSYP
	      if(((oct->y+2.*dxcur)==1.)&&(inei2==3)){
		Wp=&(cell->rfield);
	      }
#endif

#ifdef TRANSZP
	      if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
		Wp=&(cell->rfield);
	      }
#endif

	    diffR(W0,Wm,D+2*dir+0);
	    diffR(Wp,W0,D+2*dir+1);

	  }

	  // Interpolation
	  int ix,iy,iz;
	  int icell;

	  // =================================================
	  // =================================================
	  // =================================================
	  // =================================================

	  for(iz=0;iz<2;iz++){
	    for(iy=0;iy<2;iy++){
	      for(ix=0;ix<2;ix++){
		icell=ix+iy*2+iz*4;
		interpminmod_R(W0,&Wint,D+ix,D+iy+2,D+iz+4,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Wp contains the interpolation

		memcpy(Wi+icell,&Wint,sizeof(struct Rtype));

	      }
	    }
	  }

}

//================================================================================
//================================================================================
REAL Eddington(REAL fx, REAL fy, REAL fz, REAL ee, REAL c,int i,int j)
{
  REAL c2e=ee*c*c; // 2 flop
  REAL ff=0.;
  REAL arg,chi,res=0.;
  REAL n[3];
  n[0]=0.;n[1]=0.;n[2]=0.;

  if(ee>0)
    {
      ff=SQRT(fx*fx+fy*fy+fz*fz); // 6 flop
      if(ff>0)
	{
	  n[0]=fx/ff;
	  n[1]=fy/ff;
	  n[2]=fz/ff;
	}
      ff=ff/(c*ee); // 2flop
    }

  arg=FMAX(4.-3.*ff*ff,0.); // 4 flop
  chi=(3.+4.*ff*ff)/(5.+2.*SQRT(arg)); // 7 flops

  if(i==j) res=(1.-chi)/2.*c2e; // 1 flops on average
  arg=(3.*chi-1.)/2.*c2e;
  res+=arg*n[i]*n[j];

  return res;
}



// =============================================================================================================

int rad_sweepX(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i,igrp;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];
  struct Rtype Rold;

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fact;
  REAL fp,fm;
  REAL up,um;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors

    for(i=0;i<nread;i++){ // we scan the octs

      memset(FL,0,sizeof(REAL)*NVAR_R);
      memset(FR,0,sizeof(REAL)*NVAR_R);

      // Getting the original state ===========================

      curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);

      /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);
      for(iface=0;iface<2;iface++){
	memcpy(RC+iface,curcell,sizeof(struct Rtype));
      }

      // Neighbor "MUSCL-LIKE" reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface;
	memcpy(RN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].rfield),sizeof(struct Rtype));

	int condsplit;
#ifdef COARSERAD
	condsplit=1;
#else
	condsplit=(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split);
#endif
	if(condsplit){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }

      // X DIRECTION =========================================================================

      // --------- solving the Riemann Problems LEFT

      //int fg=0;
      for(igrp=0;igrp<NGRP;igrp++){

	// E
	up=RC[0].e[igrp];
	um=RN[0].e[igrp];


	fp=RC[0].fx[igrp];
	fm=RN[0].fx[igrp];

	FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);


	//FX

	up=RC[0].fx[igrp];
	um=RN[0].fx[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,0);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,0);

	FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RC[0].fy[igrp];
	um=RN[0].fy[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,0);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,0);

	FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FZ

	up=RC[0].fz[igrp];
	um=RN[0].fz[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,0);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,0);

	FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }


      // ===========================================

      // --------- solving the Riemann Problems RIGHT

      for(igrp=0;igrp<NGRP;igrp++){

	// E
	up=RN[1].e[igrp];
	um=RC[1].e[igrp];

	fp=RN[1].fx[igrp];
	fm=RC[1].fx[igrp];

	FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	/* if((RN[1].e[0]>1.0101e60)&&(RN[1].e[0]<1.0102e60)){ */
	/*   printf("XP up=%e um=%e fp=%e fm=%e FL=%e ,ffact=%d c=%e\n",up,um,fp,fm,FR[0+igrp*NVAR_R]*dt/dx,ffact[1],c); */
	/* } */

	/* if(fg){ */
	/*   printf("coucou %e %e %e\n",um,up,FR[0+igrp*NVAR_R]); */

	/* } */

	//FX

	up=RN[1].fx[igrp];
	um=RC[1].fx[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,0);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,0);

	FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RN[1].fy[igrp];
	um=RC[1].fy[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,0);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,0);

	FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RN[1].fz[igrp];
	um=RC[1].fz[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,0);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,0);

	FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }

      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours
      /*  if((ffact[0]==0)||(ffact[1]==0)){ */
      /* 	 printf("big step 0x ec=%e en=%e flx=%e s=%d %d\n",RC[0].e[0],RN[0].e[0],FL[0]*dt/dx,stencil[i].oct[6].cell[icell].split,ffact[0]); */
      /* 	printf("big step 1x ec=%e en=%e flx=%e s=%d %d\n",RC[1].e[0],RN[1].e[0],FR[0]*dt/dx,stencil[i].oct[6].cell[icell].split,ffact[1]); */



      /* } */


      for(igrp=0;igrp<NGRP;igrp++){
	for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0];
	for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1];
      }
      /* if(fg){ */
      /* 	printf("FF %e %e %d %d\n", FL[0+0*NVAR_R], FR[0+0*NVAR_R],stencil[i].oct[ioct[vnei[0]]].cell[vcell[0]].split,stencil[i].oct[ioct[vnei[1]]].cell[vcell[1]].split); */
      /* } */

      memcpy(stencil[i].New.cell[icell].rflux+0*NVAR_R*NGRP,FL,sizeof(REAL)*NVAR_R*NGRP);
      memcpy(stencil[i].New.cell[icell].rflux+1*NVAR_R*NGRP,FR,sizeof(REAL)*NVAR_R*NGRP);


      // ready for the next cell
    }
    //ready for the next oct
  }
  return 0;
}



// =============================================================================================================

int rad_sweepY(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i,igrp;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];
  struct Rtype Rold;

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fact;
  REAL fp,fm;
  REAL up,um;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors

    for(i=0;i<nread;i++){ // we scan the octs

      memset(FL,0,sizeof(REAL)*NVAR_R);
      memset(FR,0,sizeof(REAL)*NVAR_R);

      // Getting the original state ===========================

      curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);

      /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);
      for(iface=0;iface<2;iface++){
	memcpy(RC+iface,curcell,sizeof(struct Rtype));
      }

      // Neighbor "MUSCL-LIKE" reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+2;
	memcpy(RN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].rfield),sizeof(struct Rtype));

	int condsplit;
#ifdef COARSERAD
	condsplit=1;
#else
	condsplit=(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split);
#endif
	if(condsplit){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }




      // Y DIRECTION =========================================================================

      // --------- solving the Riemann Problems FRONT



      for(igrp=0;igrp<NGRP;igrp++){

	// E
	up=RC[0].e[igrp];
	um=RN[0].e[igrp];

	fp=RC[0].fy[igrp];
	fm=RN[0].fy[igrp];

	FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

 	/* if((RN[0].e[0]>1.0101e60)&&(RN[0].e[0]<1.0102e60)){ */
	/*   printf("YM up=%e um=%e fp=%e fm=%e FL=%e ,ffact=%d c=%e\n",up,um,fp,fm,FL[0+igrp*NVAR_R]*dt/dx,ffact[0],c); */
	/* } */

	//FX

	up=RC[0].fx[igrp];
	um=RN[0].fx[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,1);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,1);

	FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RC[0].fy[igrp];
	um=RN[0].fy[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,1);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,1);

	//if(up==1.) printf("FL=%e fp=%e fm=%e || ",0.5*(fp+fm)+0.5*c*(um-up),fp,fm);

	FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RC[0].fz[igrp];
	um=RN[0].fz[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,1);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,1);

	FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }


	// ===========================================




      // --------- solving the Riemann Problems BACK

      for(igrp=0;igrp<NGRP;igrp++){

	// E
	up=RN[1].e[igrp];
	um=RC[1].e[igrp];

	fp=RN[1].fy[igrp];
	fm=RC[1].fy[igrp];

	FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	/* if((RN[1].e[0]>1.0101e60)&&(RN[1].e[0]<1.0102e60)){ */
	/*   printf("YP up=%e um=%e fp=%e fm=%e FR=%e ffact=%d c=%e\n",up,um,fp,fm,FR[0+igrp*NVAR_R]*dt/dx,ffact[1],c); */
	/* } */

	//FX

	up=RN[1].fx[igrp];
	um=RC[1].fx[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,1);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,1);

	FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RN[1].fy[igrp];
	um=RC[1].fy[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,1);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,1);

	//printf("FR=%e fp=%e fm=%e up=%e um=%e\n",0.5*(fp+fm+c*(um-up)),fp,fm, up,um);
	FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RN[1].fz[igrp];
	um=RC[1].fz[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,1);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,1);

	FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }


      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours

      /*  if((ffact[0]==0)||(ffact[1]==0)){ */
      /* 	 printf("big step 0y ec=%e en=%e flx=%e s=%d %d\n",RC[0].e[0],RN[0].e[0],FL[0]*dt/dx,stencil[i].oct[6].cell[icell].split,ffact[0]); */
      /* 	printf("big step 1y ec=%e en=%e flx=%e s=%d %d\n",RC[1].e[0],RN[1].e[0],FR[0]*dt/dx,stencil[i].oct[6].cell[icell].split,ffact[1]); */
      /* } */


      for(igrp=0;igrp<NGRP;igrp++){
	for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0];
	for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1];
      }

      //if(RC[0].e[0]-FR[0]*dt/dx<0.) abort();

      memcpy(stencil[i].New.cell[icell].rflux+2*NVAR_R*NGRP,FL,sizeof(REAL)*NVAR_R*NGRP);
      memcpy(stencil[i].New.cell[icell].rflux+3*NVAR_R*NGRP,FR,sizeof(REAL)*NVAR_R*NGRP);
      //if((stencil[i].oct[6].cell[icell].rfield.e[0]==0.)) if (((FR[2]-FL[2])!=0.)&&((FR[0]-FL[0])==0.)) abort();
      //if (((FR[2]-FL[2])!=0.)&&(stencil[i].oct[6].cell[icell].rfield.fy[0]==0.)) abort();
      //if ((FR[2]==1)&&(FL[2]==1.)) abort();

      // ready for the next cell
    }
    //ready for the next oct
  }

  return 0;
}


// ===================================================================================================

int rad_sweepZ(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i,igrp;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];
  struct Rtype Rold;

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fact;
  REAL fp,fm;
  REAL up,um;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors

    for(i=0;i<nread;i++){ // we scan the octs

      memset(FL,0,sizeof(REAL)*NVAR_R);
      memset(FR,0,sizeof(REAL)*NVAR_R);

      // Getting the original state ===========================

      curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);

      /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);
      for(iface=0;iface<2;iface++){
	memcpy(RC+iface,curcell,sizeof(struct Rtype));
      }

      // Neighbor "MUSCL-LIKE" reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+4;
	memcpy(RN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].rfield),sizeof(struct Rtype));

	int condsplit;
#ifdef COARSERAD
	condsplit=1;
#else
	condsplit=(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split);
#endif
	if(condsplit){
	  ffact[iface]=1; // we consider the contriubtion of split neighbors
	}

      }


      // Z DIRECTION =========================================================================

      // --------- solving the Riemann Problems BOTTOM


      for(igrp=0;igrp<NGRP;igrp++){

	// E
	up=RC[0].e[igrp];
	um=RN[0].e[igrp];

	fp=RC[0].fz[igrp];
	fm=RN[0].fz[igrp];

	FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

 	/* if((RN[0].e[0]>1.0101e60)&&(RN[0].e[0]<1.0102e60)){ */
	/*   printf("ZM up=%e um=%e fp=%e fm=%e FL=%e ,ffact=%d c=%e\n",up,um,fp,fm,FL[0+igrp*NVAR_R]*dt/dx,ffact[0],c); */
	/* } */

	//FX

	up=RC[0].fx[igrp];
	um=RN[0].fx[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,2);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,2);

	FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RC[0].fy[igrp];
	um=RN[0].fy[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,2);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,2);

	FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FZ

	up=RC[0].fz[igrp];
	um=RN[0].fz[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,2);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,2);

	FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }


	// ===========================================




      // --------- solving the Riemann Problems TOP

      for(igrp=0;igrp<NGRP;igrp++){

	// E
	up=RN[1].e[igrp];
	um=RC[1].e[igrp];

	fp=RN[1].fz[igrp];
	fm=RC[1].fz[igrp];

	FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

 	/* if((RN[1].e[0]>1.0101e60)&&(RN[1].e[0]<1.0102e60)){ */
	/*   printf("ZP up=%e um=%e fp=%e fm=%e FL=%e ,ffact=%d c=%e\n",up,um,fp,fm,FR[0+igrp*NVAR_R]*dt/dx,ffact[1],c); */
	/* } */

	//FX

	up=RN[1].fx[igrp];
	um=RC[1].fx[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,2);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,2);

	FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RN[1].fy[igrp];
	um=RC[1].fy[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,2);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,2);

	FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RN[1].fz[igrp];
	um=RC[1].fz[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,2);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,2);

	FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }

      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours

      /*  if((ffact[0]==0)||(ffact[1]==0)){ */
      /* 	 printf("big step 0z ec=%e en=%e flx=%e s=%d %d\n",RC[0].e[0],RN[0].e[0],FL[0]*dt/dx,stencil[i].oct[6].cell[icell].split,ffact[0]); */
      /* 	printf("big step 1z ec=%e en=%e flx=%e s=%d %d\n",RC[1].e[0],RN[1].e[0],FR[0]*dt/dx,stencil[i].oct[6].cell[icell].split,ffact[1]); */
      /* } */

      for(igrp=0;igrp<NGRP;igrp++){
	for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0];
	for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1];
      }

      memcpy(stencil[i].New.cell[icell].rflux+4*NVAR_R*NGRP,FL,sizeof(REAL)*NVAR_R*NGRP);
      memcpy(stencil[i].New.cell[icell].rflux+5*NVAR_R*NGRP,FR,sizeof(REAL)*NVAR_R*NGRP);

      // ready for the next cell
    }
    //ready for the next oct
  }

  return 0;
}



//==================================================================================
//==================================================================================

#ifdef WRADTEST
void recursive_neighbor_gather_oct_rad(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct RGRID *stencil, REAL cloc){


  // =======================================
  // This function is pretty much an overkill
  // Could be simplified

  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;
  int face[8]={0,1,2,3,4,5,6,7};
  int tflag[6]={0,0,0,0,0,0};
  REAL dxcur;
  int igrp;

#ifndef ALTG
  struct Rtype Ri[8];
#else
  struct Rtype *Ri[8];
#endif

  char child[8];

  struct OCT *oct;
  struct OCT *neioct;
  struct CELL *neicell;

  ineiloc=inei;


  if(cell->child!=NULL){
    // the oct at the right level exists
    neicell=cell->child->nei[ineiloc];
    oct=cell->child;
    dxcur=POW(0.5,oct->level);

#ifdef TRANSXP
    if(ineiloc==1){
      if((oct->x+2.*dxcur)==1.){
	neicell=cell;
	face[0]=1;
	face[1]=1;
	face[2]=3;
	face[3]=3;
	face[4]=5;
	face[5]=5;
	face[6]=7;
	face[7]=7;
	tflag[1]=1;
      }
    }
#endif


#ifdef TRANSYP
    if(ineiloc==3){
      if((oct->y+2.*dxcur)==1.){
	neicell=cell;
	face[0]=2;
	face[1]=3;
	face[2]=2;
	face[3]=3;
	face[4]=7;
	face[5]=6;
	face[6]=6;
	face[7]=7;
	tflag[3]=1;
      }
    }
#endif

#ifdef TRANSZP
    if(ineiloc==5){
      if((oct->z+2.*dxcur)==1.){
	neicell=cell;
	face[0]=4;
	face[1]=5;
	face[2]=6;
	face[3]=7;
	face[4]=4;
	face[5]=5;
	face[6]=6;
	face[7]=7;
	tflag[5]=1;

      }
    }
#endif



#ifdef TRANSXM
    if(ineiloc==0){
      if(oct->x==0.){
	neicell=cell;
	face[0]=0;
	face[1]=0;
	face[2]=2;
	face[3]=2;
	face[4]=4;
	face[5]=4;
	face[6]=6;
	face[7]=6;
	tflag[0]=1;

      }
    }
#endif

#ifdef TRANSYM
    if(ineiloc==2){
      if(oct->y==0.){
	neicell=cell;
	face[0]=0;
	face[1]=1;
	face[2]=0;
	face[3]=1;
	face[4]=4;
	face[5]=5;
	face[6]=4;
	face[7]=5;
	tflag[2]=1;

      }
    }
#endif

#ifdef TRANSZM
    if(ineiloc==4){
      if(oct->z==0.){
	neicell=cell;
	face[0]=0;
	face[1]=1;
	face[2]=2;
	face[3]=3;
	face[4]=0;
	face[5]=1;
	face[6]=2;
	face[7]=3;
	tflag[4]=1;

      }
    }
#endif
  }
  else{
    getcellnei(cell->idx, vnei, vcell); // we get the neighbors
    oct=cell2oct(cell);
    dxcur=POW(0.5,oct->level);
    if(vnei[ineiloc]==6){
      neicell=&(oct->cell[vcell[ineiloc]]);
    }
    else{
      if(oct->nei[ineiloc]->child!=NULL){
	neicell=&(oct->nei[ineiloc]->child->cell[vcell[ineiloc]]);
      }
      else{
	printf("big problem\n");
	abort();
      }

#ifdef TRANSXM
      if(ineiloc==0){
	if(oct->x==0.){
	  neicell=cell;
	  face[0]=0;
	  face[1]=0;
	  face[2]=2;
	  face[3]=2;
	  face[4]=4;
	  face[5]=4;
	  face[6]=6;
	  face[7]=6;
	  tflag[0]=1;
	}
      }
#endif
#ifdef TRANSXP
      if(ineiloc==1){
	if((oct->x+2.*dxcur)==1.){
	  neicell=cell;
	  face[0]=1;
	  face[1]=1;
	  face[2]=3;
	  face[3]=3;
	  face[4]=5;
	  face[5]=5;
	  face[6]=7;
	  face[7]=7;
	  tflag[1]=1;
	}
      }
#endif


#ifdef TRANSYP
      if(ineiloc==3){
	if((oct->y+2.*dxcur)==1.){
	  neicell=cell;
	  face[0]=2;
	  face[1]=3;
	  face[2]=2;
	  face[3]=3;
	  face[4]=7;
	  face[5]=6;
	  face[6]=6;
	  face[7]=7;
	  tflag[3]=1;

	}
      }
#endif

#ifdef TRANSZP
      if(ineiloc==5){
	if((oct->z+2.*dxcur)==1.){
	  neicell=cell;
	  face[0]=4;
	  face[1]=5;
	  face[2]=6;
	  face[3]=7;
	  face[4]=4;
	  face[5]=5;
	  face[6]=6;
	  face[7]=7;
	  tflag[5]=1;

	}
      }
#endif




#ifdef TRANSYM
      if(ineiloc==2){
	if(oct->y==0.){
	  neicell=cell;
	  face[0]=0;
	  face[1]=1;
	  face[2]=0;
	  face[3]=1;
	  face[4]=4;
	  face[5]=5;
	  face[6]=4;
	  face[7]=5;
	  tflag[2]=1;
	}
      }
#endif

#ifdef TRANSZM
      if(ineiloc==4){
	if(oct->z==0.){
	  neicell=cell;
	  face[0]=0;
	  face[1]=1;
	  face[2]=2;
	  face[3]=3;
	  face[4]=0;
	  face[5]=1;
	  face[6]=2;
	  face[7]=3;
	  tflag[4]=1;

	}
      }
#endif

    }


  }


  if(neicell->child!=NULL){
    // optimal case
    for(icell=0;icell<8;icell++){
#ifndef ALTG
      memcpy(Ri+icell,&(neicell->child->cell[icell].rfield),sizeof(struct Rtype));
      child[icell]=(neicell->child->cell[icell].child!=NULL);
#else
      Ri[icell]=&(neicell->child->cell[icell].rfield);
      child[icell]=(neicell->child->cell[icell].child!=NULL);
#endif
    }

  }
  else{
#ifndef ALTG
    int il;
    for(il=0;il<8;il++) memcpy(&Ri[il],&neicell->rfield,sizeof(struct Rtype));

    for(icell=0;icell<8;icell++){
      child[icell]=0;
    }
#else
    //coarse2fine_rad2(neicell,Ri,cloc);
    int il;
    for(il=0;il<8;il++){
      Ri[il]=&(neicell->rfield);
      child[il]=0.;
    }
#endif
  }

#ifndef ALTG
  for(icell=0;icell<8;icell++){
    memcpy(&(stencil->oct[ioct].cell[icell].rfield),Ri+face[icell],sizeof(struct Rtype)); //
#ifdef TRANSXM
#ifdef REFXM
    if(tflag[0]){
      for(igrp=0;igrp<NGRP;igrp++) stencil->oct[ioct].cell[icell].rfield.fx[igrp]*=-1.;
    }
#endif
#endif


#ifdef TRANSYM
#ifdef REFYM
    if(tflag[2]){
      for(igrp=0;igrp<NGRP;igrp++) stencil->oct[ioct].cell[icell].rfield.fy[igrp]*=-1.;
    }
#endif
#endif

#ifdef TRANSZM
#ifdef REFZM
    if(tflag[4]){
      for(igrp=0;igrp<NGRP;igrp++) stencil->oct[ioct].cell[icell].rfield.fz[igrp]*=-1.;
    }
#endif
#endif
    stencil->oct[ioct].cell[icell].split=child[face[icell]];
  }
#else
  for(icell=0;icell<8;icell++){
    memcpy(&(stencil->oct[ioct].cell[icell].rfield),Ri[face[icell]],sizeof(struct Rtype)); //
    stencil->oct[ioct].cell[icell].split=child[face[icell]];
  }
#endif

}

#else
//void recursive_neighbor_gather_oct_rad(int ioct, struct CELL *cell, struct RGRID *stencil){
void recursive_neighbor_gather_oct_rad(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct RGRID *stencil, REAL cloc){
  // =======================================
  // This function is pretty much an overkill
  // Could be simplified

  int icell, igrp;
  int vnei[6],vcell[6];


  struct OCT *oct;
  struct OCT *neioct;
  struct CELL*neicell;


  if(cell->child!=NULL){
    // the oct at the right level exists
    neicell=cell->child->nei[ioct];
    oct=cell->child;
  }
  else{
    getcellnei(cell->idx, vnei, vcell); // we get the neighbors
    oct=cell2oct(cell);

    if(vnei[ioct]==6){
      neicell=&(oct->cell[vcell[ioct]]);
    }
    else{
      if(oct->nei[ioct]->child!=NULL){
	neicell=&(oct->nei[ioct]->child->cell[vcell[ioct]]);
      }
      else{
	printf("big problem\n");
	abort();
      }

    }
  }

  if(neicell->child!=NULL){
    // optimal case
    for(icell=0;icell<8;icell++){
      struct CELLLIGHT_R *c=&(stencil->oct[ioct].cell[icell]);
      struct CELL *co=&neicell->child->cell[icell];
      memcpy(&(c->rfield),&(co->rfield),sizeof(struct Rtype)); // UU
      c->split=(co->child!=NULL);
    }
  }
  else{
    struct Rtype Ri[8];
    struct Rtype *Rii[8];
    char child[8];

    //coarse2fine_rad2(neicell,Ri,cloc);
    for(icell=0;icell<8;icell++) {
      Rii  [icell] = &neicell->rfield;
      child[icell] = 0;
    }
    for(icell=0;icell<8;icell++){
      memcpy(&(stencil->oct[ioct].cell[icell].rfield),Rii[icell],sizeof(struct Rtype)); //
      stencil->oct[ioct].cell[icell].split=child[icell];
    }
  }


}
#endif


// ===================================================================================================
// ===================================================================================================

struct OCT *gatherstencilrad(struct OCT *octstart, struct RGRID *stencil, int stride, struct CPUINFO *cpu, int *nread, REAL cloc)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  //char visit[27]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int ioct;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;

      //memset(visit,0,27*sizeof(char));
      // filling the values in the central oct


      for(icell=0;icell<8;icell++){
	memcpy(&(stencil[iread].oct[6].cell[icell].rfield),&(curoct->cell[icell].rfield),sizeof(struct Rtype)); // for calculations
	memcpy(&(stencil[iread].New.cell[icell].rfieldnew),&(curoct->cell[icell].rfieldnew),sizeof(struct Rtype)); // for calculations
	stencil[iread].oct[6].cell[icell].split=(curoct->cell[icell].child!=NULL);
      }

      cell=curoct->parent;

      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  //ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  ioct=inei;
	  //visit[ioct]=1;
	  recursive_neighbor_gather_oct_rad(ioct, inei, -1, -1, 1, cell, stencil+iread,cloc);
	}
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}


// ===================================================================================================
// ===================================================================================================

void updatefieldrad(struct OCT *octstart, struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew,REAL cloc)
{
  int i,icell,igrp;
  struct Rtype R;
  struct Rtype Rupdate;
  REAL one;
  int flx;
  REAL dtsurdx=dtnew/dxcur;
  REAL F[NFLUX_R];
  REAL SRC;
  for(i=0;i<nread;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells


      int condsplit;
#ifdef COARSERAD
      condsplit=0;
#else
      condsplit=(stencil[i].oct[6].cell[icell].split);
#endif
      if(condsplit) continue;
      memcpy(F,stencil[i].New.cell[icell].rflux,sizeof(REAL)*NFLUX_R);// New fluxes from the stencil


      // ==== updating

      // actually we compute and store the delta U only
      one=1.;
      memset(&R,0,sizeof(struct Rtype)); // setting delta U
      for(flx=0;flx<6;flx++){
	for(igrp=0;igrp<NGRP;igrp++){
	  R.e[igrp]  +=F[0+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fx[igrp] +=F[1+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fy[igrp] +=F[2+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fz[igrp] +=F[3+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	}
	one*=-1.;
      }

#ifndef WCHEM
      // adding the source contribution
      for(igrp=0;igrp<NGRP;igrp++){
      	SRC=stencil[i].oct[6].cell[icell].rfield.src;
      	R.e[igrp]  +=SRC*dtnew+EMIN;
      }
#endif

      /* // scatter back the delta Uwithin the stencil */

      // TESTING FULL UPDATE IN STENCIL APPROACH


      memcpy(&Rupdate,&stencil[i].New.cell[icell].rfieldnew,sizeof(struct Rtype));


      for(igrp=0;igrp<NGRP;igrp++){
	Rupdate.e[igrp]   +=R.e[igrp];
	Rupdate.fx[igrp]  +=R.fx[igrp];
	Rupdate.fy[igrp]  +=R.fy[igrp];
	Rupdate.fz[igrp]  +=R.fz[igrp];
      }


#if 0


	// ================================ START MEGA DIAGNOSE =========================
      if((Rupdate.e[1]<0||isnan(Rupdate.e[1]))){
	  printf("ERROR Neg rad energy New=%e org=%e delta=%e srcloc=%e xion=%e eini=%e temp=%e\n",Rupdate.e[1],stencil[i].New.cell[icell].rfieldnew.e[1],R.e[1],stencil[i].oct[6].cell[icell].rfield.src,stencil[i].oct[6].cell[icell].rfield.nhplus/stencil[i].oct[6].cell[icell].rfield.nh,stencil[i].oct[6].cell[icell].rfield.e[1],stencil[i].oct[6].cell[icell].rfield.temp);
	  one=1.;
	  for(flx=0;flx<6;flx++){
	    printf("f%d %e\n",flx,F[0+0*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one);
	    one*=-1.;
	  }

	  printf("Flux %e %e %e delta %e %e %e rflux=%e %e\n",stencil[i].New.cell[icell].rfieldnew.fx[1],stencil[i].New.cell[icell].rfieldnew.fy[1],stencil[i].New.cell[icell].rfieldnew.fz[1],R.fx[1],R.fy[1],R.fz[1],SQRT(POW(stencil[i].New.cell[icell].rfieldnew.fx[1],2)+POW(stencil[i].New.cell[icell].rfieldnew.fy[1],2)+POW(stencil[i].New.cell[icell].rfieldnew.fz[1],2))/cloc/stencil[i].New.cell[icell].rfieldnew.e[1],R.e[1]);

	  int vnei[6],vcell[6];
	  getcellnei(icell, vnei, vcell); // we get the neighbors

	  for(flx=0;flx<6;flx++){

	    printf("Nei #%d e=%e src=%e xion=%e temp=%e split=%d\n",flx,stencil[i].oct[vnei[flx]].cell[vcell[flx]].rfield.e[1],stencil[i].oct[vnei[flx]].cell[vcell[flx]].rfield.src,stencil[i].oct[vnei[flx]].cell[vcell[flx]].rfield.nhplus/stencil[i].oct[vnei[flx]].cell[vcell[flx]].rfield.nh,stencil[i].oct[vnei[flx]].cell[vcell[flx]].rfield.temp,stencil[i].oct[vnei[flx]].cell[vcell[flx]].split);
	  }

	  abort();
	}
      // ================================ END MEGA DIAGNOSE =========================
#endif
      //memcpy(&(curoct->cell[icell].rfieldnew),&Rupdate,sizeof(struct Rtype));
      memcpy(&stencil[i].New.cell[icell].rfieldnew,&Rupdate,sizeof(struct Rtype));

      if(isnan(stencil[i].New.cell[icell].rfieldnew.e[1])){
	printf(" WTFFFF %e %e %e %e ||%e\n",stencil[i].New.cell[icell].rfieldnew.e[0],stencil[i].New.cell[icell].rfieldnew.fx[0],stencil[i].New.cell[icell].rfieldnew.fy[0],stencil[i].New.cell[icell].rfieldnew.fz[0],Rupdate.e[0]);
      }

    }
  }

}

// ===========================================================================================================

struct OCT *scatterstencilrad(struct OCT *octstart, struct RGRID *stencil, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew, REAL cloc)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int ioct[7]={0,1,2,3,4,5,6};
  int vnei[6],vcell[6],inei;

  nextoct=octstart;
  iread=0;

  struct Rtype R,deltaR;
  REAL dtsurdx=dtnew/dxcur;


  int one;
  REAL F[NVAR_R*NGRP];
  int igrp;

  //printf("let's scatter\n");
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	//we scatter the values in the central cell

	memcpy(&(curoct->cell[icell].rfieldnew),&(stencil[iread].New.cell[icell].rfieldnew),sizeof(struct Rtype));

	/* if(curoct->level==7) */
	/*   if(curoct->cell[icell].rfieldnew.e[0]>0.) printf("Hey N=%e S=%p\n",curoct->cell[icell].rfieldnew.e[0],curoct->cell[icell].child); */


	// let us now deal with coarser neighbors
	getcellnei(icell, vnei, vcell); // we get the neighbors

	for(inei=0;inei<6;inei++){

#ifdef TRANSXM
	  if((curoct->x==0.)&&(inei==0)){
	    continue;
	  }
#endif

#ifdef TRANSXP
	  if(((curoct->x+2.*dxcur)==1.)&&(inei==1)){
	    continue;
	  }
#endif

#ifdef TRANSYM
	  if((curoct->y==0.)&&(inei==2)){
	    continue;
	  }
#endif

#ifdef TRANSYP
	  if(((curoct->y+2.*dxcur)==1.)&&(inei==3)){
	    continue;
	  }
#endif

#ifdef TRANSZM
	  if((curoct->z==0.)&&(inei==4)){
	    continue;
	  }
#endif

#ifdef TRANSZP
	  if(((curoct->z+2.*dxcur)==1.)&&(inei==5)){
	    continue;
	  }
#endif


#ifndef COARSERAD
	  if(vnei[inei]!=6){
	    if(curoct->nei[vnei[inei]]->child==NULL){
	      // the neighbor cell is unsplit, we update its value with fluxes

	      // initial data from the new value
	      memcpy(&R,&(curoct->nei[vnei[inei]]->rfieldnew),sizeof(struct Rtype));

	      //if(curoct->cell[icell].rfield.e[0]!=0.)  if(curoct->level==6) printf("R1=%e\n",R.e[0]);
	      // getting the flux

	      memcpy(F,stencil[iread].New.cell[icell].rflux+inei*NVAR_R*NGRP,sizeof(REAL)*NVAR_R*NGRP);

	      for(igrp=0;igrp<NGRP;igrp++){
		// update
		one=POW(-1.,inei+1);
		R.e[igrp] += F[0+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fx[igrp]+= F[1+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fy[igrp]+= F[2+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fz[igrp]+= F[3+igrp*NVAR_R]*dtsurdx*one*0.125;
	      }

	      memcpy(&(curoct->nei[vnei[inei]]->rfieldnew),&R,sizeof(struct Rtype));
	    }
	  }
#endif
	}

      }

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}

// ====================================================================================================================

int advancerad(struct OCT **firstoct, int level, struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dxcur, REAL dtnew,REAL aexp, struct RUNPARAMS *param, int chemonly){

  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot,nread;
  double t[10];
  double tg=0.,th=0.,tu=0.,ts=0.,tfu=0.,ttot=0.;
  REAL cloc; // the speed of light in code units

  cloc=aexp*param->clight*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v;
  //if(cpu->rank==RANK_DISP) printf("cloc=%e aexp=%e\n",cloc,aexp);
  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  nreadtot=0;
  //printf("rank %d nextoct=%p noct=%d\n",cpu->rank, firstoct[level-1],cpu->noct[level-1]);
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next;

      t[0]=MPI_Wtime();

      // ------------ gathering the stencil value values
      nextoct= gatherstencilrad(curoct,stencil,stride,cpu, &nread,cloc);

      //printf("nreadtot=%d on CPU=%d\n",nread,cpu->rank);


      if(nread>0){

	t[2]=MPI_Wtime();
	// ------------ solving the hydro

#ifndef COARSERAD
	int condadvec=1;
#else
	int condadvec=((level==param->lcoarse)&&(!chemonly));
#endif
	if(condadvec){
	  rad_sweepX(stencil,level,cpu->rank,nread,stride,dxcur,dtnew,cloc);
	  rad_sweepY(stencil,level,cpu->rank,nread,stride,dxcur,dtnew,cloc);
	  rad_sweepZ(stencil,level,cpu->rank,nread,stride,dxcur,dtnew,cloc);
	}
	else{
	  //printf("SKIP RAD TRANSPORT on level=%d\n on cpu %d",level,cpu->rank);
	}
	// ------------ updating values within the stencil

	t[4]=MPI_Wtime();

	if(condadvec) updatefieldrad(curoct,stencil,nread,stride,cpu,dxcur,dtnew,cloc);

	// ----------- perform physical cooling and ionisation

#ifdef WCHEM
	chemrad(stencil,nread,stride,cpu,dxcur,dtnew,param,aexp,chemonly);
#endif

	// ------------ scatter back the FLUXES

	t[6]=MPI_Wtime();

	nextoct=scatterstencilrad(curoct,stencil, nread, cpu,dxcur,dtnew,cloc);


	t[8]=MPI_Wtime();

	nreadtot+=nread;
	ts+=(t[8]-t[6]);
	tu+=(t[6]-t[4]);
	th+=(t[4]-t[2]);
	tg+=(t[2]-t[0]);
      }


    }while((nextoct!=NULL)&&(nread>0));
  }

  //if(cpu->rank==RANK_DISP) printf("CPU %d| tgat=%e tcal=%e tupchem=%e tscat=%e\n",cpu->rank,tg,th,tu,ts);

  return nreadtot;
}



// =================================================================================================
// =================================================================================================


REAL RadSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dtnew, REAL aexp, int chemonly){

  int nread,nreadtot;;
  struct OCT *curoct;
  struct OCT *nextoct;

  REAL dxcur=POW(0.5,level);
  REAL one;
  struct Rtype R;
  int icell;
  int nocthydro=cpu->noct[level-1];
  double t[10];
  int igrp;

  t[0]=MPI_Wtime();
#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&nocthydro,1,MPI_INT,MPI_SUM,cpu->comm);
#endif
  //if(cpu->rank==RANK_DISP)
    //printf("Start Radiation on %d octs with dt=%e on level %d with stride=%d and aexp=%e\n",nocthydro,dtnew,level,stride,aexp);

  // ===== COMPUTING THE FLUXES

#ifndef GPUAXL
  nreadtot=advancerad(firstoct,level,cpu,stencil,stride,dxcur,dtnew,aexp,param,chemonly);
#else
  nreadtot=advanceradGPU(firstoct,level,cpu,stencil,stride,dxcur,dtnew,aexp,param,chemonly);
#endif

  // FINAL UPDATE OF THE VALUES
  if(nreadtot>0){
    nextoct=firstoct[level-1];
    do {
      curoct=nextoct;
      nextoct=curoct->next;
#ifdef WMPI
      if(curoct->cpu!=cpu->rank) continue; // bnd octs are updated by transmission hence not required
#endif
      for(icell=0;icell<8;icell++){

	int is_unsplit;

#ifdef COARSERAD
	if(chemonly){
	  is_unsplit=(curoct->cell[icell].child==NULL);
	}
	else{
	  is_unsplit=(((curoct->cell[icell].child==NULL))||(level==param->lcoarse));
	}
#else
	is_unsplit=(curoct->cell[icell].child==NULL);
#endif

	if(is_unsplit){
	  // unsplit case

	  curoct->cell[icell].rfieldnew.src=curoct->cell[icell].rfield.src;
#ifdef WCHEM
	  curoct->cell[icell].rfieldnew.nh=curoct->cell[icell].rfield.nh;
	  E2T(&curoct->cell[icell].rfieldnew,aexp,param);
#endif

	  // Update
	  memcpy(&(curoct->cell[icell].rfield),&(curoct->cell[icell].rfieldnew),sizeof(struct Rtype));

#ifdef WRADHYD
	  // inject back thermal energy into the hydro

	  curoct->cell[icell].field.p=(GAMMA-1.)*curoct->cell[icell].rfield.eint;
	  curoct->cell[icell].field.dX=curoct->cell[icell].rfield.nhplus/curoct->cell[icell].rfield.nh*curoct->cell[icell].field.d;
	  getE(&curoct->cell[icell].field);
#endif

	}
	else{
	  // split case : -> quantities are averaged
	  struct OCT *child;
	  int i;
	  child=curoct->cell[icell].child;
	  memset(&R,0,sizeof(struct Rtype));

	  for(i=0;i<8;i++){
	    for(igrp=0;igrp<NGRP;igrp++){
	      R.e[igrp] +=child->cell[i].rfield.e[igrp]*0.125;
	      R.fx[igrp]+=child->cell[i].rfield.fx[igrp]*0.125;
	      R.fy[igrp]+=child->cell[i].rfield.fy[igrp]*0.125;
	      R.fz[igrp]+=child->cell[i].rfield.fz[igrp]*0.125;
	    }
	    R.src+=child->cell[i].rfield.src*0.125;
#ifdef WCHEM
	    //nh0+=(child->cell[i].rfield.xion*child->cell[i].rfield.nh)*0.125;
	    R.nhplus+=child->cell[i].rfield.nhplus*0.125;
	    R.eint+=child->cell[i].rfield.eint*0.125;
	    R.nh+=child->cell[i].rfield.nh*0.125;
#endif
	  }

#ifdef WCHEM
	  E2T(&R,aexp,param);
#endif
	  memcpy(&curoct->cell[icell].rfield,&R,sizeof(struct Rtype));

#ifdef WRADHYD
	    // inject back thermal energy into the hydro
	    curoct->cell[icell].field.p=(GAMMA-1.)*curoct->cell[icell].rfield.eint;
	    curoct->cell[icell].field.dX=curoct->cell[icell].rfield.nhplus/curoct->cell[icell].rfield.nh*curoct->cell[icell].field.d;
	    getE(&curoct->cell[icell].field);
#endif
	}


#ifdef COARSERAD
	if(((curoct->cell[icell].child!=NULL))&&(level==param->lcoarse)){
	  //if at levelcoarse and cell is splitted we propagate the data to the leaves
	  distribE(&curoct->cell[icell],&curoct->cell[icell],param,level);
	  }
#endif

      }

    }while(nextoct!=NULL);
  }

  t[9]=MPI_Wtime();
  return (REAL)(t[9]-t[0]);

}


// ==============================================================
void set_new_rad(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu, REAL aexp){

  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif

  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next;

      for(icell=0;icell<8;icell++) memcpy(&(curoct->cell[icell].rfieldnew),&(curoct->cell[icell].rfield),sizeof(struct Rtype));

    }while(nextoct!=NULL);
  }
}
// ==============================================================
void clean_new_rad(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu, REAL aexp){

  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif

  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next;

      for(icell=0;icell<8;icell++) memcpy(&(curoct->cell[icell].rfieldnew),&(curoct->cell[icell].rfield),sizeof(struct Rtype));

    }while(nextoct!=NULL);
  }

#ifdef WMPI
  // --------------- init for coarse octs in boundaries
  if(level>param->lcoarse){
    nextoct=firstoct[level-2];
    if((nextoct!=NULL)&&(cpu->noct[level-2]!=0)){
      do {
	curoct=nextoct;
	nextoct=curoct->next;
	if(curoct->cpu!=cpu->rank){
	  // Note: only boundary coarse octs are set to zero
	  for(icell=0;icell<8;icell++) {
	    memset(&(curoct->cell[icell].rfieldnew),0,sizeof(struct Rtype));
	  }
	}
      }while(nextoct!=NULL);
    }
  }
#endif


}


// ==============================================================
void sanity_rad(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu, REAL aexp){

//  printf("sanity rad on cpu %d\n", cpu->rank);
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  REAL cloc=aexp*param->clight*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v;

  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next;
       // sanity check

      REAL E;
      REAL F;
      for(icell=0;icell<8;icell++) {
	E=curoct->cell[icell].rfield.e[0]*cloc;
	F=SQRT(POW(curoct->cell[icell].rfield.fx[0],2)+POW(curoct->cell[icell].rfield.fy[0],2)+POW(curoct->cell[icell].rfield.fz[0],2));
	if(F/E>1.0){
	  printf("E=%e F=%e cloc=%e cE=%e\n",E/cloc,F,cloc,E);
	  abort();
	}
      }

    }while(nextoct!=NULL);
  }
}




#endif
