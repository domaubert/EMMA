
#ifdef WRAD
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>

#include <mpi.h>

#ifdef WCHEM
#include "chem_utils.h"
#endif



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
#ifdef WCHEM
    WR->xion=W2->xion-W1->xion;
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
      Wr->e[igrp]=fmax(fmax(0.,fmin(beta*Wm->e[igrp],Wp->e[igrp])),fmin(Wm->e[igrp],beta*Wp->e[igrp]));
    }
    else{
      Wr->e[igrp]=fmin(fmin(0.,fmax(beta*Wm->e[igrp],Wp->e[igrp])),fmax(Wm->e[igrp],beta*Wp->e[igrp]));
    }


    if(Wp->fx[igrp]>0){
      Wr->fx[igrp]=fmax(fmax(0.,fmin(beta*Wm->fx[igrp],Wp->fx[igrp])),fmin(Wm->fx[igrp],beta*Wp->fx[igrp]));
    }
    else{
      Wr->fx[igrp]=fmin(fmin(0.,fmax(beta*Wm->fx[igrp],Wp->fx[igrp])),fmax(Wm->fx[igrp],beta*Wp->fx[igrp]));
    }


    if(Wp->fy[igrp]>0){
      Wr->fy[igrp]=fmax(fmax(0.,fmin(beta*Wm->fy[igrp],Wp->fy[igrp])),fmin(Wm->fy[igrp],beta*Wp->fy[igrp]));
    }
    else{
      Wr->fy[igrp]=fmin(fmin(0.,fmax(beta*Wm->fy[igrp],Wp->fy[igrp])),fmax(Wm->fy[igrp],beta*Wp->fy[igrp]));
    }


    if(Wp->fz[igrp]>0){
      Wr->fz[igrp]=fmax(fmax(0.,fmin(beta*Wm->fz[igrp],Wp->fz[igrp])),fmin(Wm->fz[igrp],beta*Wp->fz[igrp]));
    }
    else{
      Wr->fz[igrp]=fmin(fmin(0.,fmax(beta*Wm->fz[igrp],Wp->fz[igrp])),fmax(Wm->fz[igrp],beta*Wp->fz[igrp]));
    }

  }

  if(Wp->src>0){
      Wr->src=fmax(fmax(0.,fmin(beta*Wm->src,Wp->src)),fmin(Wm->src,beta*Wp->src));
    }
    else{
      Wr->src=fmin(fmin(0.,fmax(beta*Wm->src,Wp->src)),fmax(Wm->src,beta*Wp->src));
    }
#ifdef WCHEM
  if(Wp->xion>0){
    Wr->xion=fmax(fmax(0.,fmin(beta*Wm->xion,Wp->xion)),fmin(Wm->xion,beta*Wp->xion));
  }
  else{
    Wr->xion=fmin(fmin(0.,fmax(beta*Wm->xion,Wp->xion)),fmax(Wm->xion,beta*Wp->xion));
  }

  if(Wp->eint>0){
    Wr->eint=fmax(fmax(0.,fmin(beta*Wm->eint,Wp->eint)),fmin(Wm->eint,beta*Wp->eint));
  }
  else{
    Wr->eint=fmin(fmin(0.,fmax(beta*Wm->eint,Wp->eint)),fmax(Wm->eint,beta*Wp->eint));
  }

  if(Wp->nh>0){
    Wr->nh=fmax(fmax(0.,fmin(beta*Wm->nh,Wp->nh)),fmin(Wm->nh,beta*Wp->nh));
  }
  else{
    Wr->nh=fmin(fmin(0.,fmax(beta*Wm->nh,Wp->nh)),fmax(Wm->nh,beta*Wp->nh));
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
    Wp->xion =W0->xion +dx*Dx->xion +dy*Dy->xion +dz*Dz->xion;
    Wp->eint =W0->eint +dx*Dx->eint +dy*Dy->eint +dz*Dz->eint;
    Wp->nh =W0->nh +dx*Dx->nh +dy*Dy->nh +dz*Dz->nh;
#endif
}

//================================================================================

void coarse2fine_rad2(struct CELL *cell, struct Rtype *Wi){ 

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
	  dxcur=pow(0.5,oct->level);

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

	  for(iz=0;iz<2;iz++){
	    for(iy=0;iy<2;iy++){
	      for(ix=0;ix<2;ix++){
		icell=ix+iy*2+iz*4;
		interpminmod_R(W0,&Wint,D,D+1,D+2,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Wp contains the interpolation
		memcpy(Wi+icell,&Wint,sizeof(struct Rtype));
	      }
	    }
	  }
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
	  dxcur=pow(0.5,oct->level);

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
      ff=sqrt(fx*fx+fy*fy+fz*fz); // 6 flop
      if(ff>0)
	{
	  n[0]=fx/ff; 
	  n[1]=fy/ff;
	  n[2]=fz/ff; 
	}
      ff=ff/(c*ee); // 2flop
    }
  
  arg=fmax(4.-3.*ff*ff,0.); // 4 flop
  chi=(3.+4.*ff*ff)/(5.+2.*sqrt(arg)); // 7 flops

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

	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }

      // X DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems LEFT

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


      for(igrp=0;igrp<NGRP;igrp++){
	for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0]; 
	for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1]; 
      }

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

	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
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

	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
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

void recursive_neighbor_gather_oct_rad(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct RGRID *stencil,char *visit){


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

  struct Rtype Ri[8];
  char child[8];

  struct OCT *oct;
  struct OCT *neioct;
  struct CELL *neicell;

  ineiloc=inei;
  

  if(cell->child!=NULL){
    // the oct at the right level exists
    neicell=cell->child->nei[ineiloc];
    oct=cell->child;
    dxcur=pow(0.5,oct->level);

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
    dxcur=pow(0.5,oct->level);
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
      memcpy(Ri+icell,&(neicell->child->cell[icell].rfield),sizeof(struct Rtype));
      child[icell]=(neicell->child->cell[icell].child!=NULL);
    }

  }
  else{
    //coarse2fine_radlin(neicell,Ri);
    int il; 
    for(il=0;il<8;il++) memcpy(&Ri[il],&neicell->rfield,sizeof(struct Rtype)); 

    for(icell=0;icell<8;icell++){
      child[icell]=0;
    }
  }


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
  
}

// ===================================================================================================
// ===================================================================================================

struct OCT *gatherstencilrad(struct OCT *octstart, struct RGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  char visit[27]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
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
	  recursive_neighbor_gather_oct_rad(ioct, inei, -1, -1, 1, cell, stencil+iread,visit);
	}
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}


// ===================================================================================================
// ===================================================================================================

void updatefieldrad(struct OCT *octstart, struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew)
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
      
      if(stencil[i].oct[6].cell[icell].split) continue;
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

      /* memcpy(&(stencil[i].New.cell[icell].deltaR),&R,sizeof(struct Rtype)); */

      // TESTING FULL UPDATE IN STENCIL APPROACH
      
      //memcpy(&Rupdate,&stencil[i].oct[6].cell[icell].rfield,sizeof(struct Rtype));
      memcpy(&Rupdate,&stencil[i].New.cell[icell].rfieldnew,sizeof(struct Rtype));

      for(igrp=0;igrp<NGRP;igrp++){
	Rupdate.e[igrp]   +=R.e[igrp];
	Rupdate.fx[igrp]  +=R.fx[igrp];
	Rupdate.fy[igrp]  +=R.fy[igrp];
	Rupdate.fz[igrp]  +=R.fz[igrp];
      }
      
      if(Rupdate.e[0]<0) abort();
      //memcpy(&(curoct->cell[icell].rfieldnew),&Rupdate,sizeof(struct Rtype));
      memcpy(&stencil[i].New.cell[icell].rfieldnew,&Rupdate,sizeof(struct Rtype));

    }
  }

}

// ===========================================================================================================

struct OCT *scatterstencilrad(struct OCT *octstart, struct RGRID *stencil, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew)
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

	/* if(curoct->x<1./32.) */
	/*   if(icell==2)  */
	/*     if(curoct->cell[0].rfieldnew.xion!=curoct->cell[2].rfieldnew.xion) abort(); */
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


	  if(vnei[inei]!=6){
	    if(curoct->nei[vnei[inei]]->child==NULL){
	      // the neighbor cell is unsplit we update its value with fluxes
	   
	      // initial data from the new value
	      memcpy(&R,&(curoct->nei[vnei[inei]]->rfieldnew),sizeof(struct Rtype));

	      //if(curoct->cell[icell].rfield.e[0]!=0.)  if(curoct->level==6) printf("R1=%e\n",R.e[0]);
	      // getting the flux

	      memcpy(F,stencil[iread].New.cell[icell].rflux+inei*NVAR_R*NGRP,sizeof(REAL)*NVAR_R*NGRP);
	      
	      for(igrp=0;igrp<NGRP;igrp++){
		// update
		one=pow(-1.,inei+1);
		R.e[igrp] += F[0+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fx[igrp]+= F[1+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fy[igrp]+= F[2+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fz[igrp]+= F[3+igrp*NVAR_R]*dtsurdx*one*0.125;
	      }
 	      
	      memcpy(&(curoct->nei[vnei[inei]]->rfieldnew),&R,sizeof(struct Rtype));
	      //if(curoct->cell[icell].rfield.e[0]!=0.) if(curoct->level==6) printf("R2=%e\n",R.e[0]);
	      
	    }
	  }
	}

      }

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}

// ====================================================================================================================

int advancerad(struct OCT **firstoct, int level, struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dxcur, REAL dtnew,REAL aexp, struct RUNPARAMS *param){

  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot,nread;
  double t[10];
  double tg=0.,th=0.,tu=0.,ts=0.,tfu=0.,ttot=0.;
  REAL cloc; // the speed of light in code units

  cloc=aexp*param->clight*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v;
  if(cpu->rank==0) printf("cloc=%e aexp=%e\n",cloc,aexp);
  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  nreadtot=0;
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next; 

      t[0]=MPI_Wtime();
  
      // ------------ gathering the stencil value values
      nextoct= gatherstencilrad(curoct,stencil,stride,cpu, &nread);
      
      t[2]=MPI_Wtime();
      // ------------ solving the hydro
      rad_sweepX(stencil,level,cpu->rank,nread,stride,dxcur,dtnew,cloc);   
      rad_sweepY(stencil,level,cpu->rank,nread,stride,dxcur,dtnew,cloc); 
      rad_sweepZ(stencil,level,cpu->rank,nread,stride,dxcur,dtnew,cloc); 
      
      // ------------ updating values within the stencil
      
      t[4]=MPI_Wtime();
      
      updatefieldrad(curoct,stencil,nread,stride,cpu,dxcur,dtnew);

      // ----------- perform physical cooling and ionisation 
#ifdef WCHEM
      chemrad(curoct,stencil,nread,stride,cpu,dxcur,dtnew,param,aexp);
#endif


      // ------------ scatter back the FLUXES
      
      t[6]=MPI_Wtime();
   
      nextoct=scatterstencilrad(curoct,stencil, nread, cpu,dxcur,dtnew);


      t[8]=MPI_Wtime();

      nreadtot+=nread;


      ts+=(t[8]-t[6]);
      tu+=(t[6]-t[4]);
      th+=(t[4]-t[2]);
      tg+=(t[2]-t[0]);

    }while(nextoct!=NULL);
  }
  
  if(cpu->rank==0) printf("CPU | tgat=%e tcal=%e tup=%e tscat=%e\n",tg,th,tu,ts);

  return nreadtot;
}



// =================================================================================================
// =================================================================================================


void RadSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dtnew, REAL aexp){
  
  int nread,nreadtot;;
  struct OCT *curoct;
  struct OCT *nextoct;
  
  REAL dxcur=pow(0.5,level);
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
  if(cpu->rank==0) printf("Start Radiation on %d octs with dt=%e on level %d with stride=%d and aexp=%e\n",nocthydro,dtnew,level,stride,aexp);

  // ===== COMPUTING THE FLUXES
  
#ifndef GPUAXL
  nreadtot=advancerad(firstoct,level,cpu,stencil,stride,dxcur,dtnew,aexp,param);
#else
  nreadtot=advanceradGPU(firstoct,level,cpu,stencil,stride,dxcur,dtnew,aexp,param);
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
	if(curoct->cell[icell].child==NULL){
	  // unsplit case

	  curoct->cell[icell].rfieldnew.src=curoct->cell[icell].rfield.src;
#ifdef WCHEM
	  curoct->cell[icell].rfieldnew.nh=curoct->cell[icell].rfield.nh;
	  E2T(&curoct->cell[icell].rfieldnew,aexp,param);
#endif

	  //if(curoct->cell[icell].rfieldnew.temp!=curoct->cell[icell].rfield.temp) abort();
	  memcpy(&(curoct->cell[icell].rfield),&(curoct->cell[icell].rfieldnew),sizeof(struct Rtype));
#ifdef WRADHYD
	  // inject back thermal energy into the hydro
	  curoct->cell[icell].field.p=(GAMMA-1.)*curoct->cell[icell].rfield.eint;
	  curoct->cell[icell].field.X=curoct->cell[icell].rfield.xion;
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
	    R.xion+=child->cell[i].rfield.xion*0.125;
	    R.eint+=child->cell[i].rfield.eint*0.125;
	    R.nh+=child->cell[i].rfield.nh*0.125;
#endif

	  }

#ifdef WCHEM
	  E2T(&R,aexp,param);
#endif
	  memcpy(&curoct->cell[icell].rfield,&R,sizeof(struct Rtype));
	}
      }
    }while(nextoct!=NULL);
  }

  t[9]=MPI_Wtime();
  if(cpu->rank==0){
#ifndef GPUAXL
    printf("==== CPU RAD TOTAL TIME =%e\n",t[9]-t[0]);
#else
    printf(" === GPU RAD TOTAL TIME =%e\n",t[9]-t[0]);
#endif
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

      //for(icell=0;icell<8;icell++) memset(&(curoct->cell[icell].fieldnew),0,sizeof(struct Wtype));
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
	F=sqrt(pow(curoct->cell[icell].rfield.fx[0],2)+pow(curoct->cell[icell].rfield.fy[0],2)+pow(curoct->cell[icell].rfield.fz[0],2));
	if(F/E>1.0){
	  abort();
	}
      }
      
    }while(nextoct!=NULL);
  }
}




#endif
