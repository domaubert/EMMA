
#ifdef WHYDRO2

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include "vector.h"
#include <mpi.h>

#ifdef GPUAXL
#include "hydro_utils_gpu.h"
#endif

#define NITERMAX 10
#define ERRTOL 1e-10
#define FRACP 1e-3


// ===================================================
//
// Note: Pressure estimates are kept in double precision for better accuracy
//
// ===================================================

// ===================================================

void getE(struct Wtype *W){
  W->E=W->p/(GAMMA-1.)+0.5*W->d*(W->u*W->u+W->v*W->v+W->w*W->w);
}


// ==================== converts U -> W
void U2W(struct Utype *U, struct Wtype *W)
{
  W->d=U->d;
  W->u=U->du/U->d;
  W->v=U->dv/U->d;
  W->w=U->dw/U->d;
  
#ifdef DUAL_E
  W->p=U->eint*(GAMMA-1.);
  W->E=U->E;
#ifdef WRADHYD
  W->dX=U->dX;
#endif
#else
  W->p=(GAMMA-1.)*(U->E-((U->du)*(U->du)+(U->dv)*(U->dv)+(U->dw)*(U->dw))/(U->d)*0.5);
#endif
  W->a=SQRT(GAMMA*W->p/W->d);
}

// ==================== converts W -> U
void W2U(struct Wtype *W, struct Utype *U)
{
  U->d=W->d;
  U->du=W->d*W->u;
  U->dv=W->d*W->v;
  U->dw=W->d*W->w;

#ifdef DUAL_E
  U->eint=W->p/(GAMMA-1.);
  U->E=W->E;

#ifdef WRADHYD
  U->dX=W->dX;
#endif
#endif

}


// ================== performs the difference between two Us

void diffU(struct Utype *U2, struct Utype *U1, struct Utype *UR){
  
  UR->d =U2->d - U1->d;
  UR->du=U2->du- U1->du;
  UR->dv=U2->dv- U1->dv;
  UR->dw=U2->dw- U1->dw;
  UR->E =U2->E - U1->E;
#ifdef DUAL_E
  UR->eint=U2->eint-U1->eint;
#endif
}

// ================== performs the difference between two Ws

void diffW(struct Wtype *W2, struct Wtype *W1, struct Wtype *WR){

  WR->d=W2->d- W1->d;
  WR->u=W2->u- W1->u;
  WR->v=W2->v- W1->v;
  WR->w=W2->w- W1->w;
  WR->p=W2->p- W1->p;
#ifdef WRADHYD
  WR->dX=W2->dX- W1->dX;
#endif
}


// ================= minmod

void minmod(struct Utype *Um, struct Utype *Up, struct Utype *Ur){

  REAL beta=1.; // 1. for MINBEE 2. for SUPERBEE
  // FLUX LIMITER

  if(Up->d>0){
    Ur->d=FMAX(FMAX(0.,FMIN(beta*Um->d,Up->d)),FMIN(Um->d,beta*Up->d));
  }
  else{
    Ur->d=FMIN(FMIN(0.,FMAX(beta*Um->d,Up->d)),FMAX(Um->d,beta*Up->d));
  }


  if(Up->du>0){
    Ur->du=FMAX(FMAX(0.,FMIN(beta*Um->du,Up->du)),FMIN(Um->du,beta*Up->du));
  }
  else{
    Ur->du=FMIN(FMIN(0.,FMAX(beta*Um->du,Up->du)),FMAX(Um->du,beta*Up->du));
  }


  if(Up->dv>0){
    Ur->dv=FMAX(FMAX(0.,FMIN(beta*Um->dv,Up->dv)),FMIN(Um->dv,beta*Up->dv));
  }
  else{
    Ur->dv=FMIN(FMIN(0.,FMAX(beta*Um->dv,Up->dv)),FMAX(Um->dv,beta*Up->dv));
  }


  if(Up->dw>0){
    Ur->dw=FMAX(FMAX(0.,FMIN(beta*Um->dw,Up->dw)),FMIN(Um->dw,beta*Up->dw));
  }
  else{
    Ur->dw=FMIN(FMIN(0.,FMAX(beta*Um->dw,Up->dw)),FMAX(Um->dw,beta*Up->dw));
  }


  if(Up->E>0){
    Ur->E=FMAX(FMAX(0.,FMIN(beta*Um->E,Up->E)),FMIN(Um->E,beta*Up->E));
  }
  else{
    Ur->E=FMIN(FMIN(0.,FMAX(beta*Um->E,Up->E)),FMAX(Um->E,beta*Up->E));
  }


}



//===============================================
//===============================================

void minmod_W(struct Wtype *Wm, struct Wtype *Wp, struct Wtype *Wr){

  REAL beta=1.; // 1. for MINBEE 2. for SUPERBEE
  // FLUX LIMITER

  if(Wp->d>0){
    Wr->d=FMAX(FMAX(0.,FMIN(beta*Wm->d,Wp->d)),FMIN(Wm->d,beta*Wp->d));
  }
  else{
    Wr->d=FMIN(FMIN(0.,FMAX(beta*Wm->d,Wp->d)),FMAX(Wm->d,beta*Wp->d));
  }

#ifdef WRADHYD
  if(Wp->dX>0){
    Wr->dX=FMAX(FMAX(0.,FMIN(beta*Wm->dX,Wp->dX)),FMIN(Wm->dX,beta*Wp->dX));
  }
  else{
    Wr->dX=FMIN(FMIN(0.,FMAX(beta*Wm->dX,Wp->dX)),FMAX(Wm->dX,beta*Wp->dX));
  }
#endif

  if(Wp->u>0){
    Wr->u=FMAX(FMAX(0.,FMIN(beta*Wm->u,Wp->u)),FMIN(Wm->u,beta*Wp->u));
  }
  else{
    Wr->u=FMIN(FMIN(0.,FMAX(beta*Wm->u,Wp->u)),FMAX(Wm->u,beta*Wp->u));
  }


  if(Wp->v>0){
    Wr->v=FMAX(FMAX(0.,FMIN(beta*Wm->v,Wp->v)),FMIN(Wm->v,beta*Wp->v));
  }
  else{
    Wr->v=FMIN(FMIN(0.,FMAX(beta*Wm->v,Wp->v)),FMAX(Wm->v,beta*Wp->v));
  }


  if(Wp->w>0){
    Wr->w=FMAX(FMAX(0.,FMIN(beta*Wm->w,Wp->w)),FMIN(Wm->w,beta*Wp->w));
  }
  else{
    Wr->w=FMIN(FMIN(0.,FMAX(beta*Wm->w,Wp->w)),FMAX(Wm->w,beta*Wp->w));
  }


  if(Wp->p>0){
    Wr->p=FMAX(FMAX(0.,FMIN(beta*Wm->p,Wp->p)),FMIN(Wm->p,beta*Wp->p));
  }
  else{
    Wr->p=FMIN(FMIN(0.,FMAX(beta*Wm->p,Wp->p)),FMAX(Wm->p,beta*Wp->p));
  }


}






// ============= interp minmod

void interpminmod(struct Utype *U0, struct Utype *Up, struct Utype *Dx, struct Utype *Dy, struct Utype *Dz,REAL dx,REAL dy,REAL dz){
  
  Up->d =U0->d  + dx*Dx->d  +dy*Dy->d  +dz*Dz->d;
  Up->du=U0->du + dx*Dx->du +dy*Dy->du +dz*Dz->du;
  Up->dv=U0->dv + dx*Dx->dv +dy*Dy->dv +dz*Dz->dv;
  Up->dw=U0->dw + dx*Dx->dw +dy*Dy->dw +dz*Dz->dw;
  Up->E =U0->E  + dx*Dx->E  +dy*Dy->E  +dz*Dz->E;
#ifdef DUAL_E
  Up->eint =U0->eint  + dx*Dx->eint  +dy*Dy->eint  +dz*Dz->eint;
#endif
}


void interpminmod_W(struct Wtype *W0, struct Wtype *Wp, struct Wtype *Dx, struct Wtype *Dy, struct Wtype *Dz,REAL dx,REAL dy,REAL dz){
  
  Wp->d =W0->d +dx*Dx->d +dy*Dy->d +dz*Dz->d;
  Wp->u =W0->u +dx*Dx->u +dy*Dy->u +dz*Dz->u;
  Wp->v =W0->v +dx*Dx->v +dy*Dy->v +dz*Dz->v;
  Wp->w =W0->w +dx*Dx->w +dy*Dy->w +dz*Dz->w;
  Wp->p =W0->p +dx*Dx->p +dy*Dy->p +dz*Dz->p;
#ifdef WRADHYD
  Wp->dX =W0->dX +dx*Dx->dX +dy*Dy->dX +dz*Dz->dX;
#endif
}


// ==============================================


// ==============================================

void coarse2fine_hydro2(struct CELL *cell, struct Wtype *Wi){ 


	  struct OCT * oct;
	  
	  struct Wtype *W0;
	  struct Wtype *Wp;
	  struct Wtype *Wm;
	  struct Wtype Wint;
	  struct Wtype Dp,Dm;
	  struct Wtype D[3];
	  struct Wtype *W;
	  int inei2;
	  int vcell[6],vnei[6];
	  int dir;
	  REAL dxcur;

	  oct=cell2oct(cell);
	  getcellnei(cell->idx, vnei, vcell); // we get the neighbors
	  dxcur=POW(0.5,oct->level);

	  W0=&(cell->field);
	  // Limited Slopes
	  for(dir=0;dir<3;dir++){
	    
	    inei2=2*dir;

	    
	    if(vnei[inei2]==6){
	      Wm=&(oct->cell[vcell[inei2]].field);
	    }
	    else{
	      if(oct->nei[vnei[inei2]]->child!=NULL){
		Wm=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].field);
	      }
	    }
	    
#ifdef TRANSXM
	    if((oct->x==0.)&&(inei2==0)){
	      Wm=&(cell->field);
	    }
#endif
	    
#ifdef TRANSYM
	    if((oct->y==0.)&&(inei2==2)){
	      Wm=&(cell->field);
	    }
#endif
	    
#ifdef TRANSZM
	    if((oct->z==0.)&&(inei2==4)){
	      Wm=&(cell->field);
	    }
#endif
	    
	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      Wp=&(oct->cell[vcell[inei2]].field);
	    }
	    else{
	      if(oct->nei[vnei[inei2]]->child!=NULL){
		Wp=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].field);
	      }
	    }


#ifdef TRANSXP
	    if(((oct->x+2.*dxcur)==1.)&&(inei2==1)){
	      Wp=&(cell->field);
	    }
#endif
	    
#ifdef TRANSYP
	    if(((oct->y+2.*dxcur)==1.)&&(inei2==3)){
	      Wp=&(cell->field);
	    }
#endif
	    
#ifdef TRANSZP
	    if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
	      Wp=&(cell->field);
	    }
#endif


	    diffW(Wp,W0,&Dp); 
	    diffW(W0,Wm,&Dm); 
	    
	    
	    minmod_W(&Dm,&Dp,D+dir);
	    
		  
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
		interpminmod_W(W0,&Wint,D,D+1,D+2,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Wp contains the interpolation
		getE(&Wint);
		Wint.a=SQRT(GAMMA*Wint.p/Wint.d);
		memcpy(Wi+icell,&Wint,sizeof(struct Wtype));

	      }
	    }
	  }

}

// =====================================================================================

void coarse2fine_hydrolin(struct CELL *cell, struct Wtype *Wi){ 


	  struct OCT * oct;
	  
	  struct Wtype *W0;
	  struct Wtype *Wp;
	  struct Wtype *Wm;
	  struct Wtype Wint;
	  struct Wtype *W;


	  struct Utype U0;
	  struct Utype Up;
	  struct Utype Um;
	  struct Utype Uint;

	  /* struct Wtype Dp,Dm; */
	  /* struct Wtype D[6]; */
#ifdef LINU
	  struct Utype D[6];
#else
	  struct Wtype D[6];
#endif
	  
	  int inei2;
	  int vcell[6],vnei[6];
	  int dir;
	  REAL dxcur;

	  oct=cell2oct(cell);
	  getcellnei(cell->idx, vnei, vcell); // we get the neighbors
	  dxcur=POW(0.5,oct->level);

	  W0=&(cell->field);
	  // Limited Slopes
	  for(dir=0;dir<3;dir++){
	    
	    inei2=2*dir;
	    if(vnei[inei2]==6){
	      Wm=&(oct->cell[vcell[inei2]].field);
	    }
	    else{
	      Wm=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].field);

	    }
#ifdef TRANSXM
	      if((oct->x==0.)&&(inei2==0)){
		Wm=&(cell->field);
	      }
#endif

#ifdef TRANSYM
	      if((oct->y==0.)&&(inei2==2)){
		Wm=&(cell->field);
	      }
#endif

#ifdef TRANSZM
	      if((oct->z==0.)&&(inei2==4)){
		Wm=&(cell->field);
	      }
#endif


	    

	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      Wp=&(oct->cell[vcell[inei2]].field);
	    }
	    else{
	      Wp=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].field);

	    }
#ifdef TRANSXP
	      if(((oct->x+2.*dxcur)==1.)&&(inei2==1)){
		Wp=&(cell->field);
	      }
#endif

#ifdef TRANSYP
	      if(((oct->y+2.*dxcur)==1.)&&(inei2==3)){
		Wp=&(cell->field);
	      }
#endif

#ifdef TRANSZP
	      if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
		Wp=&(cell->field);
	      }
#endif


#ifdef LINU
	    W2U(W0,&U0);
	    W2U(Wm,&Um);
	    W2U(Wp,&Up);

	    diffU(&U0,&Um,D+2*dir+0); 
	    diffU(&Up,&U0,D+2*dir+1); 
#else
	    diffW(W0,Wm,D+2*dir+0); 
	    diffW(Wp,W0,D+2*dir+1); 
#endif	    



	    //	    minmod_W(&Dm,&Dp,D+dir);
	    
		  
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
#ifdef LINU
		interpminmod(&U0,&Uint,D+ix,D+iy+2,D+iz+4,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Wp contains the interpolation
		U2W(&Uint,&Wint);
#else
		interpminmod_W(W0,&Wint,D+ix,D+iy+2,D+iz+4,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Wp contains the interpolation
#endif


		getE(&Wint);
		Wint.a=SQRT(GAMMA*Wint.p/Wint.d);
		memcpy(Wi+icell,&Wint,sizeof(struct Wtype));

	      }
	    }
	  }

}




// ==================== pressure solver

double frootprime(double p, struct Wtype1D *WL, struct Wtype1D *WR)
{
  
  double fL,fR;
  double AL,AR,BL,BR;

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  fL=(p>WL->p?sqrt(AL/(BL+p))*(1.-(p-WL->p)/(2.*(BL+p))):pow(p/WL->p,-(GAMMA+1)/(2.*GAMMA))/(WL->d*WL->a));
  fR=(p>WR->p?sqrt(AR/(BR+p))*(1.-(p-WR->p)/(2.*(BR+p))):pow(p/WR->p,-(GAMMA+1)/(2.*GAMMA))/(WR->d*WR->a));

  return fL+fR;
}


// ------------------------------------

double froot(double p, struct Wtype1D *WL, struct Wtype1D *WR, double *u)
{
  
  double fL,fR;
  double AL,AR,BL,BR;
  double Deltau;

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  fL=(p>WL->p?(p-WL->p)*sqrt(AL/(BL+p)):2.*WL->a/(GAMMA-1.)*(pow(p/WL->p,(GAMMA-1)/(2.*GAMMA))-1.));
  fR=(p>WR->p?(p-WR->p)*sqrt(AR/(BR+p)):2.*WR->a/(GAMMA-1.)*(pow(p/WR->p,(GAMMA-1)/(2.*GAMMA))-1.));
  
  Deltau=WR->u-WL->u;
  *u=0.5*(WL->u+WR->u)+0.5*(fR-fL);

  return fL+fR+Deltau;
}



double findPressure_Hybrid(struct Wtype1D *WL, struct Wtype1D *WR, int *niter, REAL *ustar){
  double ppvrs;
  double dbar,abar;
  double pmax,pmin,pstar;
  double AL,AR,BL,BR,GL,GR;
  int cas;
  dbar=0.5*(WL->d+WR->d);
  abar=0.5*(WL->a+WR->a);
  ppvrs=0.5*((WL->p+WR->p)+(WL->u-WR->u)*dbar*abar);
  pmax=fmax(WL->p,WR->p);
  pmin=fmin(WL->p,WR->p);
  pstar=ppvrs;
  
  if(((pmax/pmin)<2.)&&((pmin<pstar)&&(pstar<pmax))){
    // PVRS CASE
    pstar=ppvrs;
    *ustar=0.5*((WL->u+WR->u)+(WL->p-WR->p)/(dbar*abar));
    cas=0;
  }
  else{
    if(pstar<pmin){
      //TRRS CASE
      double z=(GAMMA-1.)/(2.*GAMMA);
      double iz=(2.*GAMMA)/(GAMMA-1.);
      pstar=pow((WL->a+WR->a-(GAMMA-1.)/2.*(WR->u-WL->u))/(WL->a/pow(WL->p,z)+WR->a/pow(WR->p,z)),iz);
      *ustar=WL->u-2.*WL->a/(GAMMA-1.)*(pow(pstar/WL->p,z)-1.);
      cas=1;
    }
    else{
      //TSRS CASE
      double p0;
      p0=fmax(0.,ppvrs);
      
      AL=2./((GAMMA+1.)*WL->d);
      AR=2./((GAMMA+1.)*WR->d);
      
      BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
      BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

      GL=sqrt(AL/(p0+BL));
      GR=sqrt(AR/(p0+BR));

      pstar=(GL*WL->p+GR*WR->p-(WR->u-WL->u))/(GL+GR);
      *ustar=0.5*((WL->u+WR->u)+(pstar-WR->p)*GR-(pstar-WL->p)*GL);
      cas=2;
    }
  }

  return pstar;

}


// --------------------------------------

double findPressure(struct Wtype1D *WL, struct Wtype1D *WR, int *niter, REAL *u)
{

  double ptr,pts,ppv;
  double ptr0,pts0,ppv0;
  double p,porg,dp;
  int i;
  double err;
  double unsurz=(2.0*GAMMA)/(GAMMA-1.0);
  double AL,AR,BL,BR,GL,GR;
  double pmin,pmax;
  int tag;
  double u2;

  pmin=fmin(WL->p,WR->p);
  pmax=fmax(WL->p,WR->p);
  
  // EXACT SOLVER

  // hybrid guess for pressure

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  ppv0=0.5*(WL->p+WR->p)-0.125*(WR->u-WL->u)*(WR->d+WL->d)*(WR->a+WL->a);
  ptr0=pow((WL->a+WR->a-0.5*(GAMMA-1)*(WR->u-WL->u))/(WL->a/pow(WL->p,1./unsurz)+WR->a/pow(WR->p,1./unsurz)),unsurz);

  ppv=fmax(ERRTOL,ppv0);
  ptr=fmax(ERRTOL,ptr0);
  
  GL=sqrt(AL/(ppv+BL));
  GR=sqrt(AR/(ppv+BR));

  pts0=(GL*WL->p+GR*WR->p-(WR->u-WL->u))/(GL+GR);
  pts=fmax(ERRTOL,pts0);


  if(((pmax/pmin)<2.0)&&((pmin<=ppv)&&(ppv<=pmax))){
      p=ppv;
      tag=1;
    }
  else{
    if(ppv<pmin){
      p=ptr;
      tag=2;
    }
    else{
      p=pts;
      tag=3;
    }
  }


  //p=0.5*(WL->p+WR->p);
  //p=fmax(p,ERRTOL);

  REAL  p0=p;
  *niter=0;
  for(i=0;i<NITERMAX;i++)
    {
      dp=froot(p,WL,WR,&u2)/frootprime(p,WL,WR);///fmax(0.8,fmin(1.2,(1.0-0.5*froot(p,WL,WR,u)*frootprimeprime(p,WL,WR)/frootprime(p,WL,WR))));
      if((isnan(dp))){
      	printf("froot=%e frootprime=%e\n",froot(p,WL,WR,&u2),frootprime(p,WL,WR));
      	abort();
      }
      
      if(fabs(dp)<ERRTOL) break;
      while((p-dp)<0){ 
       	dp=dp*0.5; 
      } 

      porg=p;
      p=p-dp;
      //if(frootprime(p,WL,WR)==0) abort();//printf("p0=%e dp=%e p=%e fprime=%e\n",porg,dp,p,frootprime(p,WL,WR));
      err=2.*fabs(p-porg)/(fabs(p+porg));
      *niter=*niter+1;
      //if(p<=0) p=ERRTOL;
      if(err<ERRTOL) break;
      if(froot(p,WL,WR,&u2)<ERRTOL) break;
    }

  if(i==NITERMAX){
    //printf("DIVERGENCE p0=%e dp=%e p=%e fprime=%e err=%e\n",porg,dp,p,frootprime(p,WL,WR),err);
    //abort();
  }

  /* if(p>6.4e-7){ */
  /*   printf("MAX p0=%e dp=%e p=%e fprime=%e err=%e\n",porg,dp,p,frootprime(p,WL,WR),err); */
  /*   abort(); */
  /* } */
  froot(p,WL,WR,&u2); // last calculation to get u;

  *u=(REAL)u2;
  return p;
}




// ================================= Exact Riemann Solver

void getW(struct Wtype1D *W, REAL S, struct Wtype1D_double *WL, struct Wtype1D_double *WR, REAL pstar, REAL ustar)
{
  REAL SL,SR;
  REAL SHL,STL;
  REAL SHR,STR;
  REAL astar;

  if(S<ustar)
    {
      // left of contact
      if(pstar>WL->p)
	{
	  // left Shock with shock speed SL
	  SL=WL->u-WL->a*SQRT((GAMMA+1.)/(2.*GAMMA)*pstar/WL->p+(GAMMA-1.)/(2.*GAMMA)); 
	  
	  if(S<SL)
	    {
	      // left to the shock
	      W->d=WL->d;
	      W->u=WL->u;
	      W->p=WL->p;
	      W->a=WL->a;
	    }
	  else
	    {
	      W->d=WL->d*((pstar/WL->p+(GAMMA-1.)/(GAMMA+1.))/((GAMMA-1.)/(GAMMA+1.)*pstar/WL->p+1.));
	      W->u=ustar;
	      W->p=pstar;
	      W->a=SQRT(GAMMA*W->p/W->d);
	    }
	}
      else
	{
	  // left fan
	  astar=WL->a*POW(pstar/WL->p,(GAMMA-1.)/(2.*GAMMA)); // sound speed behind the fan
	  SHL=WL->u-WL->a;
	  STL=ustar-astar;

	  if(S<SHL)
	    {
	      // left to the shock
	      W->d=WL->d;
	      W->u=WL->u;
	      W->p=WL->p;
	      W->a=WL->a;
	    }
	  else
	    {
	      if(S>STL)
		{
		  W->d=WL->d*POW(pstar/WL->p,1./GAMMA);
		  W->u=ustar;
		  W->p=pstar;
		  W->a=SQRT(GAMMA*W->p/W->d);
		}
	      else
		{
		  W->d=WL->d*POW(2./(GAMMA+1.)+(GAMMA-1.)/((GAMMA+1.)*WL->a)*(WL->u-S),2./(GAMMA-1.));
		  W->u=2./(GAMMA+1.)*(WL->a+(GAMMA-1.)/2.*WL->u+S);
		  W->p=WL->p*POW(2./(GAMMA+1.)+(GAMMA-1.)/((GAMMA+1.)*WL->a)*(WL->u-S),2.*GAMMA/(GAMMA-1.));
		  W->a=SQRT(GAMMA*W->p/W->d);
		}
	    }
	}
      //if(W->p<0) abort();
    }
  else
    {
      if(pstar>WR->p)
	{
	  // Right Shock with shock speed SR
	  SR=WR->u+WR->a*SQRT((GAMMA+1.)/(2.*GAMMA)*pstar/WR->p+(GAMMA-1.)/(2.*GAMMA)); 
	  
	  if(S>SR)
	    {
	      // Right to the shock
	      W->d=WR->d;
	      W->u=WR->u;
	      W->p=WR->p;
	      W->a=WR->a;
	    }
	  else
	    {
	      W->d=WR->d*((pstar/WR->p+(GAMMA-1.)/(GAMMA+1.))/((GAMMA-1.)/(GAMMA+1.)*pstar/WR->p+1.));
	      W->u=ustar;
	      W->p=pstar;
	      W->a=SQRT(GAMMA*W->p/W->d);
	    }
	}
      else
	{
	  // Right fan
	  astar=WR->a*POW(pstar/WR->p,(GAMMA-1.)/(2.*GAMMA)); // sound speed behind the fan
	  SHR=WR->u+WR->a;
	  STR=ustar+astar;

	  if(S>SHR)
	    {
	      // Right to the fan
	      W->d=WR->d;
	      W->u=WR->u;
	      W->p=WR->p;
	      W->a=WR->a;
	    }
	  else
	    {
	      if(S<STR)
		{
		  W->d=WR->d*POW(pstar/WR->p,1./GAMMA);
		  W->u=ustar;
		  W->p=pstar;
		  W->a=SQRT(GAMMA*W->p/W->d);
		}
	      else
		{
		  W->d=WR->d*POW(2./(GAMMA+1.)-(GAMMA-1.)/((GAMMA+1.)*WR->a)*(WR->u-S),2./(GAMMA-1.));
		  W->u=2./(GAMMA+1.)*(-WR->a+(GAMMA-1.)/2.*WR->u+S);
		  W->p=WR->p*POW(2./(GAMMA+1.)-(GAMMA-1.)/((GAMMA+1.)*WR->a)*(WR->u-S),2.*GAMMA/(GAMMA-1.));
		  W->a=SQRT(GAMMA*W->p/W->d);
		}
	    }
	}

      

    }

}


// ======================================================== Flux fonctions

void getflux_X(struct Utype *U, REAL *f)
{
  f[0]=U->du;
  f[1]=0.5*(3.-GAMMA)*U->du*U->du/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->dv*U->dv+U->dw*U->dw)/U->d;
  f[2]=U->du*U->dv/U->d;
  f[3]=U->du*U->dw/U->d;
  f[4]=GAMMA*U->du/U->d*U->E-0.5*(GAMMA-1.)*U->du/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);

#ifdef WRADHYD
  f[6]=U->du*U->dX/U->d;
#endif
}

// ---------------------------------------------------------------

void getflux_Y(struct Utype *U, REAL *f)
{
  f[0]=U->dv;
  f[1]=U->dv*U->du/U->d;
  f[2]=0.5*(3.-GAMMA)*U->dv*U->dv/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dw*U->dw)/U->d;
  f[3]=U->dv*U->dw/U->d;
  f[4]=GAMMA*U->dv/U->d*U->E-0.5*(GAMMA-1.)*U->dv/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
#ifdef WRADHYD
  f[6]=U->dv*U->dX/U->d;
#endif
}

// ---------------------------------------------------------------

void getflux_Z(struct Utype *U, REAL *f)
{
  f[0]=U->dw;
  f[1]=U->dw*U->du/U->d;
  f[2]=U->dw*U->dv/U->d;
  f[3]=0.5*(3.-GAMMA)*U->dw*U->dw/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dv*U->dv)/U->d;
  f[4]=GAMMA*U->dw/U->d*U->E-0.5*(GAMMA-1.)*U->dw/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
#ifdef WRADHYD
  f[6]=U->dw*U->dX/U->d;
#endif
}




// ===============================================================================================


void  matrix_jacobian(struct Wtype *W0, REAL dt,REAL dx,struct Wtype *Dx,struct Wtype *Dy,struct Wtype *Dz, struct Wtype *Wt){


  REAL M[25];
  REAL W[6]={0.,0.,0.,0.,0.,0.};
  REAL d[5];
  int i,j;
#ifdef WRADHYD
  REAL X;
#endif

  // =====  building the A matrix

  memset(M,0,25*sizeof(REAL));
  
  // diagonal elements
  for(i=0;i<5;i++) M[i+i*5]=W0->u;
  
  // off_diagonal elements
  M[0+1*5]=W0->d;

  M[4+1*5]=W0->d*W0->a*W0->a;

  M[1+4*5]=1./W0->d;


  // ===== First Product

  d[0]=Dx->d;
  d[1]=Dx->u;
  d[2]=Dx->v;
  d[3]=Dx->w;
  d[4]=Dx->p;

  for(j=0;j<5;j++){
    for(i=0;i<5;i++){
      W[i]+=M[i+j*5]*d[j];
    }
  }

#ifdef WRADHYD
  W[5]+=W0->u*Dx->dX+W0->dX*Dx->u;
#endif

  // =====  building the B matrix

  memset(M,0,25*sizeof(REAL));
  
  // diagonal elements
  for(i=0;i<5;i++) M[i+i*5]=W0->v;
  
  // off_diagonal elements
  M[0+2*5]=W0->d;

  M[4+2*5]=W0->d*W0->a*W0->a;

  M[2+4*5]=1./W0->d;


  // ===== Second Product

  d[0]=Dy->d;
  d[1]=Dy->u;
  d[2]=Dy->v;
  d[3]=Dy->w;
  d[4]=Dy->p;

  for(j=0;j<5;j++){
    for(i=0;i<5;i++){
      W[i]+=M[i+j*5]*d[j];
      }
  }

#ifdef WRADHYD
  W[5]+=W0->v*Dx->dX+W0->dX*Dx->v;
#endif

  // =====  building the C matrix

  memset(M,0,25*sizeof(REAL));
  
  // diagonal elements
  for(i=0;i<5;i++) M[i+i*5]=W0->w;
  
  // off_diagonal elements
  M[0+3*5]=W0->d;

  M[4+3*5]=W0->d*W0->a*W0->a;

  M[3+4*5]=1./W0->d;

  d[0]=Dz->d;
  d[1]=Dz->u;
  d[2]=Dz->v;
  d[3]=Dz->w;
  d[4]=Dz->p;

  // ===== Third Product
 
  for(j=0;j<5;j++){
    for(i=0;i<5;i++){
      W[i]+=M[i+j*5]*d[j];
      }
  }

#ifdef WRADHYD
  W[5]+=W0->w*Dx->dX+W0->w*Dx->dX;
#endif
  
  // ==== Final correction
  for(i=0;i<6;i++){
    W[i]*=(-dt/dx*0.5);
  }
  
  Wt->d=W[0];
  Wt->u=W[1];
  Wt->v=W[2];
  Wt->w=W[3];
  Wt->p=W[4];
  
#ifdef WRADHYD
  Wt->dX=W[5];
#endif

}

// ======================================================================

void MUSCL_BOUND2(struct HGRID *stencil, int ioct, int icell, struct Wtype *Wi,REAL dt,REAL dx){ 

	  struct OCT * oct;
	  
	  struct Wtype *W0;
	  struct Wtype *Wp;
	  struct Wtype *Wm;
	  struct Wtype Dp,Dm;
	  struct Wtype D[3];
	  struct Wtype Wt;
	  int inei2;
	  int vcell[6],vnei[6];
	  int dir;
	  int idir;

#ifdef WGRAV
	  REAL f[3];
	  struct Utype S;
	  struct Utype U;
#endif

	  getcellnei(icell, vnei, vcell); // we get the neighbors
	  
	  W0=&(stencil->oct[ioct].cell[icell].field);
	
	  // Limited Slopes
	  for(dir=0;dir<3;dir++){
	    
	    inei2=2*dir;
	    if(vnei[inei2]==6){
	      Wm=&(stencil->oct[ioct].cell[vcell[inei2]].field);
	    }
	    else{
	      Wm=&(stencil->oct[ioct-(int)POW(3,dir)].cell[vcell[inei2]].field);
	    }

	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      Wp=&(stencil->oct[ioct].cell[vcell[inei2]].field);
	    }
	    else{
	      Wp=&(stencil->oct[ioct+(int)POW(3,dir)].cell[vcell[inei2]].field);
	    }

	    diffW(Wp,W0,&Dp); 
	    diffW(W0,Wm,&Dm); 
	    
	    minmod_W(&Dm,&Dp,D+dir);
	  }


	  // build jacobian matrix product
	  
	  matrix_jacobian(W0,dt,dx,&D[0],&D[1],&D[2],&Wt); // Here Wt contains the evolution of the state
	  
	  // READY TO EVOLVE EXTRAPOLATED VALUE

	  REAL ix[]={-0.5,0.5,0.0,0.0,0.0,0.0};
	  REAL iy[]={0.0,0.0,-0.5,0.5,0.0,0.0};
	  REAL iz[]={0.0,0.0,0.0,0.0,-0.5,0.5};

#ifdef WGRAV
#ifndef NOCOUPLE
	  memcpy(f,stencil->oct[ioct].cell[icell].f,sizeof(REAL)*3);

#ifdef CONSERVATIVE
	  S.d =0.;
	  S.du=-W0->d*f[0]*0.5*dt;
	  S.dv=-W0->d*f[1]*0.5*dt;
	  S.dw=-W0->d*f[2]*0.5*dt;
	  S.E =-(W0->d*W0->u*f[0]+W0->d*W0->v*f[1]+W0->d*W0->w*f[2])*dt*0.5;
#endif

#endif
#endif
	  for(idir=0;idir<6;idir++){
	    Wi[idir].d = W0->d+ix[idir]*D[0].d+iy[idir]*D[1].d+iz[idir]*D[2].d+Wt.d;
	    Wi[idir].u = W0->u+ix[idir]*D[0].u+iy[idir]*D[1].u+iz[idir]*D[2].u+Wt.u;
	    Wi[idir].v = W0->v+ix[idir]*D[0].v+iy[idir]*D[1].v+iz[idir]*D[2].v+Wt.v;
	    Wi[idir].w = W0->w+ix[idir]*D[0].w+iy[idir]*D[1].w+iz[idir]*D[2].w+Wt.w;
	    Wi[idir].p = FMAX(W0->p+ix[idir]*D[0].p+iy[idir]*D[1].p+iz[idir]*D[2].p+Wt.p,PMIN);
#ifdef WRADHYD
	    Wi[idir].dX = W0->dX+ix[idir]*D[0].dX+iy[idir]*D[1].dX+iz[idir]*D[2].dX+Wt.dX;
	    //printf("%e %e %e | %e %e\n",D[0].dX,D[1].dX,D[2].dX,W0->dX,Wt.dX);
#endif
	    if(Wi[idir].d<0) abort();
	    //if(Wi[idir].p==PMIN) printf("%e %e \n",W0->p,W0->p+ix[idir]*D[0].p+iy[idir]*D[1].p+iz[idir]*D[2].p+Wt.p);


#ifdef WGRAV
#ifndef NOCOUPLE

#ifdef PRIMITIVE
	    Wi[idir].u+=-f[0]*0.5*dt;
	    Wi[idir].v+=-f[1]*0.5*dt;
	    Wi[idir].w+=-f[2]*0.5*dt;
#endif

#ifdef CONSERVATIVE
 	    W2U(&Wi[idir],&U);
	    U.d  +=S.d;
	    U.du +=S.du;
	    U.dv +=S.dv;
	    U.dw +=S.dw;
	    U.E  +=S.E;
	    U2W(&U,&Wi[idir]);
#endif

#endif
#endif
	    
	    //if(Wi[idir].p<0) abort();
	    //Wi[idir].E=Wi[idir].p/(GAMMA-1.)+0.5*Wi[idir].d*(Wi[idir].u*Wi[idir].u+Wi[idir].v*Wi[idir].v+Wi[idir].w*Wi[idir].w);
	    getE(Wi+idir);
	    Wi[idir].a=SQRT(GAMMA*Wi[idir].p/Wi[idir].d);

#ifdef WRADHYD
 	    Wi[idir].dX=W0->dX; 
#endif 

	   
	  }



	  
}

//========================= HLL TOOLS =================================================================/


void speedestimateX_HLLC(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR, REAL *pstar, REAL *ustar){

  REAL qL,qR;
  struct Wtype1D WLloc;
  struct Wtype1D WRloc;
  int n;

  WLloc.d=WL->d;
  WLloc.u=WL->u;
  WLloc.p=WL->p;
  WLloc.a=SQRT(GAMMA*WLloc.p/WLloc.d);
  
  WRloc.d=WR->d;
  WRloc.u=WR->u;
  WRloc.p=WR->p;
  WRloc.a=SQRT(GAMMA*WRloc.p/WRloc.d);

  (*pstar)=(REAL)findPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=(REAL) findPressure(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) {
    printf("neg pressure ABORT\n");
    abort();
  }

  qL=(*pstar<=WL->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;
  if((*SL)>(*SR)){
    (*SL)=FMIN(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=FMAX(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
    //abort();
  }
  if((*SL)>(*SR)){
    printf("VELOCITY ORDER ABORT\n");
    abort();
  }
  if(isnan(*ustar)){
    printf("ustar isnan ABORT %e | %e %e %e | %e %e %e\n",*ustar,WLloc.d,WLloc.u,WLloc.p,WRloc.d,WRloc.u,WRloc.p);
    abort();
  }
}

void speedestimateY_HLLC(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR, REAL *pstar, REAL *ustar){

  REAL qL,qR;
  struct Wtype1D WLloc;
  struct Wtype1D WRloc;
  int n;

  WLloc.d=WL->d;
  WLloc.u=WL->v;
  WLloc.p=WL->p;
  WLloc.a=SQRT(GAMMA*WLloc.p/WLloc.d);
  
  WRloc.d=WR->d;
  WRloc.u=WR->v;
  WRloc.p=WR->p;
  WRloc.a=SQRT(GAMMA*WRloc.p/WRloc.d);

  (*pstar)=(REAL)findPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=(REAL) findPressure(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) abort();

  qL=(*pstar<=WL->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;
  if((*SL)>(*SR)){
    (*SL)=FMIN(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=FMAX(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
    //abort();
  }
  if((*SL)>(*SR)) abort();
  if(isnan(*ustar)){
    printf("ustar isnan ABORT %e | %e %e %e | %e %e %e\n",*ustar,WLloc.d,WLloc.u,WLloc.p,WRloc.d,WRloc.u,WRloc.p);
    abort();
  }

}



void speedestimateZ_HLLC(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR, REAL *pstar, REAL *ustar){

  REAL qL,qR;
  struct Wtype1D WLloc;
  struct Wtype1D WRloc;
  int n;

  WLloc.d=WL->d;
  WLloc.u=WL->w;
  WLloc.p=WL->p;
  WLloc.a=SQRT(GAMMA*WLloc.p/WLloc.d);
  
  WRloc.d=WR->d;
  WRloc.u=WR->w;
  WRloc.p=WR->p;
  WRloc.a=SQRT(GAMMA*WRloc.p/WRloc.d);

  (*pstar)=(REAL)findPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=(REAL) findPressure(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) abort();

  qL=(*pstar<=WL->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;
  if((*SL)>(*SR)){
    (*SL)=FMIN(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=FMAX(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
    //abort();
  }
  if((*SL)>(*SR)) abort();
  if(isnan(*ustar)) abort();

}




// =============================================================================================

int hydroM_sweepZ(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt){

  int inei,icell,iface;
  int i;
  int vnei[6],vcell[6];

  REAL FL[NVAR],FR[NVAR];
  struct Utype Uold;
  struct Wtype Wold;
  REAL pstar,ustar;

  struct Wtype WT[6]; // FOR MUSCL RECONSTRUCTION
  struct Wtype WC[6]; // FOR MUSCL RECONSTRUCTION

  struct Utype UC[2];
  struct Utype UN[2];
  struct Wtype WN[2];

  int ioct[7]={12,14,10,16,4,22,13};
  int idxnei[6]={1,0,3,2,5,4};

  struct Wtype *curcell;

  REAL SL,SR;
  
  int ffact[2]={0,0};
  REAL fact;

#ifdef DUAL_E
  struct Utype Us;
  REAL ebar;
  REAL ecen=0.;
  REAL divu,divuloc;
#endif

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
      
    for(i=0;i<nread;i++){ // we scan the octs
      
      memset(FL,0,sizeof(REAL)*NVAR);
      memset(FR,0,sizeof(REAL)*NVAR);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].field);

#ifdef DUAL_E
      divu=stencil[i].New.cell[icell].divu;
#endif      

      Wold.d=curcell->d;
      Wold.u=curcell->u;
      Wold.v=curcell->v;
      Wold.w=curcell->w;
      Wold.p=curcell->p;
      Wold.a=SQRT(GAMMA*Wold.p/Wold.d);

/* #ifdef WRADHYD */
/*       Wold.dX=curcell->dX; */
/* #endif */

      W2U(&Wold,&Uold); // primitive -> conservative

      REAL eold=Uold.eint;

      /* // MUSCL STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);

      MUSCL_BOUND2(stencil+i, 13, icell, WC,dt,dx);// central


      for(iface=0;iface<2;iface++){
	inei=iface+4;
	memcpy(WC+iface,WC+inei,sizeof(struct Wtype)); // moving the data towards idx=0,1
	W2U(WC+iface,UC+iface);
      }

      // Neighbor MUSCL reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+4;
	MUSCL_BOUND2(stencil+i, ioct[vnei[inei]], vcell[inei], WT,dt,dx);// 


	memcpy(WN+iface,WT+idxnei[inei],sizeof(struct Wtype)); 
	//memcpy(WN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].field),sizeof(struct Wtype)); 
	W2U(WN+iface,UN+iface);
	
	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }


      // Z DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems BOTTOM

      // Switching to Split description

      /* 	// =========================================== */

#ifdef RIEMANN_HLLC
      speedestimateZ_HLLC(&WN[0],&WC[0],&SL,&SR,&pstar,&ustar);

      if(SL>=0.){
	getflux_Z(&UN[0],FL);
#ifdef DUAL_E
	memcpy(&Us,&UN[0],sizeof(struct Utype));
#endif

      }
      else if(SR<=0.){
	getflux_Z(&UC[0],FL);
#ifdef DUAL_E
	memcpy(&Us,&UC[0],sizeof(struct Utype));
#endif
      }
      else if((SL<0.)&&(ustar>=0.)){
	getflux_Z(&UN[0],FL);
	fact=WN[0].d*(SL-WN[0].w)/(SL-ustar);
	FL[0]+=(fact*1.                                                                      -UN[0].d )*SL;
	FL[1]+=(fact*WN[0].u                                                                 -UN[0].du)*SL;
	FL[2]+=(fact*WN[0].v                                                                 -UN[0].dv)*SL;
	FL[3]+=(fact*ustar                                                                   -UN[0].dw)*SL;
	FL[4]+=(fact*(UN[0].E/UN[0].d+(ustar-WN[0].w)*(ustar+WN[0].p/(WN[0].d*(SL-WN[0].w))))-UN[0].E )*SL;

#ifdef DUAL_E
	Us.d =(fact*1.);
	Us.du=(fact*WN[0].u);
	Us.dv=(fact*WN[0].v);
	Us.dw=(fact*ustar);
	Us.E =(fact*(UN[0].E/UN[0].d+(ustar-WN[0].w)*(ustar+WN[0].p/(WN[0].d*(SL-WN[0].w)))));
#endif

#ifdef WRADHYD
	FL[6]+=(fact*WN[0].dX/WN[0].d                                                                 -UN[0].dX)*SL;
#endif
      }
      else if((ustar<=0.)&&(SR>0.)){
	getflux_Z(&UC[0],FL);
	fact=WC[0].d*(SR-WC[0].w)/(SR-ustar);
	FL[0]+=(fact*1.                                                                      -UC[0].d )*SR;
	FL[1]+=(fact*WC[0].u                                                                 -UC[0].du)*SR;
	FL[2]+=(fact*WC[0].v                                                                 -UC[0].dv)*SR;
	FL[3]+=(fact*ustar                                                                   -UC[0].dw)*SR;
	FL[4]+=(fact*(UC[0].E/UC[0].d+(ustar-WC[0].w)*(ustar+WC[0].p/(WC[0].d*(SR-WC[0].w))))-UC[0].E )*SR;

#ifdef DUAL_E
	Us.d =(fact*1.);
	Us.du=(fact*WC[0].u);
	Us.dv=(fact*WC[0].v);
	Us.dw=(fact*ustar);
	Us.E =(fact*(UC[0].E/UC[0].d+(ustar-WC[0].w)*(ustar+WC[0].p/(WC[0].d*(SR-WC[0].w)))));
#endif

#ifdef WRADHYD
	FL[6]+=(fact*WC[0].dX/WC[0].d                                                                 -UC[0].dX)*SR;
#endif
      }

#ifdef DUAL_E
      ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d); 
      divuloc=(GAMMA-1.)*(Us.dw/Us.d)*eold;
      FL[5]=(Us.dw/Us.d*ebar);
      divu+=-divuloc;
#endif

#endif
      // ===========================================



      // --------- solving the Riemann Problems TOP


      // Switching to Split description

      //=====================================================

#ifdef RIEMANN_HLLC
      speedestimateZ_HLLC(&WC[1],&WN[1],&SL,&SR,&pstar,&ustar);

      if(SL>=0.){
	getflux_Z(&UC[1],FR);
#ifdef DUAL_E
	memcpy(&Us,&UC[1],sizeof(struct Utype));
#endif

      }
      else if(SR<=0.){
	getflux_Z(&UN[1],FR);
#ifdef DUAL_E
	memcpy(&Us,&UN[1],sizeof(struct Utype));
#endif

      }
      else if((SL<0.)&&(ustar>=0.)){
	getflux_Z(&UC[1],FR);
	fact=WC[1].d*(SL-WC[1].w)/(SL-ustar);
	FR[0]+=(fact*1.                                                                      -UC[1].d )*SL;
	FR[1]+=(fact*WC[1].u                                                                 -UC[1].du)*SL;
	FR[2]+=(fact*WC[1].v                                                                 -UC[1].dv)*SL;
	FR[3]+=(fact*ustar                                                                   -UC[1].dw)*SL;
	FR[4]+=(fact*(UC[1].E/UC[1].d+(ustar-WC[1].w)*(ustar+WC[1].p/(WC[1].d*(SL-WC[1].w))))-UC[1].E )*SL;

#ifdef DUAL_E
	Us.d =(fact*1.);
	Us.du=(fact*WC[1].u);
	Us.dv=(fact*WC[1].v);
	Us.dw=(fact*ustar);
	Us.E =(fact*(UC[1].E/UC[1].d+(ustar-WC[1].w)*(ustar+WC[1].p/(WC[1].d*(SL-WC[1].w)))));
#endif

#ifdef WRADHYD
	FR[6]+=(fact*WC[1].dX/WC[1].d                                                                 -UC[1].dX)*SL;
#endif
      }
      else if((ustar<=0.)&&(SR>0.)){
	getflux_Z(&UN[1],FR);
	fact=WN[1].d*(SR-WN[1].w)/(SR-ustar);
	FR[0]+=(fact*1.                                                                      -UN[1].d )*SR;
	FR[1]+=(fact*WN[1].u                                                                 -UN[1].du)*SR;
	FR[2]+=(fact*WN[1].v                                                                 -UN[1].dv)*SR;
	FR[3]+=(fact*ustar                                                                   -UN[1].dw)*SR;
	FR[4]+=(fact*(UN[1].E/UN[1].d+(ustar-WN[1].w)*(ustar+WN[1].p/(WN[1].d*(SR-WN[1].w))))-UN[1].E )*SR;

#ifdef DUAL_E
	Us.d =(fact*1.);
	Us.du=(fact*WN[1].u);
	Us.dv=(fact*WN[1].v);
	Us.dw=(fact*ustar);
	Us.E =(fact*(UN[1].E/UN[1].d+(ustar-WN[1].w)*(ustar+WN[1].p/(WN[1].d*(SR-WN[1].w)))));
#endif

#ifdef WRADHYD
	FR[6]+=(fact*WN[1].dX/WN[1].d                                                                 -UN[1].dX)*SR;
#endif
      }

#ifdef DUAL_E
      ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d); 
      divuloc=(GAMMA-1.)*(Us.dw/Us.d)*eold;
      FR[5]=(Us.dw/Us.d*ebar);
      divu+= divuloc;
#endif


#endif


      //========================= copy the fluxes

      // Cancelling the fluxes from splitted neighbours

      for(iface=0;iface<NVAR;iface++) FL[iface]*=ffact[0]; 
      for(iface=0;iface<NVAR;iface++) FR[iface]*=ffact[1]; 


      memcpy(stencil[i].New.cell[icell].flux+4*NVAR,FL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].New.cell[icell].flux+5*NVAR,FR,sizeof(REAL)*NVAR);

      stencil[i].New.cell[icell].divu=divu;

      // ready for the next cell
    }
    //ready for the next oct
  }

  return 0;
}




//============================================================================
// =============================================================================================

int hydroM_sweepY(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt){

  int inei,icell,iface;
  int i;
  int vnei[6],vcell[6];

  REAL FL[NVAR],FR[NVAR];
  struct Utype Uold;
  struct Wtype Wold;
  REAL pstar,ustar;

  struct Wtype WT[6]; // FOR MUSCL RECONSTRUCTION
  struct Wtype WC[6]; // FOR MUSCL RECONSTRUCTION

  struct Utype UC[2];
  struct Utype UN[2];
  struct Wtype WN[2];

  int ioct[7]={12,14,10,16,4,22,13};
  int idxnei[6]={1,0,3,2,5,4};

  struct Wtype *curcell;

  REAL SL,SR;
  
  int ffact[2]={0,0};
  REAL fact;

#ifdef DUAL_E
	struct Utype Us;
	REAL ebar;
	REAL ecen=0.;
	REAL divu,divuloc;
#endif

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
      
    for(i=0;i<nread;i++){ // we scan the octs
      
      memset(FL,0,sizeof(REAL)*NVAR);
      memset(FR,0,sizeof(REAL)*NVAR);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].field);

#ifdef DUAL_E
      divu=stencil[i].New.cell[icell].divu;
#endif      

      Wold.d=curcell->d;
      Wold.u=curcell->u;
      Wold.v=curcell->v;
      Wold.w=curcell->w;
      Wold.p=curcell->p;
      Wold.a=SQRT(GAMMA*Wold.p/Wold.d);

/* #ifdef WRADHYD */
/*       Wold.X=curcell->X; */
/* #endif */

      W2U(&Wold,&Uold); // primitive -> conservative

      REAL eold=Uold.eint;

      /* // MUSCL STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);

      MUSCL_BOUND2(stencil+i, 13, icell, WC,dt,dx);// central


      for(iface=0;iface<2;iface++){
	inei=iface+2;
	memcpy(WC+iface,WC+inei,sizeof(struct Wtype)); // moving the data towards idx=0,1
	//memcpy(WC+iface,&Wold,sizeof(struct Wtype)); // moving the data towards idx=0,1

	W2U(WC+iface,UC+iface);
      }

      // Neighbor MUSCL reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+2;
	MUSCL_BOUND2(stencil+i, ioct[vnei[inei]], vcell[inei], WT,dt,dx);// 
	memcpy(WN+iface,WT+idxnei[inei],sizeof(struct Wtype)); 
	//memcpy(WN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].field),sizeof(struct Wtype)); 

       	W2U(WN+iface,UN+iface);
	
	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }




      // Y DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems FRONT

      // Switching to Split description

/* 	// =========================================== */

#ifdef RIEMANN_HLLC
	speedestimateY_HLLC(&WN[0],&WC[0],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  getflux_Y(&UN[0],FL);
#ifdef DUAL_E
	  memcpy(&Us,&UN[0],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  getflux_Y(&UC[0],FL);
#ifdef DUAL_E
	  memcpy(&Us,&UC[0],sizeof(struct Utype));
#endif
	}
	else if((SL<0.)&&(ustar>=0.)){
	  getflux_Y(&UN[0],FL);
	  fact=WN[0].d*(SL-WN[0].v)/(SL-ustar);
	  FL[0]+=(fact*1.                                                                      -UN[0].d )*SL;
	  FL[1]+=(fact*WN[0].u                                                                 -UN[0].du)*SL;
	  FL[2]+=(fact*ustar                                                                   -UN[0].dv)*SL;
	  FL[3]+=(fact*WN[0].w                                                                 -UN[0].dw)*SL;
	  FL[4]+=(fact*(UN[0].E/UN[0].d+(ustar-WN[0].v)*(ustar+WN[0].p/(WN[0].d*(SL-WN[0].v))))-UN[0].E )*SL;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WN[0].u);
	  Us.dv=(fact*ustar);
	  Us.dw=(fact*WN[0].w);
	  Us.E =(fact*(UN[0].E/UN[0].d+(ustar-WN[0].v)*(ustar+WN[0].p/(WN[0].d*(SL-WN[0].v)))));
#endif

#ifdef WRADHYD
	  FL[6]+=(fact*WN[0].dX/WN[0].d                                                                 -UN[0].dX)*SL;
#endif

	}
	else if((ustar<=0.)&&(SR>0.)){
	  getflux_Y(&UC[0],FL);
	  fact=WC[0].d*(SR-WC[0].v)/(SR-ustar);
	  FL[0]+=(fact*1.                                                                      -UC[0].d )*SR;
	  FL[1]+=(fact*WC[0].u                                                                 -UC[0].du)*SR;
	  FL[2]+=(fact*ustar                                                                   -UC[0].dv)*SR;
	  FL[3]+=(fact*WC[0].w                                                                 -UC[0].dw)*SR;
	  FL[4]+=(fact*(UC[0].E/UC[0].d+(ustar-WC[0].v)*(ustar+WC[0].p/(WC[0].d*(SR-WC[0].v))))-UC[0].E )*SR;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WC[0].u);
	  Us.dv=(fact*ustar);
	  Us.dw=(fact*WC[0].w);
	  Us.E =(fact*(UC[0].E/UC[0].d+(ustar-WC[0].v)*(ustar+WC[0].p/(WC[0].d*(SR-WC[0].v)))));
#endif

#ifdef WRADHYD
	  FL[6]+=(fact*WC[0].dX/WC[0].d                                                                 -UC[0].dX)*SR;
#endif

	}


#ifdef DUAL_E
      ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d); 
      FL[5]=(Us.dv/Us.d*ebar);
      divuloc=(GAMMA-1.)*(Us.dv/Us.d)*eold;
      divu+=-divuloc;
#endif

/* #ifdef DUAL_E */
/* 	ebar=(Us.E-0.5*(Us.du*Us.du+Us.dw*Us.dw+Us.dv*Us.dv)/Us.d); */
/* 	ecen+=ebar; */
/* 	FL[5]=(Us.dv/Us.d*ebar); */
/* 	divu[0]=(GAMMA-1.)*(Us.dv/Us.d); */
/* #endif */

#endif
	// ===========================================




      // --------- solving the Riemann Problems BACK


      // Switching to Split description

	//=====================================================

#ifdef RIEMANN_HLLC
	speedestimateY_HLLC(&WC[1],&WN[1],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  getflux_Y(&UC[1],FR);
#ifdef DUAL_E
	  memcpy(&Us,&UC[1],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  getflux_Y(&UN[1],FR);
#ifdef DUAL_E
	  memcpy(&Us,&UN[1],sizeof(struct Utype));
#endif

	}
	else if((SL<0.)&&(ustar>=0.)){
	  getflux_Y(&UC[1],FR);
	  fact=WC[1].d*(SL-WC[1].v)/(SL-ustar);
	  FR[0]+=(fact*1.                                                                      -UC[1].d )*SL;
	  FR[1]+=(fact*WC[1].u                                                                 -UC[1].du)*SL;
	  FR[2]+=(fact*ustar                                                                   -UC[1].dv)*SL;
	  FR[3]+=(fact*WC[1].w                                                                 -UC[1].dw)*SL;
	  FR[4]+=(fact*(UC[1].E/UC[1].d+(ustar-WC[1].v)*(ustar+WC[1].p/(WC[1].d*(SL-WC[1].v))))-UC[1].E )*SL;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WC[1].u);
	  Us.dv=(fact*ustar);
	  Us.dw=(fact*WC[1].w);
	  Us.E =(fact*(UC[1].E/UC[1].d+(ustar-WC[1].v)*(ustar+WC[1].p/(WC[1].d*(SL-WC[1].v)))));
#endif

#ifdef WRADHYD
	  FR[6]+=(fact*WC[1].dX/WC[1].d                                                                 -UC[1].dX)*SL;
#endif
	}
	else if((ustar<=0.)&&(SR>0.)){
	  getflux_Y(&UN[1],FR);
	  fact=WN[1].d*(SR-WN[1].v)/(SR-ustar);
	  FR[0]+=(fact*1.                                                                      -UN[1].d )*SR;
	  FR[1]+=(fact*WN[1].u                                                                 -UN[1].du)*SR;
	  FR[2]+=(fact*ustar                                                                   -UN[1].dv)*SR;
	  FR[3]+=(fact*WN[1].w                                                                 -UN[1].dw)*SR;
	  FR[4]+=(fact*(UN[1].E/UN[1].d+(ustar-WN[1].v)*(ustar+WN[1].p/(WN[1].d*(SR-WN[1].v))))-UN[1].E )*SR;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WN[1].u);
	  Us.dv=(fact*ustar);
	  Us.dw=(fact*WN[1].w);
	  Us.E =(fact*(UN[1].E/UN[1].d+(ustar-WN[1].v)*(ustar+WN[1].p/(WN[1].d*(SR-WN[1].v)))));
#endif

#ifdef WRADHYD
	  FR[6]+=(fact*WN[1].dX/WN[1].d                                                                 -UN[1].dX)*SR;
#endif
	}


#ifdef DUAL_E
      ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d); 
      divuloc=(GAMMA-1.)*(Us.dv/Us.d)*eold;
      FR[5]=(Us.dv/Us.d*ebar);
      divu+= divuloc;
#endif

#endif

      
      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours

      for(iface=0;iface<NVAR;iface++) FL[iface]*=ffact[0]; 
      for(iface=0;iface<NVAR;iface++) FR[iface]*=ffact[1]; 
      
      //      printf("u=%e x=%e\n",FR[1],FR[6]);
      
      memcpy(stencil[i].New.cell[icell].flux+2*NVAR,FL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].New.cell[icell].flux+3*NVAR,FR,sizeof(REAL)*NVAR);

      stencil[i].New.cell[icell].divu=divu;

      // ready for the next cell
    }
    //ready for the next oct
  }

  return 0;
}

//===================================================================================================
//===================================================================================================

int hydroM_sweepX(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt){

  int inei,icell,iface;
  int i;
  int vnei[6],vcell[6];

  REAL FL[NVAR],FR[NVAR];
  struct Utype Uold;
  struct Wtype Wold;
  REAL pstar,ustar;

  struct Wtype WT[6]; // FOR MUSCL RECONSTRUCTION
  struct Wtype WC[6]; // FOR MUSCL RECONSTRUCTION

  struct Utype UC[2];
  struct Utype UN[2];
  struct Wtype WN[2];

  int ioct[7]={12,14,10,16,4,22,13};
  int idxnei[6]={1,0,3,2,5,4};

  struct Wtype *curcell;

  REAL SL,SR;
  
  int ffact[2]={0,0};
  REAL fact;

#ifdef DUAL_E
  struct Utype Us;
  REAL ebar;
  REAL ecen=0.;
  REAL divu,divuloc;
#endif

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
      
    for(i=0;i<nread;i++){ // we scan the octs
      
      memset(FL,0,sizeof(REAL)*NVAR);
      memset(FR,0,sizeof(REAL)*NVAR);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].field);
      
#ifdef DUAL_E
      divu=stencil[i].New.cell[icell].divu;
#endif      

      Wold.d=curcell->d;
      Wold.u=curcell->u;
      Wold.v=curcell->v;
      Wold.w=curcell->w;
      Wold.p=curcell->p;
      Wold.a=SQRT(GAMMA*Wold.p/Wold.d);
      getE(&Wold);

/* #ifdef WRADHYD */
/*       Wold.X=curcell->X; */
/* #endif */

      W2U(&Wold,&Uold); // primitive -> conservative
      REAL eold=Uold.eint;

      /* // MUSCL STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);

      MUSCL_BOUND2(stencil+i, 13, icell, WC,dt,dx);// central



      for(iface=0;iface<2;iface++){
	W2U(WC+iface,UC+iface);
      }

      // Neighbor MUSCL reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface;
	MUSCL_BOUND2(stencil+i, ioct[vnei[inei]], vcell[inei], WT,dt,dx);// 
	memcpy(WN+iface,WT+idxnei[inei],sizeof(struct Wtype)); 
       	W2U(WN+iface,UN+iface);
	
	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}
	
	if(isnan(WN[iface].d)){
	  printf("%e %e %e\n",stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].field.d,stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].field.u,stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].field.p);
	  abort();
	}

      }

    


      // X DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems LEFT

      // Switching to Split description

/* 	// =========================================== */

#ifdef RIEMANN_HLLC
	speedestimateX_HLLC(&WN[0],&WC[0],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  getflux_X(&UN[0],FL);
#ifdef DUAL_E
	  memcpy(&Us,&UN[0],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  getflux_X(&UC[0],FL);
#ifdef DUAL_E
	  memcpy(&Us,&UC[0],sizeof(struct Utype));
#endif
	}
	else if((SL<0.)&&(ustar>=0.)){
	  REAL fbkp1;
	  getflux_X(&UN[0],FL);
	  fbkp1=FL[1];
	  fact=WN[0].d*(SL-WN[0].u)/(SL-ustar);
	  FL[0]+=(fact*1.                                                                      -UN[0].d )*SL;
	  FL[1]+=(fact*ustar                                                                   -UN[0].du)*SL;
	  FL[2]+=(fact*WN[0].v                                                                 -UN[0].dv)*SL;
	  FL[3]+=(fact*WN[0].w                                                                 -UN[0].dw)*SL;
	  FL[4]+=(fact*(UN[0].E/UN[0].d+(ustar-WN[0].u)*(ustar+WN[0].p/(WN[0].d*(SL-WN[0].u))))-UN[0].E )*SL;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*ustar);
	  Us.dv=(fact*WN[0].v);
	  Us.dw=(fact*WN[0].w);
	  Us.E =(fact*(UN[0].E/UN[0].d+(ustar-WN[0].u)*(ustar+WN[0].p/(WN[0].d*(SL-WN[0].u)))));
#endif

#ifdef WRADHYD
	 FL[6]+=(fact*WN[0].dX/WN[0].d                                                                 -UN[0].dX)*SL;
#endif
	 
	}
	else if((ustar<=0.)&&(SR>0.)){
	  getflux_X(&UC[0],FL);
	  fact=WC[0].d*(SR-WC[0].u)/(SR-ustar);
	  FL[0]+=(fact*1.                                                                      -UC[0].d )*SR;
	  FL[1]+=(fact*ustar                                                                   -UC[0].du)*SR;
	  FL[2]+=(fact*WC[0].v                                                                 -UC[0].dv)*SR;
	  FL[3]+=(fact*WC[0].w                                                                 -UC[0].dw)*SR;
	  FL[4]+=(fact*(UC[0].E/UC[0].d+(ustar-WC[0].u)*(ustar+WC[0].p/(WC[0].d*(SR-WC[0].u))))-UC[0].E )*SR;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*ustar);
	  Us.dv=(fact*WC[0].v);
	  Us.dw=(fact*WC[0].w);
	  Us.E =(fact*(UC[0].E/UC[0].d+(ustar-WC[0].u)*(ustar+WC[0].p/(WC[0].d*(SR-WC[0].u)))));
#endif

#ifdef WRADHYD
	  FL[6]+=(fact*WC[0].dX/WC[0].d                                                                 -UC[0].dX)*SR;
#endif

	}


#ifdef DUAL_E
      ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d); 
      divuloc=(GAMMA-1.)*(Us.du/Us.d)*eold;
      FL[5]=(Us.du/Us.d*ebar);
      divu+=-divuloc;
#endif


#endif

     

	// ===========================================


      

      // --------- solving the Riemann Problems RIGHT


      // Switching to Split description

	//=====================================================

#ifdef RIEMANN_HLLC
	speedestimateX_HLLC(&WC[1],&WN[1],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  getflux_X(&UC[1],FR);
#ifdef DUAL_E
	  memcpy(&Us,&UC[1],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  getflux_X(&UN[1],FR);
#ifdef DUAL_E
	  memcpy(&Us,&UN[1],sizeof(struct Utype));
#endif

	}
	else if((SL<0.)&&(ustar>=0.)){
	  
	  getflux_X(&UC[1],FR);
	  fact=WC[1].d*(SL-WC[1].u)/(SL-ustar);


	  FR[0]+=(fact*1.                                                                      -UC[1].d )*SL;
	  FR[1]+=(fact*ustar                                                                   -UC[1].du)*SL;
	  FR[2]+=(fact*WC[1].v                                                                 -UC[1].dv)*SL;
	  FR[3]+=(fact*WC[1].w                                                                 -UC[1].dw)*SL;
	  FR[4]+=(fact*(UC[1].E/UC[1].d+(ustar-WC[1].u)*(ustar+WC[1].p/(WC[1].d*(SL-WC[1].u))))-UC[1].E )*SL;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*ustar);
	  Us.dv=(fact*WC[1].v);
	  Us.dw=(fact*WC[1].w);
	  Us.E =(fact*(UC[1].E/UC[1].d+(ustar-WC[1].u)*(ustar+WC[1].p/(WC[1].d*(SL-WC[1].u)))));
#endif

#ifdef WRADHYD
	  FR[6]+=(fact*WC[1].dX/WC[1].d                                                                 -UC[1].dX)*SL;
#endif

	}
	else if((ustar<=0.)&&(SR>0.)){
	  getflux_X(&UN[1],FR);
	  fact=WN[1].d*(SR-WN[1].u)/(SR-ustar);
	  FR[0]+=(fact*1.                                                                      -UN[1].d )*SR;
	  FR[1]+=(fact*ustar                                                                   -UN[1].du)*SR;
	  FR[2]+=(fact*WN[1].v                                                                 -UN[1].dv)*SR;
	  FR[3]+=(fact*WN[1].w                                                                 -UN[1].dw)*SR;
	  FR[4]+=(fact*(UN[1].E/UN[1].d+(ustar-WN[1].u)*(ustar+WN[1].p/(WN[1].d*(SR-WN[1].u))))-UN[1].E )*SR;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*ustar);
	  Us.dv=(fact*WN[1].v);
	  Us.dw=(fact*WN[1].w);
	  Us.E =(fact*(UN[1].E/UN[1].d+(ustar-WN[1].u)*(ustar+WN[1].p/(WN[1].d*(SR-WN[1].u)))));
#endif

#ifdef WRADHYD
	  FR[6]+=(fact*WN[1].dX/WN[1].d                                                                 -UN[1].dX)*SR;
#endif
	}

#ifdef DUAL_E
      ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d); 
      divuloc=(GAMMA-1.)*(Us.du/Us.d)*eold;
      FR[5]=(Us.du/Us.d*ebar);
      divu+= divuloc;
#endif

#endif
      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours

      for(iface=0;iface<NVAR;iface++) FL[iface]*=ffact[0]; 
      for(iface=0;iface<NVAR;iface++) FR[iface]*=ffact[1]; 


      memcpy(stencil[i].New.cell[icell].flux+0*NVAR,FL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].New.cell[icell].flux+1*NVAR,FR,sizeof(REAL)*NVAR);

      stencil[i].New.cell[icell].divu=divu;

      // ready for the next cell
    }
    //ready for the next oct
  }

  return 0;
}



//============================================================================


REAL comptstep_hydro(int levelcoarse,int levelmax,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax){
  
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icell;
  REAL aa;
  REAL va,vx,vy,vz;
  REAL dt;
  REAL Smax=0.,S1;

  // Computing new timestep
  dt=tmax;
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      if(nextoct==NULL) continue; // in case the level is empty
      dxcur=POW(0.5,level); // +1 to protect level change
      do // sweeping through the octs of level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;

	  if(curoct->cpu!=cpu->rank) continue;

	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	       vx=curoct->cell[icell].field.u; 
	       vy=curoct->cell[icell].field.v; 
	       vz=curoct->cell[icell].field.w; 
	       va=SQRT(vx*vx+vy*vy+vz*vz); 
	       aa=SQRT(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d); 
	       Smax=FMAX(Smax,va+aa); 

	    }
	}while(nextoct!=NULL);

      //printf("thydro= %e Smax=%e dxcur=%e\n",dxcur*CFL/(Smax*3.),Smax,dxcur);
      //      abort();
      dt=FMIN(dxcur*CFL/(Smax*3.),dt);
    }

  // new tstep
  //printf("original dt=%f chosen dt=%f\n",tmax,dt);


/* #ifdef WMPI */
/*   // reducing by taking the smallest time step */
/*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_REAL,MPI_MIN,cpu->comm); */
/* #endif   */

  return dt;
}


//=================================================================================
//=================================================================================
//=================================================================================


//=================================================================================

#ifdef WGRAV
REAL comptstep_force(int levelcoarse,int levelmax,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax){
  
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icell;
  REAL dtloc;
  REAL dt;
  struct Utype U;
  struct Wtype W;
  REAL f[3];
  REAL DE,V,DV;
  //Smax=FMAX(Smax,SQRT(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);

  // Computing new timestep
  dt=tmax;
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      if(nextoct==NULL) continue; // in case the level is empty
      dxcur=POW(0.5,level); // +1 to protect level change
      do // sweeping through the octs of level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      memcpy(&W,&curoct->cell[icell].field,sizeof(struct Wtype));
	      W2U(&W,&U);
	      memcpy(f,curoct->cell[icell].f,sizeof(REAL)*3);
	      V=SQRT(POW(U.du,2.)+POW(U.dv,2.)+POW(U.dw,2.));
	      DV=SQRT(POW(U.d*f[0],2.)+POW(U.d*f[1],2.)+POW(U.d*f[2],2.));
	      DE=SQRT(POW(U.du*f[0],2.)+POW(U.dv*f[1],2.)+POW(U.dw*f[2],2.));
									       
	      if((DE>0.)&&(U.E>0.)){
		dtloc=0.5*U.E/DE;
		if(dt!=tmax){
		  dt=FMIN(dt,dtloc);
		}else{
		  dt=dtloc;
		}
	      }

	      if((DV>0.)&&(V>0.)){
		dtloc=0.5*V/DV;
		if(dt!=tmax){
		  dt=FMIN(dt,dtloc);
		}else{
		  dt=dtloc;
		}
	      }

	    }
	}while(nextoct!=NULL);
    }

  // new tstep
  //printf("original dt=%f chosen dt=%f\n",tmax,dt); 
  //  abort(); 

/* #ifdef WMPI */
/*   // reducing by taking the smallest time step */
/*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_REAL,MPI_MIN,cpu->comm); */
/* #endif   */

  return dt;
}
#endif

// ==============================================================================================
// ==============================================================================================
#ifdef WGRAV
void grav_correction(int level,struct RUNPARAMS *param, struct OCT ** firstoct, struct CPUINFO *cpu, REAL dt)
{
  
  struct Utype U0,U;
  struct Wtype W;
  struct Wtype Wnew;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  struct CELL *curcell;

  curoct=firstoct[level-1];
  if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
    nextoct=curoct;
    do{
      curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
      for(icell=0;icell<8;icell++){
	int ref=0;

	if(curoct->cell[icell].child==NULL){ // Leaf cell
	  curcell=&(curoct->cell[icell]);

	  memcpy(&Wnew,&(curcell->field),sizeof(struct Wtype));
	  

#ifdef PRIMITIVE
	  Wnew.u+=(-curcell->f[0]*dt);
	  Wnew.v+=(-curcell->f[1]*dt);
	  Wnew.w+=(-curcell->f[2]*dt);
#endif

#ifdef CONSERVATIVE
	  W2U(&Wnew,&U);
	  memcpy(&U0,&U,sizeof(struct Utype));
	  // grav force correction

	  U.du+=-(U0.d*curoct->cell[icell].f[0]*dt);
	  U.dv+=-(U0.d*curoct->cell[icell].f[1]*dt);
	  U.dw+=-(U0.d*curoct->cell[icell].f[2]*dt);
	  U.E+=-(U.du*curoct->cell[icell].f[0]+U.dv*curoct->cell[icell].f[1]+U.dw*curoct->cell[icell].f[2])*dt;

	  U2W(&U,&Wnew);
#endif	
	  getE(&Wnew);
	  if(Wnew.p<0) abort();
	  memcpy(&(curcell->field),&Wnew,sizeof(struct Wtype));
	}
      }
    }while(nextoct!=NULL);
  }
}
#endif
// ==============================================================================================
// ==============================================================================================

//=====================================================================================================================
//
//    REFLECHISSONS UN PEU EN FAIT LA RECHERCHE EST RECURSIVE
//

// Structure de base

#ifdef WRADTEST
void recursive_neighbor_gather_oct(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct HGRID *stencil,char *visit){

  int ix[6]={-1,1,0,0,0,0};
  int iy[6]={0,0,-1,1,0,0};
  int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;
  int face[8]={0,1,2,3,4,5,6,7};
  int tflag[6]={0,0,0,0,0,0};
  REAL dxcur;

  struct Wtype Wi[8];
  char child[8];

  struct OCT *oct;
  struct OCT *neioct;
  struct CELL *neicell;

  if(order==1){
    ineiloc=inei;
  }
  else if(order==2){
    ineiloc=inei2;
  }
  else if(order==3){
    ineiloc=inei3;

  }


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
      /* else{ */
      /* 	printf("big problem\n"); */
      /* 	abort(); */
      /* } */

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
    
  }


  //  oct=cell2oct(cell); 



  // ============================ END TRANSMISSIVE BOUNDARIES ====================

#ifdef WGRAV
    REAL floc[8*3];//local force
#endif

  if(neicell->child!=NULL){
    // optimal case
    for(icell=0;icell<8;icell++){
      memcpy(Wi+icell,&(neicell->child->cell[icell].field),sizeof(struct Wtype));
      
      child[icell]=(neicell->child->cell[icell].child!=NULL);
    }

#ifdef WGRAV
    for(icell=0;icell<8;icell++) memcpy(floc+3*icell,neicell->child->cell[icell].f,sizeof(REAL)*3);
#endif

  }
  else{
    coarse2fine_hydro2(neicell,Wi);
    //for(icell=0;icell<8;icell++) memcpy(Wi+icell,&neicell->field,sizeof(struct Wtype));
    
/* #ifdef WRADHYD */
/*     for(icell=0;icell<8;icell++) Wi[icell].X=neicell->field.X; */
/* #endif */

    for(icell=0;icell<8;icell++){
      child[icell]=0;
    }

#ifdef WGRAV
    coarse2fine_forcelin(neicell,floc);
    //for(icell=0;icell<8;icell++) memcpy(floc+3*icell,neicell->f,sizeof(REAL)*3);
#endif

  }




  for(icell=0;icell<8;icell++){

    memcpy(&(stencil->oct[ioct].cell[icell].field),Wi+face[icell],sizeof(struct Wtype)); //
    // ---------------------------

#ifdef TRANSXM
#ifdef REFXM
    if(tflag[0]){
      stencil->oct[ioct].cell[icell].field.u*=-1.;
    }
#endif
#endif


#ifdef TRANSYM
#ifdef REFYM
    if(tflag[2]){
      stencil->oct[ioct].cell[icell].field.v*=-1.;
    }
#endif
#endif

#ifdef TRANSZM
#ifdef REFZM
    if(tflag[4]){
      stencil->oct[ioct].cell[icell].field.w*=-1.;
    }
#endif
#endif


    // ---------------------------

    stencil->oct[ioct].cell[icell].split=child[face[icell]];

#ifdef WGRAV
    memcpy(stencil->oct[ioct].cell[icell].f,floc+3*face[icell],sizeof(REAL)*3); //
#endif
  }

   
  
  // next order
  if(order==1){
    for(i=0;i<6;i++){
      if((i>>1)==(inei>>1)) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      if(visit[ioct2]) continue;
      visit[ioct2]=1;
      recursive_neighbor_gather_oct(ioct2, inei, i, -1, 2, neicell, stencil,visit);
    }
  }
  else if(order==2) {
    for(i=0;i<6;i++){
      if(((i>>1)==(inei>>1))||((i>>1)==(inei2>>1))) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      if(visit[ioct2]) continue;
      visit[ioct2]=1;
      recursive_neighbor_gather_oct(ioct2, inei, inei2, i, 3, neicell, stencil,visit);
    }
  }
}

#else

void recursive_neighbor_gather_oct(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct HGRID *stencil,char *visit){

  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;


  struct OCT *oct;
  struct OCT *neioct;
  struct CELL *neicell;

  if(order==1){
    ineiloc=inei;
  }
  else if(order==2){
    ineiloc=inei2;
  }
  else if(order==3){
    ineiloc=inei3;

  }


  if(cell->child!=NULL){
    // the oct at the right level exists
    neicell=cell->child->nei[ineiloc];
    oct=cell->child;
  }
  else{
    getcellnei(cell->idx, vnei, vcell); // we get the neighbors
    oct=cell2oct(cell);

    if(vnei[ineiloc]==6){
      neicell=&(oct->cell[vcell[ineiloc]]);
    }
    else{
      if(oct->nei[ineiloc]->child!=NULL){
	neicell=&(oct->nei[ineiloc]->child->cell[vcell[ineiloc]]);
      }
    }
  }



  if(neicell->child!=NULL){
/*     // optimal case */
/*     for(icell=0;icell<8;icell++){ */
/*       memcpy(Wi+icell,&(neicell->child->cell[icell].field),sizeof(struct Wtype)); */
/*       child[icell]=(neicell->child->cell[icell].child!=NULL); */
/*     } */

/* #ifdef WGRAV */
/*     for(icell=0;icell<8;icell++) memcpy(floc+3*icell,neicell->child->cell[icell].f,sizeof(REAL)*3); */
/* #endif */

    for(icell=0.;icell<8;icell++){
      struct CELLLIGHT_H *c=&stencil->oct[ioct].cell[icell];
      struct CELL *co=&neicell->child->cell[icell];
      memcpy(&(c->field),&(co->field),sizeof(struct Wtype)); //
      c->split=(co->child!=NULL);
#ifdef WGRAV
      memcpy(c->f,co->f,sizeof(REAL)*3); //
#endif
    }
  }
  else{

    struct Wtype Wi[8];
    char child[8];
#ifdef WGRAV
    REAL floc[8*3];//local force
#endif

    coarse2fine_hydro2(neicell,Wi);
    //for(icell=0;icell<8;icell++) memcpy(Wi+icell,&neicell->field,sizeof(struct Wtype));
    
/* #ifdef WRADHYD */
/*     for(icell=0;icell<8;icell++) Wi[icell].X=neicell->field.X; */
/* #endif */

    for(icell=0;icell<8;icell++){
      child[icell]=0;
    }

#ifdef WGRAV
    coarse2fine_forcelin(neicell,floc);
    //for(icell=0;icell<8;icell++) memcpy(floc+3*icell,neicell->f,sizeof(REAL)*3);
#endif

  

  for(icell=0;icell<8;icell++){
    struct CELLLIGHT_H *c=&stencil->oct[ioct].cell[icell];
    memcpy(&(c->field),Wi+icell,sizeof(struct Wtype)); //
    c->split=child[icell];
#ifdef WGRAV
    memcpy(c->f,floc+3*icell,sizeof(REAL)*3); //
#endif
  }

  } 
  
  // next order
  if(order==1){
    for(i=0;i<6;i++){
      if((i>>1)==(inei>>1)) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      if(visit[ioct2]) continue;
      visit[ioct2]=1;
      recursive_neighbor_gather_oct(ioct2, inei, i, -1, 2, neicell, stencil,visit);
    }
  }
  else if(order==2) {
    for(i=0;i<6;i++){
      if(((i>>1)==(inei>>1))||((i>>1)==(inei2>>1))) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      if(visit[ioct2]) continue;
      visit[ioct2]=1;
      recursive_neighbor_gather_oct(ioct2, inei, inei2, i, 3, neicell, stencil,visit);
    }
  }
}

#endif
//=====================================================================================================================


struct OCT *gatherstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  char visit[27];
  int ix[6]={-1,1,0,0,0,0};
  int iy[6]={0,0,-1,1,0,0};
  int iz[6]={0,0,0,0,-1,1};
  int ioct;
  
  memset(visit,0,27*sizeof(char));

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;

      memset(visit,0,27*sizeof(char));
      // filling the values in the central oct

      for(icell=0;icell<8;icell++){
	memcpy(&(stencil[iread].oct[13].cell[icell].field),&(curoct->cell[icell].field),sizeof(struct Wtype)); // for calculations
	//memcpy(&(stencil[iread].New.cell[icell].field),&(curoct->cell[icell].fieldnew),sizeof(struct Wtype)); // for updates

	
#ifdef DUAL_E
	stencil[iread].New.cell[icell].divu=0.;
#endif
	
#ifdef WGRAV 
 	memcpy(stencil[iread].oct[13].cell[icell].f,curoct->cell[icell].f,sizeof(REAL)*3); // 
#endif 
	stencil[iread].oct[13].cell[icell].split=(curoct->cell[icell].child!=NULL);
      }
      visit[13]=1;

      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  visit[ioct]=1;
	  recursive_neighbor_gather_oct(ioct, inei, -1, -1, 1, cell, stencil+iread,visit);
	}
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}



// ================================================================
// ================================================================
// ================================================================

void updatefield(struct OCT *octstart, struct HGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew)
{
  int i,icell;
  struct Wtype W,W0,Wi;
  struct Utype U,U0,Ui;
  REAL one;
  int flx;
  REAL dtsurdx=dtnew/dxcur;
  REAL F[NFLUX];
#ifdef DUAL_E
  REAL DE,p0,p;
#endif
  
  for(i=0;i<nread;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells
      
      if(stencil[i].oct[13].cell[icell].split) continue;
      memcpy(F,stencil[i].New.cell[icell].flux,sizeof(REAL)*NFLUX);// New fluxes from the stencil

      
      // ==== updating

      // actually we compute and store the delta U only
      one=1.;
      memset(&U,0,sizeof(struct Utype)); // setting delta U
      for(flx=0;flx<6;flx++){
	U.d +=F[0+flx*NVAR]*dtsurdx*one;
	U.du+=F[1+flx*NVAR]*dtsurdx*one;
	U.dv+=F[2+flx*NVAR]*dtsurdx*one;
	U.dw+=F[3+flx*NVAR]*dtsurdx*one;
	U.E +=F[4+flx*NVAR]*dtsurdx*one;
#ifdef DUAL_E
	U.eint+=F[5+flx*NVAR]*dtsurdx*one;
#endif

#ifdef WRADHYD
#ifndef NOADX
	U.dX+=F[6+flx*NVAR]*dtsurdx*one;
#else
	U.dX+=0.;
#endif
#endif
	one*=-1.;
      }
      // scatter back the delta Uwithin the stencil
      memcpy(&(stencil[i].New.cell[icell].deltaU),&U,sizeof(struct Utype));
    }
  }

}

// ===========================================================================================================

// ================================================================



struct OCT *scatterstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int ioct[7]={12,14,10,16,4,22,13};
  int vnei[6],vcell[6],inei;
  
  nextoct=octstart;
  iread=0;

  struct Wtype W,deltaW;
  struct Utype U,deltaU;
  REAL dtsurdx=dtnew/dxcur;


#ifdef DUAL_E
  REAL DE,p0,p;
  REAL eo,e;
#endif
  int one;
  REAL F[NVAR];

  //printf("let's scatter\n");
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	//we scatter the values in the central cell
	
	//memcpy(&(curoct->cell[icell].fieldnew),&(stencil[iread].New.cell[icell].field),sizeof(struct Wtype));
	memcpy(&deltaU,&(stencil[iread].New.cell[icell].deltaU),sizeof(struct Utype)); // getting the delta U back
	W2U(&(curoct->cell[icell].fieldnew),&U);
	
	/* if(curoct->cell[icell].rfield.src>0){ */
	/*   printf("dd=%e d=%e p=%e\n",deltaU.d,curoct->cell[icell].fieldnew.d,curoct->cell[icell].field.p);  */
	/*   abort(); */
	/* } */
	/* printf("dd=%e d=%e p=%e\n",deltaU.d,curoct->cell[icell].fieldnew.d,curoct->cell[icell].field.p); */
	/* abort(); */
	//if(U.eint+deltaU.eint<0) abort();

	U.d  +=deltaU.d;
  	U.du +=deltaU.du;
	U.dv +=deltaU.dv;
	U.dw +=deltaU.dw;
	U.E  +=deltaU.E;

#ifdef WRADHYD
#ifdef NOADX
	REAL xion=curoct->cell[icell].rfield.nhplus/curoct->cell[icell].rfield.nh;
	U.dX=U.d*xion;
#else
	U.dX  +=deltaU.dX;
#endif
	//printf("d=%e deltaX=%e dX=%e X=%e\n",U.d,deltaU.dX,U.dX,U.dX/U.d);
#endif

#ifdef DUAL_E
	U.eint+=deltaU.eint;
	U.eint+=(-stencil[iread].New.cell[icell].divu*dtsurdx);

	// HERE WE FACE A CHOICE : SHOULD WE KEEP THIS INTERNAL ENERGY OR THE CONSERVATIVE ONE (E-Ekin) ?

	REAL eint=U.eint;
	REAL ekin=0.5*(U.du*U.du+U.dv*U.dv+U.dw*U.dw)/U.d;
	
	if((U.E-ekin)>0.5*POW(stencil[iread].New.cell[icell].divu,2)){ 
	  U.eint=U.E-ekin; 
	} 

	U.eint=FMAX(U.eint*(GAMMA-1.),PMIN)/(GAMMA-1.);

	if(U.eint<0) abort();

#endif
	
	/* if(curoct->level==7){ */
	/*   if((icell%2==1 && curoct->cpu==0)||(icell%2==0 && curoct->cpu==1)){ */
	/*     if(curoct->x==(0.5-1./64)){ */
	/*       printf("rank=%d icell=%d d=%e\n",curoct->cpu,icell,deltaU.d); */
	/*     } */
	/*   } */
	/* } */

	

	U2W(&U,&(curoct->cell[icell].fieldnew)); // at this stage the central cell has been updated
	getE(&(curoct->cell[icell].fieldnew));

#ifdef WRADHYD
	//Ceiling for the ionisation fraction
	if(curoct->cell[icell].fieldnew.dX<0) curoct->cell[icell].fieldnew.dX=0.;
	if(curoct->cell[icell].fieldnew.dX>curoct->cell[icell].fieldnew.d) curoct->cell[icell].fieldnew.dX=curoct->cell[icell].fieldnew.d;
#endif	
	// ==================== let us now deal with coarser neighbors
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
	      struct Wtype Wbkp;

	      // the COARSE neighbor cell is unsplit we update its value with fluxes
	   
	      // initial data from the new value
	      memcpy(&W,&(curoct->nei[vnei[inei]]->fieldnew),sizeof(struct Wtype));
	      memcpy(&Wbkp,&(curoct->nei[vnei[inei]]->fieldnew),sizeof(struct Wtype));
	      W2U(&W,&U);


	      // getting the flux
	      memcpy(F,stencil[iread].New.cell[icell].flux+inei*NVAR,sizeof(REAL)*NVAR);
	      
	      // update
	      one=POW(-1.,inei+1);
	      U.d +=F[0]*dtsurdx*one*0.125;
	      U.du+=F[1]*dtsurdx*one*0.125;
	      U.dv+=F[2]*dtsurdx*one*0.125;
	      U.dw+=F[3]*dtsurdx*one*0.125;
	      U.E +=F[4]*dtsurdx*one*0.125;

#ifdef DUAL_E
	      U.eint+=F[5]*dtsurdx*one*0.125;
#endif

#ifdef WRADHYD
	      U.dX+=F[6]*dtsurdx*one*0.125;
#endif 	      
	      //if(U.eint<0) abort();
	      //scatter back
	      U2W(&U,&W);
	      //getE(&W);// <<<<< COULD BE CRITICAL

	      memcpy(&(curoct->nei[vnei[inei]]->fieldnew),&W,sizeof(struct Wtype));
	    }
	  }
	}

      }

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}


// ================================================================
// ================================================================

//=======================================================================
//=======================================================================

int advancehydro(struct OCT **firstoct, int level, struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL dxcur, REAL dtnew){

  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot,nread;
  double t[10];
  double tg=0.,th=0.,tu=0.,ts=0.,tfu=0.,ttot=0.;
  
  
  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  nreadtot=0;
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next; 

      t[0]=MPI_Wtime();
  
      // ------------ gathering the stencil value values
      nextoct= gatherstencil(curoct,stencil,stride,cpu, &nread);
      
      if(nread>0){
	t[2]=MPI_Wtime();
      // ------------ solving the hydro
	    
      hydroM_sweepX(stencil,level,cpu->rank,nread,stride,dxcur,dtnew);   
      hydroM_sweepY(stencil,level,cpu->rank,nread,stride,dxcur,dtnew);
      hydroM_sweepZ(stencil,level,cpu->rank,nread,stride,dxcur,dtnew);

      // ------------ updating values within the stencil

      t[4]=MPI_Wtime();

      updatefield(curoct,stencil,nread,stride,cpu,dxcur,dtnew);
      
      // ------------ scatter back the FLUXES
      
      t[6]=MPI_Wtime();
   
      nextoct=scatterstencil(curoct,stencil, nread, cpu,dxcur,dtnew);


      t[8]=MPI_Wtime();

      nreadtot+=nread;
      ts+=(t[8]-t[6]);
      tu+=(t[6]-t[4]);
      th+=(t[4]-t[2]);
      tg+=(t[2]-t[0]);
      }

      
      
    }while((nextoct!=NULL)&&(nread>0));
  }
  
  //  printf("CPU %d | tgat=%e tcal=%e tup=%e tscat=%e\n",cpu->rank,tg,th,tu,ts);
  //if(cpu->rank==RANK_DISP) 

  return nreadtot;
}

// =========================================================================================
// =========================================================================================

void HydroSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL dtnew){

  int nread,nreadtot;;
  struct OCT *curoct;
  struct OCT *nextoct;
  
  REAL dxcur=POW(0.5,level);
  REAL one;
  struct Wtype W;
  struct Utype U;
  struct Wtype *Wloc;
  struct Utype *Uloc;
  int icell;
  int nocthydro=cpu->noct[level-1];
  double t[10];

  t[0]=MPI_Wtime();
#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&nocthydro,1,MPI_INT,MPI_SUM,cpu->comm);
#endif
  if(cpu->rank==RANK_DISP) printf("Start Hydro on %d octs with dt=%e on level %d with stride=%d\n",nocthydro,dtnew,level,stride);

  // ===== COMPUTING THE FLUXES
  
#ifndef GPUAXL
  nreadtot=advancehydro(firstoct,level,cpu,stencil,stride,dxcur,dtnew);
#else
  //nreadtot=advancehydro(firstoct,level,cpu,stencil,stride,dxcur,dtnew);
  nreadtot=advancehydroGPU(firstoct,level,cpu,stencil,stride,dxcur,dtnew);
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
	  if(isnan(curoct->cell[icell].fieldnew.u)){
	    printf("unsplit\n");
	    abort();
	  }
	  memcpy(&(curoct->cell[icell].field),&(curoct->cell[icell].fieldnew),sizeof(struct Wtype));
	}
	else{
	  // split case : -> quantities are averaged
	  struct OCT *child;
	  int i;
	  child=curoct->cell[icell].child;
	  memset(&W,0,sizeof(struct Wtype));
	  memset(&U,0,sizeof(struct Utype));
	  W.d=1e-13;
	  U.d=1e-13;
	  for(i=0;i<8;i++){
	    // FIX THIS WITH CONSERVATIVE QUANT
#ifdef CONSAVG
	    /// conservative average
	    Wloc=&child->cell[i].field;
	    W2U(Wloc,Uloc);
	    U.d+=Uloc->d*0.125;
	    U.du+=Uloc->du*0.125;
	    U.dv+=Uloc->dv*0.125;
	    U.dw+=Uloc->dw*0.125;
#ifdef DUAL_E
	    U.eint+=Uloc->eint*0.125;
	    U.E+=Uloc->E*0.125;
#endif

#ifdef WRADHYD
	    U.dX+=Uloc->dX*0.125;
#endif

#else
	    /// primitive average
	    W.d+=child->cell[i].field.d*0.125;
	    W.u+=child->cell[i].field.u*0.125;
	    W.v+=child->cell[i].field.v*0.125;
	    W.w+=child->cell[i].field.w*0.125;
	    W.p+=child->cell[i].field.p*0.125;
#ifdef WRADHYD
	    W.dX+=child->cell[i].field.dX*0.125;
#endif
#endif
	  }
#ifdef CONSAVG
	  U2W(&U,&W);
#endif
	  getE(&W);
	  memcpy(&curoct->cell[icell].field,&W,sizeof(struct Wtype));
	}


	//if(curoct->cell[icell].field.w <0.) abort();

      }
    }while(nextoct!=NULL);
  }

  t[9]=MPI_Wtime();
  
  if(cpu->rank==RANK_DISP){
#ifndef GPUAXL
    printf("==== CPU HYDRO TOTAL TIME =%e\n",t[9]-t[0]);
#else
    printf(" === GPU HYDRO TOTAL TIME =%e\n",t[9]-t[0]);
    //printf(" === GPU (DISABLED) HYDRO TOTAL TIME =%e\n",t[9]-t[0]);
#endif
  }
  
}


// ==============================================================
// ==============================================================
void clean_new_hydro(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu){

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

      for(icell=0;icell<8;icell++) {
	memcpy(&(curoct->cell[icell].fieldnew),&(curoct->cell[icell].field),sizeof(struct Wtype));
      }
     
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
	    memset(&(curoct->cell[icell].fieldnew),0,sizeof(struct Wtype));
	    curoct->cell[icell].fieldnew.d=1e-13;
	  }
	}
      }while(nextoct!=NULL);
    }
  }
#endif

}



#endif
