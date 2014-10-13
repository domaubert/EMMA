
#ifdef WHYDRO2

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include <mpi.h>
#include <cudpp.h>
#include "gpu_type.h"

#define NITERMAX 10
#define ERRTOL 1e-10


extern "C" struct OCT *gatherstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread);
extern "C" struct OCT *scatterstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew);
extern "C" void create_hydstencil_GPU(struct CPUINFO *cpu, int stride);
extern "C" int advancehydroGPU(struct OCT **firstoct, int level, struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL dxcur, REAL dtnew);
extern "C" void create_pinned_stencil(struct HGRID **stencil, int stride);
extern "C" void destroy_hydstencil_GPU(struct CPUINFO *cpu, int stride);
extern "C" void destroy_pinned_stencil(struct HGRID **stencil, int stride);

// ===================================================================
void create_hydstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaMalloc((void **)&(cpu->hyd_stencil),sizeof(struct HGRID)*stride);
}

// ===================================================================
void create_pinned_stencil(struct HGRID **stencil, int stride){
  cudaMallocHost( (void**)stencil, sizeof(struct HGRID)*stride );
}

// ===================================================================
void destroy_hydstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaFree(cpu->hyd_stencil);
}

// ===================================================================
void destroy_pinned_stencil(struct HGRID **stencil, int stride){
  cudaFreeHost(stencil);
}



// ==============================================================================================================

__device__ void dgetE(struct Wtype *W){
  W->E=W->p/(GAMMA-1.)+0.5*W->d*(W->u*W->u+W->v*W->v+W->w*W->w);
}

// =======================================================

__device__ void getcellnei_gpu_hydro(int cindex, int *neip, int *cell)
{
  switch(cindex){
  case 0:
    neip[0]=0;cell[0]=1;
    neip[1]=6;cell[1]=1;
    neip[2]=2;cell[2]=2;
    neip[3]=6;cell[3]=2;
    neip[4]=4;cell[4]=4;
    neip[5]=6;cell[5]=4;
    break;
  case 1:
    neip[0]=6;cell[0]=0;
    neip[1]=1;cell[1]=0;
    neip[2]=2;cell[2]=3;
    neip[3]=6;cell[3]=3;
    neip[4]=4;cell[4]=5;
    neip[5]=6;cell[5]=5;
    break;
  case 2:
    neip[0]=0;cell[0]=3;
    neip[1]=6;cell[1]=3;
    neip[2]=6;cell[2]=0;
    neip[3]=3;cell[3]=0;
    neip[4]=4;cell[4]=6;
    neip[5]=6;cell[5]=6;
    break;
  case 3:
    neip[0]=6;cell[0]=2;
    neip[1]=1;cell[1]=2;
    neip[2]=6;cell[2]=1;
    neip[3]=3;cell[3]=1;
    neip[4]=4;cell[4]=7;
    neip[5]=6;cell[5]=7;
    break;
  case 4:
    neip[0]=0;cell[0]=5;
    neip[1]=6;cell[1]=5;
    neip[2]=2;cell[2]=6;
    neip[3]=6;cell[3]=6;
    neip[4]=6;cell[4]=0;
    neip[5]=5;cell[5]=0;
    break;
  case 5:
    neip[0]=6;cell[0]=4;
    neip[1]=1;cell[1]=4;
    neip[2]=2;cell[2]=7;
    neip[3]=6;cell[3]=7;
    neip[4]=6;cell[4]=1;
    neip[5]=5;cell[5]=1;
    break;
  case 6:
    neip[0]=0;cell[0]=7;
    neip[1]=6;cell[1]=7;
    neip[2]=6;cell[2]=4;
    neip[3]=3;cell[3]=4;
    neip[4]=6;cell[4]=2;
    neip[5]=5;cell[5]=2;
    break;
  case 7:
    neip[0]=6;cell[0]=6;
    neip[1]=1;cell[1]=6;
    neip[2]=6;cell[2]=5;
    neip[3]=3;cell[3]=5;
    neip[4]=6;cell[4]=3;
    neip[5]=5;cell[5]=3;
    break;
  }

}

// ==================== converts U -> W
__device__ void dU2W(struct Utype *U, struct Wtype *W)
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
__device__ void dW2U(struct Wtype *W, struct Utype *U)
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



// ---------------------------------------------------------------
__device__ void dgetflux_X(struct Utype *U, REAL *f)
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

__device__ void dgetflux_Y(struct Utype *U, REAL *f)
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

__device__ void dgetflux_Z(struct Utype *U, REAL *f)
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



// ================== performs the difference between two Us

__device__ void ddiffU(struct Utype *U2, struct Utype *U1, struct Utype *UR){
  
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

__device__ void ddiffW(struct Wtype *W2, struct Wtype *W1, struct Wtype *WR){

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
__device__ void dminmod(struct Utype *Um, struct Utype *Up, struct Utype *Ur){

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

__device__ void dminmod_W(struct Wtype *Wm, struct Wtype *Wp, struct Wtype *Wr){

  REAL beta=1.; // 1. for MINBEE 2. for SUPERBEE
  // FLUX LIMITER

  if(Wp->d>0){
    Wr->d=FMAX(FMAX(0.,FMIN(beta*Wm->d,Wp->d)),FMIN(Wm->d,beta*Wp->d));
  }
  else{
    Wr->d=FMIN(FMIN(0.,FMAX(beta*Wm->d,Wp->d)),FMAX(Wm->d,beta*Wp->d));
  }


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


// ============= interp minmod ====================================================

__device__ void dinterpminmod(struct Utype *U0, struct Utype *Up, struct Utype *Dx, struct Utype *Dy, struct Utype *Dz,REAL dx,REAL dy,REAL dz){
  
  Up->d =U0->d  + dx*Dx->d  +dy*Dy->d  +dz*Dz->d;
  Up->du=U0->du + dx*Dx->du +dy*Dy->du +dz*Dz->du;
  Up->dv=U0->dv + dx*Dx->dv +dy*Dy->dv +dz*Dz->dv;
  Up->dw=U0->dw + dx*Dx->dw +dy*Dy->dw +dz*Dz->dw;
  Up->E =U0->E  + dx*Dx->E  +dy*Dy->E  +dz*Dz->E;
#ifdef DUAL_E
  Up->eint =U0->eint  + dx*Dx->eint  +dy*Dy->eint  +dz*Dz->eint;
#endif

}

// ============= interp minmod ====================================================

__device__ void dinterpminmod_W(struct Wtype *W0, struct Wtype *Wp, struct Wtype *Dx, struct Wtype *Dy, struct Wtype *Dz,REAL dx,REAL dy,REAL dz){
  
  Wp->d =W0->d +dx*Dx->d +dy*Dy->d +dz*Dz->d;
  Wp->u =W0->u +dx*Dx->u +dy*Dy->u +dz*Dz->u;
  Wp->v =W0->v +dx*Dx->v +dy*Dy->v +dz*Dz->v;
  Wp->w =W0->w +dx*Dx->w +dy*Dy->w +dz*Dz->w;
  Wp->p =W0->p +dx*Dx->p +dy*Dy->p +dz*Dz->p;
#ifdef WRADHYD
  Wp->dX =W0->dX +dx*Dx->dX +dy*Dy->dX +dz*Dz->dX;
#endif

}

//========================================================================================================================================

__device__ void  dmatrix_jacobian(struct Wtype *W0, REAL dt,REAL dx,struct Wtype *Dx,struct Wtype *Dy,struct Wtype *Dz, struct Wtype *Wt){


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

 
  for(j=0;j<5;j++){
    for(i=0;i<5;i++){
      W[i]+=M[i+j*5]*d[j];
      }
  }
#ifdef WRADHYD
  W[5]+=W0->w*Dx->dX+W0->w*Dx->dX;
#endif
  
  // ==== Final correction

  for(i=0;i<5;i++){
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

// ==============================================
__device__ void dMUSCL_BOUND2(struct HGRID *stencil, int ioct, int icell, struct Wtype *Wi,REAL dt,REAL dx){ 
  
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
  int shift;

#ifdef WGRAV
	  REAL f[3];
#ifdef CONSERVATIVE
	  struct Utype S;
	  struct Utype U;
#endif
#endif

	  getcellnei_gpu_hydro(icell, vnei, vcell); // we get the neighbors
	  
	  W0=&(stencil->oct[ioct].cell[icell].field);
	
	  // Limited Slopes
	  shift=1;
	  for(dir=0;dir<3;dir++){
	    
	    inei2=2*dir;
	    if(vnei[inei2]==6){
	      Wm=&(stencil->oct[ioct].cell[vcell[inei2]].field);
	    }
	    else{
	      Wm=&(stencil->oct[ioct-shift].cell[vcell[inei2]].field);
	    }

	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      Wp=&(stencil->oct[ioct].cell[vcell[inei2]].field);
	    }
	    else{
	      Wp=&(stencil->oct[ioct+shift].cell[vcell[inei2]].field);
	    }

	    ddiffW(Wp,W0,&Dp); 
	    ddiffW(W0,Wm,&Dm); 
	    
	    dminmod_W(&Dm,&Dp,D+dir);
	    shift*=3;
	  }


	  // build jacobian matrix product
	  
	  dmatrix_jacobian(W0,dt,dx,&D[0],&D[1],&D[2],&Wt); // Here Wt contains the evolution of the state

	  
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

	     /* if(Wi[idir].p<0) abort(); */
	    /* if(Wi[idir].d<0) abort(); */


#ifdef WGRAV
#ifndef NOCOUPLE

#ifdef PRIMITIVE
	    Wi[idir].u+=-f[0]*0.5*dt;
	    Wi[idir].v+=-f[1]*0.5*dt;
	    Wi[idir].w+=-f[2]*0.5*dt;
#endif

#ifdef CONSERVATIVE
 	    dW2U(&Wi[idir],&U);
	    U.d  +=S.d;
	    U.du +=S.du;
	    U.dv +=S.dv;
	    U.dw +=S.dw;
	    U.E  +=S.E;
	    dU2W(&U,&Wi[idir]);
#endif

#endif
#endif
	    dgetE(Wi+idir);
	    Wi[idir].a=SQRT(GAMMA*Wi[idir].p/Wi[idir].d);
#ifdef WRADHYD
 	    Wi[idir].dX=W0->dX; 
#endif 
	  }



	  
}

//========================================================================================
__device__ REAL dfrootprime(REAL p, struct Wtype1D *WL, struct Wtype1D *WR)
{
  
  REAL fL,fR;
  REAL AL,AR,BL,BR;

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  fL=(p>WL->p?SQRT(AL/(BL+p))*(1.-(p-WL->p)/(2.*(BL+p))):POW(p/WL->p,-(GAMMA+1)/(2.*GAMMA))/(WL->d*WL->a));
  fR=(p>WR->p?SQRT(AR/(BR+p))*(1.-(p-WR->p)/(2.*(BR+p))):POW(p/WR->p,-(GAMMA+1)/(2.*GAMMA))/(WR->d*WR->a));

  return fL+fR;
}


// ------------------------------------

__device__ REAL dfroot(REAL p, struct Wtype1D *WL, struct Wtype1D *WR, REAL *u)
{
  
  REAL fL,fR;
  REAL AL,AR,BL,BR;
  REAL Deltau;

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  fL=(p>WL->p?(p-WL->p)*SQRT(AL/(BL+p)):2.*WL->a/(GAMMA-1.)*(POW(p/WL->p,(GAMMA-1)/(2.*GAMMA))-1.));
  fR=(p>WR->p?(p-WR->p)*SQRT(AR/(BR+p)):2.*WR->a/(GAMMA-1.)*(POW(p/WR->p,(GAMMA-1)/(2.*GAMMA))-1.));
  
  Deltau=WR->u-WL->u;
  *u=0.5*(WL->u+WR->u)+0.5*(fR-fL);

  return fL+fR+Deltau;
}


//========================================================================================
//========================================================================================
__device__ REAL dfindPressure(struct Wtype1D *WL, struct Wtype1D *WR, int *niter, REAL *u)
{

  REAL ptr,pts,ppv;
  REAL ptr0,pts0,ppv0;
  REAL p,porg,dp;
  int i;
  REAL err;
  REAL unsurz=(2.0*GAMMA)/(GAMMA-1.0);
  REAL AL,AR,BL,BR,GL,GR;
  REAL pmin,pmax;
  REAL u2;

  pmin=FMIN(WL->p,WR->p);
  pmax=FMAX(WL->p,WR->p);
  
  // EXACT SOLVER

  // hybrid guess for pressure

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  ppv0=0.5*(WL->p+WR->p)-0.125*(WR->u-WL->u)*(WR->d+WL->d)*(WR->a+WL->a);
  ptr0=POW((WL->a+WR->a-0.5*(GAMMA-1)*(WR->u-WL->u))/(WL->a/POW(WL->p,1./unsurz)+WR->a/POW(WR->p,1./unsurz)),unsurz);

  ppv=FMAX(ERRTOL,ppv0);
  ptr=FMAX(ERRTOL,ptr0);
  
  GL=SQRT(AL/(ppv+BL));
  GR=SQRT(AR/(ppv+BR));

  pts0=(GL*WL->p+GR*WR->p-(WR->u-WL->u))/(GL+GR);
  pts=FMAX(ERRTOL,pts0);


  if(((pmax/pmin)<2.0)&&((pmin<=ppv)&&(ppv<=pmax))){
      p=ppv;
    }
  else{
    if(ppv<pmin){
      p=ptr;
    }
    else{
      p=pts;
    }
  }


  //p=0.5*(WL->p+WR->p);
  //p=FMAX(p,ERRTOL);

  *niter=0;
  for(i=0;i<NITERMAX;i++)
    {
      dp=dfroot(p,WL,WR,&u2)/dfrootprime(p,WL,WR);

      if(FABS(dp)<ERRTOL) break;
      while((p-dp)<0){ 
       	dp=dp*0.5; 
      } 

      porg=p;
      p=p-dp;
      err=2.*FABS(p-porg)/(FABS(p+porg));
      *niter=*niter+1;
      if(err<ERRTOL) break;
      if(dfroot(p,WL,WR,&u2)<ERRTOL) break;
    }

  dfroot(p,WL,WR,&u2); // last calculation to get u;

  *u=(REAL)u2;
  return p;
}


//========================================================================================
//========================================================================================
__device__ REAL dfindPressure_Hybrid(struct Wtype1D *WL, struct Wtype1D *WR, int *niter, REAL *ustar){
  REAL ppvrs;
  REAL dbar,abar;
  REAL pmax,pmin,pstar;
  REAL AL,AR,BL,BR,GL,GR;
  dbar=0.5*(WL->d+WR->d);
  abar=0.5*(WL->a+WR->a);
  ppvrs=0.5*((WL->p+WR->p)+(WL->u-WR->u)*dbar*abar);
  pmax=FMAX(WL->p,WR->p);
  pmin=FMIN(WL->p,WR->p);
  pstar=ppvrs;
  
  //printf("dbar=%e abar=%e ppvrs=%e pmax=%e pmin=%e pstar=%e\n",dbar,abar,ppvrs,pmax,pmin,pstar);

  if(((pmax/pmin)<2.)&&((pmin<pstar)&&(pstar<pmax))){
    // PVRS CASE
    pstar=ppvrs;
    *ustar=0.5*((WL->u+WR->u)+(WL->p-WR->p)/(dbar*abar));
  }
  else{
    if(pstar<pmin){
      //TRRS CASE
      REAL z=(GAMMA-1.)/(2.*GAMMA);
      REAL iz=(2.*GAMMA)/(GAMMA-1.);
      pstar=POW((WL->a+WR->a-(GAMMA-1.)/2.*(WR->u-WL->u))/(WL->a/POW(WL->p,z)+WR->a/POW(WR->p,z)),iz);
      *ustar=WL->u-2.*WL->a/(GAMMA-1.)*(POW(pstar/WL->p,z)-1.);
    }
    else{
      //TSRS CASE
      REAL p0;
      p0=FMAX(0.,ppvrs);
      
      AL=2./((GAMMA+1.)*WL->d);
      AR=2./((GAMMA+1.)*WR->d);
      
      BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
      BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

      GL=SQRT(AL/(p0+BL));
      GR=SQRT(AR/(p0+BR));

      pstar=(GL*WL->p+GR*WR->p-(WR->u-WL->u))/(GL+GR);
      *ustar=0.5*((WL->u+WR->u)+(pstar-WR->p)*GR-(pstar-WL->p)*GL);
    }
  }

  return pstar;

}




//====================================================================
__device__ void dspeedestimateX_HLLC(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR, REAL *pstar, REAL *ustar){

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

  //printf("%e %e %e %e ||| %e %e %e %e\n",WLloc.d,WLloc.u,WLloc.p,WLloc.a,WLloc.d,WRloc.u,WRloc.p,WRloc.a);

#if 1
  (*pstar)= dfindPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=dfindPressure(&WLloc,&WRloc,&n,ustar);
  //if((*pstar)<=0) printf("shhh pstar=%e %e %d\n",*pstar,*ustar,n);

  qL=(*pstar<=WL->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;
  if((*SL)>(*SR)){
    (*SL)=FMIN(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=FMAX(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
  }
#endif
  
  //if(isnan(*ustar)) printf("Hehey\n");
}

//====================================================================

void __device__ dspeedestimateY_HLLC(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR, REAL *pstar, REAL *ustar){

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

  (*pstar)=dfindPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=dfindPressure(&WLloc,&WRloc,&n,ustar);
  //  if((*pstar)<0) abort();

  qL=(*pstar<=WL->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;

  if((*SL)>(*SR)){
    (*SL)=FMIN(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=FMAX(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
    //abort();
  }
  //  if((*SL)>(*SR)) abort();
  //if(isnan(*ustar)) printf("Hehey y\n");

}



//====================================================================

void __device__ dspeedestimateZ_HLLC(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR, REAL *pstar, REAL *ustar){

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

  (*pstar)=dfindPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=dfindPressure(&WLloc,&WRloc,&n,ustar);
  //if((*pstar)<0) abort();

  qL=(*pstar<=WL->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:SQRT(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;
  if((*SL)>(*SR)){
    (*SL)=FMIN(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=FMAX(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
    //abort();
  }
  //if((*SL)>(*SR)) abort();
  //if(isnan(*ustar)) printf("Hehey z\n");

}


// =============================================================================================

__global__ void dhydroM_sweepZ(struct HGRID *stencil, int nread,REAL dx, REAL dt){

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
  REAL divu,divuloc;
#endif
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;

  i=bx*blockDim.x+tx;
  if(i<nread){
  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei_gpu_hydro(icell, vnei, vcell); // we get the neighbors
      
      
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

      dW2U(&Wold,&Uold); // primitive -> conservative

      REAL eold=Uold.eint;

      /* // MUSCL STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);

      dMUSCL_BOUND2(stencil+i, 13, icell, WC,dt,dx);// central

      for(iface=0;iface<2;iface++){
	inei=iface+4;
	memcpy(WC+iface,WC+inei,sizeof(struct Wtype)); // moving the data towards idx=0,1
	dW2U(WC+iface,UC+iface);
      }

      // Neighbor MUSCL reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+4;
	dMUSCL_BOUND2(stencil+i, ioct[vnei[inei]], vcell[inei], WT,dt,dx);// 
	memcpy(WN+iface,WT+idxnei[inei],sizeof(struct Wtype)); 
	//memcpy(WN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].field),sizeof(struct Wtype)); 
	dW2U(WN+iface,UN+iface);
	
	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }


      // Z DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems BOTTOM

      // Switching to Split description

      /* 	// =========================================== */

#ifdef RIEMANN_HLLC
      dspeedestimateZ_HLLC(&WN[0],&WC[0],&SL,&SR,&pstar,&ustar);

      if(SL>=0.){
	dgetflux_Z(&UN[0],FL);
#ifdef DUAL_E
	memcpy(&Us,&UN[0],sizeof(struct Utype));
#endif

      }
      else if(SR<=0.){
	dgetflux_Z(&UC[0],FL);
#ifdef DUAL_E
	memcpy(&Us,&UC[0],sizeof(struct Utype));
#endif
      }
      else if((SL<0.)&&(ustar>=0.)){
	dgetflux_Z(&UN[0],FL);
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
	dgetflux_Z(&UC[0],FL);
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
      dspeedestimateZ_HLLC(&WC[1],&WN[1],&SL,&SR,&pstar,&ustar);

      if(SL>=0.){
	dgetflux_Z(&UC[1],FR);
#ifdef DUAL_E
	memcpy(&Us,&UC[1],sizeof(struct Utype));
#endif

      }
      else if(SR<=0.){
	dgetflux_Z(&UN[1],FR);
#ifdef DUAL_E
	memcpy(&Us,&UN[1],sizeof(struct Utype));
#endif

      }
      else if((SL<0.)&&(ustar>=0.)){
	dgetflux_Z(&UC[1],FR);
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
	dgetflux_Z(&UN[1],FR);
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

      //ready for the next oct
  }
  }
}




//============================================================================
// =============================================================================================

__global__ void dhydroM_sweepY(struct HGRID *stencil,int nread,REAL dx, REAL dt){

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
  REAL divu,divuloc;
#endif
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
	
  i=bx*blockDim.x+tx;
  if(i<nread){
  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei_gpu_hydro(icell, vnei, vcell); // we get the neighbors
      
  
      
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
      dW2U(&Wold,&Uold); // primitive -> conservative

      REAL eold=Uold.eint;

      /* // MUSCL STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);

      dMUSCL_BOUND2(stencil+i, 13, icell, WC,dt,dx);// central

      for(iface=0;iface<2;iface++){
	inei=iface+2;
	memcpy(WC+iface,WC+inei,sizeof(struct Wtype)); // moving the data towards idx=0,1
	//memcpy(WC+iface,&Wold,sizeof(struct Wtype)); // moving the data towards idx=0,1

	dW2U(WC+iface,UC+iface);
      }

      // Neighbor MUSCL reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+2;
	dMUSCL_BOUND2(stencil+i, ioct[vnei[inei]], vcell[inei], WT,dt,dx);// 
	memcpy(WN+iface,WT+idxnei[inei],sizeof(struct Wtype)); 
	//memcpy(WN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].field),sizeof(struct Wtype)); 

       	dW2U(WN+iface,UN+iface);
	
	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }




      // Y DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems FRONT

      // Switching to Split description

/* 	// =========================================== */

#ifdef RIEMANN_HLLC
      dspeedestimateY_HLLC(&WN[0],&WC[0],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  dgetflux_Y(&UN[0],FL);
#ifdef DUAL_E
	  memcpy(&Us,&UN[0],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  dgetflux_Y(&UC[0],FL);
#ifdef DUAL_E
	  memcpy(&Us,&UC[0],sizeof(struct Utype));
#endif
	}
	else if((SL<0.)&&(ustar>=0.)){
	  dgetflux_Y(&UN[0],FL);
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
	  dgetflux_Y(&UC[0],FL);
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


#endif
	// ===========================================




      // --------- solving the Riemann Problems BACK


      // Switching to Split description

	//=====================================================

#ifdef RIEMANN_HLLC
	dspeedestimateY_HLLC(&WC[1],&WN[1],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  dgetflux_Y(&UC[1],FR);
#ifdef DUAL_E
	  memcpy(&Us,&UC[1],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  dgetflux_Y(&UN[1],FR);
#ifdef DUAL_E
	  memcpy(&Us,&UN[1],sizeof(struct Utype));
#endif

	}
	else if((SL<0.)&&(ustar>=0.)){
	  dgetflux_Y(&UC[1],FR);
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
	  dgetflux_Y(&UN[1],FR);
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
	
	memcpy(stencil[i].New.cell[icell].flux+2*NVAR,FL,sizeof(REAL)*NVAR);
	memcpy(stencil[i].New.cell[icell].flux+3*NVAR,FR,sizeof(REAL)*NVAR);
	
	stencil[i].New.cell[icell].divu=divu;
	
    //ready for the next oct
  }
  }
}

//===================================================================================================
//===================================================================================================

__global__ void dhydroM_sweepX(struct HGRID *stencil, int nread,REAL dx, REAL dt){

  //printf("IN\n");
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
  REAL divu,divuloc;
#endif
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  
  i=bx*blockDim.x+tx;
  if(i<nread){
  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei_gpu_hydro(icell, vnei, vcell); // we get the neighbors
      
  
      
      memset(FL,0,sizeof(REAL)*NVAR);
      memset(FR,0,sizeof(REAL)*NVAR);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].field);
      
#ifdef DUAL_E
      divu=stencil[i].New.cell[icell].divu;
#endif      

      Wold.d=curcell->d;
      Wold.u=curcell->u;;
      Wold.v=curcell->v;;
      Wold.w=curcell->w;;
      Wold.p=curcell->p;;
      Wold.a=SQRT(GAMMA*Wold.p/Wold.d);

      //printf("dt=%e dx=%e Old =%e %e %e %e %e\n",dt,dx,Wold.d,Wold.u,Wold.v,Wold.w,Wold.p);
/* #ifdef WRADHYD */
/*       Wold.X=curcell->X; */
/* #endif */
      dW2U(&Wold,&Uold); // primitive -> conservative
      REAL eold=Uold.eint;

      /* // MUSCL STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);

      dMUSCL_BOUND2(stencil+i, 13, icell, WC,dt,dx);// central
      for(iface=0;iface<2;iface++){
	dW2U(WC+iface,UC+iface);
      }

      // Neighbor MUSCL reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface;
	dMUSCL_BOUND2(stencil+i, ioct[vnei[inei]], vcell[inei], WT,dt,dx);// 
	memcpy(WN+iface,WT+idxnei[inei],sizeof(struct Wtype)); 
       	dW2U(WN+iface,UN+iface);
	
	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }




      // X DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems LEFT

      // Switching to Split description

/* 	// =========================================== */

#ifdef RIEMANN_HLLC
      //printf("Ho %e %e %e %e|| %e %e %e %e\n",WC[0].d,WC[0].u,WC[0].p,WC[0].a,WN[0].d,WN[0].u,WN[0].p,WN[0].a);
      dspeedestimateX_HLLC(&WN[0],&WC[0],&SL,&SR,&pstar,&ustar);
      //printf("Ha\n");

      if(SL>=0.){
	dgetflux_X(&UN[0],FL);
#ifdef DUAL_E
	memcpy(&Us,&UN[0],sizeof(struct Utype));
#endif

	}
      else if(SR<=0.){
	dgetflux_X(&UC[0],FL);
#ifdef DUAL_E
	memcpy(&Us,&UC[0],sizeof(struct Utype));
#endif
	}
      else if((SL<0.)&&(ustar>=0.)){
	dgetflux_X(&UN[0],FL);
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
	dgetflux_X(&UC[0],FL);
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
      dspeedestimateX_HLLC(&WC[1],&WN[1],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  dgetflux_X(&UC[1],FR);
#ifdef DUAL_E
	  memcpy(&Us,&UC[1],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  dgetflux_X(&UN[1],FR);
#ifdef DUAL_E
	  memcpy(&Us,&UN[1],sizeof(struct Utype));
#endif

	}
	else if((SL<0.)&&(ustar>=0.)){
	  dgetflux_X(&UC[1],FR);
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
	  dgetflux_X(&UN[1],FR);
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

    //ready for the next oct
  }
  }
}


// ==============================================================================================================
// ==============================================================================================================

__global__ void dupdatefield(struct HGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew)
{
  int i,icell;
  struct Utype U;
  REAL one;
  int flx;
  REAL dtsurdx=dtnew/dxcur;
  REAL F[NFLUX];


  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  i=bx*blockDim.x+tx;
  if(i<nread){
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

// =======================================================

//=======================================================================
//=======================================================================

int advancehydroGPU(struct OCT **firstoct, int level, struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL dxcur, REAL dtnew){

  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *curoct0;
  int nreadtot,nread;
  /* REAL t[10]; */
  /* REAL tg=0.,th=0.,tu=0.,ts=0.; */
  int is;
  int offset;
  CUDA_CHECK_ERROR("Hydro start");
  
  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  nreadtot=0;
  int ng;
  int nt;

  cudaStream_t stream[cpu->nstream]; // Not fully regular expression
  int vnread[cpu->nstream];
  // creating the streams
  for(is=0;is<cpu->nstream;is++){
    cudaStreamCreate(&stream[is]);
  }
  
  // Calculations

  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct0=nextoct;
      curoct=curoct0;
      
      //t[0]=MPI_Wtime();
  
#if 1
      offset=0;
      // streaming ====================
      for(is=0;is<cpu->nstream;is++){

	// ------------ gathering the stencil value values
	//printf("offser=%d\n",offset);
	curoct=nextoct;
	if(curoct!=NULL){
	  nextoct= gatherstencil(curoct,stencil+offset,stride/cpu->nstream,cpu, vnread+is);
	  if(vnread[is]!=0){
	  ng=((vnread[is]-1)/cpu->nthread)+1; // +1 is for leftovers
	
	  if(ng==1){
	    nt=vnread[is];
	  }
	  else{
	    nt=cpu->nthread;
	  }

	  dim3 gridoct(ng);
	  dim3 blockoct(nt);
	
	
#ifndef NOCOMP     
	  cudaMemcpyAsync(cpu->hyd_stencil+offset,stencil+offset,vnread[is]*sizeof(struct HGRID),cudaMemcpyHostToDevice,stream[is]);  

	  //printf("Sweep hydro dt=%e dx=%e\n",dtnew,dxcur);
	  //CUDA_CHECK_ERROR("Sweep hydro");
	  // ------------ solving the hydro
	  dhydroM_sweepX<<<gridoct,blockoct,0,stream[is]>>>(cpu->hyd_stencil+offset,vnread[is],dxcur,dtnew);
	  dhydroM_sweepY<<<gridoct,blockoct,0,stream[is]>>>(cpu->hyd_stencil+offset,vnread[is],dxcur,dtnew);  
	  dhydroM_sweepZ<<<gridoct,blockoct,0,stream[is]>>>(cpu->hyd_stencil+offset,vnread[is],dxcur,dtnew);  
	
	  //printf("Sweep hydro stop\n");
	  // ------------ updating values within the stencil

	  dupdatefield<<<gridoct,blockoct,0,stream[is]>>>(cpu->hyd_stencil+offset,vnread[is],stride,cpu,dxcur,dtnew);
	
	  cudaMemcpyAsync(stencil+offset,cpu->hyd_stencil+offset,vnread[is]*sizeof(struct HGRID),cudaMemcpyDeviceToHost,stream[is]);  
#endif
	
	  offset+=vnread[is];
	  }
	}
      }
#endif

      
	/* dev_updatefield<<<gridoct2,blockoct2>>>(cpu->hyd_stencil,nread,stride,cpu,dxcur,dtnew); */
      
	cudaDeviceSynchronize();
	// ------------ scatter back the FLUXES
	//cudaMemcpy(stencil,cpu->hyd_stencil,nread*sizeof(struct HGRID),cudaMemcpyDeviceToHost);  
	nread=offset;
	nextoct=scatterstencil(curoct0,stencil, nread, cpu,dxcur,dtnew);
      
	//t[8]=MPI_Wtime();
      
	nreadtot+=nread;
      
      
	/* ts+=(t[8]-t[6]); */
	/* tu+=(t[6]-t[4]); */
	/* th+=(t[4]-t[2]); */
	/* tg+=(t[2]-t[0]); */
	//printf("Start Error Hyd =%s nreadtot=%d\n",cudaGetErrorString(cudaGetLastError()),nreadtot);
      }while(nextoct!=NULL);
    
  }
  //printf("GPU | tgat=%e tcal=%e tup=%e tscat=%e\n",tg,th,tu,ts);

  // Destroying the streams
  for(is=0;is<cpu->nstream;is++){
    cudaStreamDestroy(stream[is]);
  }

  //printf("Start Error Hyd =%s nreadtot=%d\n",cudaGetErrorString(cudaGetLastError()),nreadtot);
  CUDA_CHECK_ERROR("Hydro Stop");

  return nreadtot;
}

#endif
