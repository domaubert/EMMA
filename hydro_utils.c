#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include "vector.h"

#include <mpi.h>

#define NITERMAX 10
#define ERRTOL 1e-10
#define CFL 0.4

#ifdef WHYDRO2

/* // ==================== converts U -> W */
/* void U2W(struct Utype *U, struct Wtype *W) */
/* { */
/*   W->d=U->d; */
/*   W->u=U->du/U->d; */
/*   W->v=U->dv/U->d; */
/*   W->w=U->dw/U->d; */
/*   W->p=(GAMMA-1.)*(U->E-(pow(W->u,2)+pow(W->v,2)+pow(W->w,2))*0.5*U->d); */
/*   W->a=sqrt(GAMMA*W->p/W->d); */
/* } */

/* // ==================== converts W -> U */
/* void W2U(struct Wtype *W, struct Utype *U) */
/* { */
/*   U->d=W->d; */
/*   U->du=W->d*W->u; */
/*   U->dv=W->d*W->v; */
/*   U->dw=W->d*W->w; */
/*   U->E=W->d*(0.5*(pow(W->u,2)+pow(W->v,2)+pow(W->w,2))+W->p/((GAMMA-1.)*W->d)); */
/* } */


// ==================== converts U -> W
void U2W(struct Utype *U, struct Wtype *W)
{
  W->d=U->d;
  W->u=U->du/U->d;
  W->v=U->dv/U->d;
  W->w=U->dw/U->d;
  W->p=(GAMMA-1.)*(U->E-((U->du)*(U->du)+(U->dv)*(U->dv)+(U->dw)*(U->dw))/(U->d)*0.5);
#ifdef DUAL_E
  if(W->p/((GAMMA-1.)*U->E)<1e-3) W->p=U->eint*(GAMMA-1.);
#endif
  W->a=sqrt(GAMMA*W->p/W->d);
}

// ==================== converts W -> U
void W2U(struct Wtype *W, struct Utype *U)
{
  U->d=W->d;
  U->du=W->d*W->u;
  U->dv=W->d*W->v;
  U->dw=W->d*W->w;
  U->E=W->d*(0.5*(W->u*W->u+W->v*W->v+W->w*W->w)+W->p/((GAMMA-1.)*W->d));

#ifdef DUAL_E
  U->eint=W->p/(GAMMA-1.);
#endif

}


// ================== performs the difference between two Us

void diffU(struct Utype *U2, struct Utype *U1, struct Utype *UR){
  
  UR->d =U2->d - U1->d;
  UR->du=U2->du- U1->du;
  UR->dv=U2->dv- U1->dv;
  UR->dw=U2->dw- U1->dw;
  UR->E =U2->E - U1->E;
}

// ================== performs the difference between two Ws

void diffW(struct Wtype *W2, struct Wtype *W1, struct Wtype *WR){

  WR->d=W2->d- W1->d;
  WR->u=W2->u- W1->u;
  WR->v=W2->v- W1->v;
  WR->w=W2->w- W1->w;
  WR->p=W2->p- W1->p;
}


// ================= minmod

void minmod(struct Utype *Um, struct Utype *Up, struct Utype *Ur){

  REAL beta=2.; // 1. for MINBEE 2. for SUPERBEE
  // FLUX LIMITER

  if(Up->d>0){
    Ur->d=fmax(fmax(0.,fmin(beta*Um->d,Up->d)),fmin(Um->d,beta*Up->d));
  }
  else{
    Ur->d=fmin(fmin(0.,fmax(beta*Um->d,Up->d)),fmax(Um->d,beta*Up->d));
  }


  if(Up->du>0){
    Ur->du=fmax(fmax(0.,fmin(beta*Um->du,Up->du)),fmin(Um->du,beta*Up->du));
  }
  else{
    Ur->du=fmin(fmin(0.,fmax(beta*Um->du,Up->du)),fmax(Um->du,beta*Up->du));
  }


  if(Up->dv>0){
    Ur->dv=fmax(fmax(0.,fmin(beta*Um->dv,Up->dv)),fmin(Um->dv,beta*Up->dv));
  }
  else{
    Ur->dv=fmin(fmin(0.,fmax(beta*Um->dv,Up->dv)),fmax(Um->dv,beta*Up->dv));
  }


  if(Up->dw>0){
    Ur->dw=fmax(fmax(0.,fmin(beta*Um->dw,Up->dw)),fmin(Um->dw,beta*Up->dw));
  }
  else{
    Ur->dw=fmin(fmin(0.,fmax(beta*Um->dw,Up->dw)),fmax(Um->dw,beta*Up->dw));
  }


  if(Up->E>0){
    Ur->E=fmax(fmax(0.,fmin(beta*Um->E,Up->E)),fmin(Um->E,beta*Up->E));
  }
  else{
    Ur->E=fmin(fmin(0.,fmax(beta*Um->E,Up->E)),fmax(Um->E,beta*Up->E));
  }


}


//===============================================
void minmod2(struct Utype *Um, struct Utype *Up, struct Utype *Ur){
  REAL r;
  REAL xi;
  REAL w=0.;
  REAL beta=1.0;
  REAL xd,xu,xv,xw,xe;
  // SLOPE LIMITER

  if(Up->d==0.){
    xi=0.;}
  else{
    r=Um->d/Up->d;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xd=xi;

  

  if(Up->du==0.){
    xi=0.;}
  else{
    r=Um->du/Up->du;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xu=xi;
  //printf(" %e ",xi);
  
  if(Up->dv==0.){
    xi=0.;}
  else{
    r=Um->dv/Up->dv;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xi=xv;
  //printf(" %e ",xi);
 
  if(Up->dw==0.){
    xi=0.;}
  else{
    r=Um->dw/Up->dw;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xw=xi;
  //printf(" %e ",xi);

  if(Up->E==0.){
    xi=0.;}
  else{
    r=Um->E/Up->E;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xe=xi;
  
  /* xi=fmin(xd,xu);  */
  /* xi=fmin(xv,xi);  */
  /* xi=fmin(xw,xi);  */
  /* xi=fmin(xe,xi);  */

  //printf(" %e %e %e %e %e\n",xd,xu,xv,xw,xe);

  //  abort();
  Ur->d =(0.5*(1.+w)*Um->d +0.5*(1.-w)*Up->d )*xd;
  Ur->du=(0.5*(1.+w)*Um->du+0.5*(1.-w)*Up->du)*xu;
  Ur->dv=(0.5*(1.+w)*Um->dv+0.5*(1.-w)*Up->dv)*xv;
  Ur->dw=(0.5*(1.+w)*Um->dw+0.5*(1.-w)*Up->dw)*xw;
  Ur->E =(0.5*(1.+w)*Um->E +0.5*(1.-w)*Up->E )*xe;

}

//===============================================
//===============================================

void minmod_W(struct Wtype *Wm, struct Wtype *Wp, struct Wtype *Wr){

  REAL beta=1.; // 1. for MINBEE 2. for SUPERBEE
  // FLUX LIMITER

  if(Wp->d>0){
    Wr->d=fmax(fmax(0.,fmin(beta*Wm->d,Wp->d)),fmin(Wm->d,beta*Wp->d));
  }
  else{
    Wr->d=fmin(fmin(0.,fmax(beta*Wm->d,Wp->d)),fmax(Wm->d,beta*Wp->d));
  }


  if(Wp->u>0){
    Wr->u=fmax(fmax(0.,fmin(beta*Wm->u,Wp->u)),fmin(Wm->u,beta*Wp->u));
  }
  else{
    Wr->u=fmin(fmin(0.,fmax(beta*Wm->u,Wp->u)),fmax(Wm->u,beta*Wp->u));
  }


  if(Wp->v>0){
    Wr->v=fmax(fmax(0.,fmin(beta*Wm->v,Wp->v)),fmin(Wm->v,beta*Wp->v));
  }
  else{
    Wr->v=fmin(fmin(0.,fmax(beta*Wm->v,Wp->v)),fmax(Wm->v,beta*Wp->v));
  }


  if(Wp->w>0){
    Wr->w=fmax(fmax(0.,fmin(beta*Wm->w,Wp->w)),fmin(Wm->w,beta*Wp->w));
  }
  else{
    Wr->w=fmin(fmin(0.,fmax(beta*Wm->w,Wp->w)),fmax(Wm->w,beta*Wp->w));
  }


  if(Wp->p>0){
    Wr->p=fmax(fmax(0.,fmin(beta*Wm->p,Wp->p)),fmin(Wm->p,beta*Wp->p));
  }
  else{
    Wr->p=fmin(fmin(0.,fmax(beta*Wm->p,Wp->p)),fmax(Wm->p,beta*Wp->p));
  }


}


//===============================================
//===============================================

void minmod2_W(struct Wtype *Wm, struct Wtype *Wp, struct Wtype *Wr){
  REAL r;
  REAL xi;
  REAL w=0.;
  REAL beta=1.0;
  REAL xd,xu,xv,xw,xe;
  // SLOPE LIMITER

  if(Wp->d==0.){
    xi=0.;}
  else{
    r=Wm->d/Wp->d;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xd=xi;

  if(Wp->u==0.){
    xi=0.;}
  else{
    r=Wm->u/Wp->u;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xu=xi;
  //printf(" %e ",xi);
  
  if(Wp->v==0.){
    xi=0.;}
  else{
    r=Wm->v/Wp->v;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xi=xv;
  //printf(" %e ",xi);
 
  if(Wp->w==0.){
    xi=0.;}
  else{
    r=Wm->w/Wp->w;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xw=xi;
  //printf(" %e ",xi);

  if(Wp->p==0.){
    xi=0.;}
  else{
    r=Wm->p/Wp->p;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xe=xi;
  
  /* xi=fmin(xd,xu);  */
  /* xi=fmin(xv,xi);  */
  /* xi=fmin(xw,xi);  */
  /* xi=fmin(xe,xi);  */

  //printf(" %e %e %e %e %e\n",xd,xu,xv,xw,xe);

  //  abort();
  Wr->d =(0.5*(1.+w)*Wm->d +0.5*(1.-w)*Wp->d)*xd;
  Wr->u =(0.5*(1.+w)*Wm->u +0.5*(1.-w)*Wp->u)*xu;
  Wr->v =(0.5*(1.+w)*Wm->v +0.5*(1.-w)*Wp->v)*xv;
  Wr->w =(0.5*(1.+w)*Wm->w +0.5*(1.-w)*Wp->w)*xw;
  Wr->p =(0.5*(1.+w)*Wm->p +0.5*(1.-w)*Wp->p)*xe*0;

}




// ============= interp minmod

void interpminmod(struct Utype *U0, struct Utype *Up, struct Utype *Dx, struct Utype *Dy, struct Utype *Dz,REAL dx,REAL dy,REAL dz){
  
  Up->d =U0->d  + dx*Dx->d  +dy*Dy->d  +dz*Dz->d;
  Up->du=U0->du + dx*Dx->du +dy*Dy->du +dz*Dz->du;
  Up->dv=U0->dv + dx*Dx->dv +dy*Dy->dv +dz*Dz->dv;
  Up->dw=U0->dw + dx*Dx->dw +dy*Dy->dw +dz*Dz->dw;
  Up->E =U0->E  + dx*Dx->E  +dy*Dy->E  +dz*Dz->E;

}


void interpminmod_W(struct Wtype *W0, struct Wtype *Wp, struct Wtype *Dx, struct Wtype *Dy, struct Wtype *Dz,REAL dx,REAL dy,REAL dz){
  
  Wp->d =W0->d +dx*Dx->d +dy*Dy->d +dz*Dz->d;
  Wp->u =W0->u +dx*Dx->u +dy*Dy->u +dz*Dz->u;
  Wp->v =W0->v +dx*Dx->v +dy*Dy->v +dz*Dz->v;
  Wp->w =W0->w +dx*Dx->w +dy*Dy->w +dz*Dz->w;
  Wp->p =W0->p +dx*Dx->p +dy*Dy->p +dz*Dz->p;

}


// ==============================================

void coarse2fine_hydro(struct CELL *cell, struct Wtype *Wi){ 


	  struct OCT * oct;
	  
	  struct Utype U0;
	  struct Utype Up;
	  struct Utype Um;
	  struct Utype Dp,Dm;
	  struct Utype D[3];
	  struct Wtype *W;
	  int inei2;
	  int vcell[6],vnei[6];
	  int dir;
	  REAL dxcur;

	  oct=cell2oct(cell);
	  getcellnei(cell->idx, vnei, vcell); // we get the neighbors
	  dxcur=pow(0.5,oct->level);

	  W=&(cell->field);
	  W2U(W,&U0);
	  // Limited Slopes
	  for(dir=0;dir<3;dir++){
	    
	    inei2=2*dir;
	    if(vnei[inei2]==6){
	      W=&(oct->cell[vcell[inei2]].field);
	    }
	    else{
	      W=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].field);

#ifdef TRANSXM
	      if((oct->nei[vnei[inei2]]->child->x-oct->x)>0.5){
		W=&(cell->field);
	      }
#endif

#ifdef TRANSYM
	      if((oct->nei[vnei[inei2]]->child->y-oct->y)>0.5){
		W=&(cell->field);
	      }
#endif

#ifdef TRANSZM
	      // if((oct->nei[vnei[inei2]]->child->z-oct->z)>0.5){
	      if((oct->z==0.)&&(inei2==4)){
		W=&(cell->field);
#ifdef REFZM
		W->w*=-1.0;
		dxcur=1./pow(2,oct->level);
		W->p=W->p+GRAV*W->d*dxcur;
#endif
	      }
#endif


	    }
	    W2U(W,&Um);

	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      W=&(oct->cell[vcell[inei2]].field);
	    }
	    else{
	      W=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].field);

#ifdef TRANSXP
	      if((oct->nei[vnei[inei2]]->child->x-oct->x)<0.){
		W=&(cell->field);
	      }
#endif

#ifdef TRANSYP
	      if((oct->nei[vnei[inei2]]->child->y-oct->y)<0.){
		W=&(cell->field);
	      }
#endif

#ifdef TRANSZP
	      // if((oct->nei[vnei[inei2]]->child->z-oct->z)<0.){
	      if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
		W=&(cell->field);
#ifdef REFZP
	    W->w*=-1.0;
	    dxcur=1./pow(2,oct->level);
	    W->p=W->p-GRAV*W->d*dxcur;
#endif

	      }
#endif

	    }
	    W2U(W,&Up);

	    diffU(&Up,&U0,&Dp); 
	    diffU(&U0,&Um,&Dm); 
	    
	    
	    minmod2(&Dm,&Dp,D+dir);
	  }

	  // Interpolation
	  int ix,iy,iz;
	  int icell;

	  for(iz=0;iz<2;iz++){
	    for(iy=0;iy<2;iy++){
	      for(ix=0;ix<2;ix++){
		icell=ix+iy*2+iz*4;
		interpminmod(&U0,&Up,D,D+1,D+2,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation
		U2W(&Up,Wi+icell);
	      }
	    }
	  }

}


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
	  dxcur=pow(0.5,oct->level);

	  W0=&(cell->field);
	  // Limited Slopes
	  for(dir=0;dir<3;dir++){
	    
	    inei2=2*dir;
	    if(vnei[inei2]==6){
	      Wm=&(oct->cell[vcell[inei2]].field);
	    }
	    else{
	      Wm=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].field);

#ifdef TRANSXM
	      if((oct->nei[vnei[inei2]]->child->x-oct->x)>0.5){
		Wm=&(cell->field);
	      }
#endif

#ifdef TRANSYM
	      if((oct->nei[vnei[inei2]]->child->y-oct->y)>0.5){
		Wm=&(cell->field);
	      }
#endif

#ifdef TRANSZM
	      //if((oct->nei[vnei[inei2]]->child->z-oct->z)>0.5){
	      if((oct->z==0.)&&(inei2==4)){
		Wm=&(cell->field);
#ifdef REFZM
		// !!!!DANGEREUX CA !
		Wm->w*=-1.0;
		dxcur=1./pow(2,oct->level);
		Wm->p=Wm->p+GRAV*Wm->d*dxcur;
#endif
	      }
#endif


	    }
	    

	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      Wp=&(oct->cell[vcell[inei2]].field);
	    }
	    else{
	      Wp=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].field);

#ifdef TRANSXP
	      if((oct->nei[vnei[inei2]]->child->x-oct->x)<0.){
		Wp=&(cell->field);
	      }
#endif

#ifdef TRANSYP
	      if((oct->nei[vnei[inei2]]->child->y-oct->y)<0.){
		Wp=&(cell->field);
	      }
#endif

#ifdef TRANSZP
	      //if((oct->nei[vnei[inei2]]->child->z-oct->z)<0.){
	      if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
		Wp=&(cell->field);
#ifdef REFZP
		Wp->w*=-1.0;
		dxcur=1./pow(2,oct->level);
		Wp->p=Wp->p-GRAV*Wp->d*dxcur;
#endif
		
	      }
#endif

	    }

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
		Wint.a=sqrt(GAMMA*Wint.p/Wint.d);
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
  double Deltau;

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
    printf("DIVERGENCE p0=%e dp=%e p=%e fprime=%e err=%e\n",porg,dp,p,frootprime(p,WL,WR),err);
    abort();
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
	  SL=WL->u-WL->a*sqrt((GAMMA+1.)/(2.*GAMMA)*pstar/WL->p+(GAMMA-1.)/(2.*GAMMA)); 
	  
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
	      W->a=sqrt(GAMMA*W->p/W->d);
	    }
	}
      else
	{
	  // left fan
	  astar=WL->a*pow(pstar/WL->p,(GAMMA-1.)/(2.*GAMMA)); // sound speed behind the fan
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
		  W->d=WL->d*pow(pstar/WL->p,1./GAMMA);
		  W->u=ustar;
		  W->p=pstar;
		  W->a=sqrt(GAMMA*W->p/W->d);
		}
	      else
		{
		  W->d=WL->d*pow(2./(GAMMA+1.)+(GAMMA-1.)/((GAMMA+1.)*WL->a)*(WL->u-S),2./(GAMMA-1.));
		  W->u=2./(GAMMA+1.)*(WL->a+(GAMMA-1.)/2.*WL->u+S);
		  W->p=WL->p*pow(2./(GAMMA+1.)+(GAMMA-1.)/((GAMMA+1.)*WL->a)*(WL->u-S),2.*GAMMA/(GAMMA-1.));
		  W->a=sqrt(GAMMA*W->p/W->d);
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
	  SR=WR->u+WR->a*sqrt((GAMMA+1.)/(2.*GAMMA)*pstar/WR->p+(GAMMA-1.)/(2.*GAMMA)); 
	  
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
	      W->a=sqrt(GAMMA*W->p/W->d);
	    }
	}
      else
	{
	  // Right fan
	  astar=WR->a*pow(pstar/WR->p,(GAMMA-1.)/(2.*GAMMA)); // sound speed behind the fan
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
		  W->d=WR->d*pow(pstar/WR->p,1./GAMMA);
		  W->u=ustar;
		  W->p=pstar;
		  W->a=sqrt(GAMMA*W->p/W->d);
		}
	      else
		{
		  W->d=WR->d*pow(2./(GAMMA+1.)-(GAMMA-1.)/((GAMMA+1.)*WR->a)*(WR->u-S),2./(GAMMA-1.));
		  W->u=2./(GAMMA+1.)*(-WR->a+(GAMMA-1.)/2.*WR->u+S);
		  W->p=WR->p*pow(2./(GAMMA+1.)-(GAMMA-1.)/((GAMMA+1.)*WR->a)*(WR->u-S),2.*GAMMA/(GAMMA-1.));
		  W->a=sqrt(GAMMA*W->p/W->d);
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

}

// ---------------------------------------------------------------

void getflux_Y(struct Utype *U, REAL *f)
{
  f[0]=U->dv;
  f[1]=U->dv*U->du/U->d;
  f[2]=0.5*(3.-GAMMA)*U->dv*U->dv/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dw*U->dw)/U->d;
  f[3]=U->dv*U->dw/U->d;
  f[4]=GAMMA*U->dv/U->d*U->E-0.5*(GAMMA-1.)*U->dv/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}

// ---------------------------------------------------------------

void getflux_Z(struct Utype *U, REAL *f)
{
  f[0]=U->dw;
  f[1]=U->dw*U->du/U->d;
  f[2]=U->dw*U->dv/U->d;
  f[3]=0.5*(3.-GAMMA)*U->dw*U->dw/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dv*U->dv)/U->d;
  f[4]=GAMMA*U->dw/U->d*U->E-0.5*(GAMMA-1.)*U->dw/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}



// ==================================================================

void MUSCL_BOUND(struct HGRID *stencil, int ioct, int icell, struct Utype *Ui,REAL dt,REAL dx){ 

	  struct OCT * oct;
	  
	  struct Utype U0;
	  struct Utype Up;
	  struct Utype Um;
	  struct Utype Dp,Dm;
	  struct Utype D[3];
	  struct Wtype *W;
	  int inei2;
	  int vcell[6],vnei[6];
	  int dir;

	  getcellnei(icell, vnei, vcell); // we get the neighbors
	  
	  W=&(stencil->oct[ioct].cell[icell].field);
	  W2U(W,&U0);
	
	  // Limited Slopes
	  for(dir=0;dir<3;dir++){
	    
	    inei2=2*dir;
	    if(vnei[inei2]==6){
	      W=&(stencil->oct[ioct].cell[vcell[inei2]].field);
	    }
	    else{
	      W=&(stencil->oct[ioct-(int)pow(3,dir)].cell[vcell[inei2]].field);
	    }
	    W2U(W,&Um);

	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      W=&(stencil->oct[ioct].cell[vcell[inei2]].field);
	    }
	    else{
	      W=&(stencil->oct[ioct+(int)pow(3,dir)].cell[vcell[inei2]].field);
	    }
	    W2U(W,&Up);

	    diffU(&Up,&U0,&Dp); 
	    diffU(&U0,&Um,&Dm); 
	    
	    minmod2(&Dm,&Dp,D+dir);
	  }


	  //Computing the Boundary Extrapolated Values

	  REAL ix[]={-0.5,0.5,0.0,0.0,0.0,0.0};
	  REAL iy[]={0.0,0.0,-0.5,0.5,0.0,0.0};
	  REAL iz[]={0.0,0.0,0.0,0.0,-0.5,0.5};
	  
	  int idir;
	  for(idir=0;idir<6;idir++){
	    interpminmod(&U0,Ui+idir,D,D+1,D+2,ix[idir],iy[idir],iz[idir]); // Up contains the interpolation
	  }



	  // READY TO EVOLVE EXTRAPOLATED VALUE

	  REAL FL[5],FR[5];
	  REAL GL[5],GR[5];
	  REAL HL[5],HR[5];
	  
	  getflux_X(Ui+0,FL);
	  getflux_X(Ui+1,FR);

	  getflux_Y(Ui+2,GL);
	  getflux_Y(Ui+3,GR);

	  getflux_Z(Ui+4,HL);
	  getflux_Z(Ui+5,HR);

	  for(idir=0;idir<6;idir++){
	    Ui[idir].d +=((FL[0]-FR[0])+(GL[0]-GR[0])+(HL[0]-HR[0]))*0.5*dt/dx;
	    Ui[idir].du+=((FL[1]-FR[1])+(GL[1]-GR[1])+(HL[1]-HR[1]))*0.5*dt/dx;
	    Ui[idir].dv+=((FL[2]-FR[2])+(GL[2]-GR[2])+(HL[2]-HR[2]))*0.5*dt/dx;
	    Ui[idir].dw+=((FL[3]-FR[3])+(GL[3]-GR[3])+(HL[3]-HR[3]))*0.5*dt/dx;
	    Ui[idir].E +=((FL[4]-FR[4])+(GL[4]-GR[4])+(HL[4]-HR[4]))*0.5*dt/dx;

	    REAL p;
	    p=Ui[idir].E-0.5*(Ui[idir].du*Ui[idir].du+Ui[idir].dv*Ui[idir].dv+Ui[idir].dw*Ui[idir].dw)/Ui[idir].d;

	    if(p<0){
	      //printf("ouch !\n");
	      //abort();
	    }
	  }

	  
}

// ===============================================================================================


void  matrix_jacobian(struct Wtype *W0, REAL dt,REAL dx,struct Wtype *Dx,struct Wtype *Dy,struct Wtype *Dz, struct Wtype *Wt){


  REAL M[25];
  REAL W[5]={0.,0.,0.,0.,0.};
  REAL d[5];
  int i,j;


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
  
  // ==== Final correction
  for(i=0;i<5;i++){
    W[i]*=(-dt/dx*0.5);
  }
  
  Wt->d=W[0];
  Wt->u=W[1];
  Wt->v=W[2];
  Wt->w=W[3];
  Wt->p=W[4];
  
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
	      Wm=&(stencil->oct[ioct-(int)pow(3,dir)].cell[vcell[inei2]].field);
	    }

	    inei2=2*dir+1;
	    if(vnei[inei2]==6){
	      Wp=&(stencil->oct[ioct].cell[vcell[inei2]].field);
	    }
	    else{
	      Wp=&(stencil->oct[ioct+(int)pow(3,dir)].cell[vcell[inei2]].field);
	    }

	    if(Wp->p<0) abort();
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
	    Wi[idir].p = W0->p+ix[idir]*D[0].p+iy[idir]*D[1].p+iz[idir]*D[2].p+Wt.p;

	    if(Wi[idir].p<0) abort();
	    if(Wi[idir].d<0) abort();


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
	    
	    if(Wi[idir].p<0) abort();
	    Wi[idir].a=sqrt(GAMMA*Wi[idir].p/Wi[idir].d);
	    //    if((idir==2)&&(Wi[idir].d>1.)) abort();
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
  WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);
  
  WRloc.d=WR->d;
  WRloc.u=WR->u;
  WRloc.p=WR->p;
  WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

  (*pstar)=(REAL)findPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=(REAL) findPressure(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) abort();

  qL=(*pstar<=WL->p?1.:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;
  if((*SL)>(*SR)){
    (*SL)=fminf(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=fmaxf(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
    //abort();
  }
  if((*SL)>(*SR)) abort();
  if(isnan(*ustar)) abort();
}

void speedestimateY_HLLC(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR, REAL *pstar, REAL *ustar){

  REAL qL,qR;
  struct Wtype1D WLloc;
  struct Wtype1D WRloc;
  int n;

  WLloc.d=WL->d;
  WLloc.u=WL->v;
  WLloc.p=WL->p;
  WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);
  
  WRloc.d=WR->d;
  WRloc.u=WR->v;
  WRloc.p=WR->p;
  WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

  (*pstar)=(REAL)findPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=(REAL) findPressure(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) abort();

  qL=(*pstar<=WL->p?1.:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;

  if((*SL)>(*SR)){
    (*SL)=fminf(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=fmaxf(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
    //abort();
  }
  if((*SL)>(*SR)) abort();

}




void speedestimateZ_HLLC(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR, REAL *pstar, REAL *ustar){

  REAL qL,qR;
  struct Wtype1D WLloc;
  struct Wtype1D WRloc;
  int n;

  WLloc.d=WL->d;
  WLloc.u=WL->w;
  WLloc.p=WL->p;
  WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);
  
  WRloc.d=WR->d;
  WRloc.u=WR->w;
  WRloc.p=WR->p;
  WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

  (*pstar)=(REAL)findPressure_Hybrid(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) (*pstar)=(REAL) findPressure(&WLloc,&WRloc,&n,ustar);
  if((*pstar)<0) abort();

  qL=(*pstar<=WL->p?1.:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WL->p-1.)));
  qR=(*pstar<=WR->p?1.:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((*pstar)/WR->p-1.)));
  
  *SL=WLloc.u-WLloc.a*qL;
  *SR=WRloc.u+WRloc.a*qR;
  if((*SL)>(*SR)){
    (*SL)=fminf(WLloc.u-WLloc.a,WRloc.u-WRloc.a);
    (*SR)=fmaxf(WLloc.u+WLloc.a,WRloc.u+WRloc.a);
    //abort();
  }
  if((*SL)>(*SR)) abort();

}


void speedestimateX(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR){//, REAL *pstar, REAL *ustar){

  REAL splus;
  splus=fmax(fabs(WL->u)+WL->a,fabs(WR->u)+WR->a);
  *SL=-splus;
  *SR= splus;


  /* REAL qL,qR; */
  /* REAL pstar,ustar; */
  
  /* pstar=0.5*(WL->p+WR->p)-0.5*(WR->u-WL->u)*(0.25*(WL->d+WR->d)*(WL->a+WR->a)); */
  /* ustar=0.5*(WL->u+WR->u)-0.5*(WR->p-WL->p)/(0.25*(WL->d+WR->d)*(WL->a+WR->a)); */

  /* qL=(pstar<=WL->p?1.:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((pstar)/WL->p-1.))); */
  /* qR=(pstar<=WR->p?1.:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((pstar)/WR->p-1.))); */
  
  /* *SL=WL->u-WL->a*qL; */
  /* *SR=WR->u+WR->a*qR; */
  
  /* if(*SL>*SR) abort(); */

}

/* void speedestimateY(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR){//, REAL *pstar, REAL *ustar){ */

/*   REAL splus; */
/*   /\* splus=fmax(fabs(WL->u)+WL->a,fabs(WR->u)+WR->a); *\/ */
/*   /\* *SL=-splus; *\/ */
/*   /\* *SR= splus; *\/ */


/*   REAL qL,qR; */
/*   REAL pstar,ustar; */
  
/*   pstar=0.5*(WL->p+WR->p)-0.5*(WR->v-WL->v)*(0.25*(WL->d+WR->d)*(WL->a+WR->a)); */
/*   ustar=0.5*(WL->v+WR->v)-0.5*(WR->p-WL->p)/(0.25*(WL->d+WR->d)*(WL->a+WR->a)); */

/*   qL=((pstar)<=WL->p?1:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((pstar)/WL->p-1.))); */
/*   qR=((pstar)<=WR->p?1:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((pstar)/WR->p-1.))); */
  
/*   *SL=WL->v-WL->a*qL; */
/*   *SR=WR->v+WR->a*qR; */
  
/*   if(*SL>*SR) abort(); */

/* } */

/* void speedestimateZ(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR){// REAL *pstar, REAL *ustar){ */

/*   REAL splus; */
/*   /\* splus=fmax(fabs(WL->u)+WL->a,fabs(WR->u)+WR->a); *\/ */
/*   /\* *SL=-splus; *\/ */
/*   /\* *SR= splus; *\/ */


/* REAL qL,qR; */
/* REAL pstar,ustar; */
  
/*   pstar=0.5*(WL->p+WR->p)-0.5*(WR->w-WL->w)*(0.25*(WL->d+WR->d)*(WL->a+WR->a)); */
/*   ustar=0.5*(WL->w+WR->w)-0.5*(WR->p-WL->p)/(0.25*(WL->d+WR->d)*(WL->a+WR->a)); */

/*   qL=((pstar)<=WL->p?1:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((pstar)/WL->p-1.))); */
/*   qR=((pstar)<=WR->p?1:sqrt(1.+(GAMMA+1.)/(2.*GAMMA)*((pstar)/WR->p-1.))); */
  
/*   *SL=WL->w-WL->a*qL; */
/*   *SR=WR->w+WR->a*qR; */
  
/*   if(*SL>*SR) abort(); */

/* } */

void W2Ustar_x(struct Wtype *W, struct Utype *U, REAL SK, REAL SS){
  
  REAL fact=W->d*(SK-W->u)/(SK-SS);
  REAL EK=0.5*W->d*(W->u*W->u+W->v*W->v+W->w*W->w)+W->p/(GAMMA-1.);

  U->d=fact;
  U->du=fact*SS;
  U->dv=fact*W->v;
  U->dw=fact*W->w;
  U->E=fact*(EK/W->d+(SS-W->u)*(SS+W->p/(W->d*(SK-W->u))));

}



void speedestimateY(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR){

  REAL splus;
  splus=fmax(fabs(WL->v)+WL->a,fabs(WR->v)+WR->a);
  *SL=-splus;
  *SR= splus;



}

void speedestimateZ(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR){

  REAL splus;
  splus=fmax(fabs(WL->w)+WL->a,fabs(WR->w)+WR->a);
  *SL=-splus;
  *SR= splus;
}



int hydroM(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt){

  int inei,icell,icellcoarse;
  int i;
  REAL temp;
  REAL res,restot=0.;
  int vnei[6],vcell[6];

  REAL FL[NVAR],FR[NVAR];
  REAL GL[NVAR],GR[NVAR];
  REAL HL[NVAR],HR[NVAR];


  REAL Smax;
  struct Wtype1D_double WRloc, WLloc;
  struct Utype Uold,Unew;
  struct Wtype Wold,Wnew;
  int idxL,idxR;
  REAL pstar,ustar;
  int n;
  struct Wtype1D Wtest;
  struct Wtype   Wtest3D;
  struct Utype   Utest3D;
  struct Utype UC[6];
  struct Utype UT[6];
  struct Wtype WT[6];
  struct Utype UN[6];
  struct Wtype WN[6];
  struct Wtype WC[6];
  int ioct[7]={12,14,10,16,4,22,13};
  int idxnei[6]={1,0,3,2,5,4};

  struct Wtype *curcell;
  struct Wtype *neicell;

  REAL SL,SR;
#ifdef RIEMANN_HLL
  
  REAL Fr[NVAR],Fl[NVAR];
#endif

  int ffact[6]={0,0,0,0,0,0};


  int tagr=0,tagl=0;
  REAL fact;

  //printf("let's do some hydro\n");
  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
      if(idxnei[0]!=1) abort();
      
    for(i=0;i<nread;i++){ // we scan the octs
      tagr=0;
      tagl=0;
      
      if(idxnei[0]!=1) abort();
      memset(FL,0,sizeof(REAL)*NVAR);
      memset(FR,0,sizeof(REAL)*NVAR);
      memset(HL,0,sizeof(REAL)*NVAR);
      memset(HR,0,sizeof(REAL)*NVAR);
      memset(GL,0,sizeof(REAL)*NVAR);
      memset(GR,0,sizeof(REAL)*NVAR);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].field);
      
      Wold.d=curcell->d;
      Wold.u=curcell->u;;
      Wold.v=curcell->v;;
      Wold.w=curcell->w;;
      Wold.p=curcell->p;;
      Wold.a=sqrt(GAMMA*Wold.p/Wold.d);

      W2U(&Wold,&Uold); // primitive -> conservative

      /* // MUSCL STATE RECONSTRUCTION */
#ifdef NOFLUX
      memset(ffact,0,sizeof(int)*6);
#endif

      MUSCL_BOUND2(stencil+i, 13, icell, WC,dt,dx);// central
      for(inei=0;inei<6;inei++){
	MUSCL_BOUND2(stencil+i, ioct[vnei[inei]], vcell[inei], WT,dt,dx);// 
	memcpy(WN+inei,WT+idxnei[inei],sizeof(struct Wtype)); 
	

       	W2U(WN+inei,UN+inei);
      	W2U(WC+inei,UC+inei);
	

#ifdef NOFLUX
	if(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].child!=NULL){
	  ffact[inei]=1; // we cancel the contriubtion of split neighbors
	}
#endif

      }

#ifdef DUAL_E
	struct Utype Us;
	REAL ebar;
	REAL ecen=0.;
	REAL divu[6]={0.,0.,0.,0.,0.,0.};
#endif


      // X DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems LEFT
      // Switching to Split description
      if(!ffact[0]){      

/* #ifdef RIEMANN_EXACT */

/* 	WLloc.d=WN[0].d; */
/* 	WLloc.u=WN[0].u; */
/* 	WLloc.p=WN[0].p; */
/* 	WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d); */

/* 	WRloc.d=WC[0].d; */
/* 	WRloc.u=WC[0].u; */
/* 	WRloc.p=WC[0].p; */
/* 	WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d); */

/* 	// Riemann Solver */
/* 	if((WLloc.p<0)||(WRloc.p<0)){ */
/* 	  printf("hajzehr\n"); */
/* 	  abort(); */
/* 	} */

/* 	pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar); */
/* 	getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
/* 	Wtest3D.d=Wtest.d; */
/* 	Wtest3D.u=Wtest.u; */
/* 	Wtest3D.p=Wtest.p; */
/* 	Wtest3D.a=Wtest.a; */
      
/* 	if(isnan(Wtest.p)){ */
/* 	  printf("ouhla NAN\n"); */
/* 	  abort(); */
/* 	} */

/* 	// Passive advection */
/* 	if(ustar>0.) */
/* 	  { */
/* 	    Wtest3D.v=WN[0].v; */
/* 	    Wtest3D.w=WN[0].w; */
/* 	  } */
/* 	else */
/* 	  { */
/* 	    Wtest3D.v=WC[0].v; */
/* 	    Wtest3D.w=WC[0].w; */
/* 	  } */
      
/* 	W2U(&Wtest3D,&Utest3D); */
      
/* 	// Getting the fluxes LEFT */
/* 	getflux_X(&Utest3D,FL); */

/* #ifdef DUAL_E */
/* 	FL[5]=Wtest3D.p*Wtest3D.u; */
/* #endif */

/* #endif */

	// =============================================


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
	  getflux_X(&UN[0],FL);
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
	  Us.E =(fact*(UC[0].E/UC[0].d+(ustar-WC[0].u)*(ustar+WC[0].p/(WC[0].d*(SL-WC[0].u)))));
#endif


	}

#ifdef DUAL_E
	ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d);
	ecen+=ebar;
	FL[5]=(Us.du/Us.d*ebar);
	divu[0]=(GAMMA-1.)*(Us.du/Us.d);
#endif

#endif




/* 	// ============================================= */
/* #ifdef RIEMANN_HLL */
/* 	speedestimateX(&WN[0],&WC[0],&SL,&SR);//,&pstar,&ustar); */

/* 	if((SL<0.)&&(SR>0)){ */
/* 	  getflux_X(&UN[0],Fl); */
/* 	  getflux_X(&UC[0],Fr); */

/* 	  FL[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UC[0].d -UN[0].d ))/(SR-SL); */
/* 	  FL[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UC[0].du-UN[0].du))/(SR-SL); */
/* 	  FL[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UC[0].dv-UN[0].dv))/(SR-SL); */
/* 	  FL[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UC[0].dw-UN[0].dw))/(SR-SL); */
/* 	  FL[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UC[0].E -UN[0].E ))/(SR-SL); */

/* #ifdef DUAL_E */

/* 	  //FL[5]=pstar*ustar; */

/* 	  Fl[5]=WN[0].p*WN[0].u; */
/* 	  Fr[5]=WC[0].p*WC[0].u; */
/* 	  FL[5]=(SR*Fl[5]-SL*Fr[5]+SL*SR*(WC[0].p -WN[0].p ))/(SR-SL); */
/* #endif */

/* 	  tagl=1; */
/* 	} */
/* 	else{ */
/* 	  if(SL>=0.){ */
/* 	    getflux_X(&UN[0],FL); */
/* #ifdef DUAL_E */
/* 	    FL[5]=WN[0].p*WN[0].u; */
/* #endif */
/* 	  } */
/* 	  else if(SR<=0.){ */
/* 	    getflux_X(&UC[0],FL); */
/* #ifdef DUAL_E */
/* 	    FL[5]=WC[0].p*WC[0].u; */
/* #endif */
/* 	  } */
/* 	} */
      
/* #ifdef DUAL_E */
/* 	// divergence contribution */
/* 	FL[5]+=0.5*Wold.p*(GAMMA-1.)/dt*stencil[i].oct[ioct[vnei[0]]].cell[vcell[0]].field.u; */
/* #endif */
      
/* #endif */
      
      } 

      // --------- solving the Riemann Problems RIGHT

      if(!ffact[1]){
	// Switching to Split description

/* #ifdef RIEMANN_EXACT */
/* 	WRloc.d=WN[1].d; */
/* 	WRloc.u=WN[1].u; */
/* 	WRloc.p=WN[1].p; */
/* 	WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d); */

/* 	WLloc.d=WC[1].d; */
/* 	WLloc.u=WC[1].u; */
/* 	WLloc.p=WC[1].p; */
/* 	WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d); */

/* 	if((WLloc.p<0)||(WRloc.p<0)){ */
/* 	  printf("hajzehr\n"); */
/* 	  abort(); */
/* 	} */

/* 	// Riemann Solver */
/* 	pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar); */
/* 	getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
/* 	Wtest3D.d=Wtest.d; */
/* 	Wtest3D.u=Wtest.u; */
/* 	Wtest3D.p=Wtest.p; */
/* 	Wtest3D.a=Wtest.a; */
      
/* 	// Passive advection */
/* 	if(ustar<0.) */
/* 	  { */
/* 	    Wtest3D.v=WN[1].v; */
/* 	    Wtest3D.w=WN[1].w; */
/* 	  } */
/* 	else */
/* 	  { */
/* 	    Wtest3D.v=WC[1].v; */
/* 	    Wtest3D.w=WC[1].w; */
/* 	  } */
      
/* 	W2U(&Wtest3D,&Utest3D); */
      
/* 	// Getting the fluxes RIGHT */
/* 	getflux_X(&Utest3D,FR); */
/* #ifdef DUAL_E */
/* 	FR[5]=Wtest3D.p*Wtest3D.u; */
/* #endif */

/* #endif */
	// =======================================================================

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
	}

#ifdef DUAL_E
	ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d);
	ecen+=ebar;
	FR[5]=(Us.du/Us.d*ebar);
	divu[1]=(GAMMA-1.)*(Us.du/Us.d);
#endif

#endif




/* 	// ======================================================================= */
/* #ifdef RIEMANN_HLL */
/* 	speedestimateX(&WC[1],&WN[1],&SL,&SR);//,&pstar,&ustar); */
/* 	/\* SL=WC[1].u-WC[1].a; *\/ */
/* 	/\* SR=WN[1].u+WN[1].a; *\/ */

/* 	if((SL<0.)&&(SR>0)){ */
/* 	  getflux_X(&UC[1],Fl); */
/* 	  getflux_X(&UN[1],Fr); */
	
/* 	  FR[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UN[1].d -UC[1].d ))/(SR-SL); */
/* 	  FR[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UN[1].du-UC[1].du))/(SR-SL); */
/* 	  FR[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UN[1].dv-UC[1].dv))/(SR-SL); */
/* 	  FR[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UN[1].dw-UC[1].dw))/(SR-SL); */
/* 	  FR[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UN[1].E -UC[1].E ))/(SR-SL); */

/* #ifdef DUAL_E */

/* 	  //FR[5]=pstar*ustar; */
/* 	  Fl[5]=WC[1].p*WC[1].u; */
/* 	  Fr[5]=WN[1].p*WN[1].u; */
/* 	  FR[5]=(SR*Fl[5]-SL*Fr[5]+SL*SR*(WN[1].p -WC[1].p))/(SR-SL); */
/* #endif */
/* 	} */
/* 	else{ */
/* 	  if(SL>=0.){ */
/* 	    getflux_X(&UC[1],FR); */
/* #ifdef DUAL_E */
/* 	    FR[5]=WC[1].p*WC[1].u; */
/* #endif */
/* 	  } */
/* 	  else if(SR<=0.){ */
/* 	    getflux_X(&UN[1],FR); */
/* #ifdef DUAL_E */
/* 	    FR[5]=WN[1].p*WN[1].u; */
/* #endif */
/* 	  } */
/* 	} */
/* 	if(isnan(FR[5])) abort(); */
      
/* #ifdef DUAL_E */
/* 	// we add the divergence component */
/* 	FR[5]+=0.5*Wold.p*(GAMMA-1.)/dt*stencil[i].oct[ioct[vnei[1]]].cell[vcell[1]].field.u; */
/* #endif */

      
/* #endif */

       } 

       // Y DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems FRONT

     
      // Switching to Split description

      if(!ffact[2]){

/* #ifdef RIEMANN_EXACT */
/* 	WLloc.d=WN[2].d; */
/* 	WLloc.u=WN[2].v; */
/* 	WLloc.p=WN[2].p; */
/* 	WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d); */

/* 	WRloc.d=WC[2].d; */
/* 	WRloc.u=WC[2].v; */
/* 	WRloc.p=WC[2].p; */
/* 	WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d); */

/* 	// Riemann Solver */
/* 	if((WLloc.p<0)||(WRloc.p<0)){ */
/* 	  printf("hajzehr\n"); */
/* 	  abort(); */
/* 	} */
/* 	pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar); */

/* 	getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
/* 	Wtest3D.d=Wtest.d; */
/* 	Wtest3D.v=Wtest.u; */
/* 	Wtest3D.p=Wtest.p; */
/* 	Wtest3D.a=Wtest.a; */
      
/* 	// Passive advection */
/* 	if(ustar>0.) */
/* 	  { */
/* 	    Wtest3D.u=WN[2].u; */
/* 	    Wtest3D.w=WN[2].w; */
/* 	  } */
/* 	else */
/* 	  { */
/* 	    Wtest3D.u=WC[2].u; */
/* 	    Wtest3D.w=WC[2].w; */
/* 	  } */
      
/* 	W2U(&Wtest3D,&Utest3D); */
      
/* 	// Getting the fluxes LEFT */
/* 	getflux_Y(&Utest3D,GL); */

/* #ifdef DUAL_E */
/* 	GL[5]=Wtest3D.p*Wtest3D.v; */
/* #endif */

/* #endif */


	// =============================================

#ifdef RIEMANN_HLLC
	speedestimateY_HLLC(&WN[2],&WC[2],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  getflux_Y(&UN[2],GL);
#ifdef DUAL_E
	  memcpy(&Us,&UN[2],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  getflux_Y(&UC[2],GL);
#ifdef DUAL_E
	  memcpy(&Us,&UC[2],sizeof(struct Utype));
#endif

	}
	else if((SL<0.)&&(ustar>=0.)){
	  getflux_Y(&UN[2],GL);
	  fact=WN[2].d*(SL-WN[2].v)/(SL-ustar);
	  GL[0]+=(fact*1.                                                                      -UN[2].d )*SL;
	  GL[1]+=(fact*WN[2].u                                                                 -UN[2].du)*SL;
	  GL[2]+=(fact*ustar                                                                   -UN[2].dv)*SL;
	  GL[3]+=(fact*WN[2].w                                                                 -UN[2].dw)*SL;
	  GL[4]+=(fact*(UN[2].E/UN[2].d+(ustar-WN[2].v)*(ustar+WN[2].p/(WN[2].d*(SL-WN[2].v))))-UN[2].E )*SL;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WN[2].u);
	  Us.dv=(fact*ustar);
	  Us.dw=(fact*WN[2].w);
	  Us.E =(fact*(UN[2].E/UN[2].d+(ustar-WN[2].v)*(ustar+WN[2].p/(WN[2].d*(SL-WN[2].v)))));
#endif

	}
	else if((ustar<=0.)&&(SR>0.)){
	  getflux_Y(&UC[2],GL);
	  fact=WC[2].d*(SR-WC[2].v)/(SR-ustar);
	  GL[0]+=(fact*1.                                                                      -UC[2].d )*SR;
	  GL[1]+=(fact*WC[2].u                                                                 -UC[2].du)*SR;
	  GL[2]+=(fact*ustar                                                                   -UC[2].dv)*SR;
	  GL[3]+=(fact*WC[2].w                                                                 -UC[2].dw)*SR;
	  GL[4]+=(fact*(UC[2].E/UC[2].d+(ustar-WC[2].v)*(ustar+WC[2].p/(WC[2].d*(SR-WC[2].v))))-UC[2].E )*SR;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WC[2].u);
	  Us.dv=(fact*ustar);
	  Us.dw=(fact*WC[2].w);
	  Us.E =(fact*(UC[2].E/UC[2].d+(ustar-WC[2].v)*(ustar+WC[2].p/(WC[2].d*(SR-WC[2].v)))));
#endif

	}

#ifdef DUAL_E
	ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d);
	ecen+=ebar;
	GL[5]=(Us.dv/Us.d*ebar);
	divu[2]=(GAMMA-1.)*(Us.dv/Us.d);
#endif

#endif

	// =============================================


/* #ifdef RIEMANN_HLL */
/* 	speedestimateY(&WN[2],&WC[2],&SL,&SR);//&pstar,&ustar); */
/* 	/\* SL=WN[2].v-WN[2].a; *\/ */
/* 	/\* SR=WC[2].v+WC[2].a; *\/ */

/* 	if((SL<0.)&&(SR>0)){ */
/* 	  getflux_Y(&UN[2],Fl); */
/* 	  getflux_Y(&UC[2],Fr); */
	
/* 	  GL[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UC[2].d -UN[2].d ))/(SR-SL); */
/* 	  GL[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UC[2].du-UN[2].du))/(SR-SL); */
/* 	  GL[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UC[2].dv-UN[2].dv))/(SR-SL); */
/* 	  GL[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UC[2].dw-UN[2].dw))/(SR-SL); */
/* 	  GL[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UC[2].E -UN[2].E ))/(SR-SL); */

/* #ifdef DUAL_E */
/* 	  Fl[5]=WN[2].p*WN[2].v; */
/* 	  Fr[5]=WC[2].p*WC[2].v; */
/* 	  GL[5]=(SR*Fl[5]-SL*Fr[5]+SL*SR*(WC[2].p -WN[2].p ))/(SR-SL); */
/* #endif */

/* 	} */
/* 	else{ */
/* 	  if(SL>=0.){ */
/* 	    getflux_Y(&UN[2],GL); */

/* #ifdef DUAL_E */
/* 	    GL[5]=WN[2].p*WN[2].v; */
/* #endif */

/* 	  } */
/* 	  else if(SR<=0.){ */
/* 	    getflux_Y(&UC[2],GL); */
/* #ifdef DUAL_E */
/* 	    GL[5]=WC[2].p*WC[2].v; */
/* #endif */

/* 	  } */
/* 	} */

/* #ifdef DUAL_E */
/* 	// divergence contribution */
/* 	GL[5]+=0.5*Wold.p*(GAMMA-1.)/dt*stencil[i].oct[ioct[vnei[2]]].cell[vcell[2]].field.v; */
/* #endif */
      
/* #endif */
      } 

      // --------- solving the Riemann Problems BACK


      // Switching to Split description
      if(!ffact[3]){

/* #ifdef RIEMANN_EXACT */
/* 	WRloc.d=WN[3].d; */
/* 	WRloc.u=WN[3].v; */
/* 	WRloc.p=WN[3].p; */
/* 	WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d); */

/* 	WLloc.d=WC[3].d; */
/* 	WLloc.u=WC[3].v; */
/* 	WLloc.p=WC[3].p; */
/* 	WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d); */

/* 	if((WLloc.p<0)||(WRloc.p<0)){ */
/* 	  printf("hajzehr\n"); */
/* 	  abort(); */
/* 	} */
/* 	// Riemann Solver */
/* 	pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar); */
/* 	getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
/* 	Wtest3D.d=Wtest.d; */
/* 	Wtest3D.v=Wtest.u; */
/* 	Wtest3D.p=Wtest.p; */
/* 	Wtest3D.a=Wtest.a; */
      
/* 	// Passive advection */
/* 	if(ustar<0.) */
/* 	  { */
/* 	    Wtest3D.u=WN[3].u; */
/* 	    Wtest3D.w=WN[3].w; */
/* 	  } */
/* 	else */
/* 	  { */
/* 	    Wtest3D.u=WC[3].u; */
/* 	    Wtest3D.w=WC[3].w; */
/* 	  } */
      
/* 	W2U(&Wtest3D,&Utest3D); */
      
/* 	// Getting the fluxes RIGHT */
/* 	getflux_Y(&Utest3D,GR); */
/* #ifdef DUAL_E */
/* 	GR[5]=Wtest3D.p*Wtest3D.v; */
/* #endif */

/* #endif */
/* 	// ============================================= */

#ifdef RIEMANN_HLLC
	speedestimateY_HLLC(&WC[3],&WN[3],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  getflux_Y(&UC[3],GR);
#ifdef DUAL_E
	  memcpy(&Us,&UC[3],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  getflux_Y(&UN[3],GR);
#ifdef DUAL_E
	  memcpy(&Us,&UN[3],sizeof(struct Utype));
#endif
	}
	else if((SL<0.)&&(ustar>=0.)){
	  getflux_Y(&UC[3],GR);
	  fact=WC[3].d*(SL-WC[3].v)/(SL-ustar);
	  GR[0]+=(fact*1.                                                                      -UC[3].d )*SL;
	  GR[1]+=(fact*WC[3].u                                                                 -UC[3].du)*SL;
	  GR[2]+=(fact*ustar                                                                   -UC[3].dv)*SL;
	  GR[3]+=(fact*WC[3].w                                                                 -UC[3].dw)*SL;
	  GR[4]+=(fact*(UC[3].E/UC[3].d+(ustar-WC[3].v)*(ustar+WC[3].p/(WC[3].d*(SL-WC[3].v))))-UC[3].E )*SL;
#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WC[3].u);
	  Us.dv=(fact*ustar);
	  Us.dw=(fact*WC[3].w);
	  Us.E =(fact*(UC[3].E/UC[3].d+(ustar-WC[3].v)*(ustar+WC[3].p/(WC[3].d*(SL-WC[3].v)))));
#endif

	}
	else if((ustar<=0.)&&(SR>0.)){
	  getflux_Y(&UN[3],GR);
	  fact=WN[3].d*(SR-WN[3].v)/(SR-ustar);
	  GR[0]+=(fact*1.                                                                      -UN[3].d )*SR;
	  GR[1]+=(fact*WN[3].u                                                                 -UN[3].du)*SR;
	  GR[2]+=(fact*ustar                                                                   -UN[3].dv)*SR;
	  GR[3]+=(fact*WN[3].w                                                                 -UN[3].dw)*SR;
	  GR[4]+=(fact*(UN[3].E/UN[3].d+(ustar-WN[3].v)*(ustar+WN[3].p/(WN[3].d*(SR-WN[3].v))))-UN[3].E )*SR;
#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WN[3].u);
	  Us.dv=(fact*ustar);
	  Us.dw=(fact*WN[3].w);
	  Us.E =(fact*(UN[3].E/UN[3].d+(ustar-WN[3].v)*(ustar+WN[3].p/(WN[3].d*(SR-WN[3].v)))));
#endif
	}

#ifdef DUAL_E
	ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d);
	ecen+=ebar;
	GR[5]=(Us.dv/Us.d*ebar);
	divu[3]=(GAMMA-1.)*(Us.dv/Us.d);
#endif

#endif
	// =============================================

/* #ifdef RIEMANN_HLL */
/* 	speedestimateY(&WC[3],&WN[3],&SL,&SR);//,&pstar,&ustar); */
/* 	/\* SL=WC[3].v-WC[3].a; *\/ */
/* 	/\* SR=WN[3].v+WN[3].a; *\/ */

/* 	if((SL<0.)&&(SR>0)){ */
/* 	  getflux_Y(&UC[3],Fl); */
/* 	  getflux_Y(&UN[3],Fr); */
	
/* 	  GR[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UN[3].d -UC[3].d ))/(SR-SL); */
/* 	  GR[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UN[3].du-UC[3].du))/(SR-SL); */
/* 	  GR[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UN[3].dv-UC[3].dv))/(SR-SL); */
/* 	  GR[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UN[3].dw-UC[3].dw))/(SR-SL); */
/* 	  GR[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UN[3].E -UC[3].E ))/(SR-SL); */

/* #ifdef DUAL_E */
/* 	  Fl[5]=WC[3].p*WC[3].v; */
/* 	  Fr[5]=WN[3].p*WN[3].v; */
/* 	  GR[5]=(SR*Fl[5]-SL*Fr[5]+SL*SR*(WN[3].p -WC[3].p))/(SR-SL); */
/* #endif */


/* 	} */
/* 	else{ */
/* 	  if(SL>=0.){ */
/* 	    getflux_Y(&UC[3],GR); */

/* #ifdef DUAL_E */
/* 	    GR[5]=WC[3].p*WC[3].v; */
/* #endif */

/* 	  } */
/* 	  else if(SR<=0.){ */
/* 	    getflux_Y(&UN[3],GR); */
/* #ifdef DUAL_E */
/* 	    GR[5]=WN[3].p*WN[3].v; */
/* #endif */
/* 	  } */
/* 	} */

/* 	// we add the divergence component */
/* 	GR[5]+=0.5*Wold.p*(GAMMA-1.)/dt*stencil[i].oct[ioct[vnei[3]]].cell[vcell[3]].field.v; */
/* #endif */

      }

      // Z DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems BOTTOM

      // Switching to Split description

      if(!ffact[4]){

/* #ifdef RIEMANN_EXACT */
/* 	WLloc.d=WN[4].d; */
/* 	WLloc.u=WN[4].w; */
/* 	WLloc.p=WN[4].p; */
/* 	WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d); */

/* 	WRloc.d=WC[4].d; */
/* 	WRloc.u=WC[4].w; */
/* 	WRloc.p=WC[4].p; */
/* 	WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d); */

/* 	if((WLloc.p<0)||(WRloc.p<0)){ */
/* 	  printf("hajzehr\n"); */
/* 	  abort(); */
/* 	} */
/* 	// Riemann Solver */
/* 	pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar); */

/* 	getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
/* 	Wtest3D.d=Wtest.d; */
/* 	Wtest3D.w=Wtest.u; */
/* 	Wtest3D.p=Wtest.p; */
/* 	Wtest3D.a=Wtest.a; */
      
/* 	// Passive advection */
/* 	if(ustar>0.) */
/* 	  { */
/* 	    Wtest3D.u=WN[4].u; */
/* 	    Wtest3D.v=WN[4].v; */
/* 	  } */
/* 	else */
/* 	  { */
/* 	    Wtest3D.u=WC[4].u; */
/* 	    Wtest3D.v=WC[4].v; */
/* 	  } */
      
/* 	W2U(&Wtest3D,&Utest3D); */
      
/* 	// Getting the fluxes LEFT */
/* 	getflux_Z(&Utest3D,HL); */

/* #ifdef DUAL_E */
/* 	HL[5]=Wtest3D.p*Wtest3D.w; */
/* #endif */

/* #endif */

/* 	// =========================================== */

#ifdef RIEMANN_HLLC
	speedestimateZ_HLLC(&WN[4],&WC[4],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  getflux_Z(&UN[4],HL);
#ifdef DUAL_E
	  memcpy(&Us,&UN[4],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  getflux_Z(&UC[4],HL);
#ifdef DUAL_E
	  memcpy(&Us,&UC[4],sizeof(struct Utype));
#endif
	}
	else if((SL<0.)&&(ustar>=0.)){
	  getflux_Z(&UN[4],HL);
	  fact=WN[4].d*(SL-WN[4].w)/(SL-ustar);
	  HL[0]+=(fact*1.                                                                      -UN[4].d )*SL;
	  HL[1]+=(fact*WN[4].u                                                                 -UN[4].du)*SL;
	  HL[2]+=(fact*WN[4].v                                                                 -UN[4].dv)*SL;
	  HL[3]+=(fact*ustar                                                                   -UN[4].dw)*SL;
	  HL[4]+=(fact*(UN[4].E/UN[4].d+(ustar-WN[4].w)*(ustar+WN[4].p/(WN[4].d*(SL-WN[4].w))))-UN[4].E )*SL;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WN[4].u);
	  Us.dv=(fact*WN[4].v);
	  Us.dw=(fact*ustar);
	  Us.E =(fact*(UN[4].E/UN[4].d+(ustar-WN[4].w)*(ustar+WN[4].p/(WN[4].d*(SL-WN[4].w)))));
#endif

	}
	else if((ustar<=0.)&&(SR>0.)){
	  getflux_Z(&UC[4],HL);
	  fact=WC[4].d*(SR-WC[4].w)/(SR-ustar);
	  HL[0]+=(fact*1.                                                                      -UC[4].d )*SR;
	  HL[1]+=(fact*WC[4].u                                                                 -UC[4].du)*SR;
	  HL[2]+=(fact*WC[4].v                                                                 -UC[4].dv)*SR;
	  HL[3]+=(fact*ustar                                                                   -UC[4].dw)*SR;
	  HL[4]+=(fact*(UC[4].E/UC[4].d+(ustar-WC[4].w)*(ustar+WC[4].p/(WC[4].d*(SR-WC[4].w))))-UC[4].E )*SR;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WC[4].u);
	  Us.dv=(fact*WC[4].v);
	  Us.dw=(fact*ustar);
	  Us.E =(fact*(UC[4].E/UC[4].d+(ustar-WC[4].w)*(ustar+WC[4].p/(WC[4].d*(SR-WC[4].w)))));
#endif

	}

#ifdef DUAL_E
	ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d);
	ecen+=ebar;
	HL[5]=(Us.dw/Us.d*ebar);
	divu[4]=(GAMMA-1.)*(Us.dw/Us.d);
#endif

#endif
	// ===========================================

/* #ifdef RIEMANN_HLL */
/* 	speedestimateZ(&WN[4],&WC[4],&SL,&SR);//,&pstar,&ustar); */
/* 	/\* SL=WN[4].w-WN[4].a; *\/ */
/* 	/\* SR=WC[4].w+WC[4].a; *\/ */

/* 	if((SL<0.)&&(SR>0.)){ */
/* 	  getflux_Z(&UN[4],Fl); */
/* 	  getflux_Z(&UC[4],Fr); */
	
/* 	  HL[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UC[4].d -UN[4].d ))/(SR-SL); */
/* 	  HL[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UC[4].du-UN[4].du))/(SR-SL); */
/* 	  HL[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UC[4].dv-UN[4].dv))/(SR-SL); */
/* 	  HL[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UC[4].dw-UN[4].dw))/(SR-SL); */
/* 	  HL[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UC[4].E -UN[4].E ))/(SR-SL); */

/* #ifdef DUAL_E */
/* 	  Fl[5]=WN[4].p*WN[4].w; */
/* 	  Fr[5]=WC[4].p*WC[4].w; */
/* 	  HL[5]=(SR*Fl[5]-SL*Fr[5]+SL*SR*(WC[4].p -WN[4].p ))/(SR-SL); */
/* #endif */

/* 	} */
/* 	else{ */
/* 	  if(SL>=0.){ */
/* 	    getflux_Z(&UN[4],HL); */
/* #ifdef DUAL_E */
/* 	    HL[5]=WN[4].p*WN[4].w; */
/* #endif */
/* 	  } */
/* 	  else if(SR<=0.){ */
/* 	    getflux_Z(&UC[4],HL); */
/* #ifdef DUAL_E */
/* 	    HL[5]=WC[4].p*WC[4].w; */
/* #endif */
/* 	  } */
/* 	} */
      

/* #ifdef DUAL_E */
/* 	// divergence contribution */
/* 	HL[5]+=0.5*Wold.p*(GAMMA-1.)/dt*stencil[i].oct[ioct[vnei[4]]].cell[vcell[4]].field.w; */
/* #endif */

/* #endif */

      }

      // --------- solving the Riemann Problems TOP


      // Switching to Split description

      if(!ffact[5]){

/* #ifdef RIEMANN_EXACT */
/* 	WRloc.d=WN[5].d; */
/* 	WRloc.u=WN[5].w; */
/* 	WRloc.p=WN[5].p; */
/* 	WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d); */

/* 	WLloc.d=WC[5].d; */
/* 	WLloc.u=WC[5].w; */
/* 	WLloc.p=WC[5].p; */
/* 	WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d); */

/* 	if((WLloc.p<0)||(WRloc.p<0)){ */
/* 	  printf("hajzehr\n"); */
/* 	  abort(); */
/* 	} */
/* 	// Riemann Solver */
/* 	pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar); */
/* 	getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
/* 	Wtest3D.d=Wtest.d; */
/* 	Wtest3D.w=Wtest.u; */
/* 	Wtest3D.p=Wtest.p; */
/* 	Wtest3D.a=Wtest.a; */
      
/* 	// Passive advection */
/* 	if(ustar<0.) */
/* 	  { */
/* 	    Wtest3D.u=WN[5].u; */
/* 	    Wtest3D.v=WN[5].v; */
/* 	  } */
/* 	else */
/* 	  { */
/* 	    Wtest3D.u=WC[5].u; */
/* 	    Wtest3D.v=WC[5].v; */
/* 	  } */
      
/* 	W2U(&Wtest3D,&Utest3D); */
      
/* 	// Getting the fluxes RIGHT */
/* 	getflux_Z(&Utest3D,HR); */

/* #ifdef DUAL_E */
/* 	HR[5]=Wtest3D.p*Wtest3D.w; */
/* #endif */

/* #endif */

	//=====================================================

#ifdef RIEMANN_HLLC
	speedestimateZ_HLLC(&WC[5],&WN[5],&SL,&SR,&pstar,&ustar);

	if(SL>=0.){
	  getflux_Z(&UC[5],HR);
#ifdef DUAL_E
	  memcpy(&Us,&UC[5],sizeof(struct Utype));
#endif

	}
	else if(SR<=0.){
	  getflux_Z(&UN[5],HR);
#ifdef DUAL_E
	  memcpy(&Us,&UN[5],sizeof(struct Utype));
#endif

	}
	else if((SL<0.)&&(ustar>=0.)){
	  getflux_Z(&UC[5],HR);
	  fact=WC[5].d*(SL-WC[5].w)/(SL-ustar);
	  HR[0]+=(fact*1.                                                                      -UC[5].d )*SL;
	  HR[1]+=(fact*WC[5].u                                                                 -UC[5].du)*SL;
	  HR[2]+=(fact*WC[5].v                                                                 -UC[5].dv)*SL;
	  HR[3]+=(fact*ustar                                                                   -UC[5].dw)*SL;
	  HR[4]+=(fact*(UC[5].E/UC[5].d+(ustar-WC[5].w)*(ustar+WC[5].p/(WC[5].d*(SL-WC[5].w))))-UC[5].E )*SL;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WC[5].u);
	  Us.dv=(fact*WC[5].v);
	  Us.dw=(fact*ustar);
	  Us.E =(fact*(UC[5].E/UC[5].d+(ustar-WC[5].w)*(ustar+WC[5].p/(WC[5].d*(SL-WC[5].w)))));
#endif

	}
	else if((ustar<=0.)&&(SR>0.)){
	  getflux_Z(&UN[5],HR);
	  fact=WN[5].d*(SR-WN[5].w)/(SR-ustar);
	  HR[0]+=(fact*1.                                                                      -UN[5].d )*SR;
	  HR[1]+=(fact*WN[5].u                                                                 -UN[5].du)*SR;
	  HR[2]+=(fact*WN[5].v                                                                 -UN[5].dv)*SR;
	  HR[3]+=(fact*ustar                                                                   -UN[5].dw)*SR;
	  HR[4]+=(fact*(UN[5].E/UN[5].d+(ustar-WN[5].w)*(ustar+WN[5].p/(WN[5].d*(SR-WN[5].w))))-UN[5].E )*SR;

#ifdef DUAL_E
	  Us.d =(fact*1.);
	  Us.du=(fact*WN[5].u);
	  Us.dv=(fact*WN[5].v);
	  Us.dw=(fact*ustar);
	  Us.E =(fact*(UN[5].E/UN[5].d+(ustar-WN[5].w)*(ustar+WN[5].p/(WN[5].d*(SR-WN[5].w)))));
#endif
	}

#ifdef DUAL_E
	ebar=(Us.E-0.5*(Us.du*Us.du+Us.dv*Us.dv+Us.dw*Us.dw)/Us.d);
	ecen+=ebar;
	HR[5]=(Us.dw/Us.d*ebar);
	divu[5]=(GAMMA-1.)*(Us.dw/Us.d);
#endif

#endif


/* #ifdef RIEMANN_HLL */
/* 	speedestimateZ(&WC[5],&WN[5],&SL,&SR);//,&pstar,&ustar); */
/* 	/\* SL=WC[5].w-WC[5].a; *\/ */
/* 	/\* SR=WN[5].w+WN[5].a; *\/ */

/* 	if((SL<0.)&&(SR>0)){ */
/* 	  getflux_Z(&UC[5],Fl); */
/* 	  getflux_Z(&UN[5],Fr); */
	
/* 	  HR[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UN[5].d -UC[5].d ))/(SR-SL); */
/* 	  HR[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UN[5].du-UC[5].du))/(SR-SL); */
/* 	  HR[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UN[5].dv-UC[5].dv))/(SR-SL); */
/* 	  HR[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UN[5].dw-UC[5].dw))/(SR-SL); */
/* 	  HR[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UN[5].E -UC[5].E ))/(SR-SL); */

/* #ifdef DUAL_E */
/* 	  Fl[5]=WC[5].p*WC[5].w; */
/* 	  Fr[5]=WN[5].p*WN[5].w; */
/* 	  HR[5]=(SR*Fl[5]-SL*Fr[5]+SL*SR*(WN[5].p -WC[5].p))/(SR-SL); */
/* #endif */
	
/* 	} */
/* 	else{ */
/* 	  if(SL>=0.){ */
/* 	    getflux_Z(&UC[5],HR); */
/* #ifdef DUAL_E */
/* 	    HR[5]=WC[5].p*WC[5].w; */
/* #endif */
/* 	  } */
/* 	  else if(SR<=0.){ */
/* 	    getflux_Z(&UN[5],HR); */
/* #ifdef DUAL_E */
/* 	    HR[5]=WN[5].p*WN[5].w; */
/* #endif */
/* 	  } */
/* 	} */

/* #endif */

      }


#ifdef DUAL_E
      FL[5]+=divu[0]*ebar/6.; // CHANGE TO /6. when in 3D
      FR[5]+=divu[1]*ebar/6.; // CHANGE TO /6. when in 3D
      GL[5]+=divu[2]*ebar/6.; // CHANGE TO /6. when in 3D */
      GR[5]+=divu[3]*ebar/6.; // CHANGE TO /6. when in 3D */
      HL[5]+=divu[4]*ebar/6.; // CHANGE TO /6. when in 3D
      HR[5]+=divu[5]*ebar/6.; // CHANGE TO /6. when in 3D
#endif

      
      //========================= copy the fluxes

      //if((ffact[4])&&(!ffact[3])) abort();

      memcpy(stencil[i].New.cell[icell].flux+0*NVAR,FL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].New.cell[icell].flux+1*NVAR,FR,sizeof(REAL)*NVAR);
      memcpy(stencil[i].New.cell[icell].flux+2*NVAR,GL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].New.cell[icell].flux+3*NVAR,GR,sizeof(REAL)*NVAR);
      memcpy(stencil[i].New.cell[icell].flux+4*NVAR,HL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].New.cell[icell].flux+5*NVAR,HR,sizeof(REAL)*NVAR);


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

  //Smax=fmax(Smax,sqrt(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);
  // Computing new timestep
  dt=tmax;
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      if(nextoct==NULL) continue; // in case the level is empty
      dxcur=pow(0.5,level); // +1 to protect level change
      do // sweeping through the octs of level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	       vx=curoct->cell[icell].field.u; 
	       vy=curoct->cell[icell].field.v; 
	       vz=curoct->cell[icell].field.w; 
	       va=sqrt(vx*vx+vy*vy+vz*vz); 
	       aa=sqrt(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d); 
	       Smax=fmax(Smax,va+aa); 

	    }
	}while(nextoct!=NULL);

      //printf("thydro= %e Smax=%e dxcur=%e\n",dxcur*CFL/(Smax*3.),Smax,dxcur);
      //      abort();
      dt=fmin(dxcur*CFL/(Smax*3.),dt);
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

REAL L_comptstep_hydro(int level, struct RUNPARAMS *param,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax){
  
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icell;
  REAL aa;
  REAL va,vx,vy,vz;
  REAL dt;
  REAL Smax=0.,S1;

  //Smax=fmax(Smax,sqrt(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);
  // Computing new timestep
  dt=tmax;
  // setting the first oct
  
  nextoct=firstoct[level-1];
  
  if(nextoct!=NULL){
    dxcur=pow(0.5,level); // +1 to protect level change
    do // sweeping through the octs of level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    vx=curoct->cell[icell].field.u; 
	    vy=curoct->cell[icell].field.v; 
	    vz=curoct->cell[icell].field.w; 
	    va=sqrt(vx*vx+vy*vy+vz*vz); 
	    aa=sqrt(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d); 
	    Smax=fmax(Smax,va+aa); 

	  }
      }while(nextoct!=NULL);
  }

  dt=fmin(dxcur*CFL/(Smax*3.),dt);

  /* #ifdef WMPI */
  /*   // reducing by taking the smallest time step */
  /*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_REAL,MPI_MIN,cpu->comm); */
  /* #endif   */

  return dt;
}

//=================================================================================
#ifdef WGRAV
REAL L_comptstep_ff(int level,struct RUNPARAMS *param,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax){
  
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icell;
  REAL dtloc;
  REAL dt;

  dt=tmax;
  // setting the first oct
      
  nextoct=firstoct[level-1];
      
  if(nextoct!=NULL){
    dxcur=pow(0.5,level); // +1 to protect level change
    do // sweeping through the octs of level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    dtloc=0.1*sqrt(2.*M_PI/(3.*curoct->cell[icell].gdata.d*aexp));
	    dt=fmin(dt,dtloc);
	  }
      }while(nextoct!=NULL);
  }
  return dt;
}


REAL comptstep_ff(int levelcoarse,int levelmax,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax){
  
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icell;
  REAL dtloc;
  REAL dt;

  //Smax=fmax(Smax,sqrt(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);

  // Computing new timestep
  dt=tmax;
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      if(nextoct==NULL) continue; // in case the level is empty
      dxcur=pow(0.5,level); // +1 to protect level change
      do // sweeping through the octs of level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      dtloc=0.1*sqrt(2.*M_PI/(3.*curoct->cell[icell].gdata.d*aexp));
	      dt=fmin(dt,dtloc);
	    }
	}while(nextoct!=NULL);
    }

  // new tstep
  // printf("original dt=%f chosen dt=%f\n",tmax,dt);


/* #ifdef WMPI */
/*   // reducing by taking the smallest time step */
/*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_REAL,MPI_MIN,cpu->comm); */
/* #endif   */

  return dt;
}

#endif

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
  //Smax=fmax(Smax,sqrt(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);

  // Computing new timestep
  dt=tmax;
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      if(nextoct==NULL) continue; // in case the level is empty
      dxcur=pow(0.5,level); // +1 to protect level change
      do // sweeping through the octs of level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      memcpy(&W,&curoct->cell[icell].field,sizeof(struct Wtype));
	      W2U(&W,&U);
	      memcpy(f,curoct->cell[icell].f,sizeof(REAL)*3);
	      V=sqrt(pow(U.du,2.)+pow(U.dv,2.)+pow(U.dw,2.));
	      DV=sqrt(pow(U.d*f[0],2.)+pow(U.d*f[1],2.)+pow(U.d*f[2],2.));
	      DE=sqrt(pow(U.du*f[0],2.)+pow(U.dv*f[1],2.)+pow(U.dw*f[2],2.));
									       
	      if((DE>0.)&&(U.E>0.)){
		dtloc=0.5*U.E/DE;
		if(dt!=tmax){
		  dt=fmin(dt,dtloc);
		}else{
		  dt=dtloc;
		}
	      }

	      if((DV>0.)&&(V>0.)){
		dtloc=0.5*V/DV;
		if(dt!=tmax){
		  dt=fmin(dt,dtloc);
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
	  memcpy(&(curcell->field),&Wnew,sizeof(struct Wtype));
	}
      }
    }while(nextoct!=NULL);
  }
}
#endif
// ==============================================================================================
// ==============================================================================================
#ifdef WHYDRO2
//=====================================================================================================================
//
//    REFLECHISSONS UN PEU EN FAIT LA RECHERCHE EST RECURSIVE
//

// Structure de base

void recursive_neighbor_gather_oct(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct HGRID *stencil,char *visit){

  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;
  static int face[8]={0,1,2,3,4,5,6,7};
  REAL dxcur;

  struct Wtype Wi[8];
#ifdef NOFLUX
  struct OCT *child[8];
#endif

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
      else{
	printf("big problem\n");
	abort();
      }
    }
  }


  /* oct=cell2oct(cell); */
  /* neioct=cell2oct(neicell); */
  /* dxcur=pow(0.5,oct->level); */

  // ============================ TRANSMISSIVE BOUNDARIES ====================
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
      }
    }
#endif

#ifdef TRANSZM
    if(ineiloc==4){
      /* dist=neioct->z-oct->z; */
      /* if(dist>0.5){ */
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
      }
    }
#endif


  // ============================ END TRANSMISSIVE BOUNDARIES ====================

#ifdef WGRAV
    REAL floc[8*3];//local force
#endif

  if(neicell->child!=NULL){
    // optimal case
    for(icell=0;icell<8;icell++){
      memcpy(Wi+icell,&(neicell->child->cell[icell].field),sizeof(struct Wtype));
      //if(icell==0) if(ioct==24) if(neicell->child->cell[icell].field.p==0.099999999999999978) abort();
      
      
#ifdef NOFLUX
      child[icell]=neicell->child->cell[icell].child;
#endif
    }

#ifdef WGRAV
    for(icell=0;icell<8;icell++) memcpy(floc+3*icell,neicell->child->cell[icell].f,sizeof(REAL)*3);
#endif

  }
  else{
    coarse2fine_hydro2(neicell,Wi);
#ifdef NOFLUX
    for(icell=0;icell<8;icell++){
      child[icell]=NULL;
    }
#endif

#ifdef WGRAV
    for(icell=0;icell<8;icell++) memcpy(floc+3*icell,neicell->f,sizeof(REAL)*3);
#endif

  }




  for(icell=0;icell<8;icell++){
    memcpy(&(stencil->oct[ioct].cell[icell].field),Wi+face[icell],sizeof(struct Wtype)); //

#ifdef NOFLUX
    stencil->oct[ioct].cell[icell].child=child[icell];
#endif
    //if(icell==0) if(ioct==24) if(stencil->oct[ioct].cell[icell].field.p==0.099999999999999978) abort();

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
#endif

#ifdef WHYDRO2
//=====================================================================================================================

#ifndef NOFLUX
struct OCT *gatherstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  
  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int ioct;
  char visit[27]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
   
  //printf("let's gather\n");
  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;


      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(&(stencil[iread].oct[13].cell[icell].field),&(curoct->cell[icell].field),sizeof(struct Wtype)); //
#ifdef WGRAV 
 	memcpy(stencil[iread].oct[13].cell[icell].f,curoct->cell[icell].f,sizeof(REAL)*3); // 
#endif 
	

      }
      visit[13]=1;
      //abort();
      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  visit[ioct]=1;
	  recursive_neighbor_gather_oct(ioct, inei, -1, -1, 1, cell, stencil+iread);
	}

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}
#else

struct OCT *gatherstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  char visit[27]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
  int ix[6]={-1,1,0,0,0,0};
  int iy[6]={0,0,-1,1,0,0};
  int iz[6]={0,0,0,0,-1,1};
  int ioct;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;


      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(&(stencil[iread].oct[13].cell[icell].field),&(curoct->cell[icell].field),sizeof(struct Wtype)); // for calculations
	memcpy(&(stencil[iread].New.cell[icell].field),&(curoct->cell[icell].fieldnew),sizeof(struct Wtype)); // for updates

#ifdef WGRAV 
 	memcpy(stencil[iread].oct[13].cell[icell].f,curoct->cell[icell].f,sizeof(REAL)*3); // 
#endif 
	stencil[iread].oct[13].cell[icell].child=curoct->cell[icell].child;
      }
      
      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  recursive_neighbor_gather_oct(ioct, inei, -1, -1, 1, cell, stencil+iread,visit);
	}
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}
#endif

#endif


// ================================================================
// ================================================================
// ================================================================

#ifdef WHYDRO2

#ifndef NOFLUX
struct OCT *scatterstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int ioct[7]={12,14,10,16,4,22,13};
  
  nextoct=octstart;
  iread=0;

  //printf("let's scatter\n");
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(&(curoct->cell[icell].flux),&(stencil[iread].New.cell[icell].flux),sizeof(REAL)*NFLUX);
      }
      

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}

#else

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

  struct Wtype W,W0,Wi;
  struct Utype U,U0,Ui;

#ifdef DUAL_E
  REAL DE,p0,p;
  REAL eo,e;
#endif
  REAL dtsurdx=dtnew/dxcur;
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
	memcpy(&(curoct->cell[icell].fieldnew),&(stencil[iread].New.cell[icell].field),sizeof(struct Wtype));
	
	
	// let us deal with coarser neighbors
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
	      memcpy(&W,&(curoct->nei[vnei[inei]]->fieldnew),sizeof(struct Wtype));
	      W2U(&W,&U);

	      memcpy(&Ui,&U,sizeof(struct Utype));
	      memcpy(&Wi,&W,sizeof(struct Wtype));
	      memcpy(&W0,&(curoct->nei[vnei[inei]]->field),sizeof(struct Wtype));
	      W2U(&W0,&U0);

	      // getting the flux
	      memcpy(F,stencil[iread].New.cell[icell].flux+inei*NVAR,sizeof(REAL)*NVAR);
	      
	      // update
	      one=pow(-1.,inei+1);
	      U.d +=F[0]*dtsurdx*one*0.125;
	      U.du+=F[1]*dtsurdx*one*0.125;
	      U.dv+=F[2]*dtsurdx*one*0.125;
	      U.dw+=F[3]*dtsurdx*one*0.125;
	      U.E +=F[4]*dtsurdx*one*0.125;
#ifdef DUAL_E
	      U.eint+=F[5]*dtsurdx*one*0.125;
#endif
	      
	      
	      //scatter back
	      U2W(&U,&W);
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

#endif
#endif

// ================================================================
// ================================================================
// ================================================================

#ifdef NOFLUX
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
      
      if(stencil[i].oct[13].cell[icell].child!=NULL) continue;
      memcpy(&W,&(stencil[i].New.cell[icell].field),sizeof(struct Wtype));// getting the original state (already contained in the new field)
      W2U(&W,&U);
      memcpy(F,stencil[i].New.cell[icell].flux,sizeof(REAL)*NFLUX);// original fluxes

      memcpy(&Ui,&U,sizeof(struct Utype));
      memcpy(&Wi,&W,sizeof(struct Wtype));
      memcpy(&W0,&(stencil[i].oct[13].cell[icell].field),sizeof(struct Wtype));// getting the original state (already contained in the new field)
      W2U(&W0,&U0);

      
      // ==== updating
      one=1.;
      for(flx=0;flx<6;flx++){
	U.d +=F[0+flx*NVAR]*dtsurdx*one;
	U.du+=F[1+flx*NVAR]*dtsurdx*one;
	U.dv+=F[2+flx*NVAR]*dtsurdx*one;
	U.dw+=F[3+flx*NVAR]*dtsurdx*one;
	U.E +=F[4+flx*NVAR]*dtsurdx*one;
#ifdef DUAL_E
	U.eint+=F[5+flx*NVAR]*dtsurdx*one;
#endif
	one*=-1.;
      }
      // scatter back in the new value within the stencil
      U2W(&U,&W);
      if(W.p<0) abort();
      //      if(W.p>2.2e-9) abort();
      memcpy(&(stencil[i].New.cell[icell].field),&W,sizeof(struct Wtype));
    }
  }

}
#endif

// =========================================================================================
// =========================================================================================

#ifdef NOFLUX
void advancehydro(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL dtnew){

  double t0,t100,t20,t80,t200,t150;
  double th=0.,tt=0.;
  int nread,nreadtot;;
  struct OCT *curoct;
  struct OCT *nextoct;
  int vnei[6],vcell[6];
  
  int flx;
  REAL F[NFLUX];
  REAL Forg[NFLUX];
  REAL dxcur=pow(0.5,level);
  REAL dtsurdx=dtnew/dxcur;
  REAL one;
  struct Utype U,UL;
  struct Utype S;
  struct Utype U0;
  struct Wtype W,WL;
  struct Wtype Wnew;
  struct CELL *neicell;
  int icell;
  
  int nocthydro=cpu->noct[level-1];

  if(cpu->rank==0) printf("Start Hydro on %d octs with dt=%e on level %d with stride=%d\n",nocthydro,dtnew,level,stride);
  
  // ===== COMPUTING THE FLUXES

  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  nreadtot=0;
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      t0=MPI_Wtime();
      curoct=nextoct;
      nextoct=curoct->next; 
      // -------------  cleaning working arrays
      
      memset(stencil,0,stride*sizeof(struct HGRID));
	    
      // ------------ gathering the stencil value values
      nextoct= gatherstencil(curoct,stencil,stride,cpu, &nread);
      
      // ------------ solving the hydro
	    
      t20=MPI_Wtime();
      hydroM(stencil,level,cpu->rank,nread,stride,dxcur,dtnew);
      t80=MPI_Wtime();
      // ------------ updating values within the stencil

      updatefield(curoct,stencil,nread,stride,cpu,dxcur,dtnew);
      
      // ------------ scatter back the FLUXES
	    
      nextoct=scatterstencil(curoct,stencil, nread, cpu,dxcur,dtnew);

      nreadtot+=nread;
      t100=MPI_Wtime();

      th+=t80-t20;
      tt+=(t100-t0);
	  

    }while(nextoct!=NULL);
  }
  t150=MPI_Wtime();

  // FINAL UPDATE OF THE VALUES
  if(nreadtot>0){
    nextoct=firstoct[level-1];
    do {
      curoct=nextoct;
      nextoct=curoct->next; 
      for(icell=0;icell<8;icell++){
	if(curoct->cell[icell].child==NULL){
	  // unsplit case
	  memcpy(&(curoct->cell[icell].field),&(curoct->cell[icell].fieldnew),sizeof(struct Wtype));
	}
	else{
	  // split case : -> quantities are averaged
	  struct OCT *child;
	  int i;
	  child=curoct->cell[icell].child;
	  memset(&W,0,sizeof(struct Wtype));
	  for(i=0;i<8;i++){
	    W.d+=child->cell[i].field.d*0.125;
	    W.u+=child->cell[i].field.u*0.125;
	    W.v+=child->cell[i].field.v*0.125;
	    W.w+=child->cell[i].field.w*0.125;
	    W.p+=child->cell[i].field.p*0.125;

	  }
	  //if(W.v!=0.) abort();
	  //U2W(&U,&W);
	  memcpy(&curoct->cell[icell].field,&W,sizeof(struct Wtype));
	}


	//if(curoct->cell[icell].field.w <0.) abort();

      }
    }while(nextoct!=NULL);
  }



  if(cpu->rank==0) printf("\n Timings per oct [total]: \n tt=%e[%e] \n th=%e[%e] \n tf=%e[%e]\n",tt/nocthydro,tt,th/nocthydro,th,(t200-t150)/nocthydro,t200-t150);
  
}

void clean_new_hydro(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu){

  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next; 

      //for(icell=0;icell<8;icell++) memset(&(curoct->cell[icell].fieldnew),0,sizeof(struct Wtype));
      for(icell=0;icell<8;icell++) memcpy(&(curoct->cell[icell].fieldnew),&(curoct->cell[icell].field),sizeof(struct Wtype));
      
    }while(nextoct!=NULL);
  }
}
#endif

#endif