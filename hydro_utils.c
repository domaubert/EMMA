#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>

#define NITERMAX 10
#define ERRTOL 1e-10
#define CFL 0.4

#ifdef WHYDRO2

// ==================== converts U -> W
void U2W(struct Utype *U, struct Wtype *W)
{
  W->d=U->d;
  W->u=U->du/U->d;
  W->v=U->dv/U->d;
  W->w=U->dw/U->d;
  W->p=(GAMMA-1.)*(U->E-((U->du)*(U->du)+(U->dv)*(U->dv)+(U->dw)*(U->dw))/(U->d)*0.5);
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
	      if((oct->nei[vnei[inei2]]->child->z-oct->z)>0.5){
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
	      if((oct->nei[vnei[inei2]]->child->z-oct->z)<0.){
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




// ==================== pressure solver

double frootprime(double p, struct Wtype1D_double *WL, struct Wtype1D_double *WR)
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

double froot(double p, struct Wtype1D_double *WL, struct Wtype1D_double *WR, double *u)
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




// --------------------------------------

double findPressure(struct Wtype1D_double *WL, struct Wtype1D_double *WR, int *niter, REAL *u)
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

  if((WL->p<0)||(WR->p<0)){
    printf("ouh!");
    abort();
    }

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

  pmin=fmin(WL->p,WR->p);
  pmax=fmax(WL->p,WR->p);

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
      
      /* if((p-dp)<0){ */
      /* 	dp=p*0.5; */
      /* } */

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

	    REAL p;
	    p=Ui[idir].E-0.5*(Ui[idir].du*Ui[idir].du+Ui[idir].dv*Ui[idir].dv+Ui[idir].dw*Ui[idir].dw)/Ui[idir].d;

	    if(p<0){
	      printf("ouch !\n");
	      abort();
	    }

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


void  matrix_jacobian(struct Wtype *W0, REAL dt,REAL dx,struct Wtype *D, struct Wtype *Wt){

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

  d[0]=D[0].d;
  d[1]=D[0].u;
  d[2]=D[0].v;
  d[3]=D[0].w;
  d[4]=D[0].p;

  for(j=0;j<5;j++){
    for(i=0;i<5;i++){
      W[i]+=M[i+j*5]*d[j];
      }
  }
  
  if(isnan(W[4])) abort();
  // =====  building the B matrix

  memset(M,0,25*sizeof(REAL));
  
  // diagonal elements
  for(i=0;i<5;i++) M[i+i*5]=W0->v;
  
  // off_diagonal elements
  M[0+2*5]=W0->d;

  M[4+2*5]=W0->d*W0->a*W0->a;

  M[2+4*5]=1./W0->d;


  // ===== Second Product

  d[0]=D[1].d;
  d[1]=D[1].u;
  d[2]=D[1].v;
  d[3]=D[1].w;
  d[4]=D[1].p;

  for(j=0;j<5;j++){
    for(i=0;i<5;i++){
      W[i]+=M[i+j*5]*d[j];
      }
  }

  if(isnan(W[4])) abort();
  // =====  building the C matrix

  memset(M,0,25*sizeof(REAL));
  
  // diagonal elements
  for(i=0;i<5;i++) M[i+i*5]=W0->w;
  
  // off_diagonal elements
  M[0+3*5]=W0->d;

  M[4+3*5]=W0->d*W0->a*W0->a;

  M[3+4*5]=1./W0->d;

  d[0]=D[2].d;
  d[1]=D[2].u;
  d[2]=D[2].v;
  d[3]=D[2].w;
  d[4]=D[2].p;

  
  for(j=0;j<5;j++){
    for(i=0;i<5;i++){
      W[i]+=M[i+j*5]*d[j];
      }
  }
  
  if(isnan(W[4])) abort();
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

	    diffW(Wp,W0,&Dp); 
	    diffW(W0,Wm,&Dm); 
	    
	    minmod_W(&Dm,&Dp,D+dir);

	    /* if(dir==0){ */
	    /*   if(W0->x==(0.015625)){ */
	    /* 	printf("LEFT DU=%e Dm=%e Dp=%e [%e]\n",D[0].u,Dm.u,Dp.u,Dm.u/Dp.u); */
	    /*   } */
	      
	    /*   if(W0->x==(1.-0.015625)){ */
	    /* 	printf("RIGHT DU=%e Dm=%e Dp=%e [%e]\n",D[0].u,Dm.u,Dp.u,Dm.u/Dp.u); */
	    /* 	abort(); */
	    /*   } */
	    /* } */
	  }


	  // build jacobian matrix product

	  
	  matrix_jacobian(W0,dt,dx,D,&Wt); // Here Wt contains the evolution of the state
	  

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
	    Wi[idir].u = W0->u+ix[idir]*D[0].u+iy[idir]*D[1].u+iz[idir]*D[2].u+Wt.u;//-f[0]*0.5*dt;
	    Wi[idir].v = W0->v+ix[idir]*D[0].v+iy[idir]*D[1].v+iz[idir]*D[2].v+Wt.v;//-f[1]*0.5*dt;
	    Wi[idir].w = W0->w+ix[idir]*D[0].w+iy[idir]*D[1].w+iz[idir]*D[2].w+Wt.w;//-f[2]*0.5*dt;
	    Wi[idir].p = W0->p+ix[idir]*D[0].p+iy[idir]*D[1].p+iz[idir]*D[2].p+Wt.p;
	 
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
	    if(isnan(Wi[idir].p)) abort();
	    
	    Wi[idir].a=sqrt(GAMMA*Wi[idir].p/Wi[idir].d);
	    if(Wi[idir].p<0) abort();
	  }
	  


	  
}


//============================================================================

void speedestimate(struct Wtype *WL,struct Wtype *WR, REAL *SL, REAL *SR){

  REAL splus;
  splus=fmax(fabs(WL->u)+WL->a,fabs(WR->u)+WR->a);
  *SL=-splus;
  *SR= splus;

}



int hydroM(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt){

  int inei,icell,icellcoarse;
  int i;
  REAL temp;
  REAL res,restot=0.;
  int vnei[6],vcell[6];
  int vneic[6],vcellc[6];

  REAL FL[NVAR],FR[NVAR];
  REAL GL[NVAR],GR[NVAR];
  REAL HL[NVAR],HR[NVAR];

  memset(FL,0,sizeof(REAL)*NVAR);
  memset(FR,0,sizeof(REAL)*NVAR);
  memset(HL,0,sizeof(REAL)*NVAR);
  memset(HR,0,sizeof(REAL)*NVAR);
  memset(GL,0,sizeof(REAL)*NVAR);
  memset(GR,0,sizeof(REAL)*NVAR);

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

#ifdef RIEMANN_HLL
  REAL SL,SR;
  REAL Fr[NVAR],Fl[NVAR];
#endif
  int tagr=0,tagl=0;

  //printf("let's do some hydro\n");
  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
    for(i=0;i<nread;i++){ // we scan the octs
      tagr=0;
      tagl=0;
      
      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].field);
      
      Wold.d=curcell->d;
      Wold.u=curcell->u;;
      Wold.v=curcell->v;;
      Wold.w=curcell->w;;
      Wold.p=curcell->p;;
      Wold.a=sqrt(GAMMA*Wold.p/Wold.d);

      W2U(&Wold,&Uold); // primitive -> conservative
      

      /* for(inei=0;inei<6;inei++){ */
      /* 	memcpy(WC+inei,&Wold,sizeof(struct Wtype)); */
      /* 	memcpy(WN+inei,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].field),sizeof(struct Wtype)); */
      /* 	W2U(WN+inei,UN+inei); */
      /* 	W2U(WC+inei,UC+inei); */
      /* } */

      /* // MUSCL STATE RECONSTRUCTION */

      MUSCL_BOUND2(stencil+i, 13, icell, WC,dt,dx);// central
      for(inei=0;inei<6;inei++){
      	MUSCL_BOUND2(stencil+i, ioct[vnei[inei]], vcell[inei], WT,dt,dx);//
      	memcpy(WN+inei,WT+idxnei[inei],sizeof(struct Wtype));
       	W2U(WN+inei,UN+inei);
      	W2U(WC+inei,UC+inei);
      }

      
      // X DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems LEFT

      // Switching to Split description
      
#ifdef RIEMANN_EXACT

      WLloc.d=WN[0].d;
      WLloc.u=WN[0].u;
      WLloc.p=WN[0].p;
      WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);

      WRloc.d=WC[0].d;
      WRloc.u=WC[0].u;
      WRloc.p=WC[0].p;
      WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

      // Riemann Solver
      if((WLloc.p<0)||(WRloc.p<0)){
	printf("hajzehr\n");
	abort();
      }

      pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      Wtest3D.d=Wtest.d;
      Wtest3D.u=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      if(isnan(Wtest.p)){
	printf("ouhla NAN\n");
	abort();
      }

      // Passive advection
      if(ustar>0.)
	{
	  Wtest3D.v=WN[0].v;
	  Wtest3D.w=WN[0].w;
	}
      else
	{
	  Wtest3D.v=WC[0].v;
	  Wtest3D.w=WC[0].w;
	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes LEFT
      getflux_X(&Utest3D,FL);

#endif

#ifdef RIEMANN_HLL
      speedestimate(&WN[0],&WC[0],&SL,&SR);
      /* SL=WN[0].u-WN[0].a; */
      /* SR=WC[0].u+WC[0].a; */

      if((SL<0.)&&(SR>0)){
	getflux_X(&UN[0],Fl);
	getflux_X(&UC[0],Fr);

	
	FL[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UC[0].d -UN[0].d ))/(SR-SL);
	FL[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UC[0].du-UN[0].du))/(SR-SL);
	FL[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UC[0].dv-UN[0].dv))/(SR-SL);
	FL[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UC[0].dw-UN[0].dw))/(SR-SL);
	FL[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UC[0].E -UN[0].E ))/(SR-SL);

#ifdef DUAL_E
	Fl[5]=WN[0].p*WN[0].u;
	Fr[5]=WC[0].p*WC[0].u;
	FL[5]=(SR*Fl[5]-SL*Fr[5]+SL*SR*(WC[0].p -WN[0].p ))/(SR-SL);
#endif

	tagl=1;
      }
      else{
	if(SL>=0.){
	  getflux_X(&UN[0],FL);
#ifdef DUAL_E
	  FL[5]=WN[0].p*WN[0].u;
#endif
	}
	else if(SR<=0.){
	  getflux_X(&UC[0],FL);
#ifdef DUAL_E
	  FL[5]=WC[0].p*WC[0].u;
#endif
	}
      }
      

      
#endif



      // --------- solving the Riemann Problems RIGHT


      // Switching to Split description

#ifdef RIEMANN_EXACT
      WRloc.d=WN[1].d;
      WRloc.u=WN[1].u;
      WRloc.p=WN[1].p;
      WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

      WLloc.d=WC[1].d;
      WLloc.u=WC[1].u;
      WLloc.p=WC[1].p;
      WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);

      if((WLloc.p<0)||(WRloc.p<0)){
	printf("hajzehr\n");
	abort();
      }

      // Riemann Solver
      pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.u=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar<0.)
	{
	  Wtest3D.v=WN[1].v;
	  Wtest3D.w=WN[1].w;
	}
      else
	{
	  Wtest3D.v=WC[1].v;
	  Wtest3D.w=WC[1].w;
	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes RIGHT
      getflux_X(&Utest3D,FR);

      if(isnan(FL[0])){
	printf("NAN GL\n");
	abort();
      }

#endif


#ifdef RIEMANN_HLL
      speedestimate(&WC[1],&WN[1],&SL,&SR);
      /* SL=WC[1].u-WC[1].a; */
      /* SR=WN[1].u+WN[1].a; */

      if((SL<0.)&&(SR>0)){
	getflux_X(&UC[1],Fl);
	getflux_X(&UN[1],Fr);
	
	FR[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UN[1].d -UC[1].d ))/(SR-SL);
	FR[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UN[1].du-UC[1].du))/(SR-SL);
	FR[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UN[1].dv-UC[1].dv))/(SR-SL);
	FR[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UN[1].dw-UC[1].dw))/(SR-SL);
	FR[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UN[1].E -UC[1].E ))/(SR-SL);

#ifdef DUAL_E
	Fl[5]=WC[1].p*WC[1].u;
	Fr[5]=WN[1].p*WN[1].u;
	FR[5]=(SR*Fl[5]-SL*Fr[5]+SL*SR*(WN[1].p -WC[1].p))/(SR-SL);
#endif

	tagr=1;
      }
      else{
	if(SL>=0.){
	  getflux_X(&UC[1],FR);
#ifdef DUAL_E
	  FR[5]=WC[1].p*WC[1].u;
#endif
	}
	else if(SR<=0.){
	  getflux_X(&UN[1],FR);
#ifdef DUAL_E
	  FR[5]=WN[1].p*WN[1].u;
#endif
	}
      }
      if(isnan(FR[5])) abort();
      

      
#endif


      /* if(curcell->x==0.015625){ */
      /* 	printf("LEFT %e %d || %e %e %e -- %e %e %e ==> %e %e <== %e\n",curcell->x,tagl,WC[1].d,WC[1].u,WC[1].p,WN[1].d,WN[1].u,WN[1].p,FL[4],FR[4],Wold.u); */
      /* } */

      /* if(curcell->x==(1.-0.015625)){ */
      /* 	printf("RIGHT %e %d || %e %e %e -- %e %e %e ==> %e %e <== %e \n",curcell->x,tagr,WN[0].d,WN[0].u,WN[0].p,WC[0].d,WC[0].u,WC[0].p,FL[4],FR[4],Wold.u); */
      /* 	abort(); */
      /* } */


       // Y DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems FRONT

      // Switching to Split description


#ifdef RIEMANN_EXACT
      WLloc.d=WN[2].d;
      WLloc.u=WN[2].v;
      WLloc.p=WN[2].p;
      WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);

      WRloc.d=WC[2].d;
      WRloc.u=WC[2].v;
      WRloc.p=WC[2].p;
      WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

      // Riemann Solver
      if((WLloc.p<0)||(WRloc.p<0)){
	printf("hajzehr\n");
	abort();
      }
      pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar);

      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.v=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar>0.)
	{
	  Wtest3D.u=WN[2].u;
	  Wtest3D.w=WN[2].w;
	}
      else
	{
	  Wtest3D.u=WC[2].u;
	  Wtest3D.w=WC[2].w;
	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes LEFT
      getflux_Y(&Utest3D,GL);


#endif

#ifdef RIEMANN_HLL
      speedestimate(&WN[2],&WC[2],&SL,&SR);
      /* SL=WN[2].v-WN[2].a; */
      /* SR=WC[2].v+WC[2].a; */

      if((SL<0.)&&(SR>0)){
	getflux_Y(&UN[2],Fl);
	getflux_Y(&UC[2],Fr);
	
	GL[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UC[2].d -UN[2].d ))/(SR-SL);
	GL[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UC[2].du-UN[2].du))/(SR-SL);
	GL[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UC[2].dv-UN[2].dv))/(SR-SL);
	GL[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UC[2].dw-UN[2].dw))/(SR-SL);
	GL[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UC[2].E -UN[2].E ))/(SR-SL);

      }
      else{
	if(SL>=0.){
	  getflux_Y(&UN[2],GL);
	}
	else if(SR<=0.){
	  getflux_Y(&UC[2],GL);
	}
      }
      
#endif

      // --------- solving the Riemann Problems BACK


      // Switching to Split description

#ifdef RIEMANN_EXACT
      WRloc.d=WN[3].d;
      WRloc.u=WN[3].v;
      WRloc.p=WN[3].p;
      WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

      WLloc.d=WC[3].d;
      WLloc.u=WC[3].v;
      WLloc.p=WC[3].p;
      WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);

      if((WLloc.p<0)||(WRloc.p<0)){
	printf("hajzehr\n");
	abort();
      }
      // Riemann Solver
      pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.v=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar<0.)
	{
	  Wtest3D.u=WN[3].u;
	  Wtest3D.w=WN[3].w;
	}
      else
	{
	  Wtest3D.u=WC[3].u;
	  Wtest3D.w=WC[3].w;
	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes RIGHT
      getflux_Y(&Utest3D,GR);
#endif

#ifdef RIEMANN_HLL
      speedestimate(&WC[3],&WN[3],&SL,&SR);
      /* SL=WC[3].v-WC[3].a; */
      /* SR=WN[3].v+WN[3].a; */

      if((SL<0.)&&(SR>0)){
	getflux_Y(&UC[3],Fl);
	getflux_Y(&UN[3],Fr);
	
	GR[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UN[3].d -UC[3].d ))/(SR-SL);
	GR[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UN[3].du-UC[3].du))/(SR-SL);
	GR[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UN[3].dv-UC[3].dv))/(SR-SL);
	GR[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UN[3].dw-UC[3].dw))/(SR-SL);
	GR[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UN[3].E -UC[3].E ))/(SR-SL);

      }
      else{
	if(SL>=0.){
	  getflux_Y(&UC[3],GR);
	}
	else if(SR<=0.){
	  getflux_Y(&UN[3],GR);
	}
      }
      
#endif


      // Z DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems BOTTOM

      // Switching to Split description


#ifdef RIEMANN_EXACT
      WLloc.d=WN[4].d;
      WLloc.u=WN[4].w;
      WLloc.p=WN[4].p;
      WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);

      WRloc.d=WC[4].d;
      WRloc.u=WC[4].w;
      WRloc.p=WC[4].p;
      WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

      if((WLloc.p<0)||(WRloc.p<0)){
	printf("hajzehr\n");
	abort();
      }
      // Riemann Solver
      pstar=(REAL)findPressure(&WLloc,&WRloc,&n,&ustar);

      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.w=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar>0.)
	{
	  Wtest3D.u=WN[4].u;
	  Wtest3D.v=WN[4].v;
	}
      else
	{
	  Wtest3D.u=WC[4].u;
	  Wtest3D.v=WC[4].v;
	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes LEFT
      getflux_Z(&Utest3D,HL);
#endif

#ifdef RIEMANN_HLL
      speedestimate(&WN[4],&WC[4],&SL,&SR);
      /* SL=WN[4].w-WN[4].a; */
      /* SR=WC[4].w+WC[4].a; */

      if((SL<0.)&&(SR>0.)){
	getflux_Z(&UN[4],Fl);
	getflux_Z(&UC[4],Fr);
	
	HL[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UC[4].d -UN[4].d ))/(SR-SL);
	HL[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UC[4].du-UN[4].du))/(SR-SL);
	HL[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UC[4].dv-UN[4].dv))/(SR-SL);
	HL[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UC[4].dw-UN[4].dw))/(SR-SL);
	HL[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UC[4].E -UN[4].E ))/(SR-SL);

      }
      else{
	if(SL>=0.){
	  getflux_Z(&UN[4],HL);
	}
	else if(SR<=0.){
	  getflux_Z(&UC[4],HL);
	}
      }
      


#endif


      // --------- solving the Riemann Problems Top


      // Switching to Split description

#ifdef RIEMANN_EXACT
      WRloc.d=WN[5].d;
      WRloc.u=WN[5].w;
      WRloc.p=WN[5].p;
      WRloc.a=sqrt(GAMMA*WRloc.p/WRloc.d);

      WLloc.d=WC[5].d;
      WLloc.u=WC[5].w;
      WLloc.p=WC[5].p;
      WLloc.a=sqrt(GAMMA*WLloc.p/WLloc.d);

      if((WLloc.p<0)||(WRloc.p<0)){
	printf("hajzehr\n");
	abort();
      }
      // Riemann Solver
      pstar=findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.w=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar<0.)
	{
	  Wtest3D.u=WN[5].u;
	  Wtest3D.v=WN[5].v;
	}
      else
	{
	  Wtest3D.u=WC[5].u;
	  Wtest3D.v=WC[5].v;
	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes RIGHT
      getflux_Z(&Utest3D,HR);
#endif

#ifdef RIEMANN_HLL
      speedestimate(&WC[5],&WN[5],&SL,&SR);
      /* SL=WC[5].w-WC[5].a; */
      /* SR=WN[5].w+WN[5].a; */

      if((SL<0.)&&(SR>0)){
	getflux_Z(&UC[5],Fl);
	getflux_Z(&UN[5],Fr);
	
	HR[0]=(SR*Fl[0]-SL*Fr[0]+SL*SR*(UN[5].d -UC[5].d ))/(SR-SL);
	HR[1]=(SR*Fl[1]-SL*Fr[1]+SL*SR*(UN[5].du-UC[5].du))/(SR-SL);
	HR[2]=(SR*Fl[2]-SL*Fr[2]+SL*SR*(UN[5].dv-UC[5].dv))/(SR-SL);
	HR[3]=(SR*Fl[3]-SL*Fr[3]+SL*SR*(UN[5].dw-UC[5].dw))/(SR-SL);
	HR[4]=(SR*Fl[4]-SL*Fr[4]+SL*SR*(UN[5].E -UC[5].E ))/(SR-SL);
	
      }
      else{
	if(SL>=0.){
	  getflux_Z(&UC[5],HR);
	}
	else if(SR<=0.){
	  getflux_Z(&UN[5],HR);
	}
      }
      
      if(HR[3]==0.) abort();
      
#endif

      
      //========================= copy the fluxes

#ifdef DUAL_E
      GL[5]=0.;
      GR[5]=0.;
      HL[5]=0.;
      HR[5]=0.;
#endif

      memcpy(stencil[i].new.cell[icell].flux+0*NVAR,FL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].new.cell[icell].flux+1*NVAR,FR,sizeof(REAL)*NVAR);
      memcpy(stencil[i].new.cell[icell].flux+2*NVAR,GL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].new.cell[icell].flux+3*NVAR,GR,sizeof(REAL)*NVAR);
      memcpy(stencil[i].new.cell[icell].flux+4*NVAR,HL,sizeof(REAL)*NVAR);
      memcpy(stencil[i].new.cell[icell].flux+5*NVAR,HR,sizeof(REAL)*NVAR);


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
	      dtloc=0.1*sqrt(2.*M_PI/(3.*curoct->cell[icell].density*aexp));
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


#endif
