#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>

#define NITERMAX 10
#define ERRTOL 1e-6
#define CFL 0.8

#ifdef WHYDRO

// ==================== converts U -> W
void U2W(struct Utype *U, struct Wtype *W)
{
  W->d=U->d;
  W->u=U->du/U->d;
  W->v=U->dv/U->d;
  W->w=U->dw/U->d;
  W->p=(GAMMA-1.)*(U->E-((U->du)*(U->du)+(U->dv)*(U->dv)+(U->dw)*(U->dw))/(U->d)*0.5);
  W->a=sqrtf(GAMMA*W->p/W->d);
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

// ==================== pressure solver

float frootprime(float p, struct Wtype1D *WL, struct Wtype1D *WR)
{
  
  float fL,fR;
  float AL,AR,BL,BR;
  float Deltau;

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  fL=(p>WL->p?sqrtf(AL/(BL+p))*(1.-(p-WL->p)/(2.*(BL+p))):powf(p/WL->p,-(GAMMA+1)/(2.*GAMMA))/(WL->d*WL->a));
  fR=(p>WR->p?sqrtf(AR/(BR+p))*(1.-(p-WR->p)/(2.*(BR+p))):powf(p/WR->p,-(GAMMA+1)/(2.*GAMMA))/(WR->d*WR->a));

  return fL+fR;
}

// ------------------------------------

float froot(float p, struct Wtype1D *WL, struct Wtype1D *WR, float *u)
{
  
  float fL,fR;
  float AL,AR,BL,BR;
  float Deltau;

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  fL=(p>WL->p?(p-WL->p)*sqrtf(AL/(BL+p)):2.*WL->a/(GAMMA-1.)*(powf(p/WL->p,(GAMMA-1)/(2.*GAMMA))-1.));
  fR=(p>WR->p?(p-WR->p)*sqrtf(AR/(BR+p)):2.*WR->a/(GAMMA-1.)*(powf(p/WR->p,(GAMMA-1)/(2.*GAMMA))-1.));
  
  Deltau=WR->u-WL->u;
  *u=0.5*(WL->u+WR->u)+0.5*(fR-fL);

  return fL+fR+Deltau;
}

// --------------------------------------

float findPressure(struct Wtype1D *WL, struct Wtype1D *WR, int *niter, float *u)
{
  float p,porg,dp;
  int i;
  float err;

  p=0.5*(WL->p+WR->p);
  //p=0.5*(WL->p+WR->p)-0.125*(WR->u-WL->u)*(WR->d+WL->d)*(WR->a+WL->a);
  p=fmaxf(p,ERRTOL);

  *niter=0;
  for(i=0;i<NITERMAX;i++)
    {
      dp=froot(p,WL,WR,u)/frootprime(p,WL,WR);
      porg=p;
      p=p-dp;
      //printf("p0=%e dp=%e p=%e fprime=%e\n",porg,dp,p,frootprime(p,WL,WR));
      err=2.*fabsf(p-porg)/(fabsf(p+porg));
      *niter=*niter+1;
      if(err<ERRTOL) break;
    }

  froot(p,WL,WR,u); // last calculation to get u;

  return p;
}


// ================================= Riemann Solver

void getW(struct Wtype1D *W, float S, struct Wtype1D *WL, struct Wtype1D *WR, float pstar, float ustar)
{
  float SL,SR;
  float SHL,STL;
  float SHR,STR;
  float astar;

  if(S<ustar)
    {
      // left of contact
      if(pstar>WL->p)
	{
	  // left Shock with shock speed SL
	  SL=WL->u-WL->a*sqrtf((GAMMA+1)/(2.*GAMMA)*pstar/WL->p+(GAMMA-1.)/(2.*GAMMA)); 
	  
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
	      W->d=WL->d*((pstar/WL->p+(GAMMA-1.)/(GAMMA+1.))/((GAMMA-1)/(GAMMA+1)*pstar/WL->p+1.));
	      W->u=ustar;
	      W->p=pstar;
	      W->a=sqrtf(GAMMA*W->p/W->d);
	    }
	}
      else
	{
	  // left fan
	  astar=WL->a*powf(pstar/WL->p,(GAMMA-1)/(2.*GAMMA)); // sound speed behind the fan
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
		  W->d=WL->d*powf(pstar/WL->p,1./GAMMA);
		  W->u=ustar;
		  W->p=pstar;
		  W->a=sqrtf(GAMMA*W->p/W->d);
		}
	      else
		{
		  W->d=WL->d*powf(2./(GAMMA+1.)+(GAMMA-1.)/((GAMMA+1.)*WL->a)*(WL->u-S),2./(GAMMA-1.));
		  W->u=2./(GAMMA+1.)*(WL->a+(GAMMA-1.)/2.*WL->u+S);
		  W->p=WL->p*powf(2./(GAMMA+1.)+(GAMMA-1.)/((GAMMA+1.)*WL->a)*(WL->u-S),2.*GAMMA/(GAMMA-1.));
		  W->a=sqrtf(GAMMA*W->p/W->d);
		}
	    }
	}

    }
  else
    {
      if(pstar>WR->p)
	{
	  // Right Shock with shock speed SR
	  SR=WR->u+WR->a*sqrtf((GAMMA+1)/(2.*GAMMA)*pstar/WR->p+(GAMMA-1.)/(2.*GAMMA)); 
	  
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
	      W->d=WR->d*((pstar/WR->p+(GAMMA-1.)/(GAMMA+1.))/((GAMMA-1)/(GAMMA+1)*pstar/WR->p+1.));
	      W->u=ustar;
	      W->p=pstar;
	      W->a=sqrtf(GAMMA*W->p/W->d);
	    }
	}
      else
	{
	  // Right fan
	  astar=WR->a*powf(pstar/WR->p,(GAMMA-1)/(2.*GAMMA)); // sound speed behind the fan
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
		  W->d=WR->d*powf(pstar/WR->p,1./GAMMA);
		  W->u=ustar;
		  W->p=pstar;
		  W->a=sqrtf(GAMMA*W->p/W->d);
		}
	      else
		{
		  W->d=WR->d*powf(2./(GAMMA+1.)-(GAMMA-1.)/((GAMMA+1.)*WR->a)*(WR->u-S),2./(GAMMA-1.));
		  W->u=2./(GAMMA+1.)*(-WR->a+(GAMMA-1.)/2.*WR->u+S);
		  W->p=WR->p*powf(2./(GAMMA+1.)-(GAMMA-1.)/((GAMMA+1.)*WR->a)*(WR->u-S),2.*GAMMA/(GAMMA-1.));
		  W->a=sqrtf(GAMMA*W->p/W->d);
		}
	    }
	}

      

    }

}


// ======================================================== Flux fonctions

void getflux_X(struct Utype *U, float *f)
{
  f[0]=U->du;
  f[1]=0.5*(3.-GAMMA)*U->du*U->du/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->dv*U->dv+U->dw*U->dw)/U->d;
  f[2]=U->du*U->dv/U->d;
  f[3]=U->du*U->dw/U->d;
  f[4]=GAMMA*U->du/U->d*U->E-0.5*(GAMMA-1.)*U->du/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}

// ---------------------------------------------------------------

void getflux_Y(struct Utype *U, float *f)
{
  f[0]=U->dv;
  f[1]=U->dv*U->du/U->d;
  f[2]=0.5*(3.-GAMMA)*U->dv*U->dv/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dw*U->dw)/U->d;
  f[3]=U->dv*U->dw/U->d;
  f[4]=GAMMA*U->dv/U->d*U->E-0.5*(GAMMA-1.)*U->dv/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}

// ---------------------------------------------------------------

void getflux_Z(struct Utype *U, float *f)
{
  f[0]=U->dw;
  f[1]=U->dw*U->du/U->d;
  f[2]=U->dw*U->dv/U->d;
  f[3]=0.5*(3.-GAMMA)*U->dw*U->dw/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dv*U->dv)/U->d;
  f[4]=GAMMA*U->dw/U->d*U->E-0.5*(GAMMA-1.)*U->dw/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}


//============================================================================
int hydrosolve(struct MULTIVECT *data, int level, int curcpu, int nread,int stride,float dx, float dt){

  int inei,icell,icellcoarse;
  int i;
  float temp;
  float res,restot=0.;
  int vnei[6],vcell[6];
  int vneic[6],vcellc[6];
  int idxnei;

  float FL[5],FR[5];
  float GL[5],GR[5];
  float HL[5],HR[5];
  float Smax;
  struct Wtype1D WRloc, WLloc;
  struct Utype Uold,Unew;
  struct Wtype Wold,Wnew;
  int idxL,idxR;
  float pstar,ustar;
  int n;
  struct Wtype1D Wtest;
  struct Wtype   Wtest3D;
  struct Utype   Utest3D;
  int flag;

  for(i=0;i<stride;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells
      flag=0;

      // we skip octs which do not belong to the current cpu
/* #ifdef WMPI */
/*       if(data->veccpu[i]!=curcpu){ */
/*       	data->vec_dnew[i+icell*stride]=data->vec_d[i+icell*stride]; */
/*       	data->vec_unew[i+icell*stride]=data->vec_u[i+icell*stride]; */
/*       	data->vec_vnew[i+icell*stride]=data->vec_v[i+icell*stride]; */
/*       	data->vec_wnew[i+icell*stride]=data->vec_w[i+icell*stride]; */
/*       	data->vec_pnew[i+icell*stride]=data->vec_p[i+icell*stride]; */
/*       	continue; */
/*       } */
/* #endif */

      temp=0.;
      getcellnei(icell, vnei, vcell); // we get the neighbors

      // Getting the original state ===========================
      idxR=i+icell*stride;

      Wold.d=data->vec_d[idxR];
      Wold.u=data->vec_u[idxR];
      Wold.v=data->vec_v[idxR];
      Wold.w=data->vec_w[idxR];
      Wold.p=data->vec_p[idxR];
      Wold.a=sqrtf(GAMMA*Wold.p/Wold.d);

      W2U(&Wold,&Uold); // primitive -> conservative

      // X DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems LEFT
      int vc;
      inei=0; // we go to the left
      if(vnei[inei]==6){
	idxL=i+vcell[inei]*stride;
      }
      else{
	idxL=data->vecnei[i+vnei[inei]*stride];
	vc=vcell[inei];
#ifdef TRANSXM
	if(idxL==i){
	  // we found a transmissive boundary
	  vc=vc+((vc&1)==0?1:-1); // flipping the cells
	}
#endif
	idxL=idxL+vc*stride;
      }
      idxR=i+icell*stride;


      // Switching to Split description
      WLloc.d=data->vec_d[idxL];
      WLloc.u=data->vec_u[idxL];
      WLloc.p=data->vec_p[idxL];
      WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d);

      WRloc.d=data->vec_d[idxR];
      WRloc.u=data->vec_u[idxR];
      WRloc.p=data->vec_p[idxR];
      WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d);
      //if(WLloc.d!=WRloc.d) printf("R== %e %e %e p=%e u=%e niter=%d\n", WLloc.d,WRloc.d,Wtest3D.d,pstar,ustar,n);

      // Riemann Solver
      pstar=findPressure(&WLloc,&WRloc,&n,&ustar);

      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.u=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      //if(WLloc.d!=WRloc.d) printf("L== %e %e %e p=%e u=%e niter=%d\n", WLloc.d,WRloc.d,Wtest3D.d,pstar,ustar,n);
      /* if((WLloc.d==1.0)&&(WRloc.d==0.125)){ */
      /* 	printf("L== %e %e %e p=%e u=%e niter=%d\n", WLloc.d,WRloc.d,Wtest3D.d,pstar,ustar,n); */
      /* 	flag=1; */
      /* } */

      // Passive advection
      if(ustar>0.)
	{
	  Wtest3D.v=data->vec_v[idxL];
	  Wtest3D.w=data->vec_w[idxL];
	}
      else
	{
	  Wtest3D.v=data->vec_v[idxR];
	  Wtest3D.w=data->vec_w[idxR];
	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes LEFT
      getflux_X(&Utest3D,FL);

      // --------- solving the Riemann Problems RIGHT

      inei=1; // we go to the right
      idxL=i+icell*stride;
      if(vnei[inei]==6){
	idxR=i+vcell[inei]*stride;
      }
      else{
	idxR=data->vecnei[i+vnei[inei]*stride];
	vc=vcell[inei];
#ifdef TRANSXP	
	if(idxR==i){
	  // we found a transmissive boundary
	  vc=vc+((vc&1)==0?1:-1); // flipping the cells
	}
#endif
	idxR=idxR+vc*stride;
      }

      // Switching to Split description
      WLloc.d=data->vec_d[idxL];
      WLloc.u=data->vec_u[idxL];
      WLloc.p=data->vec_p[idxL];
      WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d);

      WRloc.d=data->vec_d[idxR];
      WRloc.u=data->vec_u[idxR];
      WRloc.p=data->vec_p[idxR];
      WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d);

      // Riemann Solver
      pstar=findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);

      
      Wtest3D.d=Wtest.d;
      Wtest3D.u=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;

      //if(WLloc.d!=WRloc.d) printf("R== %e %e %e p=%e u=%e niter=%d\n", WLloc.d,WRloc.d,Wtest3D.d,pstar,ustar,n);
      
      // Passive advection
      if(ustar>0.)
	{
	  Wtest3D.v=data->vec_v[idxL];
	  Wtest3D.w=data->vec_w[idxL];
	}
      else
	{
	  Wtest3D.v=data->vec_v[idxR];
	  Wtest3D.w=data->vec_w[idxR];
	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes RIGHT
      getflux_X(&Utest3D,FR);

      /* if(((WLloc.d==1.)&&(WRloc.d==0.125))){ */
      /* 	//printf("L== %e %e %e p=%e u=%e niter=%d FL=%e FR=%e\n", WLloc.d,WRloc.d,Wtest3D.d,pstar,ustar,n,FL[0],FR[0]); */
      /* 	printf("L== %e %e %e p=%e u=%e niter=%d FL=%e FR=%e exp=%e\n dt=%e dx=%e", WLloc.d,WRloc.d,Wtest3D.d,pstar,ustar,n,FL[0],FR[0],1.+dt/dx*(FL[0]-FR[0]),dt,dx); */

      /* 	//flag=1; */
      /* 	abort(); */
      /* } */

      /* // Y DIRECTION ========================================================================= */
      
      /* // --------- solving the Riemann Problems FRONT */

      /* inei=2; // we go to the left */
      /* if(vnei[inei]==6){ */
      /* 	idxL=i+vcell[inei]*stride; */
      /* } */
      /* else{ */
      /* 	idxL=data->vecnei[i+vnei[inei]*stride]+vcell[inei]*stride; */
      /* } */
      /* idxR=i+icell*stride; */

      /* // Switching to Split description */
      /* WLloc.d=data->vec_d[idxL]; */
      /* WLloc.u=data->vec_v[idxL]; */
      /* WLloc.p=data->vec_p[idxL]; */
      /* WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d); */

      /* WRloc.d=data->vec_d[idxR]; */
      /* WRloc.u=data->vec_v[idxR]; */
      /* WRloc.p=data->vec_p[idxR]; */
      /* WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d); */

      /* // Riemann Solver */
      /* pstar=findPressure(&WLloc,&WRloc,&n,&ustar); */
      /* getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
      /* Wtest3D.d=Wtest.d; */
      /* Wtest3D.u=Wtest.u; */
      /* Wtest3D.p=Wtest.p; */
      /* Wtest3D.a=Wtest.a; */
      
      /* // Passive advection */
      /* if(ustar>0.) */
      /* 	{ */
      /* 	  Wtest3D.v=data->vec_u[idxL]; */
      /* 	  Wtest3D.w=data->vec_w[idxL]; */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  Wtest3D.v=data->vec_u[idxR]; */
      /* 	  Wtest3D.w=data->vec_w[idxR]; */
      /* 	} */
      
      /* W2U(&Wtest3D,&Utest3D); */
      
      /* // Getting the fluxes LEFT */
      /* getflux_Y(&Utest3D,GL); */

      /* // --------- solving the Riemann Problems BACK */

      /* inei=3; // we go to the right */
      /* idxL=i+icell*stride; */
      /* if(vnei[inei]==6){ */
      /* 	idxR=i+vcell[inei]*stride; */
      /* } */
      /* else{ */
      /* 	idxR=data->vecnei[i+vnei[inei]*stride]+vcell[inei]*stride; */
      /* } */

      /* // Switching to Split description */
      /* WLloc.d=data->vec_d[idxL]; */
      /* WLloc.u=data->vec_v[idxL]; */
      /* WLloc.p=data->vec_p[idxL]; */
      /* WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d); */

      /* WRloc.d=data->vec_d[idxR]; */
      /* WRloc.u=data->vec_v[idxR]; */
      /* WRloc.p=data->vec_p[idxR]; */
      /* WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d); */

      /* // Riemann Solver */
      /* pstar=findPressure(&WLloc,&WRloc,&n,&ustar); */
      /* getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
      /* Wtest3D.d=Wtest.d; */
      /* Wtest3D.u=Wtest.u; */
      /* Wtest3D.p=Wtest.p; */
      /* Wtest3D.a=Wtest.a; */
      
      /* // Passive advection */
      /* if(ustar>0.) */
      /* 	{ */
      /* 	  Wtest3D.v=data->vec_u[idxL]; */
      /* 	  Wtest3D.w=data->vec_w[idxL]; */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  Wtest3D.v=data->vec_u[idxR]; */
      /* 	  Wtest3D.w=data->vec_w[idxR]; */
      /* 	} */
      
      /* W2U(&Wtest3D,&Utest3D); */
      
      /* // Getting the fluxes RIGHT */
      /* getflux_Y(&Utest3D,GR); */

      /* // Z DIRECTION ========================================================================= */
      
      /* // --------- solving the Riemann Problems BOTTOM */

      /* inei=4; // we go to the left */
      /* if(vnei[inei]==6){ */
      /* 	idxL=i+vcell[inei]*stride; */
      /* } */
      /* else{ */
      /* 	idxL=data->vecnei[i+vnei[inei]*stride]+vcell[inei]*stride; */
      /* } */
      /* idxR=i+icell*stride; */

      /* // Switching to Split description */
      /* WLloc.d=data->vec_d[idxL]; */
      /* WLloc.u=data->vec_w[idxL]; */
      /* WLloc.p=data->vec_p[idxL]; */
      /* WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d); */

      /* WRloc.d=data->vec_d[idxR]; */
      /* WRloc.u=data->vec_w[idxR]; */
      /* WRloc.p=data->vec_p[idxR]; */
      /* WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d); */

      /* // Riemann Solver */
      /* pstar=findPressure(&WLloc,&WRloc,&n,&ustar); */
      /* getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
      /* Wtest3D.d=Wtest.d; */
      /* Wtest3D.u=Wtest.u; */
      /* Wtest3D.p=Wtest.p; */
      /* Wtest3D.a=Wtest.a; */
      
      /* // Passive advection */
      /* if(ustar>0.) */
      /* 	{ */
      /* 	  Wtest3D.v=data->vec_u[idxL]; */
      /* 	  Wtest3D.w=data->vec_v[idxL]; */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  Wtest3D.v=data->vec_u[idxR]; */
      /* 	  Wtest3D.w=data->vec_v[idxR]; */
      /* 	} */
      
      /* W2U(&Wtest3D,&Utest3D); */
      
      /* // Getting the fluxes LEFT */
      /* getflux_Z(&Utest3D,HL); */

      /* // --------- solving the Riemann Problems BACK */

      /* inei=5; // we go to the right */
      /* idxL=i+icell*stride; */
      /* if(vnei[inei]==6){ */
      /* 	idxR=i+vcell[inei]*stride; */
      /* } */
      /* else{ */
      /* 	idxR=data->vecnei[i+vnei[inei]*stride]+vcell[inei]*stride; */
      /* } */

      /* // Switching to Split description */
      /* WLloc.d=data->vec_d[idxL]; */
      /* WLloc.u=data->vec_w[idxL]; */
      /* WLloc.p=data->vec_p[idxL]; */
      /* WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d); */

      /* WRloc.d=data->vec_d[idxR]; */
      /* WRloc.u=data->vec_w[idxR]; */
      /* WRloc.p=data->vec_p[idxR]; */
      /* WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d); */

      /* // Riemann Solver */
      /* pstar=findPressure(&WLloc,&WRloc,&n,&ustar); */
      /* getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar); */
      
      /* Wtest3D.d=Wtest.d; */
      /* Wtest3D.u=Wtest.u; */
      /* Wtest3D.p=Wtest.p; */
      /* Wtest3D.a=Wtest.a; */
      
      /* // Passive advection */
      /* if(ustar>0.) */
      /* 	{ */
      /* 	  Wtest3D.v=data->vec_u[idxL]; */
      /* 	  Wtest3D.w=data->vec_v[idxL]; */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  Wtest3D.v=data->vec_u[idxR]; */
      /* 	  Wtest3D.w=data->vec_v[idxR]; */
      /* 	} */
      
      /* W2U(&Wtest3D,&Utest3D); */
      
      /* // Getting the fluxes RIGHT */
      /* getflux_Z(&Utest3D,HR); */

      // Updating the data ====================================================================
      // Unsplit scheme

      Unew.d  =Uold.d  +dt/dx*((FL[0]-FR[0]));//+(GL[0]-GR[0])+(HL[0]-HR[0]));
      Unew.du =Uold.du +dt/dx*((FL[1]-FR[1]));//+(GL[1]-GR[1])+(HL[1]-HR[1]));
      Unew.dv =Uold.dv +dt/dx*((FL[2]-FR[2]));//+(GL[2]-GR[2])+(HL[2]-HR[2]));
      Unew.dw =Uold.dw +dt/dx*((FL[3]-FR[3]));//+(GL[3]-GR[3])+(HL[3]-HR[3]));
      Unew.E  =Uold.E  +dt/dx*((FL[4]-FR[4]));//+(GL[4]-GR[4])+(HL[4]-HR[4]));
		  
      U2W(&Unew,&Wnew);
      
      //if(Wold.d==1.) printf("R== %e %e %e p=%e u=%e niter=%d FL=%e FR=%e dt/dx=%e\n", WLloc.d,WRloc.d,Wnew.d,pstar,ustar,n,FL[0],FR[0],dt/dx);

      idxR=i+icell*stride;
      
      // sending the data back to the vector
      data->vec_dnew[idxR]=Wnew.d;
      data->vec_unew[idxR]=Wnew.u;
      data->vec_vnew[idxR]=Wnew.v;
      data->vec_wnew[idxR]=Wnew.w;
      data->vec_pnew[idxR]=Wnew.p;

      // ready for the next cell
    }
    //ready for the next oct
  }

  // full update
  memcpy(data->vec_d,data->vec_dnew,sizeof(float)*stride*8);
  memcpy(data->vec_u,data->vec_unew,sizeof(float)*stride*8);
  memcpy(data->vec_v,data->vec_vnew,sizeof(float)*stride*8);
  memcpy(data->vec_w,data->vec_wnew,sizeof(float)*stride*8);
  memcpy(data->vec_p,data->vec_pnew,sizeof(float)*stride*8);

  return 0;
}


//============================================================================


float comptstep_hydro(int levelcoarse,int levelmax,struct OCT** firstoct, float fa, float fa2, struct CPUINFO* cpu, float tmax){
  
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  float dxcur;
  int icell;
  float aa;
  float va,vx,vy,vz;
  float dt;
  float Smax=0.;

  //Smax=fmaxf(Smax,sqrtf(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);

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
	      vx=curoct->cell[icell].u;
	      vy=curoct->cell[icell].v;
	      vz=curoct->cell[icell].w;
	      va=sqrtf(vx*vx+vy*vy+vz*vz);

	      aa=sqrtf(GAMMA*curoct->cell[icell].p/curoct->cell[icell].d);
	      Smax=fmaxf(Smax,va+aa);
	    }
	}while(nextoct!=NULL);

      dt=fminf(dxcur*CFL/(Smax*3.),dt);
    }

  // new tstep
 printf("original dt=%f chosen dt=%f\n",tmax,dt);


/* #ifdef WMPI */
/*   // reducing by taking the smallest time step */
/*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_FLOAT,MPI_MIN,cpu->comm); */
/* #endif   */

  return dt;
}


#endif
