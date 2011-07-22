#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>

#define NITERMAX 100
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

  float ptr,pts,ppv;
  float p,porg,dp;
  int i;
  float err;
  float unsurz=(2.0*GAMMA)/(GAMMA-1.0);
  float AL,AR,BL,BR,GL,GR;
  float pmin,pmax;
  int tag;

  // hybrid guess for pressure

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;


  ppv=fmaxf(ERRTOL,0.5*(WL->p+WR->p)-0.125*(WR->u-WL->u)*(WR->d+WL->d)*(WR->a+WL->a));
  ptr=fmaxf(ERRTOL,pow((WL->a+WR->a-0.5*(GAMMA-1)*(WR->u-WL->u))/(WL->a/pow(WL->p,1./unsurz)+WR->a/pow(WR->p,1./unsurz)),unsurz));
  
  GL=sqrt(AL/(ppv+BL));
  GR=sqrt(AR/(ppv+BR));

  pts=fmaxf(ERRTOL,(GL*WL->p+GR*WR->p-(WR->u-WL->u))/(GL+GR));

  pmin=fminf(WL->p,WR->p);
  pmax=fmaxf(WL->p,WR->p);

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
  //p=fmaxf(p,ERRTOL);

  *niter=0;
  for(i=0;i<NITERMAX;i++)
    {
      dp=froot(p,WL,WR,u)/frootprime(p,WL,WR);
      if((isnan(dp))||(p-dp)<0){
      	printf("froot=%e frootprime=%e\n",froot(p,WL,WR,u),frootprime(p,WL,WR));
      	abort();
      }

      porg=p;
      p=p-dp;
      //if(frootprime(p,WL,WR)==0) abort();//printf("p0=%e dp=%e p=%e fprime=%e\n",porg,dp,p,frootprime(p,WL,WR));
      err=2.*fabsf(p-porg)/(fabsf(p+porg));
      *niter=*niter+1;
      //if(p<=0) p=ERRTOL;
      if(err<ERRTOL) break;
    }

  if(i==NITERMAX){
    printf("DIVERGENCE p0=%e dp=%e p=%e fprime=%e err=%e\n",porg,dp,p,frootprime(p,WL,WR),err);
    abort();
  }

  if(p>6.4e-7){
    printf("MAX p0=%e dp=%e p=%e fprime=%e err=%e\n",porg,dp,p,frootprime(p,WL,WR),err);
    abort();
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


// =================================================================================================================

#ifdef SELFGRAV
void correct_grav_hydro(struct OCT *octstart, struct CPUINFO *cpu, float dt)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int icell;

  struct Utype Uold,Unew;
  struct Wtype Wold,Wnew;
  float fx,fy,fz;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      // filling the values
      for(icell=0;icell<8;icell++){

	Wold.d=curoct->cell[icell].d;
	Wold.u=curoct->cell[icell].u;
	Wold.v=curoct->cell[icell].v;
	Wold.w=curoct->cell[icell].w;
	Wold.p=curoct->cell[icell].p;
	Wold.a=sqrtf(GAMMA*Wold.p/Wold.d);

	W2U(&Wold,&Uold); // primitive -> conservative

	// we store the gravitational force in the new fields
	fx=curoct->cell[icell].fx;
	fy=curoct->cell[icell].fy;
	fz=curoct->cell[icell].fz;

	
	// implicit update
	Unew.d  =Uold.d  ;
	Unew.du =Uold.du +Unew.d*fx*dt*0.5;
	Unew.dv =Uold.dv +Unew.d*fy*dt*0.5;
	Unew.dw =Uold.dw +Unew.d*fz*dt*0.5;
	Unew.E  =Uold.E  +(Unew.du*fx+Unew.dv*fy+Unew.dw*fz)*dt*0.5;
      
	U2W(&Unew,&Wnew);

	curoct->cell[icell].u=Wnew.u;
	curoct->cell[icell].v=Wnew.v;
	curoct->cell[icell].w=Wnew.w;
	curoct->cell[icell].p=Wnew.p;

      }
    }while(nextoct!=NULL);
  }
}
#endif

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

  for(icell=0;icell<8;icell++){ // we scan the cells
      getcellnei(icell, vnei, vcell); // we get the neighbors
      for(i=0;i<stride;i++){ // we scan the octs

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

      /* WRloc.d=data->vec_d[idxR]; */
      /* WRloc.u=data->vec_u[idxR]; */
      /* WRloc.p=data->vec_p[idxR]; */
      /* WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d); */

      WRloc.d=Wold.d;
      WRloc.u=Wold.u;//data->vec_u[idxR];
      WRloc.p=Wold.p;//data->vec_p[idxR];
      WRloc.a=Wold.a;//sqrtf(GAMMA*WRloc.p/WRloc.d);

      // Riemann Solver
      pstar=findPressure(&WLloc,&WRloc,&n,&ustar);

      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.u=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
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
      /* WLloc.d=data->vec_d[idxL]; */
      /* WLloc.u=data->vec_u[idxL]; */
      /* WLloc.p=data->vec_p[idxL]; */
      /* WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d); */

      WLloc.d=Wold.d;
      WLloc.u=Wold.u;
      WLloc.p=Wold.p;
      WLloc.a=Wold.a;

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

      // Y DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems FRONT

      inei=2; // we go to the left
      if(vnei[inei]==6){
      	idxL=i+vcell[inei]*stride;
      }
      else{
      	idxL=data->vecnei[i+vnei[inei]*stride];
	vc=vcell[inei];
#ifdef TRANSYM
	if(idxL==i){
	  // we found a transmissive boundary
	  vc=icell;
	}
#endif
	idxL=idxL+vc*stride;
      }
      idxR=i+icell*stride;

      // Switching to Split description
      WLloc.d=data->vec_d[idxL];
      WLloc.u=data->vec_v[idxL];
      WLloc.p=data->vec_p[idxL];
      WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d);

      /* WRloc.d=data->vec_d[idxR]; */
      /* WRloc.u=data->vec_v[idxR]; */
      /* WRloc.p=data->vec_p[idxR]; */
      /* WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d); */

      WRloc.d=Wold.d;
      WRloc.u=Wold.v;
      WRloc.p=Wold.p;
      WRloc.a=Wold.a;

      // Riemann Solver
      pstar=findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.v=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar>0.)
      	{
      	  Wtest3D.u=data->vec_u[idxL];
      	  Wtest3D.w=data->vec_w[idxL];
      	}
      else
      	{
      	  Wtest3D.u=data->vec_u[idxR];
      	  Wtest3D.w=data->vec_w[idxR];
      	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes LEFT
      getflux_Y(&Utest3D,GL);

      //if(WLloc.d!=WRloc.d) printf("%e %e %e %e %e %e %e\n",WLloc.d,WRloc.d,GL[0],WLloc.u,WRloc.u,pstar,ustar);

      // --------- solving the Riemann Problems BACK

      inei=3; // we go to the right
      idxL=i+icell*stride;
      if(vnei[inei]==6){
      	idxR=i+vcell[inei]*stride;
      }
      else{
      	idxR=data->vecnei[i+vnei[inei]*stride];
	vc=vcell[inei];
#ifdef TRANSYP
	if(idxR==i) vc=icell;
#endif
	idxR=idxR+vc*stride;
	
      }

      // Switching to Split description
      /* WLloc.d=data->vec_d[idxL]; */
      /* WLloc.u=data->vec_v[idxL]; */
      /* WLloc.p=data->vec_p[idxL]; */
      /* WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d); */


      WLloc.d=Wold.d;
      WLloc.u=Wold.v;
      WLloc.p=Wold.p;
      WLloc.a=Wold.a;

      WRloc.d=data->vec_d[idxR];
      WRloc.u=data->vec_v[idxR];
      WRloc.p=data->vec_p[idxR];
      WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d);

      // Riemann Solver
      pstar=findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.v=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar>0.)
      	{
      	  Wtest3D.u=data->vec_u[idxL];
      	  Wtest3D.w=data->vec_w[idxL];
      	}
      else
      	{
      	  Wtest3D.u=data->vec_u[idxR];
      	  Wtest3D.w=data->vec_w[idxR];
      	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes RIGHT
      getflux_Y(&Utest3D,GR);

      // Z DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems BOTTOM

      inei=4; // we go to the left
      if(vnei[inei]==6){
      	idxL=i+vcell[inei]*stride;
      }
      else{
      	idxL=data->vecnei[i+vnei[inei]*stride];
	vc=vcell[inei];
#ifdef TRANSZM
	if(idxL==i) vc=icell;
#endif
	idxL=idxL+vc*stride;
      }
      idxR=i+icell*stride;

      // Switching to Split description
      WLloc.d=data->vec_d[idxL];
      WLloc.u=data->vec_w[idxL];
      WLloc.p=data->vec_p[idxL];
      WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d);

      WRloc.d=Wold.d;
      WRloc.u=Wold.w;
      WRloc.p=Wold.p;
      WRloc.a=Wold.a;

      /* WRloc.d=data->vec_d[idxR]; */
      /* WRloc.u=data->vec_w[idxR]; */
      /* WRloc.p=data->vec_p[idxR]; */
      /* WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d); */

      // Riemann Solver
      pstar=findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.w=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar>0.)
      	{
      	  Wtest3D.u=data->vec_u[idxL];
      	  Wtest3D.v=data->vec_v[idxL];
      	}
      else
      	{
      	  Wtest3D.u=data->vec_u[idxR];
      	  Wtest3D.v=data->vec_v[idxR];
      	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes LEFT
      getflux_Z(&Utest3D,HL);

      // --------- solving the Riemann Problems BACK

      inei=5; // we go to the right
      idxL=i+icell*stride;
      if(vnei[inei]==6){
      	idxR=i+vcell[inei]*stride;
      }
      else{
      	idxR=data->vecnei[i+vnei[inei]*stride];
	vc=vcell[inei];
#ifdef TRANSZP
	if(idxR==i) vc=icell;
#endif
	idxR=idxR+vc*stride;
      }

      // Switching to Split description
      /* WLloc.d=data->vec_d[idxL]; */
      /* WLloc.u=data->vec_w[idxL]; */
      /* WLloc.p=data->vec_p[idxL]; */
      /* WLloc.a=sqrtf(GAMMA*WLloc.p/WLloc.d); */

      WLloc.d=Wold.d;
      WLloc.u=Wold.w;
      WLloc.p=Wold.p;
      WLloc.a=Wold.a;

      WRloc.d=data->vec_d[idxR];
      WRloc.u=data->vec_w[idxR];
      WRloc.p=data->vec_p[idxR];
      WRloc.a=sqrtf(GAMMA*WRloc.p/WRloc.d);

      // Riemann Solver
      pstar=findPressure(&WLloc,&WRloc,&n,&ustar);
      getW(&Wtest,0., &WLloc, &WRloc, pstar, ustar);
      
      Wtest3D.d=Wtest.d;
      Wtest3D.w=Wtest.u;
      Wtest3D.p=Wtest.p;
      Wtest3D.a=Wtest.a;
      
      // Passive advection
      if(ustar>0.)
      	{
      	  Wtest3D.u=data->vec_u[idxL];
      	  Wtest3D.v=data->vec_v[idxL];
      	}
      else
      	{
      	  Wtest3D.u=data->vec_u[idxR];
      	  Wtest3D.v=data->vec_v[idxR];
      	}
      
      W2U(&Wtest3D,&Utest3D);
      
      // Getting the fluxes RIGHT
      getflux_Z(&Utest3D,HR);

      // Updating the data ====================================================================
      // Unsplit scheme

      idxR=i+icell*stride;

      float fx=0.,fy=0.,fz=0.; // @ this stage the forces are stored in the updated vectors
     
#ifdef SELFGRAV
      fx=-data->vec_unew[idxR];
      fy=-data->vec_vnew[idxR];
      fz=-data->vec_wnew[idxR];
#endif
      //if(Uold.d>0.5) printf("g=%e org=%e\n",fx*Uold.d*dt,Uold.du);
      Unew.d  =Uold.d  +dt/dx*((FL[0]-FR[0])+(GL[0]-GR[0])+(HL[0]-HR[0]));
      Unew.du =Uold.du +dt/dx*((FL[1]-FR[1])+(GL[1]-GR[1])+(HL[1]-HR[1]))+Uold.d*fx*dt*0.5;
      Unew.dv =Uold.dv +dt/dx*((FL[2]-FR[2])+(GL[2]-GR[2])+(HL[2]-HR[2]))+Uold.d*fy*dt*0.5;
      Unew.dw =Uold.dw +dt/dx*((FL[3]-FR[3])+(GL[3]-GR[3])+(HL[3]-HR[3]))+Uold.d*fz*dt*0.5;
      Unew.E  =Uold.E  +dt/dx*((FL[4]-FR[4])+(GL[4]-GR[4])+(HL[4]-HR[4]))+(Uold.du*fx+Uold.dv*fy+Uold.dw*fz)*dt*0.5;
      
      U2W(&Unew,&Wnew);
      
      
      // sending the data back to the vector
      data->vec_dnew[idxR]=Wnew.d;
      data->vec_unew[idxR]=Wnew.u;
      data->vec_vnew[idxR]=Wnew.v;
      data->vec_wnew[idxR]=Wnew.w;
      data->vec_pnew[idxR]=(Wnew.p<0?Wold.p:Wnew.p);

      /* if(Wnew.p<0){ */
      /* 	printf("negative pressure in hydro\n"); */
      /* 	abort(); */
      /* } */

      // ready for the next cell
    }
    //ready for the next oct
  }

  /* // full update */
  /* memcpy(data->vec_d,data->vec_dnew,sizeof(float)*stride*8); */
  /* memcpy(data->vec_u,data->vec_unew,sizeof(float)*stride*8); */
  /* memcpy(data->vec_v,data->vec_vnew,sizeof(float)*stride*8); */
  /* memcpy(data->vec_w,data->vec_wnew,sizeof(float)*stride*8); */
  /* memcpy(data->vec_p,data->vec_pnew,sizeof(float)*stride*8); */

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
  // printf("original dt=%f chosen dt=%f\n",tmax,dt);


/* #ifdef WMPI */
/*   // reducing by taking the smallest time step */
/*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_FLOAT,MPI_MIN,cpu->comm); */
/* #endif   */

  return dt;
}

float comptstep_ff(int levelcoarse,int levelmax,struct OCT** firstoct, float aexp, struct CPUINFO* cpu, float tmax){
  
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  float dxcur;
  int icell;
  float dtloc;
  float dt;

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
	      dtloc=0.1*sqrt(2.*M_PI/(3.*curoct->cell[icell].density*aexp));
	      dt=fminf(dt,dtloc);
	    }
	}while(nextoct!=NULL);
    }

  // new tstep
  // printf("original dt=%f chosen dt=%f\n",tmax,dt);


/* #ifdef WMPI */
/*   // reducing by taking the smallest time step */
/*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_FLOAT,MPI_MIN,cpu->comm); */
/* #endif   */

  return dt;
}


#endif
