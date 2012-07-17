#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define NITERMAX 30
#define ERRTOL 1e-10
#define GAMMA (1.4)
#define NX 128
#define CFL 0.8
#define STEPMAX 250
#define REAL double
#define REAL2 double

struct Wtype{
  REAL d;   // density
  REAL u;   // velocity
  REAL v;   // velocity
  REAL w;   // velocity
  REAL p;   // pressure
  REAL a;   // sound speed
};

struct Utype{
  REAL d;    // density
  REAL du;   // momentum
  REAL dv;   // momentum
  REAL dw;   // momentum
  REAL E;    // Energy
};


struct Wtype1D{
  REAL2 d;   // density
  REAL2 u;   // velocity
  REAL2 p;   // pressure
  REAL2 a;   // sound speed
};

struct Utype1D{
  REAL d;    // density
  REAL du;   // momentum
  REAL E;    // Energy
};


void U2W(struct Utype *U, struct Wtype *W)
{
  W->d=U->d;
  W->u=U->du/U->d;
  W->v=U->dv/U->d;
  W->w=U->dw/U->d;
  W->p=(GAMMA-1.)*(U->E-((U->du)*(U->du)+(U->dv)*(U->dv)+(U->dw)*(U->dw))/(U->d)*0.5);
  W->a=sqrtf(GAMMA*W->p/W->d);
}

void W2U(struct Wtype *W, struct Utype *U)
{
  U->d=W->d;
  U->du=W->d*W->u;
  U->dv=W->d*W->v;
  U->dw=W->d*W->w;
  U->E=W->d*(0.5*(W->u*W->u+W->v*W->v+W->w*W->w)+W->p/((GAMMA-1.)*W->d));
}

REAL2 frootprime(REAL2 p, struct Wtype1D *WL, struct Wtype1D *WR)
{
  
  REAL fL,fR;
  REAL AL,AR,BL,BR;
  REAL Deltau;

  AL=2./((GAMMA+1.)*WL->d);
  AR=2./((GAMMA+1.)*WR->d);
  
  BL=(GAMMA-1.)/(GAMMA+1.)*WL->p;
  BR=(GAMMA-1.)/(GAMMA+1.)*WR->p;

  fL=(p>WL->p?sqrtf(AL/(BL+p))*(1.-(p-WL->p)/(2.*(BL+p))):powf(p/WL->p,-(GAMMA+1)/(2.*GAMMA))/(WL->d*WL->a));
  fR=(p>WR->p?sqrtf(AR/(BR+p))*(1.-(p-WR->p)/(2.*(BR+p))):powf(p/WR->p,-(GAMMA+1)/(2.*GAMMA))/(WR->d*WR->a));

  return fL+fR;
}

REAL2 froot(REAL2 p, struct Wtype1D *WL, struct Wtype1D *WR, REAL2 *u)
{
  
  REAL2 fL,fR;
  REAL2 AL,AR,BL,BR;
  REAL2 Deltau;

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


REAL2 findPressure(struct Wtype1D *WL, struct Wtype1D *WR, int *niter, REAL *u)
{
  REAL2 p,porg,dp;
  int i;
  REAL2 err;
  REAL2 u2;

  p=0.5*(WL->p+WR->p);
  p=0.5*(WL->p+WR->p)-0.125*(WR->u-WL->u)*(WR->d+WL->d)*(WR->a+WL->a);
  //p=fmaxf(p,ERRTOL);
  //printf("p=%e\n",p);
  *niter=0;
  for(i=0;i<NITERMAX;i++)
    {
      dp=froot(p,WL,WR,&u2)/frootprime(p,WL,WR);
      porg=p;
      if(dp>p) dp=p*0.5;
      p=p-dp;
      err=2.*fabs(p-porg)/(fabs(p+porg));
      //printf("p0=%e dp=%e p=%e fprime=%e err=%e froot=%e\n",porg,dp,p,frootprime(p,WL,WR),err,froot(p,WL,WR,&u2));
      *niter=*niter+1;
      
      if(err<ERRTOL) break;
    }
  //printf("p0=%e dp=%e p=%e fprime=%e err=%e\n",porg,dp,p,frootprime(p,WL,WR),err);
  froot(p,WL,WR,&u2); // last calculation to get u;
  *u=(float)u2;
  
  return p;
}


void getW(struct Wtype1D *W, REAL S, struct Wtype1D *WL, struct Wtype1D *WR, REAL pstar, REAL ustar)
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
      	  	  if(W->p<0) abort();

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


void getflux_X(struct Utype *U, REAL *f)
{
  f[0]=U->du;
  f[1]=0.5*(3.-GAMMA)*U->du*U->du/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->dv*U->dv+U->dw*U->dw)/U->d;
  f[2]=U->du*U->dv/U->d;
  f[3]=U->du*U->dw/U->d;
  f[4]=GAMMA*U->du/U->d*U->E-0.5*(GAMMA-1.)*U->du/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}


void getflux_Y(struct Utype *U, REAL *f)
{
  f[0]=U->dv;
  f[1]=U->dv*U->du/U->d;
  f[2]=0.5*(3.-GAMMA)*U->dv*U->dv/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dw*U->dw)/U->d;
  f[3]=U->dv*U->dw/U->d;
  f[4]=GAMMA*U->dv/U->d*U->E-0.5*(GAMMA-1.)*U->dv/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}


void getflux_Z(struct Utype *U, REAL *f)
{
  f[0]=U->dw;
  f[1]=U->dw*U->du/U->d;
  f[2]=U->dw*U->dv/U->d;
  f[3]=0.5*(3.-GAMMA)*U->dw*U->dw/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dv*U->dv)/U->d;
  f[4]=GAMMA*U->dw/U->d*U->E-0.5*(GAMMA-1.)*U->dw/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}

//==============================================================
//==============================================================
//==============================================================
//==============================================================
//==============================================================

int main()
{
  
  struct Wtype1D WL, WR;
  REAL ustar;
  REAL pstar;
  int n;
  int i,j,k;
  int idx,idxL,idxR;
  REAL dx=1./NX;
  REAL t;
  struct Wtype *Warray;
  struct Utype *Uarray;
  

  struct Wtype1D Wtest;
  struct Wtype   Wtest3D;
  struct Utype   Utest3D;

  FILE *fp;
  float X0;
  float TMAX;

  //TEST 1
  REAL x[NX+2];
  for(i=0;i<NX;i++)
    {
      x[i+1]=i*dx;
    }
  
  x[0]=x[1];
  x[NX+1]=x[NX];

  //***************************************************************************
  // THEORETICAL RESULT
  //***************************************************************************

  float dl,ul,pl,dr,ur,pr;

  fp=fopen("Input.riemann","r");
  fscanf(fp,"%f %f %f",&dl,&ul,&pl); 
  fscanf(fp,"%f %f %f",&dr,&ur,&pr); 
  fscanf(fp,"%f %f",&X0,&TMAX); 
  fclose(fp);

  WL.d=dl;
  WL.u=ul;
  WL.p=pl;

  WR.d=dr;
  WR.u=ur;
  WR.p=pr;


  printf("dL=%e uL=%e pL=%e\n",WL.d,WL.u,WL.p);
  printf("dR=%e uR=%e pR=%e\n",WR.d,WR.u,WR.p);
  printf("X0=%e TMAX=%e\n",X0,TMAX);


  WL.a=sqrt(GAMMA*WL.p/WL.d);
  WR.a=sqrt(GAMMA*WR.p/WR.d);
  
  pstar=(float)findPressure(&WL,&WR,&n,&ustar);
  printf("TEST 1 niter=%d pstart=%e ustar=%e\n",n,pstar,ustar);

  fp=fopen("test1.dat","w");
  t=TMAX;
  for(i=0;i<NX;i++)
    {
      //printf("i=%d\n",i);
      getW(&Wtest, (x[i]-X0)/t, &WL, &WR, pstar, ustar);
      //printf("%e %e %e %e\n",x[i],Wtest.d,Wtest.u,Wtest.p);
      fprintf(fp,"%e %e %e %e\n",x[i],Wtest.d,Wtest.u,Wtest.p);
    }
  fclose(fp);

  getW(&Wtest, 0., &WL, &WR, pstar, ustar);
  printf("%e %e %e %e\n",0.,Wtest.d,Wtest.u,Wtest.p);

  return 0;
}



  /* /\* // TEST 1 *\/ */
  /* WL.d=1.; */
  /* WL.u  =0.75; */
  /* WL.p  =1.; */

  /* WR.d=0.125; */
  /* WR.u  =0.; */
  /* WR.p  =0.1; */
  /* X0=0.3; */
  /* TMAX=0.2; */

  /* // TEST 2 */
  /* WL.d=1.; */
  /* WL.u  =-2.; */
  /* WL.p  =0.4; */

  /* WR.d=1.; */
  /* WR.u  =2.; */
  /* WR.p  =0.4; */
  /* X0=0.5; */
  /* TMAX=0.15: */

  // TEST 3
  /* WL.d=1.; */
  /* WL.u  =0.; */
  /* WL.v  =0.; */
  /* WL.w  =0.; */
  /* WL.p  =1000.; */

  /* WR.d=1.; */
  /* WR.u  =0.; */
  /* WL.v  =0.; */
  /* WL.w  =0.; */
  /* WR.p  =0.01; */
  /* X0=0.5; */
  /* TMAX=0.012; */
