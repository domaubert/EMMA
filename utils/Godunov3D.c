#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define NITERMAX 10
#define ERRTOL 1e-6
#define GAMMA (1.4)
#define NX 32
#define CFL 0.8
#define STEPMAX 5


struct Wtype{
  float d;   // density
  float u;   // velocity
  float v;   // velocity
  float w;   // velocity
  float p;   // pressure
  float a;   // sound speed
};

struct Utype{
  float d;    // density
  float du;   // momentum
  float dv;   // momentum
  float dw;   // momentum
  float E;    // Energy
};


struct Wtype1D{
  float d;   // density
  float u;   // velocity
  float p;   // pressure
  float a;   // sound speed
};

struct Utype1D{
  float d;    // density
  float du;   // momentum
  float E;    // Energy
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


void getflux_X(struct Utype *U, float *f)
{
  f[0]=U->du;
  f[1]=0.5*(3.-GAMMA)*U->du*U->du/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->dv*U->dv+U->dw*U->dw)/U->d;
  f[2]=U->du*U->dv/U->d;
  f[3]=U->du*U->dw/U->d;
  f[4]=GAMMA*U->du/U->d*U->E-0.5*(GAMMA-1.)*U->du/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}


void getflux_Y(struct Utype *U, float *f)
{
  f[0]=U->dv;
  f[1]=U->dv*U->du/U->d;
  f[2]=0.5*(3.-GAMMA)*U->dv*U->dv/U->d+(GAMMA-1.)*U->E-0.5*(GAMMA-1.)*(U->du*U->du+U->dw*U->dw)/U->d;
  f[3]=U->dv*U->dw/U->d;
  f[4]=GAMMA*U->dv/U->d*U->E-0.5*(GAMMA-1.)*U->dv/(U->d*U->d)*(U->du*U->du+U->dv*U->dv+U->dw*U->dw);
}


void getflux_Z(struct Utype *U, float *f)
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

void BCZMtrans(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(j=0;j<NX;j++)
    {
      for(i=0;i<NX;i++)
	{
	  idx= (i+1)+(j+1)*(NX+2)+(  1)*(NX+2)*(NX+2);
	  idxb=(i+1)+(j+1)*(NX+2)+(  0)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCZMper(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(j=0;j<NX;j++)
    {
      for(i=0;i<NX;i++)
	{
	  idx= (i+1)+(j+1)*(NX+2)+(  NX)*(NX+2)*(NX+2);
	  idxb=(i+1)+(j+1)*(NX+2)+(  0)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}



void BCZPtrans(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(j=0;j<NX;j++)
    {
      for(i=0;i<NX;i++)
	{
	  idx= (i+1)+(j+1)*(NX+2)+(NX  )*(NX+2)*(NX+2);
	  idxb=(i+1)+(j+1)*(NX+2)+(NX+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCZPper(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(j=0;j<NX;j++)
    {
      for(i=0;i<NX;i++)
	{
	  idx= (i+1)+(j+1)*(NX+2)+(   1)*(NX+2)*(NX+2);
	  idxb=(i+1)+(j+1)*(NX+2)+(NX+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCYMtrans(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(k=0;k<NX;k++)
    {
      for(i=0;i<NX;i++)
	{
	  idx= (i+1)+(  1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  idxb=(i+1)+(  0)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCYMper(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(k=0;k<NX;k++)
    {
      for(i=0;i<NX;i++)
	{
	  idx= (i+1)+( NX)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  idxb=(i+1)+(  0)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCYPtrans(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(k=0;k<NX;k++)
    {
      for(i=0;i<NX;i++)
	{
	  idx= (i+1)+(NX  )*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  idxb=(i+1)+(NX+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCYPper(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(k=0;k<NX;k++)
    {
      for(i=0;i<NX;i++)
	{
	  idx= (i+1)+(   1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  idxb=(i+1)+(NX+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}


void BCXMtrans(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(k=0;k<NX;k++)
    {
      for(j=0;j<NX;j++)
	{
	  idx= (1  )+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  idxb=(0  )+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCXMper(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(k=0;k<NX;k++)
    {
      for(j=0;j<NX;j++)
	{
	  idx= (NX )+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  idxb=(0  )+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCXPtrans(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(k=0;k<NX;k++)
    {
      for(j=0;j<NX;j++)
	{
	  idx= (NX  )+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  idxb=(NX+1)+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

void BCXPper(struct Wtype *W, struct Utype *U)
{
  int i,j,k;
  int idx,idxb;
  for(k=0;k<NX;k++)
    {
      for(j=0;j<NX;j++)
	{
	  idx= (1   )+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  idxb=(NX+1)+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	  memcpy(&(W[idxb]),&(W[idx]),sizeof(struct Wtype));
	  memcpy(&(U[idxb]),&(U[idx]),sizeof(struct Utype));
	}
    }
}

//==============================================================
//==============================================================
//==============================================================

int main()
{
  
  struct Wtype WL, WR;
  float pstar,ustar;
  int n;
  int i,j,k;
  int idx,idxL,idxR;
  float dx=1./NX;
  float t;
  struct Wtype *Warray;
  struct Utype *Uarray;
  
  Warray=calloc((NX+2)*(NX+2)*(NX+2),sizeof(struct Wtype));
  Uarray=calloc((NX+2)*(NX+2)*(NX+2),sizeof(struct Utype));

  struct Wtype1D Wtest;
  struct Wtype   Wtest3D;
  struct Utype   Utest3D;

  FILE *fp;
  float X0;
  float TMAX;

  //TEST 1
  float x[NX+2];
  for(i=0;i<NX;i++)
    {
      x[i+1]=i*dx;
    }
  
  //***************************************************************************
  // THEORETICAL RESULT
  //***************************************************************************

  /* // TEST 1 */
  WL.d=1.;
  WL.u  =0.;
  WL.v=0.;
  WL.w=0.;
  WL.p  =1.;

  WR.d=0.125;
  WR.u  =0.;
  WR.v=0.;
  WR.w=0.;
  WR.p  =0.1;
  X0=0.3;
  TMAX=0.25;

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

  WL.a=sqrtf(GAMMA*WL.p/WL.d);
  WR.a=sqrtf(GAMMA*WR.p/WR.d);
  
  /* pstar=findPressure(&WL,&WR,&n,&ustar); */
  /* printf("TEST 1 niter=%d pstart=%e ustar=%e\n",n,pstar,ustar); */

  /* fp=fopen("test1.dat","w"); */
  /* t=TMAX; */
  /* for(i=0;i<NX;i++) */
  /*   { */
  /*     //printf("i=%d\n",i); */
  /*     getW(&Wtest, (x[i]-X0)/t, &WL, &WR, pstar, ustar); */
  /*     //      printf("%e %e %e %e\n",x[i],Wtest.d,Wtest.u,Wtest.p); */
  /*     fprintf(fp,"%e %e %e %e\n",x[i],Wtest.d,Wtest.u,Wtest.p); */
  /*   } */
  /* fclose(fp); */

  //***************************************************************************
  // NUMERICAL RESULT
  //***************************************************************************
  

  // ICs   ****************************************************

  for(k=0;k<NX;k++)
    {
      for(j=0;j<NX;j++)
	{
	  for(i=0;i<NX;i++)
	    {
	      idx =i  +1+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);

	      if(x[i+1]<=X0)
		{
		  Warray[idx].d=WL.d;
		  Warray[idx].u=WL.u;
		  Warray[idx].v=WL.v;
		  Warray[idx].w=WL.w;
		  Warray[idx].p=WL.p;
		  Warray[idx].a=WL.a;
		  W2U(&Warray[idx],&Uarray[idx]);
		}
	      else
		{
		  Warray[idx].d=WR.d;
		  Warray[idx].u=WR.u;
		  Warray[idx].v=WR.v;
		  Warray[idx].w=WR.w;
		  Warray[idx].p=WR.p;
		  Warray[idx].a=WR.a;
		  W2U(&Warray[idx],&Uarray[idx]);
		}
	    }
	}
    }

  // boundaries


  BCXMtrans(Warray,Uarray);
  BCXPtrans(Warray,Uarray);

  BCYMper(Warray,Uarray);
  BCYPper(Warray,Uarray);

  BCZMper(Warray,Uarray);
  BCZPper(Warray,Uarray);

  float dt;
  int count;
  float FL[5],FR[5];
  float GL[5],GR[5];
  float HL[5],HR[5];
  float Smax;
  struct Wtype1D WRloc, WLloc,dtest;
  int d;
  t=0.;
  count=1;
  while(t<TMAX)
    {
      // Computing timestep
      Smax=0.;
      for(i=0;i<(NX+2)*(NX+2)*(NX+2);i++)
	{
	  Smax=fmaxf(Smax,sqrtf(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);
	}

      dt=dx*CFL/Smax/3.;

      printf("Timestep #%d t=%e dt=%e\n",count,t,dt);

      // setting the time step
      for(k=0;k<NX;k++)
	{
	  for(j=0;j<NX;j++)
	    {
	      for(i=0;i<NX;i++)
		{

		  d=0;
		  idx =i  +1+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
		  //printf("%d\n",idx);
		  
		  // X DIRECTION =========================================================================
		  // solving the Riemann Problems LEFT

		  idxL=i-1+1+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
		  idxR=i  +1+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);

		  // Switching to Split description
		  WLloc.d=Warray[idxL].d;
		  WLloc.u=Warray[idxL].u;
		  WLloc.p=Warray[idxL].p;
		  WLloc.a=Warray[idxL].a;
		  
		  WRloc.d=Warray[idxR].d;
		  WRloc.u=Warray[idxR].u;
		  WRloc.p=Warray[idxR].p;
		  WRloc.a=Warray[idxR].a;

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
		      Wtest3D.v=Warray[idxL].v;
		      Wtest3D.w=Warray[idxL].w;
		    }
		  else
		    {
		      Wtest3D.v=Warray[idxR].v;
		      Wtest3D.w=Warray[idxR].w;
		    }

		  W2U(&Wtest3D,&Utest3D);
		  if((WLloc.d==1.)&&(WRloc.d==1.)){
		    d=1.;
		  }

		  // Getting the fluxes LEFT
		  getflux_X(&Utest3D,FL);

		  // solving the Riemann Problems RIGHT

		  idxL=i  +1+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
		  idxR=i+1+1+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);

		  // Switching to Split description
		  WLloc.d=Warray[idxL].d;
		  WLloc.u=Warray[idxL].u;
		  WLloc.p=Warray[idxL].p;
		  WLloc.a=Warray[idxL].a;
		  
		  WRloc.d=Warray[idxR].d;
		  WRloc.u=Warray[idxR].u;
		  WRloc.p=Warray[idxR].p;
		  WRloc.a=Warray[idxR].a;

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
		      Wtest3D.v=Warray[idxL].v;
		      Wtest3D.w=Warray[idxL].w;
		    }
		  else
		    {
		      Wtest3D.v=Warray[idxR].v;
		      Wtest3D.w=Warray[idxR].w;
		    }


		  W2U(&Wtest3D,&Utest3D);
		  
		  // Getting the fluxes RIGHT
		  getflux_X(&Utest3D,FR);

		  /* if(((WLloc.d==1.)&&(WRloc.d==0.125))&&(d==1)){ */
		  /*   //printf("L== %e %e %e p=%e u=%e niter=%d FL=%e FR=%e exp=%e\n", WLloc.d,WRloc.d,Wtest3D.d,pstar,ustar,n,FL[0],FR[0],1.+dt/dx*(FL[0]-FR[0])); */
		  /*   printf("L== %e %e %e p=%e u=%e niter=%d FL=%e FR=%e exp=%e\n dt=%e dx=%e", WLloc.d,WRloc.d,Wtest3D.d,pstar,ustar,n,FL[0],FR[0],1.+dt/dx*(FL[0]-FR[0]),dt,dx); */
		  /*   //flag=1; */
		  /*   abort(); */
		  /* } */


		  // Y DIRECTION =========================================================================
		  // solving the Riemann Problems LEFT

		  idxL=i  +1+(j-1+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
		  idxR=i  +1+(j  +1)*(NX+2)+(k+1)*(NX+2)*(NX+2);

		  // Switching to Split description
		  WLloc.d=Warray[idxL].d;
		  WLloc.u=Warray[idxL].v;
		  WLloc.p=Warray[idxL].p;
		  WLloc.a=Warray[idxL].a;
		  
		  WRloc.d=Warray[idxR].d;
		  WRloc.u=Warray[idxR].v;
		  WRloc.p=Warray[idxR].p;
		  WRloc.a=Warray[idxR].a;

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
		      Wtest3D.u=Warray[idxL].u;
		      Wtest3D.w=Warray[idxL].w;
		    }
		  else
		    {
		      Wtest3D.u=Warray[idxR].u;
		      Wtest3D.w=Warray[idxR].w;
		    }


		  W2U(&Wtest3D,&Utest3D);
		  
		  // Getting the fluxes LEFT
		  getflux_Y(&Utest3D,GL);

		  // solving the Riemann Problems RIGHT

		  idxL=i  +1+(j  +1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
		  idxR=i  +1+(j+1+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);

		  // Switching to Split description
		  WLloc.d=Warray[idxL].d;
		  WLloc.u=Warray[idxL].v;
		  WLloc.p=Warray[idxL].p;
		  WLloc.a=Warray[idxL].a;
		  
		  WRloc.d=Warray[idxR].d;
		  WRloc.u=Warray[idxR].v;
		  WRloc.p=Warray[idxR].p;
		  WRloc.a=Warray[idxR].a;

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
		      Wtest3D.u=Warray[idxL].u;
		      Wtest3D.w=Warray[idxL].w;
		    }
		  else
		    {
		      Wtest3D.u=Warray[idxR].u;
		      Wtest3D.w=Warray[idxR].w;
		    }


		  W2U(&Wtest3D,&Utest3D);
		  
		  // Getting the fluxes RIGHT
		  getflux_Y(&Utest3D,GR);
		  

		  // Z DIRECTION =========================================================================
		  // solving the Riemann Problems LEFT

		  idxL=i  +1+(j  +1)*(NX+2)+(k-1+1)*(NX+2)*(NX+2);
		  idxR=i  +1+(j  +1)*(NX+2)+(k  +1)*(NX+2)*(NX+2);

		  // Switching to Split description
		  WLloc.d=Warray[idxL].d;
		  WLloc.u=Warray[idxL].w;
		  WLloc.p=Warray[idxL].p;
		  WLloc.a=Warray[idxL].a;
		  
		  WRloc.d=Warray[idxR].d;
		  WRloc.u=Warray[idxR].w;
		  WRloc.p=Warray[idxR].p;
		  WRloc.a=Warray[idxR].a;

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
		      Wtest3D.u=Warray[idxL].u;
		      Wtest3D.v=Warray[idxL].v;
		    }
		  else
		    {
		      Wtest3D.u=Warray[idxR].u;
		      Wtest3D.v=Warray[idxR].v;
		    }


		  W2U(&Wtest3D,&Utest3D);
		  
		  // Getting the fluxes LEFT
		  getflux_Z(&Utest3D,HL);

		  // solving the Riemann Problems RIGHT

		  idxL=i  +1+(j  +1)*(NX+2)+(k  +1)*(NX+2)*(NX+2);
		  idxR=i  +1+(j  +1)*(NX+2)+(k+1+1)*(NX+2)*(NX+2);

		  // Switching to Split description
		  WLloc.d=Warray[idxL].d;
		  WLloc.u=Warray[idxL].w;
		  WLloc.p=Warray[idxL].p;
		  WLloc.a=Warray[idxL].a;
		  
		  WRloc.d=Warray[idxR].d;
		  WRloc.u=Warray[idxR].w;
		  WRloc.p=Warray[idxR].p;
		  WRloc.a=Warray[idxR].a;

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
		      Wtest3D.u=Warray[idxL].u;
		      Wtest3D.v=Warray[idxL].v;
		    }
		  else
		    {
		      Wtest3D.u=Warray[idxR].u;
		      Wtest3D.v=Warray[idxR].v;
		    }


		  W2U(&Wtest3D,&Utest3D);
		  
		  // Getting the fluxes RIGHT
		  getflux_Z(&Utest3D,HR);
		  

		  // Updating the data ====================================================================
		  
 		  Uarray[idx].d  =Uarray[idx].d  +dt/dx*(FL[0]-FR[0])+dt/dx*(GL[0]-GR[0])+dt/dx*(HL[0]-HR[0]);
		  Uarray[idx].du =Uarray[idx].du +dt/dx*(FL[1]-FR[1])+dt/dx*(GL[1]-GR[1])+dt/dx*(HL[1]-HR[1]);
		  Uarray[idx].dv =Uarray[idx].dv +dt/dx*(FL[2]-FR[2])+dt/dx*(GL[2]-GR[2])+dt/dx*(HL[2]-HR[2]);
		  Uarray[idx].dw =Uarray[idx].dw +dt/dx*(FL[3]-FR[3])+dt/dx*(GL[3]-GR[3])+dt/dx*(HL[3]-HR[3]);
		  Uarray[idx].E  =Uarray[idx].E  +dt/dx*(FL[4]-FR[4])+dt/dx*(GL[4]-GR[4])+dt/dx*(HL[4]-HR[4]);
		  
		}
	    }
	}

      printf("cells ok\n");

      //memcpy(&Uarray[0],&Uarray[1],3*sizeof(float));
      //      memcpy(&Uarray[NX],&Uarray[NX-1],3*sizeof(float));

      // setting boundaries 
      BCXMtrans(Warray,Uarray);
      BCXPtrans(Warray,Uarray);
      
      BCYMper(Warray,Uarray);
      BCYPper(Warray,Uarray);
      
      BCZMper(Warray,Uarray);
      BCZPper(Warray,Uarray);

      printf("Bound ok\n");

      // preparing next data (primitive variables) + timestep
      for(i=0;i<(NX+2)*(NX+2)*(NX+2);i++)
	{
	  U2W(&Uarray[i],&Warray[i]);
	}

      // update time
      t+=dt;
      count++;
      if(count>STEPMAX) break;
    }


  puts("DONE");

  fp=fopen("test1.num","wb");
  int nn=NX;
  float buff[(NX+2)*(NX+2)*(NX+2)];
  fwrite(&nn,sizeof(int),1,fp);
  for(k=0;k<NX;k++)
    {
      for(j=0;j<NX;j++)
	{
	  for(i=0;i<NX;i++)
	    {
	      
	      idx =i  +1+(j+1)*(NX+2)+(k+1)*(NX+2)*(NX+2);
	      buff[idx]=Warray[idx].d;
	    }
	}
    }

  fwrite(buff,sizeof(float),(NX+2)*(NX+2)*(NX+2),fp);
  fclose(fp);

  return 0;
}
