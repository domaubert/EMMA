#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXITER 16
#define MAXJ 5

#include "prototypes.h"



/*
  =========================================================================================
  =========================================================================================

  Friedmann Models from Martel & Shapiro 1997 Supercomoving coordinates

  t_tilde(aexp=1.0)-t_tilde(aexp)=sqrt(omegam)/2.*integ(faexp_tilde,aexp,1.0,omegam,omegav);
  t      (aexp=1.0)-t      (aexp)=(1./H0        )*integ(faexp      ,aexp,1.0,omegam,omegav);

  =========================================================================================
  =========================================================================================


// NOTE ALL the calculations are performed in double precision here

 */



// =====================================================================================
double faexp_tilde(double a, double om, double ov){

  // For conversion between aexp and t_tilde

  return 1./sqrt((1.0-om-ov)*a*a*a*a+om*a*a*a+ov*a*a*a*a*a*a);
}

// =====================================================================================

double faexp(double a, double om, double ov){

  // For conversion between aexp and t

  return a*a/sqrt((1.0-om-ov)*a*a*a*a+om*a*a*a+ov*a*a*a*a*a*a);
}


// =====================================================================================

double integ_da_dt_tilde(double a, double b, double om, double ov,double tol)
{
  /* Romberg integration */

  double h,g0,fourj,gmax,error,g1;
  double g[MAXJ+2];
  int nint;
  int i,j,jmax,k;
  double res;
  h=0.5*(b-a);
  gmax=h*(faexp_tilde(a,om,ov)+faexp_tilde(b,om,ov));
  g[1]=gmax;
  nint=1;
  error=1e20;
  i=0;
  do{
    i++;
    if((i>MAXITER)||((i>5)&&(fabs(error)<tol))) break;
    g0=0.;
    for(k=1;k<=nint;k++) g0+=faexp_tilde(a+(k+k-1)*h,om,ov);
    g0=0.5*g[1]+h*g0;
    h=0.5*h;
    nint=nint+nint;
    jmax=(i<MAXJ?i:MAXJ);
    fourj=1.0;
    for(j=1;j<=jmax;j++){
      fourj=4.0*fourj;
      g1=g0+(g0-g[j])/(fourj-1.0);
      g[j]=g0;
      g0=g1;
    }

    if(fabs(g0)>tol){
      error=1.0-gmax/g0;
    }
    else{
      error=gmax;
    }
    gmax=g0;
    g[jmax+1]=g0;
  }while(1);

  res=g0;

  if((i>MAXITER)&&(fabs(error)>tol)){
    printf("ROMBINT FAILED TO CONVERGE integ=%e error=%e\n",res,error);
  }

  return res;

}

// =====================================================================================

double integ_da_dt(double a, double b, double om, double ov,double tol)
{
  /* Romberg integration */

  double h,g0,fourj,gmax,error,g1;
  double g[MAXJ+2];
  int nint;
  int i,j,jmax,k;
  double res;
  h=0.5*(b-a);
  gmax=h*(faexp(a,om,ov)+faexp(b,om,ov));
  g[1]=gmax;
  nint=1;
  error=1e20;
  i=0;
  do{
    i++;
    if((i>MAXITER)||((i>5)&&(fabs(error)<tol))) break;
    g0=0.;
    for(k=1;k<=nint;k++) g0+=faexp(a+(k+k-1)*h,om,ov);
    g0=0.5*g[1]+h*g0;
    h=0.5*h;
    nint=nint+nint;
    jmax=(i<MAXJ?i:MAXJ);
    fourj=1.0;
    for(j=1;j<=jmax;j++){
      fourj=4.0*fourj;
      g1=g0+(g0-g[j])/(fourj-1.0);
      g[j]=g0;
      g0=g1;
    }

    if(fabs(g0)>tol){
      error=1.0-gmax/g0;
    }
    else{
      error=gmax;
    }
    gmax=g0;
    g[jmax+1]=g0;
  }while(1);

  res=g0;

  if((i>MAXITER)&&(fabs(error)>tol)){
    printf("ROMBINT FAILED TO CONVERGE integ=%e error=%e\n",res,error);
  }

  return res;

}

// =====================================================================================

double interp_aexp(double ttilde,double *tab_aexp,double *tab_ttilde){

  int iloc=1;
  while(ttilde>tab_ttilde[iloc]){
    iloc++;
  }

  return tab_aexp[iloc-1]+(ttilde-tab_ttilde[iloc-1])*(tab_aexp[iloc]-tab_aexp[iloc-1])/(tab_ttilde[iloc]-tab_ttilde[iloc-1]);

}

// =====================================================================================


void compute_friedmann(double amin, double amax, int ntab, double omegam, double omegav, double *tab_aexp, double *tab_ttilde, double *tab_t)
{

  //FILE *fp;
  //  fp=fopen("data/fried.txt","w");
  int i;
  for(i=0;i<ntab;i++){
    //printf("%d\n",i);
    tab_aexp[i]=amin+(amax-amin)/(ntab-1)*i;
    tab_t[i]=-integ_da_dt(tab_aexp[i],1.0,omegam,omegav,1e-8); // in hubble time units
    tab_ttilde[i]=-0.5*sqrt(omegam)*integ_da_dt_tilde(tab_aexp[i],1.0,omegam,omegav,1e-8);
    //fprintf(fp,"%lf %lf %f\n",tab_aexp[i],tab_t[i],tab_ttilde[i]);
  }
  //fclose(fp);

}



// =====================================================================================
double dladt(double a, double omegam, double omegav)
{
  double eta;
  eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav);
  return a*eta;
}



// =====================================================================================
double ddplus(double a,double omegam, double omegav)
{
  double eta;
  double res;

  if(a==0.) return 0.;
  eta=sqrt(omegam/a+omegav*a*a+1.-omegam-omegav);
  res=2.5/(eta*eta*eta);
  return res;

}

// =====================================================================================

double integ_ddplus(double a, double b, double om, double ov,double tol)
{
  /* Romberg integration */

  double h,g0,fourj,gmax,error,g1;
  double g[MAXJ+2];
  int nint;
  int i,j,jmax,k;
  double res;
  h=0.5*(b-a);
  gmax=h*(ddplus(a,om,ov)+ddplus(b,om,ov));
  g[1]=gmax;
  nint=1;
  error=1e20;
  i=0;
  do{
    i++;
    if((i>MAXITER)||((i>5)&&(fabs(error)<tol))) break;
    g0=0.;
    for(k=1;k<=nint;k++) g0+=ddplus(a+(k+k-1)*h,om,ov);
    g0=0.5*g[1]+h*g0;
    h=0.5*h;
    nint=nint+nint;
    jmax=(i<MAXJ?i:MAXJ);
    fourj=1.0;
    for(j=1;j<=jmax;j++){
      fourj=4.0*fourj;
      g1=g0+(g0-g[j])/(fourj-1.0);
      g[j]=g0;
      g0=g1;
    }

    if(fabs(g0)>tol){
      error=1.0-gmax/g0;
    }
    else{
      error=gmax;
    }
    gmax=g0;
    g[jmax+1]=g0;
  }while(1);

  res=g0;

  if((i>MAXITER)&&(fabs(error)>tol)){
    printf("ROMBINT FAILED TO CONVERGE integ=%e error=%e\n",res,error);
  }

  return res;

}

// =====================================================================================
double dplus(double a, double omegam, double omegav)
{
  double eta;
  eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav);
  return eta*integ_ddplus(0.,a,omegam,omegav,1e-8)/a;
}

// =====================================================================================

double fomega(double a, double omegam, double omegav)
{
  double omegak,eta,res;

  if((omegam==1.0)&&(omegav==0.)){
    res=1.0;
    return res;
  }

  omegak=1.-omegam-omegav;
  eta=sqrt(omegam/a+omegav*a*a+omegak);
  res=(2.5/dplus(a,omegam,omegav)-1.5*omegam/a-omegak)/(eta*eta);

  return res;

}

// ========================================================================================

REAL a2t(struct RUNPARAMS *param, REAL a){
#ifdef TESTCOSMO
  REAL t=integ_da_dt(1e-8,a,param->cosmo->om,param->cosmo->ov,1e-8);
  return t/(param->cosmo->H0*1e3/1e6/PARSEC *365*24*3600);
#endif // TESTCOSMO
}
