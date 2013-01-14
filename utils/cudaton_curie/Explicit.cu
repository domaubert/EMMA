#include "params.h"
#include "common.h"
#include "bnd.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Explicit.h"
#include "GPU.h"
#include "bnd.h"
#include "Atomic.h"

#define ONE 0.9999999f
#define ZERO 0.0000001f
#define NCELL4 ((NCELLZ+NBOUND2)*(NCELLY+NBOUND2)*(NCELLX+NBOUND2))
#define TINY 1e-26
#define FLAGSOURCE 5.
#define KBOLTZ 1.3806e-23
#define EPSCOOL 0.0001

#define REAL float

//**********************************************************************************
//**********************************************************************************

__device__ float cuCompute_FaceLF(float fl, float fr, float ul, float ur)
{
  return (fl+fr-ur+ul)*0.5f;
}


//**********************************************************************************

#ifdef SDISCRETE

__global__ void cuAddSource(float *cuegy_new, float *cuflx_new, float *cusrc0, int *cusrc0pos,float dt, float dx, int nsource, float aexp, float c)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  int i,j,k;
  
  i=cusrc0pos[tx+bx*blockDim.x+0*nsource];
  j=cusrc0pos[tx+bx*blockDim.x+1*nsource];
  k=cusrc0pos[tx+bx*blockDim.x+2*nsource];
  
  float s=cusrc0[tx+bx*blockDim.x]; // WARNING : Sources should be considered as densities !!!!
  //float s=cusrc0[tx+bx*blockDim.x]/dx/dx/dx;

#ifdef SSHARP // Sharp point-like sources
  int idx;
  idx=(i+NBOUND)+(j+NBOUND)*(NCELLX+NBOUND2)+(k+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  cuegy_new[idx]+=s*dt; // 2 FLOPS
#else
  int idx;
  float sloc;
  float eloc;
  // center 1
  sloc=s*dt*0.125;

  idx=(i+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;

  // face centers 6
  sloc=s*dt*0.0625;

  idx=(i+1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=c*eloc;
  cuflx_new[idx+1*NCELL4]=0.;
  cuflx_new[idx+2*NCELL4]=0.;

  idx=(i+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=0.;
  cuflx_new[idx+1*NCELL4]=c*eloc;
  cuflx_new[idx+2*NCELL4]=0.;


  idx=(i+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=0.;
  cuflx_new[idx+1*NCELL4]=0.;
  cuflx_new[idx+2*NCELL4]=c*eloc;


  idx=(i-1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc;
  cuflx_new[idx+1*NCELL4]=0.;
  cuflx_new[idx+2*NCELL4]=0.;


  idx=(i+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=0.;
  cuflx_new[idx+1*NCELL4]=-c*eloc;
  cuflx_new[idx+2*NCELL4]=0.;

  idx=(i+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;

  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=0.;
  cuflx_new[idx+1*NCELL4]=0.;
  cuflx_new[idx+2*NCELL4]=-c*eloc;

  // edge centers 12
  sloc=s*dt*0.03125;

  idx=(i-1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc*rsqrtf(2.);
  cuflx_new[idx+1*NCELL4]=-c*eloc*rsqrtf(2.);
  cuflx_new[idx+2*NCELL4]=0.;

  idx=(i-1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc*rsqrtf(2.);
  cuflx_new[idx+1*NCELL4]=+c*eloc*rsqrtf(2.);
  cuflx_new[idx+2*NCELL4]=0.;

  idx=(i+1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=+c*eloc*rsqrtf(2.);
  cuflx_new[idx+1*NCELL4]=-c*eloc*rsqrtf(2.);
  cuflx_new[idx+2*NCELL4]=0.;

  idx=(i+1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=+c*eloc*rsqrtf(2.);
  cuflx_new[idx+1*NCELL4]=+c*eloc*rsqrtf(2.);
  cuflx_new[idx+2*NCELL4]=0.;

  idx=(i-1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc*rsqrtf(2.);
  cuflx_new[idx+1*NCELL4]=0.;
  cuflx_new[idx+2*NCELL4]=-c*eloc*rsqrtf(2.);

  idx=(i-1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc*rsqrtf(2.);
  cuflx_new[idx+1*NCELL4]=0.;
  cuflx_new[idx+2*NCELL4]=+c*eloc*rsqrtf(2.);

  idx=(i+1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=+c*eloc*rsqrtf(2.);
  cuflx_new[idx+1*NCELL4]=0.;
  cuflx_new[idx+2*NCELL4]=-c*eloc*rsqrtf(2.);

  idx=(i+1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=+c*eloc*rsqrtf(2.);
  cuflx_new[idx+1*NCELL4]=0.;
  cuflx_new[idx+2*NCELL4]=+c*eloc*rsqrtf(2.);

  idx=(i+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=0.;
  cuflx_new[idx+1*NCELL4]=-c*eloc*rsqrtf(2.);
  cuflx_new[idx+2*NCELL4]=-c*eloc*rsqrtf(2.);

  idx=(i+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=0.;
  cuflx_new[idx+1*NCELL4]=+c*eloc*rsqrtf(2.);
  cuflx_new[idx+2*NCELL4]=-c*eloc*rsqrtf(2.);

  idx=(i+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=0.;
  cuflx_new[idx+1*NCELL4]=-c*eloc*rsqrtf(2.);
  cuflx_new[idx+2*NCELL4]=+c*eloc*rsqrtf(2.);

  idx=(i+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=0.;
  cuflx_new[idx+1*NCELL4]=+c*eloc*rsqrtf(2.);
  cuflx_new[idx+2*NCELL4]=+c*eloc*rsqrtf(2.);


  // corners
  sloc=s*dt*0.015625;

  idx=(i+1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=c*eloc*rsqrtf(3.);
  cuflx_new[idx+1*NCELL4]=c*eloc*rsqrtf(3.);
  cuflx_new[idx+2*NCELL4]=c*eloc*rsqrtf(3.)
;
  idx=(i-1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc*rsqrtf(3.);
  cuflx_new[idx+1*NCELL4]=c*eloc*rsqrtf(3.);
  cuflx_new[idx+2*NCELL4]=c*eloc*rsqrtf(3.);

  idx=(i+1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=c*eloc*rsqrtf(3.);
  cuflx_new[idx+1*NCELL4]=-c*eloc*rsqrtf(3.);
  cuflx_new[idx+2*NCELL4]=c*eloc*rsqrtf(3.);

  idx=(i+1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=c*eloc*rsqrtf(3.);
  cuflx_new[idx+1*NCELL4]=c*eloc*rsqrtf(3.);
  cuflx_new[idx+2*NCELL4]=-c*eloc*rsqrtf(3.);

  idx=(i-1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc*rsqrtf(3.);
  cuflx_new[idx+1*NCELL4]=-c*eloc*rsqrtf(3.);
  cuflx_new[idx+2*NCELL4]=c*eloc*rsqrtf(3.);

  idx=(i+1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=c*eloc*rsqrtf(3.);
  cuflx_new[idx+1*NCELL4]=-c*eloc*rsqrtf(3.);
  cuflx_new[idx+2*NCELL4]=-c*eloc*rsqrtf(3.);

  idx=(i-1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc*rsqrtf(3.);
  cuflx_new[idx+1*NCELL4]=c*eloc*rsqrtf(3.);
  cuflx_new[idx+2*NCELL4]=-c*eloc*rsqrtf(3.);

  idx=(i-1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  eloc=cuegy_new[idx]+sloc;
  cuegy_new[idx]=eloc;
  cuflx_new[idx+0*NCELL4]=-c*eloc*rsqrtf(3.);
  cuflx_new[idx+1*NCELL4]=-c*eloc*rsqrtf(3.);
  cuflx_new[idx+2*NCELL4]=-c*eloc*rsqrtf(3.);


#endif


}


__global__ void cuSpotSource(float *cuxion, int *cusrc0pos, int nsource)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  int i,j,k;
  
  i=cusrc0pos[tx+bx*blockDim.x+0*nsource];
  j=cusrc0pos[tx+bx*blockDim.x+1*nsource];
  k=cusrc0pos[tx+bx*blockDim.x+2*nsource];
  

#ifdef SSHARP // Sharp point-like sources
  int idx;
  idx=(i+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  cuxion[idx]=FLAGSOURCE;
#else
  int idx;

  idx=(i+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  cuxion[idx]=FLAGSOURCE;
  
  // face centers 6

  idx=(i+1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
  cuxion[idx]=FLAGSOURCE;

  idx=(i+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;


  idx=(i+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;


  idx=(i-1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;


  idx=(i+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 

  cuxion[idx]=FLAGSOURCE;

  // edge centers 12

  idx=(i-1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i-1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i-1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i-1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+1+NBOUND)+(j+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  // corners
  idx=(i+1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;
;
  idx=(i-1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i-1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k+1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i+1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i-1+NBOUND)+(j+1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
 
  cuxion[idx]=FLAGSOURCE;

  idx=(i-1+NBOUND)+(j-1+NBOUND)*(NCELLS+NBOUND2)+(k-1+NBOUND)*(NCELLS+NBOUND2)*(NCELLS+NBOUND2);
   
  cuxion[idx]=FLAGSOURCE;


#endif


}

#endif

//**********************************************************************************
__global__ void cuComputeELF(float *cuegy, float *cuflx, float *cuegy_new, float c, float dx, float dt, int iter, float aexp, float egy_min)
{
  int tx=threadIdx.x+NBOUND;
  int bx=blockIdx.x +NBOUND;
  int by=blockIdx.y +NBOUND;

  REAL um1,up1,fm1,fp1,u0;

  // REMINDER LF flux : (fl+fr-ur+ul)*0.5f;
  //  f_icjcks_p =cuCompute_FaceLF(f[2+idx*3],f[2+idxp*3],c*e[idx],c*e[idxp]);

  REAL res;
  REAL dtsurdx=dt/dx;

  // Divergence along Z

  int baseidu=(tx)+bx*(NCELLX+NBOUND2);
  um1=cuegy[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];
  u0 =cuegy[baseidu+(by  )*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];
  up1=cuegy[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];
  
  int baseidf=2*NCELL4+(tx)+bx*(NCELLX+NBOUND2);
  fm1=cuflx[baseidf+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];
  fp1=cuflx[baseidf+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];

  res=u0-0.5*((fp1-fm1)+c*(2*u0-um1-up1))*dtsurdx; // 10 FLOP

  // Divergence along Y

  baseidu=(tx)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  um1=cuegy[baseidu+(bx-1)*(NCELLX+NBOUND2)];
  u0 =cuegy[baseidu+(bx  )*(NCELLX+NBOUND2)];
  up1=cuegy[baseidu+(bx+1)*(NCELLX+NBOUND2)];
  
  baseidf=1*NCELL4+(tx)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  fm1=cuflx[baseidf+(bx-1)*(NCELLX+NBOUND2)];
  fp1=cuflx[baseidf+(bx+1)*(NCELLX+NBOUND2)];
  
  res=res-0.5*((fp1-fm1)+c*(2*u0-um1-up1))*dtsurdx; // 10 FLOP
  

  //Divergence along X

  __shared__ float u[NCELLX+NBOUND2],f[NCELLX+NBOUND2];

  baseidu=bx*(NCELLX+NBOUND2)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  baseidf=0*NCELL4+bx*(NCELLX+NBOUND2)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  u[tx]=cuegy [baseidu+tx];
  f[tx]=cuflx [baseidf+tx];


  
  if(tx-NBOUND==0) 
    {
      u[NBOUND-1]=cuegy[baseidu+tx-1];
      f[NBOUND-1]=cuflx[baseidf+tx-1];
    }

  if(tx-NBOUND==blockDim.x-1) 
    {
      u[NCELLX+NBOUND]=cuegy[baseidu+tx+1];
      f[NCELLX+NBOUND]=cuflx[baseidf+tx+1];
    }

  __syncthreads();
  
  res=res-0.5*((f[tx+1]-f[tx-1])+c*(2*u[tx]-u[tx+1]-u[tx-1]))*dtsurdx; // 10 FLOP

  //if(res<0) printf("res=%e fp1=%e fm1=%e on idf=%d ide=%d\n",res,f[tx+1],f[tx-1],baseidf+tx,baseidu+tx);

#ifndef SDISCRETE
  cuegy_new[baseidu+tx]=fmaxf(res,egy_min);
#else
  cuegy_new[baseidu+tx]=fmaxf(res,egy_min);
#endif
}





//**********************************************************************************

__device__ float Eddington(float fx, float fy, float fz, float ee, float c,int i,int j)
{

  float c2e=ee*c*c; // 2 flop
  float ff=0.;
  float arg,chi,res=0.;
  float n[3];
#ifdef ISOTROP
  
  if(i==j) res=1./3.;

#else
  n[0]=0.;n[1]=0.;n[2]=0.;

  if(ee>0)
    {
      ff=sqrtf(fx*fx+fy*fy+fz*fz); // 6 flop
      if(ff>0)
	{
	  n[0]=fx/ff; 
	  n[1]=fy/ff;
	  n[2]=fz/ff; 
	}
      ff=ff/(c*ee); // 2flop
    }
  
  arg=fmaxf(4.-3.*ff*ff,0.); // 4 flop
  chi=(3.+4.*ff*ff)/(5.+2.*sqrtf(arg)); // 7 flops

  if(i==j) res=(1.-chi)/2.*c2e; // 1 flops on average
  arg=(3.*chi-1.)/2.*c2e;
  res+=arg*n[i]*n[j];
#endif

  return res;
}




//**********************************************************************************
__global__ void cuComputeF_TOTAL_LF(float *cuflx, float *cuflx_new, float c, float dx, float dt, int iter, float *cuegy, float aexp)
{
  int tx=threadIdx.x+NBOUND;
  int bx=blockIdx.x +NBOUND;
  int by=blockIdx.y +NBOUND;

  float fm1,fp1;

  // REMINDER LF flux : (fl+fr-ur+ul)*0.5f;
  //  f_icjcks_p =cuCompute_FaceLF(f[2+idx*3],f[2+idxp*3],c*e[idx],c*e[idxp]);

  float resfx, resfy, resfz;

  __shared__ float u[(NCELLX+NBOUND2)*3],fp[(NCELLX+NBOUND2)*3],fm[(NCELLX+NBOUND2)*3],ep[(NCELLX+NBOUND2)],em[(NCELLX+NBOUND2)];

  //================================================ Z DIRECTION =============================================
  
  int baseidu=0*NCELL4+(tx)+bx*(NCELLX+NBOUND2);

  u[0*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // FX local cell
  u[1*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4]; // FX local cell
  u[2*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+2*NCELL4]; // FX local cell

  ep[tx]=cuegy[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // E Cell+1
  em[tx]=cuegy[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // E Cell+1

  fm[0*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // E Cell+1
  fm[1*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4]; // E Cell+1
  fm[2*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4*2]; // E Cell+1


  fp[0*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // E Cell+1
  fp[1*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4]; // E Cell+1
  fp[2*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4*2]; // E Cell+1


  __syncthreads();

  // FX Divergence along Z
  
  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,0,2);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,0,2);

  resfx=u[tx+0*(NCELLX+NBOUND2)]-0.5*((fp1-fm1)+c*(2*u[tx+0*(NCELLX+NBOUND2)]-fm[tx+0*(NCELLX+NBOUND2)]-fp[tx+0*(NCELLX+NBOUND2)]))/dx*dt; 
  

 // FY Divergence along Z


  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,1,2);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,1,2);
  resfy=u[tx+1*(NCELLX+NBOUND2)]-0.5*((fp1-fm1)+c*(2*u[tx+1*(NCELLX+NBOUND2)]-fm[tx+1*(NCELLX+NBOUND2)]-fp[tx+1*(NCELLX+NBOUND2)]))/dx*dt; 
  

  // FZ Divergence along Z

  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,2,2);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,2,2);
  resfz=u[tx+2*(NCELLX+NBOUND2)]-0.5*((fp1-fm1)+c*(2*u[tx+2*(NCELLX+NBOUND2)]-fm[tx+2*(NCELLX+NBOUND2)]-fp[tx+2*(NCELLX+NBOUND2)]))/dx*dt; 
  
  __syncthreads();


  //================================================ Y DIRECTION =============================================
  
  baseidu=0*NCELL4+(tx)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  u[0*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx)*(NCELLX+NBOUND2)]; // FX local cell
  u[1*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx)*(NCELLX+NBOUND2)+NCELL4]; // FX local cell
  u[2*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx)*(NCELLX+NBOUND2)+2*NCELL4]; // FX local cell

  ep[tx]=cuegy[baseidu+(bx+1)*(NCELLX+NBOUND2)]; // E Cell+1
  em[tx]=cuegy[baseidu+(bx-1)*(NCELLX+NBOUND2)]; // E Cell+1

  fm[0*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx-1)*(NCELLX+NBOUND2)]; // E Cell+1
  fm[1*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx-1)*(NCELLX+NBOUND2)+NCELL4]; // E Cell+1
  fm[2*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx-1)*(NCELLX+NBOUND2)+NCELL4*2]; // E Cell+1

  fp[0*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx+1)*(NCELLX+NBOUND2)]; // E Cell+1
  fp[1*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx+1)*(NCELLX+NBOUND2)+NCELL4]; // E Cell+1
  fp[2*(NCELLX+NBOUND2)+tx]=cuflx[baseidu+(bx+1)*(NCELLX+NBOUND2)+NCELL4*2]; // E Cell+1

  __syncthreads();

  // FX Divergence along Y
  
  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,0,1);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,0,1);
  resfx=resfx-0.5*((fp1-fm1)+c*(2*u[tx+0*(NCELLX+NBOUND2)]-fm[tx+0*(NCELLX+NBOUND2)]-fp[tx+0*(NCELLX+NBOUND2)]))/dx*dt; 

 // FY Divergence along Y

  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,1,1);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,1,1);
  resfy=resfy-0.5*((fp1-fm1)+c*(2*u[tx+1*(NCELLX+NBOUND2)]-fm[tx+1*(NCELLX+NBOUND2)]-fp[tx+1*(NCELLX+NBOUND2)]))/dx*dt; 
  
  // FZ Divergence along Y

  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,2,1);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,2,1);
  resfz=resfz-0.5*((fp1-fm1)+c*(2*u[tx+2*(NCELLX+NBOUND2)]-fm[tx+2*(NCELLX+NBOUND2)]-fp[tx+2*(NCELLX+NBOUND2)]))/dx*dt; 
  


  __syncthreads();


  //================================================ X DIRECTION =============================================

  baseidu=0*NCELL4+bx*(NCELLX+NBOUND2)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  
  u[0*(NCELLX+NBOUND2)+tx]=cuflx[tx+baseidu]; // FX local cell
  u[1*(NCELLX+NBOUND2)+tx]=cuflx[tx+baseidu+NCELL4]; // FX local cell
  u[2*(NCELLX+NBOUND2)+tx]=cuflx[tx+baseidu+2*NCELL4]; // FX local cell
  ep[tx]=cuegy[tx+baseidu]; // E Cell+1

  if(tx-NBOUND==0) 
    {
      u[NBOUND-1+0*(NCELLX+NBOUND2)]=cuflx[baseidu+tx-1];
      u[NBOUND-1+1*(NCELLX+NBOUND2)]=cuflx[baseidu+tx-1+NCELL4];
      u[NBOUND-1+2*(NCELLX+NBOUND2)]=cuflx[baseidu+tx-1+2*NCELL4];
      ep[NBOUND-1]=cuegy[tx-1+baseidu]; 
    }

  if(tx-NBOUND==blockDim.x-1) 
    {
      u[NCELLX+NBOUND+0*(NCELLX+NBOUND2)]=cuflx[baseidu+tx+1];
      u[NCELLX+NBOUND+1*(NCELLX+NBOUND2)]=cuflx[baseidu+tx+1+NCELL4];
      u[NCELLX+NBOUND+2*(NCELLX+NBOUND2)]=cuflx[baseidu+tx+1+2*NCELL4];
      ep[NCELLX+NBOUND]=cuegy[tx+1+baseidu]; 
    }

  __syncthreads();


  // FX Divergence along X

  fp1=Eddington(u[0*(NCELLX+NBOUND2)+tx+1],u[1*(NCELLX+NBOUND2)+tx+1],u[2*(NCELLX+NBOUND2)+tx+1],ep[tx+1],c,0,0);
  fm1=Eddington(u[0*(NCELLX+NBOUND2)+tx-1],u[1*(NCELLX+NBOUND2)+tx-1],u[2*(NCELLX+NBOUND2)+tx-1],ep[tx-1],c,0,0);
  resfx=resfx-0.5*((fp1-fm1)+c*(2*u[tx+0*(NCELLX+NBOUND2)]-u[tx+1+0*(NCELLX+NBOUND2)]-u[tx-1+0*(NCELLX+NBOUND2)]))/dx*dt;
  

  // FY Divergence along X

  fp1=Eddington(u[0*(NCELLX+NBOUND2)+tx+1],u[1*(NCELLX+NBOUND2)+tx+1],u[2*(NCELLX+NBOUND2)+tx+1],ep[tx+1],c,1,0);
  fm1=Eddington(u[0*(NCELLX+NBOUND2)+tx-1],u[1*(NCELLX+NBOUND2)+tx-1],u[2*(NCELLX+NBOUND2)+tx-1],ep[tx-1],c,1,0);
  resfy=resfy-0.5*((fp1-fm1)+c*(2*u[tx+1*(NCELLX+NBOUND2)]-u[tx+1+1*(NCELLX+NBOUND2)]-u[tx-1+1*(NCELLX+NBOUND2)]))/dx*dt;
  

  // FZ Divergence along X

  fp1=Eddington(u[0*(NCELLX+NBOUND2)+tx+1],u[1*(NCELLX+NBOUND2)+tx+1],u[2*(NCELLX+NBOUND2)+tx+1],ep[tx+1],c,2,0);
  fm1=Eddington(u[0*(NCELLX+NBOUND2)+tx-1],u[1*(NCELLX+NBOUND2)+tx-1],u[2*(NCELLX+NBOUND2)+tx-1],ep[tx-1],c,2,0);
  resfz=resfz-0.5*((fp1-fm1)+c*(2*u[tx+2*(NCELLX+NBOUND2)]-u[tx+1+2*(NCELLX+NBOUND2)]-u[tx-1+2*(NCELLX+NBOUND2)]))/dx*dt;
  

  cuflx_new[baseidu+tx]=resfx;
  cuflx_new[baseidu+tx+NCELL4]=resfy;
  cuflx_new[baseidu+tx+2*NCELL4]=resfz;

}



//**********************************************************************************
__device__ float cufindroot3_2(float a,float b,float c,float d,float xorg)
{
  
  int i;
  float f,df,x;
  x=xorg;

  for(i=0;i<=10;i++)
    {
      f=a*x*x*x+b*x*x+c*x+d;
      df=3*a*x*x+2*b*x+c;
      if(fabsf(f/(df*x))<0.00001) break;
      x=x-f/df;
    }
  
  if(x>ONE) x=ONE;
  if(x<ZERO) x=ZERO;
  return x;
  //  if((x<=ONE)&&(x>=ZERO)) return x;
  //return xorg;

}
__device__ float cufindroot3(float a,float b,float c,float d,float xorg)
{
  float Q,R;
  float A,B,th;
  float x1,x2,x3;
  float x0;

  Q=((b/a)*(b/a)-3.*(c/a))/9.;
  R=(2*(b/a)*(b/a)*(b/a)-9*b*c/a/a+27*d/a)/54.;
  
  if(R*R<Q*Q*Q)
    {
      th=acosf(R/sqrtf(Q*Q*Q));
      x1=-2.*sqrtf(Q)*cosf(th/3.)-b/3./a;
      x2=-2.*sqrtf(Q)*cosf((th+2*M_PI)/3.)-b/3./a;
      x3=-2.*sqrtf(Q)*cosf((th-2*M_PI)/3.)-b/3./a;
      if((x1>=0)&&(x1<=1))
	{
	  x0=x1;
	}
      else if((x2>=0)&&(x2<=1))
	{
	  x0=x2;
	}
      else if((x3>=0)&&(x3<=1))
	{
	  x0=x3;
	}
      else
	{
	  x0=xorg;
	}
      
    }
  else
    {

      A=-copysignf(1.,R)*powf(fabsf(R)+sqrtf(R*R-Q*Q*Q),1./3.);
      B=0.;
      if(A!=0) B=Q/A;
      x0=(A+B)-b/3./a;
      if((x0>1)||(x0<0))
	{
#ifdef DEVICE_EMULATION
	  puts("ERROR in root finder, anormal ion fraction");
	  abort();
#endif
	}
    }
  
  //  if(idx==154869){printf("a=%e b=%e c=%e d=%e\n",a,b,c,d); printf("R=%e Q=%e x0=%e\n",R,Q,x0);}

  //x0=x3;
  return x0;
}


//**********************************************************************************
//**********************************************************************************

__device__ float cufindrootacc(float a,float b,float c)
{
  float q=-0.5*(b+copysign(1.,b)*sqrt(b*b-4*a*c));
  float x1=q/a;
  float x2=c/q;
  float x0=x1;
  if((x2<=ONE)&&(x2>=ZERO)) x0=x2;
  return x0;
}


//=========================================================
//=========================================================

__device__ float cucompute_alpha_b(float temp, float unit_number, float aexp)
{
  // CASE B recombination rate m**3 s*-1
  // temperature should be given in Kelvin
  
  float alpha_b,lambda;
  lambda=2e0*157807e0/temp;
  alpha_b=2.753e-14*powf(lambda,1.5)/powf(1e0+powf(lambda/2.740,0.407),2.242); //cm3/s
#ifdef COSMO
  alpha_b=alpha_b*1e-6*unit_number/(aexp*aexp*aexp); //m3/s
#else
  alpha_b=alpha_b*1e-6*unit_number; //m3/s
#endif
  return alpha_b;
}

//=========================================================
//=========================================================

__device__ float cucompute_alpha_a(float temp, float unit_number, float aexp)
{
  // CASE A recombination rate m**3 s*-1
  // temperature should be given in Kelvin
  
  float alpha_a,lambda;
  lambda=2e0*157807e0/temp;
  alpha_a=1.269e-13*powf(lambda,1.503)/powf(1e0+powf(lambda/0.522,0.470),1.923); //cm3/s
#ifdef COSMO
  alpha_a=alpha_a*1e-6*unit_number/(aexp*aexp*aexp); //m3/s
#else
  alpha_a=alpha_a*1e-6*unit_number; //m3/s
#endif
  return alpha_a;
}

//=========================================================
//=========================================================

__device__ float cucompute_beta(float temp, float unit_number, float aexp)
{
  // Collizional ionization rate m**3 s*-1
  // temperature in Kelvin
  float beta,T5;
  T5=temp/1e5;
  beta=5.85e-11*sqrtf(temp)/(1+sqrtf(T5))*expf(-(157809e0/temp)); //cm3/s
#ifdef COSMO
  beta=beta*1e-6*unit_number/(aexp*aexp*aexp); // !m3/s
#else
  beta=beta*1e-6*unit_number; // !m3/s
#endif
  return beta;
}

//**********************************************************************************
//**********************************************************************************
__device__ void cuCompCooling(float temp, float x, float nH, float *lambda, float *tcool, float aexp,float CLUMPF)
{

 
  float c1,c2,c3,c4,c5,c6;
  float unsurtc;
  float nh2;

  nh2=nH*1e-6;// ! m-3 ==> cm-3
  

  // Collisional Ionization Cooling

  c1=expf(-157809.1e0/temp)*1.27e-21*sqrtf(temp)/(1.f+sqrtf(temp/1e5))*x*(1.f-x)*nh2*nh2*CLUMPF;
  

  // Case A Recombination Cooling

  c2=1.778e-29*temp*powf(2e0*157807e0/temp,1.965e0)/powf(1.f+powf(2e0*157807e0/temp/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2*CLUMPF;
  
  
  // Case B Recombination Cooling

  c6=3.435e-30*temp*powf(2e0*157807e0/temp,1.970e0)/powf(1.f+(powf(2e0*157807e0/temp/2.250e0,0.376e0)),3.720e0)*x*x*nh2*nh2*CLUMPF;
  

  // Collisional excitation cooling

  c3=expf(-118348e0/temp)*7.5e-19/(1+sqrtf(temp/1e5))*x*(1.f-x)*nh2*nh2*CLUMPF;
  
  
  // Bremmsstrahlung

  c4=1.42e-27*1.5e0*sqrtf(temp)*x*x*nh2*nh2*CLUMPF;
  
  // Compton Cooling
  
  c5=1.017e-37*powf(2.727/aexp,4)*(temp-2.727/aexp)*nh2*x;
  
  // Overall Cooling
  
  *lambda=c1+c2+c3+c4+c5+c6;// ! erg*cm-3*s-1
  

  // Unit Conversion

  *lambda=(*lambda)*1e-7*1e6;// ! J*m-3*s-1

  // cooling times

  unsurtc=fmaxf(c1,c2);
  unsurtc=fmaxf(unsurtc,c3);
  unsurtc=fmaxf(unsurtc,c4);
  unsurtc=fmaxf(unsurtc,fabs(c5));
  unsurtc=fmaxf(unsurtc,c6)*1e-7;// ==> J/cm3/s

  *tcool=1.5e0*nh2*(1+x)*1.3806e-23*temp/unsurtc; //Myr
}



//**********************************************************************************
//**********************************************************************************
__device__ void cuCompCoolingDerivate(float temp, float x, float nH, float *lambda, float *tcool, float aexp)
{
  float c1,c2,c3,c4,c5,c6;
  float nh2;

  nh2=nH*1e-6;// ! m-3 ==> cm-3

  // Collisional Ionization Cooling

  //  c1=1.27e-21*sqrtf(temp)/(1e0+sqrtf(temp/1e5))*expf(-157809.1e0/temp)*x*(1-x)*nh2*nh2;

  c1=expf(-157809.1e0/temp)*(2.00418e-10+6.337776e-13*sqrtf(temp)+6.35e-16*temp)/powf(1000.+3.16228*sqrtf(temp),2)/powf(temp,1.5)*x*(1-x)*nh2*nh2;
  
  // Case A Recombination Cooling

  //c2=1.778e-29*temp*powf(2e0*157807e0/temp,1.965e0)/powf(1e0+powf(2e0*157807e0/temp/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2;

  c2=(-1.09724e-18/powf(temp,1.965)/powf(1.+784.853/powf(temp,0.502),2.697)+1.20745e-15/powf(temp,2.467)/powf(1.+784.853/powf(temp,0.502),3.697))*x*x*nh2*nh2;

  // Case B Recombination Cooling

  //  c6=3.435e-30*temp*powf(2e0*157807e0/temp,1.970e0)/powf(1e0+(powf(2e0*157807e0/temp/2.250e0,0.376e0)),3.720e0)*x*x*nh2*nh2;

  c6=(-2.27005e-19/powf(temp,1.97)/powf(1.+86.1514/powf(temp,0.376),3.72)+2.82005e-17/powf(temp,2.346)/powf(1.+86.1514/powf(temp,0.376),4.72))*x*x*nh2*nh2;


  // Collisional excitation cooling

  //  c3=7.5e-19/(1+sqrtf(temp/1e5))*expf(-118348e0/temp)*x*(1-x)*nh2*nh2;

  c3=expf(-118348e0/temp)*(8.8761e-9+(2.80687e-11-1.18585e-16*temp)*sqrtf(temp))/powf(316.228*temp+powf(temp,1.5),2)*x*(1-x)*nh2*nh2;
  
  // Bremmsstrahlung

  //  c4=1.42e-27*1.5e0*sqrtf(temp)*x*x*nh2*nh2;

  c4=1.065e-27/sqrtf(temp)*x*x*nh2*nh2;

  // Compton Cooling
  
  //c5=1.017e-37*powf(2.727/aexp,4)*(temp-2.727/aexp)*nh2*x;
  c5=1.017e-37*powf(2.727/aexp,4)*nh2*x;

  // Overall Cooling
  
  *lambda=c1+c2+c3+c4+c5+c6;// ! erg*cm-3*s-1*K-1

  // Unit Conversion

  *lambda=(*lambda)*1e-1;// ! J*m-3*s-1*K-1

}

//**********************************************************************************

#ifndef TESTCOOL
__global__ void cuComputeTemp(float *cuxion, float *cudensity, float *cutemperature, float *cuegy_new, float fudgecool, float c, float dt,float unit_number, int ncvgcool, float aexp, float hubblet, float *cuflx_new, float CLUMPF, float egy_min, float fesc, float boost, float *cusrc)
{
  int 	tx=threadIdx.x,
	bx=blockIdx.x,
	by=blockIdx.y,
	idx1=tx+bx*blockDim.x+by*gridDim.x*blockDim.x,
	k=idx1/(NCELLX*NCELLY),
	j=(idx1-k*(NCELLX*NCELLY))/NCELLX,
	i=idx1-k*(NCELLX*NCELLY)-j*(NCELLX),
	idx=(i+NBOUND)+(j+NBOUND)*(NCELLX+NBOUND2)+(k+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2), // following a convention a[k,j,i] where i varies the first,
	idloc=tx,
	idloc3=3*idloc,
	igrp,
 	nitcool=0;

  float hnu0=13.6*1.6022e-19,
	Cool,
	tcool,
	dtcool,
	tcool1,
	currentcool_t=0.f,
	alpha,
	alphab,
	beta,
	eint,
	xt,
	ai_tmp1=0.,
	hnu[NGRP],		// ! Average Photon Energy (J)
	factgrp[NGRP],		
	alphae[NGRP],
	alphai[NGRP],		
	et[NGRP],
	p[NGRP];

  __shared__ float
	egyloc[BLOCKCOOL*NGRP],
    floc[3*BLOCKCOOL*NGRP],
	x0[BLOCKCOOL],
	nH[BLOCKCOOL],
    tloc[BLOCKCOOL],
	srcloc[BLOCKCOOL];

#ifdef S_X
  float N2[NGRP];
  float F2[NGRP];
  float E0overI[NGRP];
#endif

  c=c*aexp; 			// switch back to physical velocity
  SECTION_EFFICACE;
  FACTGRP;

  x0[idloc]=cuxion[idx];
  nH[idloc]=cudensity[idx]*unit_number/(aexp*aexp*aexp);
  tloc[idloc]=cutemperature[idx]; 
  srcloc[idloc]=cusrc[idx]*unit_number/(aexp*aexp*aexp); 




  for (igrp=0;igrp<NGRP;igrp++)
	{			// switch to physical units, chemistry remains unchanged with and without cosmo
	  //egyloc[idloc+igrp*BLOCKCOOL]=cuegy_new[idx+igrp*NCELL4];
	  egyloc[idloc+igrp*BLOCKCOOL]=cuegy_new[idx+igrp*NCELL4]*unit_number/(aexp*aexp*aexp); 
	  floc[0+idloc3+igrp*BLOCKCOOL*3]=cuflx_new[0*NCELL4+idx+igrp*NCELL4*3]/(aexp*aexp);
	  floc[1+idloc3+igrp*BLOCKCOOL*3]=cuflx_new[1*NCELL4+idx+igrp*NCELL4*3]/(aexp*aexp);
	  floc[2+idloc3+igrp*BLOCKCOOL*3]=cuflx_new[2*NCELL4+idx+igrp*NCELL4*3]/(aexp*aexp);
 	 }

  __syncthreads();

//  float ncomov=nH[idloc]*(aexp*aexp*aexp);
//  float CLUMPF2=fminf(powf(ncomov/0.25,1.)*(ncomov>0.25)+(ncomov<=0.25)*1.,100.);
//  float CLUMPI=fminf(powf(ncomov/0.5,1.2)*(ncomov>0.5)+(ncomov<=0.5)*1.,10.);
//  float CLUMPF2=100.*powf(ncomov/(ncomov+5.),2)+1.f;
//  float CLUMPI=100.*powf(ncomov/(ncomov+8.),2)+1.f; 
  
//  float CLUMPF2=fmaxf(fminf(1000.f*powf(ncomov,2.5),70.),1.); 
//  float CLUMPI=fmaxf(fminf(20.f*powf(ncomov,2.5),20.),1.); 

//  float CLUMPF2=fmaxf(fminf(1000.f*powf(ncomov,2.5),200.),1.); 
//  float CLUMPI=fmaxf(fminf(20.f*powf(ncomov,2.5),100.),1.); 

//  float CLUMPI=fminf(expf(ncomov/0.5),20.);
//  float CLUMPF2=fminf(expf(ncomov/0.5),20.);
//  float CLUMPF2=fminf(expf(ncomov/0.25),120.);

#ifdef WCLUMP
  float CLUMPF2=fminf(fmaxf(powf(nH[idloc]/6.,0.7),1.),40.);
  float CLUMPI=1.;
#else
  float CLUMPF2=1.;
  float CLUMPI=1.;
#endif

//  float CLUMPF2=fmaxf(fminf(powf(nH[idloc]/400.)^2.5+powf(nH[idloc]/8.)^0.7*(1-x0[idloc]),100.),1.); 
//  float CLUMPI==fmaxf(fminf(powf(nH[idloc]/400.)^2.5+powf(nH[idloc]/8.)^0.2*(1-x0[idloc]),100.),1.); 

  for(igrp=0;igrp<NGRP;igrp++)
	{alphai[igrp] *= CLUMPI;
	 alphae[igrp] *= CLUMPI;}

  while(currentcool_t<dt)
    {
      nitcool++;

      eint=1.5*nH[idloc]*KBOLTZ*(1.f+x0[idloc])*tloc[idloc];

      //== Getting a timestep
      cuCompCooling(tloc[idloc],x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
      
      ai_tmp1=0.;
      for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += ((alphae[igrp])*hnu[igrp]-(alphai[igrp])*hnu0)*egyloc[idloc+igrp*BLOCKCOOL];

      tcool=fabsf(eint/(nH[idloc]*(1.0f-x0[idloc])*ai_tmp1-Cool));
      ai_tmp1=0.;
      dtcool=fminf(fudgecool*tcool,dt-currentcool_t);
      
      alpha=cucompute_alpha_a(tloc[idloc],1.,1.)*CLUMPF2;
      alphab=cucompute_alpha_b(tloc[idloc],1.,1.)*CLUMPF2;
      beta=cucompute_beta(tloc[idloc],1.,1.)*CLUMPF2;
      
      //== Update
      //      egyloc[idloc]=((alpha-alphab)*x0[idloc]*x0[idloc]*nH[idloc]*nH[idloc]*dtcool+egyloc[idloc])/(1.f+dtcool*(alphai*(1.f-x0[idloc])*nH[idloc]+3*hubblet));

      // ABSORPTION
  bool test = 0;
  for(igrp=0;igrp<NGRP;igrp++)
    {
      ai_tmp1 = alphai[igrp];
      et[igrp]=((alpha-alphab)*x0[idloc]*x0[idloc]*nH[idloc]*nH[idloc]*dtcool*factgrp[igrp]+egyloc[idloc+igrp*BLOCKCOOL]+srcloc[idloc]*dtcool*fesc*boost*factgrp[igrp])/(1.f+dtcool*(ai_tmp1*(1.f-x0[idloc])*nH[idloc]+3*hubblet));
      //et[igrp]=egyloc[idloc+igrp*BLOCKCOOL];
      
      if(et[igrp]<0) 	{test=1;}
      
      p[igrp]=(1.f+(alphai[igrp]*nH[idloc]*(1-x0[idloc])+2*hubblet)*dtcool);
    }
  ai_tmp1=0.;

  if (test) 
    {fudgecool/=10.f; 
      continue;	} 

  // IONISATION
#ifndef S_X
  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*et[igrp];}
#else
  N2[0]=1.0f;
  float pp=(1.f-powf(x0[idloc],0.4092f)); 
  if(pp<0.f) pp=0.f; 
    
    N2[1]=1.0f+0.3908f*powf(pp,1.7592f)*E0overI[1]; 
    if(N2[1]<1.0f) N2[1]=1.0f; 
    
    //N2[1]=1.0f;

  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*et[igrp]*N2[igrp];}
#endif

  xt=1.f-(alpha*x0[idloc]*x0[idloc]*nH[idloc]*dtcool+(1.f -x0[idloc]))/(1.f+dtcool*(beta*x0[idloc]*nH[idloc]+ai_tmp1));
  ai_tmp1=0.;

  if((xt>1.f)||(xt<0.f)) 
    {fudgecool/=10.f; 
      continue;	} 

  cuCompCooling(tloc[idloc],xt,nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);

#ifdef COOLING
  // HEATING
#ifndef S_X
  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0));}
#else
  float pp2;
  F2[0]=1.0f;
  F2[1]=1.0f;
  pp2=1.0f-powf(xt,0.2663f); 
  if(pp2<0.f) pp2=0.f; 
  F2[1]=0.9971f*(1.0f-powf(pp2,1.3163f)); 
  
  if(F2[1]>1.0f) F2[1]=1.0f; 
  if(F2[1]<0.0f) F2[1]=0.0f; 
  
  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp];}
#endif

  eint=(eint+dtcool*(nH[idloc]*(1.f-xt)*(ai_tmp1)-Cool))/(1.f+3*hubblet*dtcool);
  //if(eint==0.) printf("NULL TEMP dtcool=%e Cool=%e ai=%e\n",dtcool,Cool,ai_tmp1);
  ai_tmp1=0;

  if(eint<0.f) 
    {
      fudgecool/=10.f; 
      continue;
    } 
#endif

  for(igrp =0;igrp<NGRP;igrp++)
	{egyloc[idloc+igrp*BLOCKCOOL]=et[igrp];
	 floc[0+idloc3+igrp*BLOCKCOOL*3]=floc[0+idloc3+igrp*BLOCKCOOL*3]/p[igrp];
	 floc[1+idloc3+igrp*BLOCKCOOL*3]=floc[1+idloc3+igrp*BLOCKCOOL*3]/p[igrp];
	 floc[2+idloc3+igrp*BLOCKCOOL*3]=floc[2+idloc3+igrp*BLOCKCOOL*3]/p[igrp];	
	}
  
  x0[idloc]=xt;

#ifdef COOLING
      tloc[idloc]=eint/(1.5f*nH[idloc]*KBOLTZ*(1.f+x0[idloc]));
#endif
      currentcool_t+=dtcool;
      if((nitcool==ncvgcool)&&(ncvgcool!=0)) break;
    }

for(igrp=0;igrp<NGRP;igrp++)
  {
    cuegy_new[idx+igrp*NCELL4]=fmax(egyloc[idloc+igrp*BLOCKCOOL]*aexp*aexp*aexp,egy_min*factgrp[igrp]);
    cuflx_new[0*NCELL4+idx+igrp*NCELL4*3]=floc[0+idloc3+igrp*BLOCKCOOL*3]*aexp*aexp;
    cuflx_new[1*NCELL4+idx+igrp*NCELL4*3]=floc[1+idloc3+igrp*BLOCKCOOL*3]*aexp*aexp;
    cuflx_new[2*NCELL4+idx+igrp*NCELL4*3]=floc[2+idloc3+igrp*BLOCKCOOL*3]*aexp*aexp;
  }
 
//if(idx==2395605) printf("cuegy_new in cool=%e\n",cuegy_new[idx]);

  cutemperature[idx]=tloc[idloc];
  cuxion[idx]=x0[idloc];
  __syncthreads();
}


#else
__global__ void cuComputeTemp(float *cuxion, float *cudensity, float *cutemperature, float *cuegy_new, float fudgecool, float c, float dt,float unit_number, int ncvgcool, float aexp, float hubblet,float *cuflx, float CLUMP,float egy_min)
{
  float hnu;
  float Heat,Cool,Cool2,DCool;
  float tcool,dtcool,tcool1;
  float dxdt,tempnew=0.;
  float currentcool_t;
  int nitcool;

  
  //float sige=1.0970e-22;				// ! ED average cross-section (m2) for black body T=10^5 K
  //float sign=1.6279e-22;				// ! ND average cross-section (m2) for black body T=10^5 K
  //hnu=29.64e0;//     ! Average Photon Energy (eV)
  hnu=AVG_EGY;//     ! Average Photon Energy (eV)
  float sige=AVG_CSE;
  float sign=AVG_CSN;

  //float alpha,beta,alphai=1.6279e-22*c,alphae=1.0970e-22*c;
  float alpha,beta,alphai=sign*c,alphae=sige*c;

  int idx,idloc;


  int tx=threadIdx.x;
  int bx=blockIdx.x;
  int by=blockIdx.y;

  int idx1=tx+bx*blockDim.x+by*gridDim.x*blockDim.x;
  int k=idx1/(NCELLX*NCELLY);
  int j=(idx1-k*(NCELLX*NCELLY))/NCELLX;
  int i=idx1-k*(NCELLX*NCELLY)-j*(NCELLX);
  idx=(i+NBOUND)+(j+NBOUND)*(NCELLX+NBOUND2)+(k+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2); // following a convention a[k,j,i] where i varies the first
  
  
  idloc=tx;

  __shared__ float egyloc[BLOCKCOOL];
  __shared__ float x0[BLOCKCOOL];
  __shared__ float nH[BLOCKCOOL];
  __shared__ float tinit[BLOCKCOOL];
  
  
  x0[idloc]=cuxion[idx];
  nH[idloc]=cudensity[idx]*unit_number/(aexp*aexp*aexp);
  egyloc[idloc]=cuegy_new[idx]*unit_number/(aexp*aexp*aexp); // switch to physical units, chemistry remains unchanged with and without cosmo
  tinit[idloc]=cutemperature[idx];

  __syncthreads();

  //Heat=nH[idloc]*(1e0-x0[idloc])*egyloc[idloc]*c*(sige*hnu-sign*13.6e0)*1.6022e-19;
  Heat=nH[idloc]*(1e0-x0[idloc])*egyloc[idloc]*(alphae*hnu-alphai*13.6e0)*1.6022e-19;

  //cuCompCoolingDerivate(tinit[idloc],x0[idloc],nH[idloc],&DCool,&tcool,aexp);
  cuCompCooling(tinit[idloc]*(1.+EPSCOOL),x0[idloc],nH[idloc],&Cool2,&tcool1,aexp);
  cuCompCooling(tinit[idloc],x0[idloc],nH[idloc],&Cool,&tcool1,aexp);
  DCool=(Cool2-Cool)/(EPSCOOL*tinit[idloc]);

  float fboltz=1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc]);

  alpha=cucompute_alpha_a(tinit[idloc],1.,1.);
  beta=cucompute_beta(tinit[idloc],1.,1.);
  dxdt=-alpha*nH[idloc]*x0[idloc]*x0[idloc]+beta*nH[idloc]*x0[idloc]*(1e0-x0[idloc])+alphai*egyloc[idloc]*(1-x0[idloc]);
  tcool=fboltz*tinit[idloc]/fabsf(Heat-Cool+DCool*tinit[idloc]);
  tcool=fminf(tcool,1./fabsf(DCool/fboltz+3*hubblet+dxdt/(1+x0[idloc])));
  dtcool=tcool*fudgecool;


  if(dt<=dtcool)
    {
      //tempnew=tinit[idloc]+(Heat-Cool)*dt/fboltz/(1.+DCool*dt/fboltz);
      tempnew=tinit[idloc]*(1.+(Heat-Cool+DCool*tinit[idloc])*dt/(fboltz*tinit[idloc]))/(1.+dt*(DCool/fboltz+3*hubblet+dxdt/(1.+x0[idloc])));
      tinit[idloc]=tempnew;
    }
  else
    {
      currentcool_t=0e0;
      nitcool=0;
      while(currentcool_t<dt)
	{
	  //cuCompCoolingDerivate(tinit[idloc],x0[idloc],nH[idloc],&DCool,&tcool1,aexp);
	  cuCompCooling(tinit[idloc]*(1.+EPSCOOL),x0[idloc],nH[idloc],&Cool2,&tcool1,aexp);
	  cuCompCooling(tinit[idloc],x0[idloc],nH[idloc],&Cool,&tcool1,aexp);
	  DCool=(Cool2-Cool)/(EPSCOOL*tinit[idloc]);
	  
	  alpha=cucompute_alpha_a(tinit[idloc],1.,1.);
	  beta=cucompute_beta(tinit[idloc],1.,1.);
	  dxdt=-alpha*nH[idloc]*x0[idloc]*x0[idloc]+beta*nH[idloc]*x0[idloc]*(1e0-x0[idloc])+alphai*egyloc[idloc]*(1-x0[idloc]);

	  tcool=fboltz*tinit[idloc]/fabsf(Heat-Cool+DCool*tinit[idloc]);
	  tcool=fminf(tcool,1./fabsf(DCool/fboltz+3*hubblet+dxdt/(1+x0[idloc])));
	  dtcool=fminf(tcool*fudgecool,dt-currentcool_t);
	  
	  //tempnew=tinit[idloc]+(Heat-Cool)*dtcool/fboltz/(1.+DCool*dtcool/fboltz);
	  tempnew=tinit[idloc]*(1.+(Heat-Cool+DCool*tinit[idloc])*dtcool/(fboltz*tinit[idloc]))/(1.+dtcool*(DCool/fboltz+3*hubblet+dxdt/(1.+x0[idloc])));
	  
	  tinit[idloc]=tempnew;
	  
	  currentcool_t=currentcool_t+dtcool;
	  nitcool++;
	  //if(((nitcool>=ncvgcool)&&(ncvgcool!=0))&&(fabs(tempnew-tinit[idloc])/tinit[idloc]<1e-2)) break;
	}
    }

  cutemperature[idx]=tinit[idloc];
  __syncthreads();
}

//**********************************************************************************
//**********************************************************************************
__global__ void cuComputeIon(float *cuegy_new, float *cuflx_new, float *cuxion, float *cudensity, float *cutemperature, float dt, float c, float egy_min, float unit_number, float aexp)
{
  int idx,idloc,idloc3;
  float x0;//,xorg,nH;//,enew;
  //float sign=1.6279e-22; //ND average cross-section (m2) for black bodyT=10^5 K
  float sign=AVG_CSN;
#ifdef COSMO
  float alpha,alphab,beta,alphai=sign*c*unit_number/(aexp*aexp);
#else
  float alpha,alphab,beta,alphai=sign*c*unit_number;
#endif

  float n,p,q;

  int tx=threadIdx.x;
  int bx=blockIdx.x;
  int by=blockIdx.y;
  
  idx=(tx+NBOUND)+(bx+NBOUND)*(NCELLX+NBOUND2)+(by+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  
  idloc=tx;
  idloc3=3*idloc;

  __shared__ float enew[NCELLX];
  __shared__ float floc[3*NCELLX];
  __shared__ float xorg[NCELLX];
  __shared__ float nH[NCELLX];
  __shared__ float temploc[NCELLX];
  
  xorg[idloc]=cuxion[idx];
  x0=xorg[idloc];
  nH[idloc]=cudensity[idx];
  enew[idloc]=cuegy_new[idx];
  floc[0+idloc3]=cuflx_new[0*NCELL4+idx];
  floc[1+idloc3]=cuflx_new[1*NCELL4+idx];
  floc[2+idloc3]=cuflx_new[2*NCELL4+idx];
  temploc[idloc]=cutemperature[idx];

  alpha=cucompute_alpha_a(temploc[idloc],unit_number,aexp);
  alphab=cucompute_alpha_b(temploc[idloc],unit_number,aexp);
  beta=cucompute_beta(temploc[idloc],unit_number,aexp);
  
  
  float m = (beta+alphab)*nH[idloc]*nH[idloc]*dt;
  n = nH[idloc]-(alpha+beta)*nH[idloc]/alphai-alphab*nH[idloc]*nH[idloc]*dt-2*beta*nH[idloc]*nH[idloc]*dt;
  p = -nH[idloc]*(1e0+x0)-enew[idloc]-1e0/(alphai*dt)+beta*nH[idloc]/alphai+beta*nH[idloc]*nH[idloc]*dt;
  q = enew[idloc]+nH[idloc]*x0+x0/(alphai*dt);
      
  
  m*=unit_number;
  n*=unit_number;
  p*=unit_number;
  q*=unit_number;
  
  
  x0=cufindroot3_2(m,n,p,q,x0);

  enew[idloc]=enew[idloc]-nH[idloc]*(x0-xorg[idloc])-alphab*nH[idloc]*nH[idloc]*x0*x0*dt+beta*nH[idloc]*nH[idloc]*x0*(1-x0)*dt;
  enew[idloc]=fmaxf(egy_min,enew[idloc]);
  
  p=(1.f+alphai*nH[idloc]*(1-x0)*dt);
  floc[0+idloc3]=floc[0+idloc3]/p;
  floc[1+idloc3]=floc[1+idloc3]/p;
  floc[2+idloc3]=floc[2+idloc3]/p;
  
  
  
  
  q=floc[0+idloc3]*floc[0+idloc3]+floc[1+idloc3]*floc[1+idloc3]+floc[2+idloc3]*floc[2+idloc3];
  q=sqrtf(q)/(c*enew[idloc]);
  
  if(q>1.)
    {
      floc[0+idloc3]/=q;
      floc[1+idloc3]/=q;
      floc[2+idloc3]/=q;
    }

  cuxion[idx]=x0;
  cuflx_new[0*NCELL4+idx]=floc[0+idloc3];
  cuflx_new[1*NCELL4+idx]=floc[1+idloc3];
  cuflx_new[2*NCELL4+idx]=floc[2+idloc3];
  cuegy_new[idx]=enew[idloc];
}

#endif
