#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "params.h"
#include "bnd.h"
#include "GPU.h"

#define NCELL4 ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))

//-------------------------------------------------------------------------------------------------------

__global__ void cusetboundaryref_xm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=NBOUND+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=NBOUND-1+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];

  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]=-fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundaryref_xp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=NCELLX+NBOUND-1+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=NCELLX+NBOUND+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];

  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]=-fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundaryref_ym(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(NBOUND)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(NBOUND-1)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]=-fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundaryref_yp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(NCELLY+NBOUND-1)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(NCELLY+NBOUND)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]=-fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}
__global__ void cusetboundaryref_zm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(by)*(NCELLX+NBOUND2)+(NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(by)*(NCELLX+NBOUND2)+(NBOUND-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]=-fr[2];

}

__global__ void cusetboundaryref_zp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(by)*(NCELLX+NBOUND2)+(NCELLZ+NBOUND-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(by)*(NCELLX+NBOUND2)+(NCELLZ+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]=-fr[2];

}


//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------

__global__ void cusetboundarytrans_xm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=NBOUND+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=NBOUND-1+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];

  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundarytrans_xp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=NCELLX+NBOUND-1+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=NCELLX+NBOUND+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];

  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundarytrans_ym(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(NBOUND)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(NBOUND-1)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundarytrans_yp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(NCELLY+NBOUND-1)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(NCELLY+NBOUND)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}
__global__ void cusetboundarytrans_zm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(by)*(NCELLX+NBOUND2)+(NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(by)*(NCELLX+NBOUND2)+(NBOUND-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundarytrans_zp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(by)*(NCELLX+NBOUND2)+(NCELLZ+NBOUND-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(by)*(NCELLX+NBOUND2)+(NCELLZ+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}


//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------

__global__ void cusetboundaryper_xm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=NCELLX+NBOUND-1+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=NBOUND-1+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];

  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;
  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundaryper_xp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=NBOUND+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=NCELLX+NBOUND+(tx)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];

  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundaryper_ym(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(NCELLY+NBOUND-1)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(NBOUND-1)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundaryper_yp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(NBOUND)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(NCELLY+NBOUND)*(NCELLX+NBOUND2)+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}
__global__ void cusetboundaryper_zm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(by)*(NCELLX+NBOUND2)+(NCELLZ+NBOUND-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(by)*(NCELLX+NBOUND2)+(NBOUND-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}

__global__ void cusetboundaryper_zp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx)
{
  int tx=threadIdx.x+NBOUND;
  int by=blockIdx.x+NBOUND;
  float er,xr,dr,tr,fr[3];

  int idx=(tx)+(by)*(NCELLX+NBOUND2)+(NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxm1=(tx)+(by)*(NCELLX+NBOUND2)+(NCELLZ+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  er=cuegy[idx];
  xr=cuxion[idx];
  dr=cudensity[idx];
  tr=cutemperature[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];
  
  cuegy[idxm1]=er;
  cudensity[idxm1]=dr;
  cuxion[idxm1]=xr;
  cutemperature[idxm1]=tr;

  cuflx[idxm1+0*NCELL4]= fr[0];
  cuflx[idxm1+1*NCELL4]= fr[1];
  cuflx[idxm1+2*NCELL4]= fr[2];

}


