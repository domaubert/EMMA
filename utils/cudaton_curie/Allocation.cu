#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "params.h"
#include "bnd.h"
#include "GPU.h"
#include "Atomic.h"
#ifdef WMPI
#include "communication.h"
#endif


extern "C" void Allocation(void);
extern "C" void cuAllocation(void);


//*************************************************************************
//*************************************************************************

#define CUERR() printf("\n %s \n",cudaGetErrorString(cudaGetLastError()))

//*************************************************************************
//*************************************************************************
void Allocation(void)
{
  egy=  (float *) calloc(NGRP*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2),sizeof(float));
  flx=  (float *) calloc(NGRP*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)*3,sizeof(float));
//  dedd=  (float *) calloc((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)*3*3,sizeof(float));

//  egy_new=  (float *) calloc(NGRP*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2),sizeof(float));
//  flx_new=  (float *) calloc(NGRP*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)*3,sizeof(float));

#ifndef SDISCRETE
  src0=  (float *) calloc((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2),sizeof(float));
#else
  src0=  (float *) calloc(nsource,sizeof(float));
  src0pos=  (int *) calloc(nsource*3,sizeof(int));
#endif

  xion=  (float *) calloc((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2),sizeof(float));
  density=  (float *) calloc((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2),sizeof(float));
  temperature=  (float *) calloc((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2),sizeof(float));

#ifdef WMPI
  buff=(float*) calloc(NBUFF,sizeof(float)); // 4 for EGY + FLUX // HERE WE ASSUME THAT NCELLZ >=NCELLX,Y
#endif
}

void cuAllocation(void)
{
  cudaMalloc((void**)&cuegy,NGRP*((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(float)); //1
  cudaMalloc((void**)&cuflx,NGRP*((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)*3)*sizeof(float)); //3

  cudaMalloc((void**)&cuegy_new,NGRP*((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(float));//1 
  cudaMalloc((void**)&cuflx_new,NGRP*((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)*3)*sizeof(float)); //3
 
#ifdef SDISCRETE
  cudaMalloc((void**)&cusrc0,nsource*sizeof(float));//1
  cudaMalloc((void**)&cusrc0pos,3*nsource*sizeof(int));//1
#else
  cudaMalloc((void**)&cusrc0,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(float));//1
#endif

  cudaMalloc((void**)&cuxion,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(float));// 1
  cudaMalloc((void**)&cudensity,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(float));//1
  cudaMalloc((void**)&cutemperature,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(float));//1

#ifdef WMPI
  cudaMalloc((void**)&cubuff,NBUFF*sizeof(float)); 
#endif

  cudaThreadSynchronize();
  
  //  CUERR();
  
}
