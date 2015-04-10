
#include <stdio.h>
#include <stdlib.h>
#include <cudpp.h>
#include "prototypes.h"
#include "gpu_type.h"

extern "C" void initlocaldevice(int,int);
extern "C" int countdevices(int);
extern "C" void checkdevice(int);
extern "C" void CPU2GPU(float *gpupt, float *cpupt, int noctet);
extern "C" void GPU2CPU(float *cpupt, float *gpupt, int noctet);
extern "C" void GPU2GPU(float *cpupt, float *gpupt, int noctet);
extern "C" void CPU2GPU_INT(int *gpupt, int *cpupt, int noctet);
extern "C" void GPU2CPU_INT(int *cpupt, int *gpupt, int noctet);
extern "C" void CPU2GPU_UINT(unsigned int *gpupt, unsigned int *cpupt, int noctet);
extern "C" void GPU2CPU_UINT(unsigned int *gpupt, unsigned int *cpupt, int noctet);
extern "C" REAL * GPUallocREAL(int nmem);
extern "C" unsigned long int GPUallocScanPlan(int stride);

//#define CUERR() printf("\n %s \n",cudaGetErrorString(cudaGetLastError()))

//************************************************************************************************************************
//************************************************************************************************************************

void CPU2GPU(float *gpupt, float *cpupt, int noctet)
{
  cudaMemcpy(gpupt,cpupt,noctet,cudaMemcpyHostToDevice);  
}

void GPU2CPU(float *cpupt, float *gpupt, int noctet)
{
  cudaMemcpy(cpupt,gpupt,noctet,cudaMemcpyDeviceToHost);  
}

void GPU2GPU(float *cpupt, float *gpupt, int noctet)
{
  cudaMemcpy(cpupt,gpupt,noctet,cudaMemcpyDeviceToDevice);  
}

void CPU2GPU_INT(int *gpupt, int *cpupt, int noctet)
{
  cudaMemcpy(gpupt,cpupt,noctet,cudaMemcpyHostToDevice);  
}

void GPU2CPU_INT(int *cpupt, int *gpupt, int noctet)
{
  cudaMemcpy(cpupt,gpupt,noctet,cudaMemcpyDeviceToHost);  
}

void CPU2GPU_UINT(unsigned int *gpupt, unsigned int *cpupt, int noctet)
{
  cudaMemcpy(gpupt,cpupt,noctet,cudaMemcpyHostToDevice);  
}

void GPU2CPU_UINT(unsigned int *cpupt, unsigned int *gpupt, int noctet)
{
  cudaMemcpy(cpupt,gpupt,noctet,cudaMemcpyDeviceToHost);  
}

int countdevices(int rank)
{
  int count;
  cudaGetDeviceCount(&count);
  printf("%d device(s) found locally by proc %d\n",count,rank);
  return count;
}

void initlocaldevice(int rank, int count)
{
  int devicenum=rank%count;
  struct cudaDeviceProp prop;
  cudaSetDevice(0);
  cudaGetDeviceProperties(&prop,0);
  printf("Local Device #%d (%s)initialized for proc #%d with %d B of VRAM\n",devicenum,prop.name,rank,(int)prop.totalGlobalMem);
}


void checkdevice(int rank)
{
  int idevice;
  cudaGetDevice(&idevice);
  printf("Local Device #%d dedicated to proc #%d\n",idevice,rank);
}


REAL * GPUallocREAL(int nmem){
  REAL *resA;
  cudaMalloc((void **)&resA,sizeof(REAL)*nmem);
  printf("POINT=%p\n",resA);
  return resA;
}

unsigned long int GPUallocScanPlan(int nmem){
  struct CUPARAM *cuparam;

  cudaMallocHost( (void**)&cuparam, sizeof(struct CUPARAM));

  cudppCreate(&(cuparam->theCudpp));
    
  cuparam->config.algorithm = CUDPP_SCAN;
  cuparam->config.datatype = CUDPP_REAL; // homemade type
  cuparam->config.options=CUDPP_OPTION_BACKWARD | CUDPP_OPTION_INCLUSIVE;
  cuparam->config.op=CUDPP_ADD;

  cuparam->scanplan =0;
  cudppPlan(cuparam->theCudpp,&(cuparam->scanplan), cuparam->config, nmem, 1, 0);

  return (unsigned long int) cuparam;
}
