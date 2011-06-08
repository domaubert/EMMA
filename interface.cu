
#include <stdio.h>
#include <stdlib.h>

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
  //printf("%d device(s) found locally by proc %d\n",count,rank);
  return count;
}

void initlocaldevice(int rank, int count)
{
  int devicenum=rank%count;
  struct cudaDeviceProp prop;
  cudaSetDevice(devicenum);
  cudaGetDeviceProperties(&prop,devicenum);
  printf("Local Device #%d initialized for proc #%d with %d B of VRAM\n",devicenum,rank,(int)prop.totalGlobalMem);
}


void checkdevice(int rank)
{
  int idevice;
  cudaGetDevice(&idevice);
  printf("Local Device #%d dedicated to proc #%d\n",idevice,rank);
}

