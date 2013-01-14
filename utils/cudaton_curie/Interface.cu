#include <stdio.h>
#include <stdlib.h>

#ifndef WMPI
extern "C" void initlocaldevice(int,int);
void initlocaldevice(int rank, int count)
{
  int devicenum=rank%count;
  struct cudaDeviceProp prop;
  cudaSetDevice(devicenum);
  cudaGetDeviceProperties(&prop,devicenum);
  printf("Local Device #%d initialized for proc #%d with %d Mb of VRAM\n",devicenum,rank,prop.totalGlobalMem);
}

#else

#include "GPU.h"
#include "bnd.h"
#include "communication.h"

extern "C" void initlocaldevice(int,int);
extern "C" void initlocaldeviceTITANE(int rank, int count, int *hostnum, int nhost);
extern "C" int countdevices(int);
extern "C" void checkdevice(int);

#define NCELL4 ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))
#define CUERR() printf("\n %s \n",cudaGetErrorString(cudaGetLastError()))
//************************************************************************************************************************
//************************************************************************************************************************

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
  countdevices(rank);
  printf("%d device(s) found locally by proc %d\n",count,rank);
  cudaSetDevice(devicenum);
  cudaGetDeviceProperties(&prop,devicenum);
  printf("Local Device #%d initialized for proc #%d with %d Mb of VRAM\n",devicenum,rank,prop.totalGlobalMem);
}

void initlocaldeviceTITANE(int rank, int count, int *hostnum, int nhost)
{
  int devicenum=0;
  int i;

  for(i=0;i<nhost;i++)
    {
      if((hostnum[i]==hostnum[rank])&&(i!=rank))
	{
	  //shared host found
	  if(rank>i)
	    {
	      devicenum=1;
	    }
	}
    }
  struct cudaDeviceProp prop;
  cudaSetDevice(devicenum);
  cudaGetDeviceProperties(&prop,devicenum);
}

void checkdevice(int rank)
{
  int idevice;
  cudaGetDevice(&idevice);
  printf("Local Device #%d dedicated to proc #%d\n",idevice,rank);
}

//************************************************************************************************************************
//************************************************************************************************************************
//************************************************************************************************************************
//************************************************************************************************************************

__global__ void bufgather_zm(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(tx+NBOUND)+(bx+NBOUND)*(NCELLX+NBOUND2)+(NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLX);

  er=cuegy[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];

  cubuff[idxbuf+0*NCELLX*NCELLY]=er;
  cubuff[idxbuf+1*NCELLX*NCELLY]=fr[0];
  cubuff[idxbuf+2*NCELLX*NCELLY]=fr[1];
  cubuff[idxbuf+3*NCELLX*NCELLY]=fr[2];
}

__global__ void bufgather_zp(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(tx+NBOUND)+(bx+NBOUND)*(NCELLX+NBOUND2)+(NCELLZ+NBOUND-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLX);

  er=cuegy[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];

  cubuff[idxbuf+0*NCELLX*NCELLY]=er;
  cubuff[idxbuf+1*NCELLX*NCELLY]=fr[0];
  cubuff[idxbuf+2*NCELLX*NCELLY]=fr[1];
  cubuff[idxbuf+3*NCELLX*NCELLY]=fr[2];
}

//************************************************************************************************************************

__global__ void bufgather_ym(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(tx+NBOUND)+(NBOUND)*(NCELLX+NBOUND2)+(bx+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLX);

  er=cuegy[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];

  cubuff[idxbuf+0*NCELLX*NCELLZ]=er;
  cubuff[idxbuf+1*NCELLX*NCELLZ]=fr[0];
  cubuff[idxbuf+2*NCELLX*NCELLZ]=fr[1];
  cubuff[idxbuf+3*NCELLX*NCELLZ]=fr[2];
}

__global__ void bufgather_yp(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(tx+NBOUND)+(NCELLY+NBOUND-1)*(NCELLX+NBOUND2)+(bx+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLX);

  er=cuegy[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];

  cubuff[idxbuf+0*NCELLX*NCELLZ]=er;
  cubuff[idxbuf+1*NCELLX*NCELLZ]=fr[0];
  cubuff[idxbuf+2*NCELLX*NCELLZ]=fr[1];
  cubuff[idxbuf+3*NCELLX*NCELLZ]=fr[2];
}

//************************************************************************************************************************

__global__ void bufgather_xm(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(NBOUND)+(tx+NBOUND)*(NCELLX+NBOUND2)+(bx+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLY);

  er=cuegy[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];

  cubuff[idxbuf+0*NCELLY*NCELLZ]=er;
  cubuff[idxbuf+1*NCELLY*NCELLZ]=fr[0];
  cubuff[idxbuf+2*NCELLY*NCELLZ]=fr[1];
  cubuff[idxbuf+3*NCELLY*NCELLZ]=fr[2];
}

__global__ void bufgather_xp(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(NCELLX+NBOUND-1)+(tx+NBOUND)*(NCELLX+NBOUND2)+(bx+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLY);

  er=cuegy[idx];
  fr[0]=cuflx[idx+0*NCELL4];
  fr[1]=cuflx[idx+1*NCELL4];
  fr[2]=cuflx[idx+2*NCELL4];

  cubuff[idxbuf+0*NCELLY*NCELLZ]=er;
  cubuff[idxbuf+1*NCELLY*NCELLZ]=fr[0];
  cubuff[idxbuf+2*NCELLY*NCELLZ]=fr[1];
  cubuff[idxbuf+3*NCELLY*NCELLZ]=fr[2];
}


//************************************************************************************************************************
//************************************************************************************************************************
//************************************************************************************************************************
//************************************************************************************************************************

__global__ void bufscatter_zm(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(tx+NBOUND)+(bx+NBOUND)*(NCELLX+NBOUND2)+(NBOUND-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLX);

  er=cubuff[idxbuf+0*NCELLX*NCELLY];
  fr[0]=cubuff[idxbuf+1*NCELLX*NCELLY];
  fr[1]=cubuff[idxbuf+2*NCELLX*NCELLY];
  fr[2]=cubuff[idxbuf+3*NCELLX*NCELLY];

  cuegy[idx]=er;
  cuflx[idx+0*NCELL4]=fr[0];
  cuflx[idx+1*NCELL4]=fr[1];
  cuflx[idx+2*NCELL4]=fr[2];

}

__global__ void bufscatter_zp(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(tx+NBOUND)+(bx+NBOUND)*(NCELLX+NBOUND2)+(NCELLZ+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLX);
          
  er=cubuff[idxbuf+0*NCELLX*NCELLY];
  fr[0]=cubuff[idxbuf+1*NCELLX*NCELLY];
  fr[1]=cubuff[idxbuf+2*NCELLX*NCELLY];
  fr[2]=cubuff[idxbuf+3*NCELLX*NCELLY];

  cuegy[idx]=er;
  cuflx[idx+0*NCELL4]=fr[0];
  cuflx[idx+1*NCELL4]=fr[1];
  cuflx[idx+2*NCELL4]=fr[2];

}


//************************************************************************************************************************

__global__ void bufscatter_ym(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(tx+NBOUND)+(NBOUND-1)*(NCELLX+NBOUND2)+(bx+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLX);

  er=cubuff[idxbuf+0*NCELLX*NCELLZ];
  fr[0]=cubuff[idxbuf+1*NCELLX*NCELLZ];
  fr[1]=cubuff[idxbuf+2*NCELLX*NCELLZ];
  fr[2]=cubuff[idxbuf+3*NCELLX*NCELLZ];

  cuegy[idx]=er;
  cuflx[idx+0*NCELL4]=fr[0];
  cuflx[idx+1*NCELL4]=fr[1];
  cuflx[idx+2*NCELL4]=fr[2];
}


__global__ void bufscatter_yp(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(tx+NBOUND)+(NCELLY+NBOUND)*(NCELLX+NBOUND2)+(bx+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLX);
          
  er=cubuff[idxbuf+0*NCELLX*NCELLZ];
  fr[0]=cubuff[idxbuf+1*NCELLX*NCELLZ];
  fr[1]=cubuff[idxbuf+2*NCELLX*NCELLZ];
  fr[2]=cubuff[idxbuf+3*NCELLX*NCELLZ];

  cuegy[idx]=er;
  cuflx[idx+0*NCELL4]=fr[0];
  cuflx[idx+1*NCELL4]=fr[1];
  cuflx[idx+2*NCELL4]=fr[2];

}

//************************************************************************************************************************

__global__ void bufscatter_xm(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(NBOUND-1)+(tx+NBOUND)*(NCELLX+NBOUND2)+(bx+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLY);

  er=cubuff[idxbuf+0*NCELLY*NCELLZ];
  fr[0]=cubuff[idxbuf+1*NCELLY*NCELLZ];
  fr[1]=cubuff[idxbuf+2*NCELLY*NCELLZ];
  fr[2]=cubuff[idxbuf+3*NCELLY*NCELLZ];

  cuegy[idx]=er;
  cuflx[idx+0*NCELL4]=fr[0];
  cuflx[idx+1*NCELL4]=fr[1];
  cuflx[idx+2*NCELL4]=fr[2];
}

__global__ void bufscatter_xp(float *cuegy, float *cuflx, float *cubuff)
{
  int tx=threadIdx.x;
  int bx=blockIdx.x;
  float er,fr[3];
  int idx=(NCELLX+NBOUND)+(tx+NBOUND)*(NCELLX+NBOUND2)+(bx+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  int idxbuf=tx+bx*(NCELLY);
          
  er=cubuff[idxbuf+0*NCELLY*NCELLZ];
  fr[0]=cubuff[idxbuf+1*NCELLY*NCELLZ];
  fr[1]=cubuff[idxbuf+2*NCELLY*NCELLZ];
  fr[2]=cubuff[idxbuf+3*NCELLY*NCELLZ];

  cuegy[idx]=er;
  cuflx[idx+0*NCELL4]=fr[0];
  cuflx[idx+1*NCELL4]=fr[1];
  cuflx[idx+2*NCELL4]=fr[2];
}

//************************************************************************************************************************

void exchange_zp(float *cuegy, float *cuflx, float *cubuff, float *buff, int *neighbor,int even)
{
  dim3 gridexchange(NCELLY);
  dim3 blockexchange(NCELLX);

  //  if(neighbor[5]==6) printf("1005 spotted even=%d\n",even);
  if(even)
    {
      //  if(neighbor[5]==6) printf("1005 spotted before gahther");
      bufgather_zp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      //      cudaThreadSynchronize();
      //      if(neighbor[5]==6) printf("1005 spotted before smemcpy");
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      //      if(neighbor[5]==6) printf("1005 spotted before siwthc");
      switchbuff(buff,neighbor[5],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      //      CUERR();
      bufscatter_zp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
  else
    {
      bufgather_zm<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[4],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_zm<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
}


void exchange_zm(float *cuegy, float *cuflx, float *cubuff, float *buff, int *neighbor,int even)
{
  dim3 gridexchange(NCELLY);
  dim3 blockexchange(NCELLX);
  
  if(!even)
    {
      bufgather_zp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[5],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_zp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
  else
    {
      bufgather_zm<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[4],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_zm<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
}

//************************************************************************************************************************

void exchange_yp(float *cuegy, float *cuflx, float *cubuff, float *buff, int *neighbor,int even)
{
  dim3 gridexchange(NCELLZ);
  dim3 blockexchange(NCELLX);

  if(even)
    {
      bufgather_yp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[3],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_yp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
  else
    {
      bufgather_ym<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[2],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_ym<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }

}


void exchange_ym(float *cuegy, float *cuflx, float *cubuff, float *buff, int *neighbor,int even)
{
  dim3 gridexchange(NCELLZ);
  dim3 blockexchange(NCELLX);
  
  if(!even)
    {
      bufgather_yp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[3],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_yp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
  else
    {
      bufgather_ym<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[2],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_ym<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
}

//************************************************************************************************************************

void exchange_xp(float *cuegy, float *cuflx, float *cubuff, float *buff, int *neighbor,int even)
{
  dim3 gridexchange(NCELLZ);
  dim3 blockexchange(NCELLY);

  if(even)
    {
      bufgather_xp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[1],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_xp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
  else
    {
      bufgather_xm<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[0],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_xm<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }

}


void exchange_xm(float *cuegy, float *cuflx, float *cubuff, float *buff, int *neighbor,int even)
{
  dim3 gridexchange(NCELLZ);
  dim3 blockexchange(NCELLY);
  
  if(!even)
    {
      bufgather_xp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[1],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_xp<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
  else
    {
      bufgather_xm<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
      cudaMemcpy(buff,cubuff,NBUFF*sizeof(float),cudaMemcpyDeviceToHost);
      switchbuff(buff,neighbor[0],NBUFF);
      cudaMemcpy(cubuff,buff,NBUFF*sizeof(float),cudaMemcpyHostToDevice);
      bufscatter_xm<<<gridexchange,blockexchange>>>(cuegy,cuflx,cubuff);
    }
}


//************************************************************************************************************************
//************************************************************************************************************************
//************************************************************************************************************************
//************************************************************************************************************************

#endif
