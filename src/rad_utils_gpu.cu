				       
#ifdef WRAD
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>

#include <mpi.h>

#ifdef WCHEM
#include "chem_utils.cuh"
#endif

#include "gpu_type.h"

extern "C" struct OCT *gatherstencilrad(struct OCT *octstart, struct RGRID *stencil, int stride, struct CPUINFO *cpu, int *nread, REAL cloc);
extern "C" struct OCT *scatterstencilrad(struct OCT *octstart, struct RGRID *stencil, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew, REAL cloc);
extern "C" int advanceradGPU (struct OCT **firstoct, int level, struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dxcur, REAL dtnew,REAL aexp, struct RUNPARAMS *param, int chemonly);


extern "C" void create_radstencil_GPU(struct CPUINFO *cpu, int stride);
extern "C" void create_pinned_stencil_rad(struct RGRID **stencil, int stride);
extern "C" void destroy_radstencil_GPU(struct CPUINFO *cpu, int stride);
extern "C" void destroy_pinned_stencil_rad(struct RGRID **stencil, int stride);
extern "C" void create_param_GPU(struct RUNPARAMS *param, struct CPUINFO *cpu);

// ===================================================================
void create_radstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaMalloc((void **)&(cpu->rad_stencil),sizeof(struct RGRID)*stride);
  cudaDeviceSynchronize(); 
  //  printf("Start Error =%s\n",cudaGetErrorString(cudaGetLastError()));

  
}

// ===================================================================
void create_pinned_stencil_rad(struct RGRID **stencil, int stride){
  cudaMallocHost( (void**)stencil, sizeof(struct RGRID)*stride );
  cudaDeviceSynchronize(); 
  //  printf("Start2 Error =%s\n",cudaGetErrorString(cudaGetLastError()));
}

// ===================================================================
void destroy_radstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaFree(cpu->rad_stencil);
}

// ===================================================================
void destroy_pinned_stencil_rad(struct RGRID **stencil, int stride){
  cudaFreeHost(stencil);
}


// ===================================================================
void create_param_GPU(struct RUNPARAMS *param, struct CPUINFO *cpu){

  struct RUNPARAMS paramgpu;
  memcpy(&paramgpu,param,sizeof(struct RUNPARAMS));

  struct COSMOPARAM CG;
  struct COSMOPARAM *dCG;
  memcpy(&CG,param->cosmo,sizeof(struct COSMOPARAM));

  cudaMalloc((void **)&(dCG),sizeof(struct COSMOPARAM)); // allocate cosmo struct on GPU
  cudaMemcpy(dCG,&CG,sizeof(struct COSMOPARAM),cudaMemcpyHostToDevice);  
  paramgpu.cosmo=dCG; // we save the pointer location

#ifdef WRAD
  REAL *space_bound;
  cudaMalloc((void **)&(space_bound),(param->atomic.ngrp_space+1)*sizeof(REAL)); // allocate cosmo struct on GPU
  cudaMemcpy(space_bound,param->atomic.space_bound,(param->atomic.ngrp_space+1)*sizeof(REAL),cudaMemcpyHostToDevice);
  paramgpu.atomic.space_bound=space_bound;

  REAL *time_bound;
  cudaMalloc((void **)&(time_bound),(param->atomic.ngrp_time+1)*sizeof(REAL)); // allocate cosmo struct on GPU
  cudaMemcpy(time_bound,param->atomic.time_bound,(param->atomic.ngrp_time+1)*sizeof(REAL),cudaMemcpyHostToDevice);
  paramgpu.atomic.time_bound=time_bound;

  REAL *hnu;
  cudaMalloc((void **)&(hnu),(param->atomic.n)*sizeof(REAL)); // allocate cosmo struct on GPU
  cudaMemcpy(hnu,param->atomic.hnu,(param->atomic.n)*sizeof(REAL),cudaMemcpyHostToDevice);
  paramgpu.atomic.hnu=hnu;

  REAL *alphae;
  cudaMalloc((void **)&(alphae),(param->atomic.n)*sizeof(REAL)); // allocate cosmo struct on GPU
  cudaMemcpy(alphae,param->atomic.alphae,(param->atomic.n)*sizeof(REAL),cudaMemcpyHostToDevice);
  paramgpu.atomic.alphae=alphae;

  REAL *alphai;
  cudaMalloc((void **)&(alphai),(param->atomic.n)*sizeof(REAL)); // allocate cosmo struct on GPU
  cudaMemcpy(alphai,param->atomic.alphai,(param->atomic.n)*sizeof(REAL),cudaMemcpyHostToDevice);
  paramgpu.atomic.alphai=alphai;

  REAL *factgrp;
  cudaMalloc((void **)&(factgrp),(param->atomic.n)*sizeof(REAL)); // allocate cosmo struct on GPU
  cudaMemcpy(factgrp,param->atomic.factgrp,(param->atomic.n)*sizeof(REAL),cudaMemcpyHostToDevice);
  paramgpu.atomic.factgrp=factgrp;
#endif


  cudaMalloc((void **)&(cpu->dparam),sizeof(struct RUNPARAMS));
  cudaMemcpy(cpu->dparam,&paramgpu,sizeof(struct RUNPARAMS),cudaMemcpyHostToDevice);  
  cudaDeviceSynchronize(); 
}


// =======================================================
__device__ void getcellnei_gpu_rad(int cindex, int *neip, int *cell)
{
  switch(cindex){
  case 0:
    neip[0]=0;cell[0]=1;
    neip[1]=6;cell[1]=1;
    neip[2]=2;cell[2]=2;
    neip[3]=6;cell[3]=2;
    neip[4]=4;cell[4]=4;
    neip[5]=6;cell[5]=4;
    break;
  case 1:
    neip[0]=6;cell[0]=0;
    neip[1]=1;cell[1]=0;
    neip[2]=2;cell[2]=3;
    neip[3]=6;cell[3]=3;
    neip[4]=4;cell[4]=5;
    neip[5]=6;cell[5]=5;
    break;
  case 2:
    neip[0]=0;cell[0]=3;
    neip[1]=6;cell[1]=3;
    neip[2]=6;cell[2]=0;
    neip[3]=3;cell[3]=0;
    neip[4]=4;cell[4]=6;
    neip[5]=6;cell[5]=6;
    break;
  case 3:
    neip[0]=6;cell[0]=2;
    neip[1]=1;cell[1]=2;
    neip[2]=6;cell[2]=1;
    neip[3]=3;cell[3]=1;
    neip[4]=4;cell[4]=7;
    neip[5]=6;cell[5]=7;
    break;
  case 4:
    neip[0]=0;cell[0]=5;
    neip[1]=6;cell[1]=5;
    neip[2]=2;cell[2]=6;
    neip[3]=6;cell[3]=6;
    neip[4]=6;cell[4]=0;
    neip[5]=5;cell[5]=0;
    break;
  case 5:
    neip[0]=6;cell[0]=4;
    neip[1]=1;cell[1]=4;
    neip[2]=2;cell[2]=7;
    neip[3]=6;cell[3]=7;
    neip[4]=6;cell[4]=1;
    neip[5]=5;cell[5]=1;
    break;
  case 6:
    neip[0]=0;cell[0]=7;
    neip[1]=6;cell[1]=7;
    neip[2]=6;cell[2]=4;
    neip[3]=3;cell[3]=4;
    neip[4]=6;cell[4]=2;
    neip[5]=5;cell[5]=2;
    break;
  case 7:
    neip[0]=6;cell[0]=6;
    neip[1]=1;cell[1]=6;
    neip[2]=6;cell[2]=5;
    neip[3]=3;cell[3]=5;
    neip[4]=6;cell[4]=3;
    neip[5]=5;cell[5]=3;
    break;
  }

}


//================================================================================
__device__ void ddiffR(struct Rtype *W2, struct Rtype *W1, struct Rtype *WR){
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    WR->e[igrp]=W2->e[igrp]- W1->e[igrp];
    WR->fx[igrp]=W2->fx[igrp]- W1->fx[igrp];
    WR->fy[igrp]=W2->fy[igrp]- W1->fy[igrp];
    WR->fz[igrp]=W2->fz[igrp]- W1->fz[igrp];
    WR->src[igrp]=W2->src[igrp]- W1->src[igrp];
  }

#ifdef SUPERNOVAE
    WR->snfb=W2->snfb- W1->snfb;
#endif
#ifdef WCHEM
    WR->nhplus=W2->nhplus-W1->nhplus;
#ifdef HELIUM
    WR->nheplus=W2->nheplus-W1->nheplus;
    WR->nhepplus=W2->nhepplus-W1->nhepplus;
#endif
    WR->eint=W2->eint-W1->eint;
    WR->nh=W2->nh-W1->nh;
#endif
}




//================================================================================
__device__ void dminmod_R(struct Rtype *Wm, struct Rtype *Wp, struct Rtype *Wr){

  REAL beta=1.; // 1. for MINBEE 2. for SUPERBEE
  // FLUX LIMITER
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){

    if(Wp->e[igrp]>0){
      Wr->e[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->e[igrp],Wp->e[igrp])),FMIN(Wm->e[igrp],beta*Wp->e[igrp]));
    }
    else{
      Wr->e[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->e[igrp],Wp->e[igrp])),FMAX(Wm->e[igrp],beta*Wp->e[igrp]));
    }


    if(Wp->fx[igrp]>0){
      Wr->fx[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->fx[igrp],Wp->fx[igrp])),FMIN(Wm->fx[igrp],beta*Wp->fx[igrp]));
    }
    else{
      Wr->fx[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->fx[igrp],Wp->fx[igrp])),FMAX(Wm->fx[igrp],beta*Wp->fx[igrp]));
    }


    if(Wp->fy[igrp]>0){
      Wr->fy[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->fy[igrp],Wp->fy[igrp])),FMIN(Wm->fy[igrp],beta*Wp->fy[igrp]));
    }
    else{
      Wr->fy[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->fy[igrp],Wp->fy[igrp])),FMAX(Wm->fy[igrp],beta*Wp->fy[igrp]));
    }


    if(Wp->fz[igrp]>0){
      Wr->fz[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->fz[igrp],Wp->fz[igrp])),FMIN(Wm->fz[igrp],beta*Wp->fz[igrp]));
    }
    else{
      Wr->fz[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->fz[igrp],Wp->fz[igrp])),FMAX(Wm->fz[igrp],beta*Wp->fz[igrp]));
    }

    if(Wp->src[igrp]>0){
      Wr->src[igrp]=FMAX(FMAX(0.,FMIN(beta*Wm->src[igrp],Wp->src[igrp])),FMIN(Wm->src[igrp],beta*Wp->src[igrp]));
    }
    else{
      Wr->src[igrp]=FMIN(FMIN(0.,FMAX(beta*Wm->src[igrp],Wp->src[igrp])),FMAX(Wm->src[igrp],beta*Wp->src[igrp]));
    }

  }




#ifdef SUPERNOVAE
  if(Wp->snfb>0){
      Wr->snfb=FMAX(FMAX(0.,FMIN(beta*Wm->snfb,Wp->snfb)),FMIN(Wm->snfb,beta*Wp->snfb));
    }
    else{
      Wr->snfb=FMIN(FMIN(0.,FMAX(beta*Wm->snfb,Wp->snfb)),FMAX(Wm->snfb,beta*Wp->snfb));
    }
#endif

#ifdef WCHEM
  if(Wp->nhplus>0){
    Wr->nhplus=FMAX(FMAX(0.,FMIN(beta*Wm->nhplus,Wp->nhplus)),FMIN(Wm->nhplus,beta*Wp->nhplus));
  }
  else{
    Wr->nhplus=FMIN(FMIN(0.,FMAX(beta*Wm->nhplus,Wp->nhplus)),FMAX(Wm->nhplus,beta*Wp->nhplus));
  }

#ifdef HELIUM
  if(Wp->nheplus>0){
    Wr->nheplus=FMAX(FMAX(0.,FMIN(beta*Wm->nheplus,Wp->nheplus)),FMIN(Wm->nheplus,beta*Wp->nheplus));
  }
  else{
    Wr->nheplus=FMIN(FMIN(0.,FMAX(beta*Wm->nheplus,Wp->nheplus)),FMAX(Wm->nheplus,beta*Wp->nheplus));
  }

  if(Wp->nhepplus>0){
    Wr->nhepplus=FMAX(FMAX(0.,FMIN(beta*Wm->nhepplus,Wp->nhepplus)),FMIN(Wm->nhepplus,beta*Wp->nhepplus));
  }
  else{
    Wr->nhepplus=FMIN(FMIN(0.,FMAX(beta*Wm->nhepplus,Wp->nhepplus)),FMAX(Wm->nhepplus,beta*Wp->nhepplus));
  }

#endif

  if(Wp->eint>0){
    Wr->eint=FMAX(FMAX(0.,FMIN(beta*Wm->eint,Wp->eint)),FMIN(Wm->eint,beta*Wp->eint));
  }
  else{
    Wr->eint=FMIN(FMIN(0.,FMAX(beta*Wm->eint,Wp->eint)),FMAX(Wm->eint,beta*Wp->eint));
  }

  if(Wp->nh>0){
    Wr->nh=FMAX(FMAX(0.,FMIN(beta*Wm->nh,Wp->nh)),FMIN(Wm->nh,beta*Wp->nh));
  }
  else{
    Wr->nh=FMIN(FMIN(0.,FMAX(beta*Wm->nh,Wp->nh)),FMAX(Wm->nh,beta*Wp->nh));
  }
#endif

}



//================================================================================
__device__ void dinterpminmod_R(struct Rtype *W0, struct Rtype *Wp, struct Rtype *Dx, struct Rtype *Dy, struct Rtype *Dz,REAL dx,REAL dy,REAL dz){
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    Wp->e[igrp] =W0->e[igrp] +dx*Dx->e[igrp] +dy*Dy->e[igrp] +dz*Dz->e[igrp];
    Wp->fx[igrp] =W0->fx[igrp] +dx*Dx->fx[igrp] +dy*Dy->fx[igrp] +dz*Dz->fx[igrp];
    Wp->fy[igrp] =W0->fy[igrp] +dx*Dx->fy[igrp] +dy*Dy->fy[igrp] +dz*Dz->fy[igrp];
    Wp->fz[igrp] =W0->fz[igrp] +dx*Dx->fz[igrp] +dy*Dy->fz[igrp] +dz*Dz->fz[igrp];
    Wp->src[igrp] =W0->src[igrp] +dx*Dx->src[igrp] +dy*Dy->src[igrp] +dz*Dz->src[igrp];
  }


#ifdef WCHEM
    Wp->nhplus =W0->nhplus + dx*Dx->nhplus + dy*Dy->nhplus + dz*Dz->nhplus;
#ifdef HELIUM
    Wp->nheplus =W0->nheplus + dx*Dx->nheplus + dy*Dy->nheplus + dz*Dz->nheplus;
    Wp->nhepplus =W0->nhepplus + dx*Dx->nhepplus + dy*Dy->nhepplus + dz*Dz->nhepplus;
#endif
    Wp->eint =W0->eint +dx*Dx->eint +dy*Dy->eint +dz*Dz->eint;
    Wp->nh =W0->nh +dx*Dx->nh +dy*Dy->nh +dz*Dz->nh;
#endif
#ifdef SUPERNOVAE
    Wp->snfb =W0->snfb +dx*Dx->snfb +dy*Dy->snfb +dz*Dz->snfb;
#endif
}



//================================================================================
//================================================================================
__device__ REAL dEddington(REAL fx, REAL fy, REAL fz, REAL ee, REAL c,int i,int j)
{
  REAL c2e=ee*c*c; // 2 flop
  REAL ff=0.;
  REAL arg,chi,res=0.;
  REAL n[3];
  n[0]=0.;n[1]=0.;n[2]=0.;

  if(ee>0)
    {
      ff=SQRT(fx*fx+fy*fy+fz*fz); // 6 flop
      if(ff>0)
	{
	  n[0]=fx/ff; 
	  n[1]=fy/ff;
	  n[2]=fz/ff; 
	}
      ff=ff/(c*ee); // 2flop
    }
  
  arg=FMAX(4.-3.*ff*ff,0.); // 4 flop
  chi=(3.+4.*ff*ff)/(5.+2.*SQRT(arg)); // 7 flops

  if(i==j) res=(1.-chi)/2.*c2e; // 1 flops on average
  arg=(3.*chi-1.)/2.*c2e;
  res+=arg*n[i]*n[j];

  return res;
}



// =============================================================================================================

__global__ void drad_sweepX(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i,igrp;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fp,fm;
  REAL up,um;

  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  
  i=bx*blockDim.x+tx;
  if(i<nread){
  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei_gpu_rad(icell, vnei, vcell); // we get the neighbors
      
    memset(FL,0,sizeof(REAL)*NVAR_R);
    memset(FR,0,sizeof(REAL)*NVAR_R);

    // Getting the original state ===========================
      
    curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);
      
    /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
    memset(ffact,0,sizeof(int)*2);
    for(iface=0;iface<2;iface++){
      memcpy(RC+iface,curcell,sizeof(struct Rtype));
    }

    // Neighbor "MUSCL-LIKE" reconstruction
    for(iface=0;iface<2;iface++){
      inei=iface;
      memcpy(RN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].rfield),sizeof(struct Rtype));

      int condsplit;
#ifdef COARSERAD
      condsplit=1;
#else
      condsplit=(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split);
#endif
      if(condsplit){
	ffact[iface]=1; // we cancel the contriubtion of split neighbors
      }
    }

    // X DIRECTION =========================================================================
      
    // --------- solving the Riemann Problems LEFT

    for(igrp=0;igrp<NGRP;igrp++){
	
      // E
      up=RC[0].e[igrp];
      um=RN[0].e[igrp];

      fp=RC[0].fx[igrp];
      fm=RN[0].fx[igrp];

      FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
      //FX

      up=RC[0].fx[igrp];
      um=RN[0].fx[igrp];

      fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,0);
      fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,0);

      FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      //FY

      up=RC[0].fy[igrp];
      um=RN[0].fy[igrp];

      fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,0);
      fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,0);

      FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      //FZ

      up=RC[0].fz[igrp];
      um=RN[0].fz[igrp];

      fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,0);
      fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,0);

      FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

    }
      

    // ===========================================

    // --------- solving the Riemann Problems RIGHT

    for(igrp=0;igrp<NGRP;igrp++){
	
      // E
      up=RN[1].e[igrp];
      um=RC[1].e[igrp];

      fp=RN[1].fx[igrp];
      fm=RC[1].fx[igrp];

      FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
      //FX

      up=RN[1].fx[igrp];
      um=RC[1].fx[igrp];

      fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,0);
      fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,0);

      FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      //FY

      up=RN[1].fy[igrp];
      um=RC[1].fy[igrp];

      fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,0);
      fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,0);

      FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      //FX

      up=RN[1].fz[igrp];
      um=RC[1].fz[igrp];

      fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,0);
      fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,0);

      FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

    }
      
    //========================= copy the fluxes
    // Cancelling the fluxes from splitted neighbours


    for(igrp=0;igrp<NGRP;igrp++){
      for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0]; 
      for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1]; 
    }

    memcpy(stencil[i].New.cell[icell].rflux+0*NVAR_R*NGRP,FL,sizeof(REAL)*NVAR_R*NGRP);
    memcpy(stencil[i].New.cell[icell].rflux+1*NVAR_R*NGRP,FR,sizeof(REAL)*NVAR_R*NGRP);

    //ready for the next oct
  }
}
}



// =============================================================================================================

__global__ void drad_sweepY(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i,igrp;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fp,fm;
  REAL up,um;
  
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  
  i=bx*blockDim.x+tx;
  if(i<nread){
  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei_gpu_rad(icell, vnei, vcell); // we get the neighbors
      
  
      
      memset(FL,0,sizeof(REAL)*NVAR_R);
      memset(FR,0,sizeof(REAL)*NVAR_R);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);
      
      /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);
      for(iface=0;iface<2;iface++){
	memcpy(RC+iface,curcell,sizeof(struct Rtype));
      }

      // Neighbor "MUSCL-LIKE" reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+2;
	memcpy(RN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].rfield),sizeof(struct Rtype));

	int condsplit;
#ifdef COARSERAD
	condsplit=1;
#else
	condsplit=(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split);
#endif
	if(condsplit){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}
	
      }




      // Y DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems FRONT


      
      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RC[0].e[igrp];
	um=RN[0].e[igrp];

	fp=RC[0].fy[igrp];
	fm=RN[0].fy[igrp];

	FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RC[0].fx[igrp];
	um=RN[0].fx[igrp];

	fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,1);
	fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,1);

	FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RC[0].fy[igrp];
	um=RN[0].fy[igrp];

	fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,1);
	fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,1);

	//if(up==1.) printf("FL=%e fp=%e fm=%e || ",0.5*(fp+fm)+0.5*c*(um-up),fp,fm);

	FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RC[0].fz[igrp];
	um=RN[0].fz[igrp];

	fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,1);
	fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,1);

	FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      

	// ===========================================


      

      // --------- solving the Riemann Problems BACK

      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RN[1].e[igrp];
	um=RC[1].e[igrp];

	fp=RN[1].fy[igrp];
	fm=RC[1].fy[igrp];

	FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RN[1].fx[igrp];
	um=RC[1].fx[igrp];

	fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,1);
	fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,1);

	FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RN[1].fy[igrp];
	um=RC[1].fy[igrp];

	fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,1);
	fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,1);

	//printf("FR=%e fp=%e fm=%e up=%e um=%e\n",0.5*(fp+fm+c*(um-up)),fp,fm, up,um);
	FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RN[1].fz[igrp];
	um=RC[1].fz[igrp];

	fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,1);
	fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,1);

	FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      
      
      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours

      for(igrp=0;igrp<NGRP;igrp++){
	for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0]; 
	for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1]; 
      }
      
      memcpy(stencil[i].New.cell[icell].rflux+2*NVAR_R*NGRP,FL,sizeof(REAL)*NVAR_R*NGRP);
      memcpy(stencil[i].New.cell[icell].rflux+3*NVAR_R*NGRP,FR,sizeof(REAL)*NVAR_R*NGRP);

    //ready for the next oct
  }
}
}


// ===================================================================================================

__global__ void drad_sweepZ(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i,igrp;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fp,fm;
  REAL up,um;
  
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  
  i=bx*blockDim.x+tx;
  if(i<nread){
    for(icell=0;icell<8;icell++){ // we scan the cells
      getcellnei_gpu_rad(icell, vnei, vcell); // we get the neighbors
      
      
      memset(FL,0,sizeof(REAL)*NVAR_R);
      memset(FR,0,sizeof(REAL)*NVAR_R);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);
      
      /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);
      for(iface=0;iface<2;iface++){
	memcpy(RC+iface,curcell,sizeof(struct Rtype));
      }

      // Neighbor "MUSCL-LIKE" reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+4;
	memcpy(RN+iface,&(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].rfield),sizeof(struct Rtype));

	int condsplit;
#ifdef COARSERAD
	condsplit=1;
#else
	condsplit=(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split);
#endif
	if(condsplit){
	  ffact[iface]=1; // we consider the contriubtion of split neighbors
	}
      }



      
      // Z DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems BOTTOM


      
      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RC[0].e[igrp];
	um=RN[0].e[igrp];

	fp=RC[0].fz[igrp];
	fm=RN[0].fz[igrp];

	FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	//FX

	up=RC[0].fx[igrp];
	um=RN[0].fx[igrp];

	fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,2);
	fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,2);

	FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RC[0].fy[igrp];
	um=RN[0].fy[igrp];

	fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,2);
	fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,2);

	FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FZ

	up=RC[0].fz[igrp];
	um=RN[0].fz[igrp];

	fp=dEddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,2);
	fm=dEddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,2);

	FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      

      // ===========================================


      

      // --------- solving the Riemann Problems TOP

      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RN[1].e[igrp];
	um=RC[1].e[igrp];

	fp=RN[1].fz[igrp];
	fm=RC[1].fz[igrp];

	FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RN[1].fx[igrp];
	um=RC[1].fx[igrp];

	fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,2);
	fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,2);

	FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RN[1].fy[igrp];
	um=RC[1].fy[igrp];

	fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,2);
	fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,2);

	FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RN[1].fz[igrp];
	um=RC[1].fz[igrp];

	fp=dEddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,2);
	fm=dEddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,2);

	FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      
      
      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours

      for(igrp=0;igrp<NGRP;igrp++){
	for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0]; 
	for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1]; 
      }



      memcpy(stencil[i].New.cell[icell].rflux+4*NVAR_R*NGRP,FL,sizeof(REAL)*NVAR_R*NGRP);
      memcpy(stencil[i].New.cell[icell].rflux+5*NVAR_R*NGRP,FR,sizeof(REAL)*NVAR_R*NGRP);

      //ready for the next oct
    }
}
}

// ===================================================================================================
// ===================================================================================================

__global__ void dupdatefieldrad(struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew, REAL cloc)
{
  int i,icell,igrp;
  struct Rtype R;
  struct Rtype Rupdate;
  REAL one;
  int flx;
  REAL dtsurdx=dtnew/dxcur;
  REAL F[NFLUX_R];

  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  i=bx*blockDim.x+tx;
  if(i<nread){
    for(icell=0;icell<8;icell++){ // we scan the cells
      int condsplit;
#ifdef COARSERAD
      condsplit=0;
#else
      condsplit=(stencil[i].oct[6].cell[icell].split);
#endif
      if(condsplit) continue;
      memcpy(F,stencil[i].New.cell[icell].rflux,sizeof(REAL)*NFLUX_R);// New fluxes from the stencil
    
    // ==== updating
    
    // actually we compute and store the delta U only
      one=1.;
      memset(&R,0,sizeof(struct Rtype)); // setting delta U
      for(flx=0;flx<6;flx++){
	for(igrp=0;igrp<NGRP;igrp++){
	  R.e[igrp]  +=F[0+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fx[igrp] +=F[1+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fy[igrp] +=F[2+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fz[igrp] +=F[3+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	}
	one*=-1.;
      }
    
#ifndef WCHEM
      // adding the source contribution
      REAL SRC;
      for(igrp=0;igrp<NGRP;igrp++){
	SRC=stencil[i].oct[6].cell[icell].rfield.src;
	R.e[igrp]  +=SRC*dtnew+EMIN;
      }
#endif
    
    // scatter back the delta Uwithin the stencil
    
    //memcpy(&(stencil[i].New.cell[icell].deltaR),&R,sizeof(struct Rtype));
    
    // TESTING FULL UPDATE IN STENCIL APPROACH
    
    memcpy(&Rupdate,&stencil[i].New.cell[icell].rfieldnew,sizeof(struct Rtype));
    
    for(igrp=0;igrp<NGRP;igrp++){ 
      Rupdate.e[igrp]   +=R.e[igrp];
      Rupdate.fx[igrp]  +=R.fx[igrp];
      Rupdate.fy[igrp]  +=R.fy[igrp];
      Rupdate.fz[igrp]  +=R.fz[igrp];
    }
      
    //memcpy(&(curoct->cell[icell].rfieldnew),&Rupdate,sizeof(struct Rtype));
    memcpy(&stencil[i].New.cell[icell].rfieldnew,&Rupdate,sizeof(struct Rtype));
    }
  }
}

// ====================================================================================================================

int advanceradGPU (struct OCT **firstoct, int level, struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dxcur, REAL dtnew,REAL aexp, struct RUNPARAMS *param, int chemonly){

  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *curoct0;
  int nreadtot,nread;
  //double t[10];
  //double tg=0.,th=0.,tu=0.,ts=0.;//,tfu=0.,ttot=0.;
  REAL cloc; // the speed of light in code units
  CUDA_CHECK_ERROR("Rad Start");

  cloc=aexp*param->clight*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v;
  //printf("cloc=%e aexp=%e\n",cloc,aexp);


  cudaStream_t stream[cpu->nstream];
  int vnread[cpu->nstream];
  int is,offset;
  // creating the streams
  for(is=0;is<cpu->nstream;is++){
    cudaStreamCreate(&stream[is]);
  }

  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  nreadtot=0;
  int ng;
  int nt;

  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do{
      curoct0=nextoct;
      curoct=curoct0;

      //t[0]=MPI_Wtime();

      // streaming ====================
      offset=0;
      for(is=0;is<cpu->nstream;is++){
	// ------------ gathering the stencil value values
	curoct=nextoct;
	if(curoct!=NULL){
	  //printf("Start Error  -2fg=%s is=%d vnread=%d offset=%d\n",cudaGetErrorString(cudaGetLastError()),is,vnread[is],offset);
	  nextoct= gatherstencilrad(curoct,stencil+offset,stride/cpu->nstream,cpu, vnread+is,cloc);
	  //printf("Start Error  -1=%s is=%d vnread=%d offset=%d\n",cudaGetErrorString(cudaGetLastError()),is,vnread[is],offset);
	  if(vnread[is]!=0){

	    ng=((vnread[is]-1)/cpu->nthread)+1; // +1 to treat leftovers
	    if(ng==1){
	      nt=vnread[is];
	    }
	    else{
	      nt=cpu->nthread;
	    }
	  
	    dim3 gridoct(ng);
	    dim3 blockoct(nt);
	  
#ifdef WCHEM
	    dim3 gridoct_chem(ng);
	    dim3 blockoct_chem(nt);
#endif      

	    //t[2]=MPI_Wtime();
	  
	    //printf("Start Error  0=%s is=%d vnread=%d offset=%d\n",cudaGetErrorString(cudaGetLastError()),is,vnread[is],offset);

	    cudaMemcpyAsync(cpu->rad_stencil+offset,stencil+offset,vnread[is]*sizeof(struct RGRID),cudaMemcpyHostToDevice,stream[is]);  
	  
	    //printf("Start Error  1=%s is=%d vnread=%d offset=%d\n",cudaGetErrorString(cudaGetLastError()),is,vnread[is],offset);
	  
#ifndef COARSERAD
	    int condadvec=1;
#else
	    int condadvec=((level==param->lcoarse)&&(!chemonly));
#endif
	    

#ifndef NOCOMP
	    if(condadvec){
/* // ------------ solving the hydro */
	      drad_sweepX<<<gridoct,blockoct,0,stream[is]>>>(cpu->rad_stencil+offset,level,cpu->rank,vnread[is],stride,dxcur,dtnew,cloc);   
	      drad_sweepY<<<gridoct,blockoct,0,stream[is]>>>(cpu->rad_stencil+offset,level,cpu->rank,vnread[is],stride,dxcur,dtnew,cloc);  
	      drad_sweepZ<<<gridoct,blockoct,0,stream[is]>>>(cpu->rad_stencil+offset,level,cpu->rank,vnread[is],stride,dxcur,dtnew,cloc);  
	    }
	  
	    //printf("Start Error  2=%s\n",cudaGetErrorString(cudaGetLastError()));
	    // ------------ updating values within the stencil
	  
	    //t[4]=MPI_Wtime();
	  
	    if(condadvec) dupdatefieldrad<<<gridoct,blockoct,0,stream[is]>>>(cpu->rad_stencil+offset,vnread[is],stride,cpu,dxcur,dtnew,cloc); 
	    //printf("Start Error  3=%s\n",cudaGetErrorString(cudaGetLastError()));
	    // ----------- perform physical cooling and ionisation 
#ifdef WCHEM
	    dchemrad<<<gridoct_chem,blockoct_chem,0,stream[is]>>>(cpu->rad_stencil+offset,vnread[is],stride,cpu,dxcur,dtnew,cpu->dparam,aexp,chemonly); 
#endif
	    
#endif
	    cudaMemcpyAsync(stencil+offset,cpu->rad_stencil+offset,vnread[is]*sizeof(struct RGRID),cudaMemcpyDeviceToHost,stream[is]);
	    //printf("Start Error  4=%s\n",cudaGetErrorString(cudaGetLastError()));

	    offset+=vnread[is];
	  }
	}
      }
      
      // ------------ scatter back the FLUXES
      //CUDA_CHECK_ERROR("Bef Sync");
      cudaDeviceSynchronize();
      //CUDA_CHECK_ERROR("Af Sync");
      //t[6]=MPI_Wtime();
   
      nread=offset;
      //printf("Start Error  5=%s\n",cudaGetErrorString(cudaGetLastError()));
      nextoct=scatterstencilrad(curoct0,stencil, nread, cpu,dxcur,dtnew,cloc);

      //printf("Start Error  6=%s\n",cudaGetErrorString(cudaGetLastError()));

      //t[8]=MPI_Wtime();

      nreadtot+=nread;
    }while(nextoct!=NULL);
  }

  // Destroying the streams
  for(is=0;is<cpu->nstream;is++){
    cudaStreamDestroy(stream[is]);
  }
  //printf("Start Error Hyd =%s nreadtot=%d\n",cudaGetErrorString(cudaGetLastError()),nreadtot);

  //printf("GPU | tgat=%e tcal=%e tup=%e tscat=%e\n",tg,th,tu,ts);
  CUDA_CHECK_ERROR("Rad Stop");

  return nreadtot;
}

#endif
