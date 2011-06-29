
#include "prototypes.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WGPU
#include "interface.h"
#include "vector_gpu.h"
#endif

float poisson_jacob(int level,int levelcoarse,int levelmax, struct OCT **firstoct,struct MULTIVECT* vectors,int stride,struct CPUINFO *cpu,float omegam,float tsim, struct PACKET **sendbuffer,struct PACKET **recvbuffer, int niter, float acc){
  
  int i;
  float dx;
  struct OCT *curoct;
  struct OCT *nextoct;
  int nread;
  float norm_d=0.,normd_org=0.;
  float res=0.;
  double tit[2];
  int iter;
  float restot=0.;
  float factdens;
  float epsilon;
  double t[10];
    
  // ------------- cleaning the vector positions
  clean_vec(levelmax,firstoct);

  // --------------- setting the first oct of the level
  curoct=firstoct[level-1];


#ifdef WMPI
  int nolevel=(curoct==NULL);
  // --------------- do we need to calculate the level ?
  MPI_Allreduce(MPI_IN_PLACE,&nolevel,1,MPI_INT,MPI_SUM,cpu->comm);
  if(nolevel==cpu->nproc){
    //we skip the calculation
    epsilon=0.;
    return epsilon;
  }
#endif

  if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
  // ------------- setting the factor of the density (cosmo vs newtonian gravity)

  if(level>=levelcoarse){
#ifndef TESTCOSMO
    factdens=4.0*M_PI;
#else
    factdens=1.5*omegam/tsim;
#endif
    dx=pow(0.5,level);
  }
  else{
    factdens=1.;
    dx=pow(0.5,level);
  }
  
  // ------------- some checks

  if(((stride*8)<pow(2,3*levelcoarse))){
    printf("still debugging please increase stride : stride=%d required=%d\n",stride,(int)pow(2,3*levelcoarse));
    abort();
  }

  // -------------  cleaning working arrays

  memset(vectors->vecpot,0,sizeof(float)*stride*8);
  memset(vectors->vecpotnew,0,sizeof(float)*stride*8);
  memset(vectors->vecden,0,sizeof(float)*stride*8);
  memset(vectors->vecnei,0,sizeof(int)*stride*6);
  memset(vectors->vecl,0,sizeof(int)*stride);
  memset(vectors->veccpu,0,sizeof(int)*stride);
  memset(vectors->vecicoarse,0,sizeof(int)*stride);

#ifdef WGPU
  cudaMemset(vectors->vecpot_d,0,sizeof(float)*stride*8);
  cudaMemset(vectors->vecpotnew_d,0,sizeof(float)*stride*8);
  cudaMemset(vectors->vecden_d,0,sizeof(float)*stride*8);
  cudaMemset(vectors->vecnei_d,0,sizeof(int)*stride*6);
  cudaMemset(vectors->vecl_d,0,sizeof(int)*stride);
  cudaMemset(vectors->veccpu_d,0,sizeof(int)*stride);
  cudaMemset(vectors->vecicoarse_d,0,sizeof(int)*stride);
#endif


  // ------------------- setting neighbors to -1
  for(i=0;i<stride*6;i++){vectors->vecnei[i]=-1;}

  

  // ---------------- gathering the neighbors
  
#ifdef NEWJACK2
  int nocttotal;
  nocttotal=countvecocts(curoct,stride,cpu,&nread); // we count the actual numbers of octs involved
  stride=(nocttotal%64==0?nocttotal:(nocttotal/64+1)*64); // let us be wild ! we change stride for nread
  //printf("new stride=%d\n",stride);
  nextoct=gathervecnei2(curoct,vectors->vecnei,stride,cpu,&nread); // actual gather
#else
  nextoct=gathervecnei(curoct,vectors->vecnei,vectors->vecpot,1,vectors->vecl,stride,cpu,&nread);
#endif
  
  
  // ------------ gathering the values
#ifdef NEWJACK2
  if(level>levelcoarse){
    //gathering the data from coarse levels
    nextoct=gathervec2(firstoct[level-2],vectors->vecden,0,vectors->vecl,vectors->vecicoarse,vectors->veccpu,stride,cpu,&nread); // density
    nextoct=gathervec2(firstoct[level-2],vectors->vecpot,1,vectors->vecl,vectors->vecicoarse,vectors->veccpu,stride,cpu,&nread); // potential
  }

  nextoct=gathervec2(curoct,vectors->vecden,0,vectors->vecl,vectors->vecicoarse,vectors->veccpu,stride,cpu,&nread); // density (only for coarse level and finer)
  nextoct=gathervec2(curoct,vectors->vecpot,1,vectors->vecl,vectors->vecicoarse,vectors->veccpu,stride,cpu,&nread); // potential

  //printf("nread=%d stride=%d\n",nread,stride);
#else
  nextoct=gathervec(curoct,vectors->vecden,0,vectors->vecl,stride,cpu,&nread); // density
  nextoct=gathervec(curoct,vectors->vecpot,1,vectors->vecl,stride,cpu,&nread); // potential
#endif

     
  // -------------- we contrast the density by removing the average value
#ifndef WGPU      
  if(level>=levelcoarse) remove_valvec(vectors->vecden,nread,stride,1.,level,vectors->vecl); // restricted to coarse + fine levels
#else
  CPU2GPU(vectors->vecden_d,vectors->vecden,sizeof(float)*stride*8);
  CPU2GPU_INT(vectors->vecl_d  ,vectors->vecl,sizeof(int)*stride);
  CPU2GPU_INT(vectors->veccpu_d  ,vectors->veccpu,sizeof(int)*stride);
  if(level>=levelcoarse) remove_valvec_GPU(vectors->vecden_d,nread,stride,1.,level,vectors->vecl_d); // restricted to coarse + fine levels
  GPU2CPU(vectors->vecden,vectors->vecden_d,sizeof(float)*stride*8);
#endif
  
  // --------------- we square and sum the density
#ifndef WGPU
    // CPU CASE
  norm_d+=square_vec(vectors->vecden,nread,stride,level,cpu->rank,vectors->vecl,vectors->veccpu)*pow(factdens,2);
#else
    // GPU CASE
  norm_d+=square_vec_GPU(vectors->vecden_d,nread,stride,level,cpu->rank,vectors->vecl_d,vectors->veccpu_d,vectors->vec2_d,vectors->vecsum_d)*pow(factdens,2);
#endif


#ifdef WGPU
  // ---------------  sending data to GPU
  CPU2GPU(vectors->vecpotnew_d,vectors->vecpotnew,sizeof(float)*stride*8);
  CPU2GPU(vectors->vecpot_d,vectors->vecpot,sizeof(float)*stride*8);
  CPU2GPU_INT(vectors->vecnei_d,vectors->vecnei,sizeof(int)*stride*6);
  CPU2GPU_INT(vectors->vecicoarse_d,vectors->vecicoarse,sizeof(int)*stride);
#endif

  }
  
#ifdef WMPI
  // ---------------  reducing the norm of the density
  normd_org=norm_d;
  MPI_Allreduce(MPI_IN_PLACE,&norm_d,1,MPI_FLOAT,MPI_SUM,cpu->comm);
#endif

  //======================== looping over iterations
  
  for(iter=0;iter<niter;iter++){
    res=0.;
    if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){

    //---------------   computing the laplacian
#ifdef NEWJACK2
#ifdef WGPU
    cudaMemset(vectors->vec2_d,0,sizeof(float)*stride*8);
    cudaMemset(vectors->vecsum_d,0,sizeof(float)*stride*8);
    res=laplacian_vec2_GPU(vectors->vecden_d,vectors->vecpot_d,vectors->vecpotnew_d,vectors->vecnei_d,vectors->vecl_d,vectors-> vecicoarse_d,vectors->veccpu_d,level,cpu->rank,nread,stride,dx,factdens,vectors->vec2_d,vectors->vecsum_d);
    //if(level==7) printf("res=%e\n",res);
#else
    res=laplacian_vec2(vectors->vecden,vectors->vecpot,vectors->vecpotnew,vectors->vecnei,vectors->vecl,vectors->vecicoarse,vectors->veccpu,level,cpu->rank,nread,stride,dx,factdens);
#endif
#else
    res=laplacian_vec(vectors->vecden,vectors->vecpot,vectors->vecpotnew,vectors->vecnei,vectors->vecl,level,nread,stride,dx,omegam,tsim);
#endif
    
    // ---------------  exchange potential evaluations
#ifdef WGPU
    GPU2GPU(vectors->vecpot_d,vectors->vecpotnew_d,sizeof(float)*stride*8);
#else
    memcpy(vectors->vecpot,vectors->vecpotnew,sizeof(float)*stride*8);
#endif
	
    
    // exchanging boundary data through the network ======================
#ifdef WMPI
#ifdef WGPU
    // ---------------  GPU 2 HOST
    GPU2CPU(vectors->vecpot,vectors->vecpot_d,sizeof(float)*8*stride);
#endif
    // ---------------  scatter back the result to the tree
    t[1]=MPI_Wtime();
    nextoct=scattervec_light(curoct,vectors->vecpot,1,stride,cpu,nread,level);  // border only
    t[2]=MPI_Wtime();
#endif
    }

    // ---------------  sending tree data through the network
#ifdef WMPI
    mpi_exchange_level(cpu,sendbuffer,recvbuffer,2,(iter==0),level);
#endif
    //mpi_exchange(cpu,sendbuffer,recvbuffer,2,(iter==0));
    
#ifdef WMPI
    if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
      t[3]=MPI_Wtime();
      nextoct=gathervec2_light(curoct,vectors->vecpot,1,stride,cpu,&nread,level); // potential
#ifdef WGPU
      // ---------------  HOST 2 GPU
      t[4]=MPI_Wtime();
      CPU2GPU(vectors->vecpot_d,vectors->vecpot,sizeof(float)*8*stride);
#endif
      t[5]=MPI_Wtime();
    }
    
    float restot;
    // ---------------  Reducing the residuals
    MPI_Allreduce(&res,&restot,1,MPI_FLOAT,MPI_SUM,cpu->comm);
    res=restot;
#endif
      // END MPI EXCHANGE =================================
    
    // ---------------  breaking condition
    if(level>=levelcoarse){
      epsilon=sqrt(res/norm_d);
    }
    else{
      epsilon=sqrt(res);
    }
    if(epsilon<acc) {
      break;
    }
    // ---------------------------


  } // next iteration ready

  if((cpu->rank==0)&&(norm_d!=0.)) printf("level=%2d iter=%4d relative residual=%e res=%e den=%e \n ",level,iter,sqrt(res/norm_d),sqrt(res),sqrt(norm_d));

  if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
#ifdef MULTIGRID
    // if multigrid enabled, we recompute the residual
    // to spare some memory we copy the residual in vecpotnew
    if(level<=levelcoarse){
#ifndef WGPU
      // CPU
      residual_vec2(vectors->vecden,vectors->vecpot,vectors->vecpotnew,vectors->vecnei,vectors->vecl,vectors->vecicoarse,vectors->veccpu,level,cpu->rank,nread,stride,dx,factdens);
#else
      residual_vec2_GPU(vectors->vecden_d,vectors->vecpot_d,vectors->vecpotnew_d,vectors->vecnei_d,vectors->vecl_d, vectors->vecicoarse_d, vectors->veccpu_d,level,cpu->rank,nread,stride,dx,factdens);
#endif
    }

#endif


    // --------------- full transfer of data to the tree
#ifdef WGPU
    GPU2CPU(vectors->vecpot,vectors->vecpot_d,sizeof(float)*8*stride);
#endif

    //---------------  scatter back the result
    nextoct=scattervec(curoct,vectors->vecpot,1,stride,cpu,nread); 

#ifdef MULTIGRID
#ifdef WGPU
    // -------------- bringing back the residual in the potnew field
    GPU2CPU(vectors->vecpotnew,vectors->vecpotnew_d,sizeof(float)*8*stride);
#endif
    // -------------- we scatter the residual in the temp field of the tree
    nextoct=scattervec(curoct,vectors->vecpotnew,2,stride,cpu,nread); 
#endif
  }

#ifdef MULTIGRID
#ifdef WMPI
  mpi_exchange_level(cpu,sendbuffer,recvbuffer,4,(iter==0),level);
#endif
#endif
  
  return epsilon;
}

//=====================================================

float  poisson_mgrid(int level,int levelcoarse,int levelmax,int levelmin, struct OCT **firstoct,struct MULTIVECT *vectors,int stride, struct CPUINFO *cpu, float omegam, float tsim, struct PACKET** sendbuffer, struct PACKET **recvbuffer, int niter, float acc){

  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  float res;

  // pre relaxation
  res=poisson_jacob(level,levelcoarse,levelmax,firstoct,vectors,stride,cpu,omegam,tsim,sendbuffer,recvbuffer,15,acc);

  if(!((level==levelcoarse)&&(res<acc))){
  // reduction
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
    do{ 
      curoct=nextoct;
      nextoct=curoct->next;
      curoct->parent->density=0.;
      for(icell=0;icell<8;icell++){
	curoct->parent->density+=curoct->cell[icell].temp*0.125; // we average the residual and store it as the new density
      }
    }while(nextoct!=NULL);
  }

  // full relaxation at coarse level or recursive call to mgrid
  if((level-1)==levelmin){
    poisson_jacob(level-1,levelcoarse,levelmax,firstoct,vectors,stride,cpu,omegam,tsim,sendbuffer,recvbuffer,niter,acc);
  }
  else{
    poisson_mgrid(level-1,levelcoarse,levelmax,levelmin,firstoct,vectors,stride,cpu,omegam,tsim,sendbuffer,recvbuffer,niter,acc);
  }

  // prolongation + correction
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
    do{ 
      curoct=nextoct;
      nextoct=curoct->next;
      for(icell=0;icell<8;icell++) {
	curoct->cell[icell].pot-=curoct->parent->pot; // we propagate the error and correct the evaluation
      }
    }while(nextoct!=NULL);
  }

  // post relaxation
  res=poisson_jacob(level,levelcoarse,levelmax,firstoct,vectors,stride,cpu,omegam,tsim,sendbuffer,recvbuffer,15,acc);
  }
  
  return res;
}
