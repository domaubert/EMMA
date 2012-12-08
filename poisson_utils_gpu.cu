#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include <mpi.h>
#include <cudpp.h>

#define NTHREAD 16

extern "C" struct OCT *gatherstencilgrav(struct OCT *octstart, struct STENGRAV *stencil, int stride, struct CPUINFO *cpu, int *nread);
extern "C" struct OCT *scatterstencilgrav(struct OCT *octstart, struct STENGRAV *stencil, int nread, int stride, struct CPUINFO *cpu);
extern "C" void clean_vecpos(int level,struct OCT **firstoct);
extern "C" struct OCT *gatherstencilgrav_nei(struct OCT *octstart, struct STENGRAV *gstencil, int stride, struct CPUINFO *cpu, int *nread);
extern "C" void update_pot_in_tree(int level,struct OCT ** firstoct,  struct CPUINFO *cpu, struct RUNPARAMS *param);
extern "C" REAL PoissonJacobiGPU(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim);
extern "C" void create_gravstencil_GPU(struct CPUINFO *cpu, int stride);
extern "C" void create_pinned_gravstencil(struct STENGRAV *gstencil, int stride);
extern "C" void destroy_pinned_gravstencil(struct STENGRAV *gstencil, int stride);
extern "C" void destroy_gravstencil_GPU(struct CPUINFO *cpu, int stride);


// ====================== structure for CUDPP =======


struct CUPARAM{
  CUDPPHandle theCudpp;
  CUDPPConfiguration config;
  CUDPPHandle scanplan;
};


// =======================================================
// =======================================================

#ifndef FASTGRAV
void create_pinned_gravstencil(struct STENGRAV *gstencil, int stride){
  struct GGRID *grav_stencil;
  cudaMallocHost( (void**)&grav_stencil, sizeof(struct GGRID)*stride );
  gstencil->stencil=grav_stencil;
  REAL *pnew;
  REAL *res;
  REAL *resLR;
  cudaMallocHost( (void**)&pnew, sizeof(REAL)*stride*8 );
  cudaMallocHost( (void**)&res, sizeof(REAL)*stride*8 );
  gstencil->pnew=pnew;
  gstencil->res=res;
#ifdef ONFLYRED
  cudaMallocHost( (void**)&resLR, sizeof(REAL)*stride);
  gstencil->resLR=resLR;
#endif
  abort();
}
#else
void create_pinned_gravstencil(struct STENGRAV *gstencil, int stride){
  REAL *d;
  REAL *p;
  REAL *pnew;
  int *nei;
  int *level;
  int *cpu;
  char *valid;
  REAL *res;
  REAL *res2;
  REAL *resLR;

  cudaMallocHost( (void**)&d, sizeof(REAL)*stride*8 );
  cudaMallocHost( (void**)&p, sizeof(REAL)*stride*8 );
  cudaMallocHost( (void**)&pnew, sizeof(REAL)*stride*8 );
  
  cudaMallocHost( (void**)&nei, sizeof(int)*stride*7 );
  cudaMallocHost( (void**)&level, sizeof(int)*stride);
  cudaMallocHost( (void**)&cpu, sizeof(int)*stride);
  cudaMallocHost( (void**)&valid, sizeof(char)*stride);
  cudaMallocHost( (void**)&res, sizeof(REAL)*stride*8);
  cudaMallocHost( (void**)&res2, sizeof(REAL)*stride*8);
  cudaMallocHost( (void**)&resLR, sizeof(REAL)*stride);

  printf("Pinned memory on Host ok\n"); 

  gstencil->d=d;
  gstencil->p=p;
  gstencil->pnew=pnew;
  gstencil->nei=nei;
  gstencil->level=level;
  gstencil->cpu=cpu;
  gstencil->valid=valid;
  gstencil->res=res;
  gstencil->res2=res2;
  gstencil->resLR=resLR;
}
#endif


#ifndef FASTGRAV
void destroy_pinned_gravstencil(struct STENGRAV *gstencil, int stride){

  cudaFreeHost(gstencil->stencil);
  cudaFreeHost(gstencil->pnew);
  cudaFreeHost(gstencil->res);
#ifdef ONFLYRED
  cudaFreeHost(gstencil->resLR);
#endif
}
#else
void destroy_pinned_gravstencil(struct STENGRAV *gstencil, int stride){
  cudaFreeHost(gstencil->d);
  cudaFreeHost(gstencil->p);
  cudaFreeHost(gstencil->pnew);
  cudaFreeHost(gstencil->level);
  cudaFreeHost(gstencil->cpu);
  cudaFreeHost(gstencil->nei);
  cudaFreeHost(gstencil->valid);
  cudaFreeHost(gstencil->res);
  cudaFreeHost(gstencil->res2);
  cudaFreeHost(gstencil->resLR);
}
#endif



#ifndef FASTGRAV
void create_gravstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaMalloc((void **)&(cpu->dev_stencil),sizeof(struct GGRID)*stride);
  cudaMalloc((void **)&(cpu->res),sizeof(REAL)*stride*8);
  cudaMalloc((void **)&(cpu->pnew),sizeof(REAL)*stride*8);
  cudaMalloc((void **)&(cpu->resLR),sizeof(REAL)*stride);
}
#else
void create_gravstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaMalloc( (void**)&(cpu->dev_stencil->d), sizeof(REAL)*stride*8 );
  cudaMalloc( (void**)&(cpu->dev_stencil->p), sizeof(REAL)*stride*8 );
  cudaMalloc( (void**)&(cpu->dev_stencil->pnew), sizeof(REAL)*stride*8 );
  
  cudaMalloc( (void**)&(cpu->dev_stencil->nei), sizeof(int)*stride*7 );
  cudaMalloc( (void**)&(cpu->dev_stencil->level), sizeof(int)*stride);
  cudaMalloc( (void**)&(cpu->dev_stencil->cpu), sizeof(int)*stride);
  cudaMalloc( (void**)&(cpu->dev_stencil->valid), sizeof(char)*stride);
  cudaMalloc( (void**)&(cpu->dev_stencil->res), sizeof(REAL)*stride*8);
  cudaMalloc( (void**)&(cpu->dev_stencil->res2), sizeof(REAL)*stride*8);
  cudaMalloc( (void**)&(cpu->dev_stencil->resLR), sizeof(REAL)*stride);

  cudaMalloc( (void**)&(cpu->gpu_stencil), sizeof(struct STENGRAV));
  cudaMemcpy(cpu->gpu_stencil,cpu->dev_stencil,sizeof(struct STENGRAV),cudaMemcpyHostToDevice);
}
#endif

#ifndef FASTGRAV
void destroy_gravstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaFree(cpu->dev_stencil);
  cudaFree(cpu->res);
  cudaFree(cpu->pnew);
  cudaFree(cpu->resLR);
  
}
#else
void destroy_gravstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaFree(cpu->dev_stencil->d);
  cudaFree(cpu->dev_stencil->p);
  cudaFree(cpu->dev_stencil->pnew);
  cudaFree(cpu->dev_stencil->level);
  cudaFree(cpu->dev_stencil->cpu);
  cudaFree(cpu->dev_stencil->nei);
  cudaFree(cpu->dev_stencil->valid);
  cudaFree(cpu->dev_stencil->res);
  cudaFree(cpu->dev_stencil->res2);
  cudaFree(cpu->dev_stencil->resLR);
}
#endif

// =======================================================
__device__ void getcellnei_gpu(int cindex, int *neip, int *cell)
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


//========================================================================
//========================================================================

#ifndef FASTGRAV
__global__ void dev_PoissonJacobi_single(struct GGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, int flag, REAL factdens, REAL *vres, REAL *stockres, REAL *stockpnew, REAL *stockresLR){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int inei,icell;
  int i;
  REAL temp;
  REAL res;
  int vnei[6],vcell[6];
  //int ioct[7]={12,14,10,16,4,22,13};
  int ioct[7]={0,1,2,3,4,5,6};
  struct Gtype *curcell;
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  i=bx*blockDim.x+tx;


  stockresLR[i]=0.;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei_gpu(icell, vnei, vcell); // we get the neighbors
      
    temp=0.;
    res=0.;
      
    curcell=&(stencil[i].oct[ioct[6]].cell[icell].gdata);
      
    // Computing the laplacian ===========================
    
    for(inei=0;inei<6;inei++){
      temp+=stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].gdata.p;
    }
 
    // setting up the residual
    res=temp;
      
    // we finish the laplacian
    temp=temp/6.0;
    temp=temp-dx*dx*curcell->d/6.0*factdens;

    // we finsih the residual
    res=res-6.0*curcell->p;
    res=res/(dx*dx)-factdens*curcell->d;

    // we store the new value of the potential
    //    stencil[i].pnew[icell]=temp;
    //stockpnew[icell*nread+i]=temp;
    stockpnew[icell+i*8]=temp;

    // we store the local residual
    if(flag) {
      vres[icell+i*8]=factdens*curcell->d; 
      stockres[icell+i*8]=factdens*curcell->d;
    }
    else{
      vres[icell+i*8]=res*res; 
      stockres[icell+i*8]=res;
    }

    stockresLR[i]+=res*0.125;
  }


}
#else
__global__ void dev_PoissonJacobi_single(struct STENGRAV *gstencil, int level, int curcpu, int nread,int stride,REAL dx, int flag, REAL factdens){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int inei,icell;
  int i;
  REAL temp;
  REAL res;
  int vnei[6],vcell[6];
  //int ioct[7]={12,14,10,16,4,22,13};
  int ival=0;
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  i=bx*blockDim.x+tx;

  struct Gtype ct;
  struct Gtype *curcell;

  if(gstencil->valid[i]!=1){
   
  }
  else{
  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei_gpu(icell, vnei, vcell); // we get the neighbors
    // we skip octs that are non valid
    
    
    if(icell==0) gstencil->resLR[i]=0.;
      
    temp=0.;
    res=0.;

    ct.p=gstencil->p[i+icell*stride];
    ct.d=gstencil->d[i+icell*stride];
    curcell=&ct;
    
    // Computing the laplacian ===========================
 
    for(inei=0;inei<6;inei++){
      int idxnei=gstencil->nei[i+vnei[inei]*stride];
      temp+=gstencil->p[idxnei+vcell[inei]*stride];
    }
 
    // setting up the residual
    res=temp;
    
    // we finish the laplacian
    temp=temp/6.0;
    temp=temp-dx*dx*curcell->d/6.0*factdens;
      
    

    // we finsih the residual
    res=res-6.0*curcell->p;
    res=res/(dx*dx)-factdens*curcell->d;

    // we store the new value of the potential
    gstencil->pnew[i+icell*stride]=temp;

    // we store the local residual
    if(flag) {
      gstencil->res[i+icell*stride]=factdens*curcell->d*factdens*curcell->d;
      //gstencil->res[i+icell*stride]=1.;
    }
    else{
      gstencil->res[icell*stride+i]=res*res;
    }

    // the low resolution residual
    gstencil->resLR[i]+=res*0.125;
    
    // ready for the next cell
    if(icell==0) ival++;
    
    //ready for the next oct
  }
  }
}


//===============================================================================
__global__ void dev_update_pot_in_stencil(struct STENGRAV *gstencil, int level, int curcpu, int nread,int stride, int flag){

  int icell;

  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  int i=bx*blockDim.x+tx;


      
  if(gstencil->valid[i]!=1){
    
  }
  else{
    for(icell=0;icell<8;icell++){ // update
      gstencil->p[icell*stride+i]=gstencil->pnew[icell*stride+i];
    }
  }

}

#endif

//=========================================================================================================


#ifndef FASTGRAV
REAL PoissonJacobiGPU(int level,struct RUNPARAMS *param, struct OCT ** firstoct, struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim)
{
  REAL dxcur;
  int iter;
  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot;
  int nread;
  REAL residual,dres;
  int icell;
  int nitmax;
  REAL factdens;

  REAL *rloc;
  cudaMallocHost((void**)&rloc,sizeof(REAL)*cpu->nstream);

  /* double t[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; */
  /* double tgat=0.,tcal=0.,tscat=0.,tall=0.; */
  /* double tstart,tstop,tglob=0.,tup=0.; */
  /* double tt; */

  // ======================= some stuf for CUDPP =======================/
  REAL *resA;
  REAL *resB;


  cudaMalloc((void **)&resA,sizeof(REAL)*stride*8);
  cudaMalloc((void **)&resB,sizeof(REAL)*stride*8);

  
  struct CUPARAM cuparam;
  
  cudppCreate(&(cuparam.theCudpp));
    
  cuparam.config.algorithm = CUDPP_SCAN;
  cuparam.config.op = CUDPP_MAX;
  cuparam.config.datatype = CUDPP_DOUBLE;
  cuparam.config.options=CUDPP_OPTION_BACKWARD | CUDPP_OPTION_INCLUSIVE;
  cuparam.config.op=CUDPP_MAX;

  cuparam.scanplan =0;
  cudppPlan(cuparam.theCudpp,&(cuparam.scanplan), cuparam.config, stride*8, 1, 0);


  //======================== END CUDPP STUFF ========================/

  cudaStream_t stream[cpu->nstream];
  int is,offset;
  // creating the streams
  for(is=0;is<cpu->nstream;is++){
    cudaStreamCreate(&stream[is]);
  }



  // Computing the factor of the density
  if(level>=param->lcoarse){
#ifndef TESTCOSMO
    factdens=4.0*M_PI;
#else
    //    factdens=6.0*tsim; WARNING JUST TESTING WITHOUT TSIM!!!
    factdens=6.0;
#endif
  }
  else{
    factdens=1.;
  }

  // Computing the max number for iteration

  if((level==param->mgridlmin)||(level>param->lcoarse)){
    nitmax=param->niter;
  }
  else{
    nitmax=param->nrelax;
  }

  dxcur=pow(0.5,level);


  // Scanning the Octs

  for(iter=0;iter<nitmax;iter++){
    // --------------- some inits for iterative solver
    residual=0.;

    // --------------- setting the first oct of the level
    nextoct=firstoct[level-1];
    nreadtot=0;
    if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
      do {
	curoct=nextoct;
	nextoct=curoct->next; 
	
	// ------------ gathering the stencil value values


	// find a way to stick data on the device
	
	nextoct=gatherstencilgrav(curoct,stencil,stride,cpu, &nread);
	dim3 gridoct((nread/(cpu->nthread*cpu->nstream)>1?(nread/(cpu->nthread*cpu->nstream)):1));
	dim3 blockoct(cpu->nthread);
		     
	// streaming data
	for(is=0;is<cpu->nstream;is++){
	  offset=is*(nread/cpu->nstream);
	  
	  cudaMemcpyAsync(cpu->dev_stencil+offset,stencil->stencil+offset,nread*sizeof(struct GGRID)/cpu->nstream,cudaMemcpyHostToDevice,stream[is]);
	  
	  // ------------ solving the hydro
	  dev_PoissonJacobi_single<<<gridoct,blockoct,0,stream[is]>>>(cpu->dev_stencil+offset,level,cpu->rank,nread/cpu->nstream,stride,dxcur,(iter==0),factdens,resA+offset*8,cpu->res+offset*8,cpu->pnew+offset*8,cpu->resLR+offset);
	  // ------------ computing the residuals

	  cudaStreamSynchronize(stream[is]);
	  cudppScan(cuparam.scanplan, resB+offset*8, resA+offset*8, nread*8/(cpu->nstream));
	  cudaMemcpyAsync(rloc+is,resB+offset*8,sizeof(REAL),cudaMemcpyDeviceToHost,stream[is]);

	  cudaMemcpyAsync(stencil->res+offset*8,cpu->res+offset*8,(nread/cpu->nstream)*sizeof(REAL)*8,cudaMemcpyDeviceToHost,stream[is]);   
	  cudaMemcpyAsync(stencil->pnew+offset*8,cpu->pnew+offset*8,(nread/cpu->nstream)*sizeof(REAL)*8,cudaMemcpyDeviceToHost,stream[is]);   
	  cudaMemcpyAsync(stencil->resLR+offset,cpu->resLR+offset,(nread/cpu->nstream)*sizeof(REAL),cudaMemcpyDeviceToHost,stream[is]);   

	  cudaStreamSynchronize(stream[is]);
	  residual=(residual>rloc[is]?residual:rloc[is])*(iter!=0);
	 
	}



	// ------------ scatter back the data

	nextoct=scatterstencilgrav(curoct,stencil, nread, stride,cpu);

	nreadtot+=nread;
	
      }while(nextoct!=NULL);
    }
    
    
    // at this stage an iteration has been completed : let's update the potential and compute the residual
    if(nreadtot>0){
      curoct=firstoct[level-1];
      if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
    	nextoct=curoct;
    	do{
    	  curoct=nextoct;
    	  nextoct=curoct->next;
    	  if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
    	  for(icell=0;icell<8;icell++){
    	    curoct->cell[icell].gdata.p=curoct->cell[icell].pnew;
    	  }
    	}while(nextoct!=NULL);
      }
    }

    if(iter>0){
      dres=sqrt(residual);
      if((dres)<param->poissonacc){
	if(level>=param->lcoarse) break;
      }
    }
    
  }
  //  printf("GPU | level=%d iter=%d res=%e tgat=%e tcal=%e tscat=%e tall=%e tup=%e tglob=%e\n",level,iter,dres,tgat/iter,tcal/iter,tscat/iter,tall/iter,tup/iter,tglob/iter);
  printf("GPU | level=%d iter=%d res=%e \n",level,iter,dres);

  cudppDestroyPlan(cuparam.scanplan);
  cudppDestroy(cuparam.theCudpp);
  cudaFree(resA);
  cudaFree(resB);
  cudaFreeHost(rloc);

  // Destroying the streams
  for(is=0;is<cpu->nstream;is++){
    cudaStreamDestroy(stream[is]);
  }
  return dres;
}

#else



//===============================================================================
REAL PoissonJacobiGPU(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim)
{
  REAL dxcur;
  int iter;
  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot;
  int nread;
  REAL fnorm,residual,dres;
  int nitmax;
  REAL factdens;
  REAL rloc;
  REAL res0;



  struct CUPARAM cuparam;
  
  cudppCreate(&(cuparam.theCudpp));
  
  cuparam.config.algorithm = CUDPP_SCAN;
  cuparam.config.datatype = CUDPP_DOUBLE;
  cuparam.config.options=CUDPP_OPTION_BACKWARD | CUDPP_OPTION_INCLUSIVE;
  cuparam.config.op=CUDPP_ADD;

  cuparam.scanplan =0;
  cudppPlan(cuparam.theCudpp,&(cuparam.scanplan), cuparam.config, stride*8, 1, 0);

  cudaStream_t stream;
  cudaStreamCreate(&stream);

  // Computing the factor of the density
  if(level>=param->lcoarse){
#ifndef TESTCOSMO
    factdens=4.0*M_PI;
#else
    factdens=6.0;
#endif
  }
  else{
    factdens=1.;
  }

  // Computing the max number for iteration

  if((level==param->mgridlmin)||(level>param->lcoarse)){
    nitmax=param->niter;
  }
  else{
    nitmax=param->nrelax;
  }

  dxcur=pow(0.5,level);
  
  double tall=0.,tcal=0.,tscat=0.,tgat=0.;
  double tglob=0.,tup=0.;

  double tstart,tstop,tt;
  //if(level==6) nitmax=10;

  int ngrid=(stride/(cpu->nthread)>1?(stride/(cpu->nthread)):1);
  dim3 gridoct(ngrid);
  dim3 blockoct(cpu->nthread);

  //printf("g=%d b=%d\n",stride/cpu->nthread,cpu->nthread);
  //printf("Start Error =%s\n",cudaGetErrorString(cudaGetLastError()));
  for(iter=0;iter<nitmax;iter++){
    tstart=MPI_Wtime();
    // --------------- setting the first oct of the level
    nextoct=firstoct[level-1];
    nreadtot=0;

    double temps[10];

    // --------------- some inits for iterative solver
    if(iter==0){
      fnorm=0.;
      residual=0.;
      nread=0;
    }
    else{
      residual=0.;
    }

    // Scanning the octs

    if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
      do {
	curoct=nextoct;
		
	
	if((iter==0)||(cpu->noct[level-1]>nread)){
	  // ------------ some cleaning
	  clean_vecpos(level,firstoct);
	  if(level>param->lcoarse) clean_vecpos(level-1,firstoct); // for two level interactions
	  
	  //------------- gathering neighbour structure
	  nextoct=gatherstencilgrav_nei(curoct,stencil,stride,cpu, &nread);
	  
	  // ------------ gathering the stencil value 
	  gatherstencilgrav(firstoct[level-1],stencil,stride,cpu, &nread); //  the whole level must be parsed to recover the values
	  // Super Big-ass copy toward the device
	  cudaMemcpy(cpu->dev_stencil->d,stencil->d,stride*sizeof(REAL)*8,cudaMemcpyHostToDevice);
	  cudaMemcpy(cpu->dev_stencil->p,stencil->p,stride*sizeof(REAL)*8,cudaMemcpyHostToDevice);
	  cudaMemcpy(cpu->dev_stencil->nei,stencil->nei,stride*sizeof(int)*7,cudaMemcpyHostToDevice);
	  cudaMemcpy(cpu->dev_stencil->level,stencil->level,stride*sizeof(int)*1,cudaMemcpyHostToDevice);
	  cudaMemcpy(cpu->dev_stencil->cpu,stencil->cpu,stride*sizeof(int)*1,cudaMemcpyHostToDevice);
	  cudaMemcpy(cpu->dev_stencil->valid,stencil->valid,stride*sizeof(char)*1,cudaMemcpyHostToDevice);

	  cudaMemset(cpu->dev_stencil->res,0,stride*sizeof(REAL)*8);
	  cudaMemset(cpu->dev_stencil->resLR,0,stride*sizeof(REAL));
	  cudaMemset(cpu->dev_stencil->pnew,0,stride*sizeof(REAL)*8);

	  
	}
	else{ 
	  nextoct=NULL;
	} 
	
	// Jacobi iteration
	dev_PoissonJacobi_single<<<gridoct,blockoct>>>(cpu->gpu_stencil,level,cpu->rank,nread,stride,dxcur,(iter==0),factdens);

	// ------------ computing the residuals


	// ==== FIX COMPUTING RESIDUAL
	
	cudaStreamSynchronize(stream);
	cudppScan(cuparam.scanplan, cpu->dev_stencil->res2, cpu->dev_stencil->res, stride*8);

	if(iter==0){
	  cudaMemcpy(&rloc,cpu->dev_stencil->res2,sizeof(REAL),cudaMemcpyDeviceToHost);
	  fnorm+=rloc;
	}
	else{
	  //residual=(residual>rloc?residual:rloc);
	  if((iter==(nitmax-1))||(level>=param->lcoarse))
	    {
	      cudaMemcpy(&rloc,cpu->dev_stencil->res2,sizeof(REAL),cudaMemcpyDeviceToHost);
	      residual+=rloc;
	    }
	}
	
	// ------------ scatter back the data in the tree


	// we scatter back in the tree b/c of a small stencil
	if(cpu->noct[level-1]>nread){

	  cudaMemcpy(stencil->pnew,cpu->dev_stencil->pnew,stride*sizeof(REAL)*8,cudaMemcpyDeviceToHost);
	  cudaMemcpy(stencil->res,cpu->dev_stencil->res,stride*sizeof(REAL)*8,cudaMemcpyDeviceToHost);
	  cudaMemcpy(stencil->resLR,cpu->dev_stencil->resLR,stride*sizeof(REAL)*1,cudaMemcpyDeviceToHost);
	  
	  scatterstencilgrav(curoct,stencil,nread,stride, cpu);
	  
  

	} 
	else{ 
	  //large stencil
	  dev_update_pot_in_stencil<<<gridoct,blockoct>>>(cpu->gpu_stencil,level,cpu->rank,nread,stride,0); 
	} 

	tall+=temps[9]-temps[0];
	tcal+=temps[7]-temps[3];
	tscat+=temps[9]-temps[7];
	tgat+=temps[3]-temps[0];
	nreadtot+=nread;
      }while(nextoct!=NULL);
    }


    // at this stage an iteration has been completed : let's update the potential in the tree (small stencil only)
    if(cpu->noct[level-1]>nread){
      if(nreadtot>0) update_pot_in_tree(level,firstoct,cpu,param);
    }

    //printf("iter=%d nreadtot=%d\n",iter,nreadtot);
    if((iter==1)&&(level>=param->lcoarse)) res0=residual;
     
    tstop=MPI_Wtime();
    tup+=(tstop-tt);
    tglob+=(tstop-tstart);
    
    if(iter>0){
      if(level<param->lcoarse){
	dres=sqrt(residual);
      }
      else{
	dres=sqrt(residual/fnorm);
      }
      if((dres)<param->poissonacc){
	if(level>=param->lcoarse) break;
      }
    }
  }


  if(cpu->noct[level-1]==nread){ // large stencil
    cudaMemcpy(stencil->pnew,cpu->dev_stencil->pnew,stride*sizeof(REAL)*8,cudaMemcpyDeviceToHost);
    cudaMemcpy(stencil->res,cpu->dev_stencil->res,stride*sizeof(REAL)*8,cudaMemcpyDeviceToHost);
    cudaMemcpy(stencil->resLR,cpu->dev_stencil->resLR,stride*sizeof(REAL)*1,cudaMemcpyDeviceToHost);
    
    scatterstencilgrav(curoct,stencil,nread,stride, cpu);
    
    if(nreadtot>0) update_pot_in_tree(level,firstoct,cpu,param);
    tt=MPI_Wtime();
  }

  if(level>param->lcoarse){
    printf("GPU | level=%d iter=%d res=%e fnorm=%e\n",level,iter,dres,fnorm);
  }
  else{
    printf("GPU === | level=%d iter=%d res=%e fnorm=%e resraw=%e res0=%e\n",level,iter,dres,fnorm,sqrt(residual),sqrt(res0));
  }

  cudppDestroyPlan(cuparam.scanplan);
  cudppDestroy(cuparam.theCudpp);
  cudaStreamDestroy(stream);


  return dres;
}

#endif
