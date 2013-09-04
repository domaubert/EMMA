#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include <mpi.h>
#include <cudpp.h>

#ifdef WGRAV

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
extern "C" void mpi_exchange_level(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field, int cmp_keys, int level);


// ====================== structure for CUDPP =======


struct CUPARAM{
  CUDPPHandle theCudpp;
  CUDPPConfiguration config;
  CUDPPHandle scanplan;
};


// =======================================================
// =======================================================

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
  //abort();
}

// ===============================================================
// ===============================================================

void destroy_pinned_gravstencil(struct STENGRAV *gstencil, int stride){

  cudaFreeHost(gstencil->stencil);
  cudaFreeHost(gstencil->pnew);
  cudaFreeHost(gstencil->res);
#ifdef ONFLYRED
  cudaFreeHost(gstencil->resLR);
#endif
}


// ===============================================================
// ===============================================================


void create_gravstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaMalloc((void **)&(cpu->dev_stencil),sizeof(struct GGRID)*stride);
  cudaMalloc((void **)&(cpu->res),sizeof(REAL)*stride*8);
  cudaMalloc((void **)&(cpu->pnew),sizeof(REAL)*stride*8);
  cudaMalloc((void **)&(cpu->resLR),sizeof(REAL)*stride);
}

void destroy_gravstencil_GPU(struct CPUINFO *cpu, int stride){
  cudaFree(cpu->dev_stencil);
  cudaFree(cpu->res);
  cudaFree(cpu->pnew);
  cudaFree(cpu->resLR);
  
}

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
      vres[icell+i*8]=factdens*curcell->d*factdens*curcell->d; 
      stockres[icell+i*8]=factdens*curcell->d;
    }
    else{
      vres[icell+i*8]=res*res; 
      stockres[icell+i*8]=res;
    }

    stockresLR[i]+=res*0.125;
  }


}

//=========================================================================================================

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
  REAL fnorm,res0=0.;
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
  cuparam.config.datatype = CUDPP_DOUBLE;
  cuparam.config.options=CUDPP_OPTION_BACKWARD | CUDPP_OPTION_INCLUSIVE;
  cuparam.config.op=CUDPP_ADD;

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
  
  fnorm=0.;
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
	  //residual=(residual>rloc[is]?residual:rloc[is])*(iter!=0);
	  if(iter==0){
	    fnorm+=rloc[is];
	  }
	  else{
	    residual+=rloc[is];
	  }
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

#ifdef WMPI
    //printf("iter=%d\n",iter);
    if((iter<=param->niter)||(iter%10==0)){
      mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,(iter==0),level); // potential field exchange
      if(iter==0){
	//if(level==7) printf("rank=%d fnorm=%e\n",cpu->rank,fnorm);
	MPI_Allreduce(MPI_IN_PLACE,&fnorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
      else{
	MPI_Allreduce(MPI_IN_PLACE,&residual,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
    }
#endif


    if((iter==1)&&(level>=param->lcoarse)) res0=residual;
    
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


  if(level>param->lcoarse){
    printf("GPU | level=%d iter=%d res=%e fnorm=%e\n",level,iter,dres,fnorm);
  }
  else{
    printf("GPU === | level=%d iter=%d res=%e fnorm=%e resraw=%e res0=%e\n",level,iter,dres,fnorm,sqrt(residual),sqrt(res0));
  }


  //  printf("GPU | level=%d iter=%d res=%e tgat=%e tcal=%e tscat=%e tall=%e tup=%e tglob=%e\n",level,iter,dres,tgat/iter,tcal/iter,tscat/iter,tall/iter,tup/iter,tglob/iter);
  //printf("GPU | level=%d iter=%d res=%e \n",level,iter,dres);

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


#endif
