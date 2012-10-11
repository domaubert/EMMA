#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include <mpi.h>
#include <cudpp.h>

#define NTHREAD 16

extern "C" struct OCT *gatherstencilgrav(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread);
extern "C" struct OCT *scatterstencilgrav(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu);
extern "C" REAL PoissonJacobiGPU(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL tsim);
extern "C" void create_stencil_GPU(struct CPUINFO *cpu, int stride);
extern "C" REAL comp_residual(struct HGRID *stencil, int level, int curcpu, int nread,int stride, int flag);


// ====================== structure for CUDPP =======


struct CUPARAM{
  CUDPPHandle theCudpp;
  CUDPPConfiguration config;
  CUDPPHandle scanplan;
};


// =======================================================
void create_stencil_GPU(struct CPUINFO *cpu, int stride){
  cudaMalloc((void **)&(cpu->dev_stencil),sizeof(struct HGRID)*stride);
  printf("%p\n",cpu->dev_stencil);
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

__global__ void dev_PoissonJacobi_single(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, int flag, REAL factdens, REAL *vres){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int inei,icell;
  int i;
  REAL temp;
  REAL res;
  int vnei[6],vcell[6];
  int ioct[7]={12,14,10,16,4,22,13};
  struct Gtype *curcell;
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;

  i=bx*blockDim.x+tx;


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
    stencil[i].New.cell[icell].pnew=temp;

    // we store the local residual
    if(flag) {
      //stencil[i].New.cell[icell].res=factdens*curcell->d;
      vres[icell*stride+i]=factdens*curcell->d;
    }
    else{
      vres[icell*stride+i]=res*res;
      //stencil[i].New.cell[icell].res=res;
    }
  }


}


//=========================================================================================================

REAL PoissonJacobiGPU(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL tsim)
{
  REAL dxcur;
  int iter;
  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot;
  int nread;
  REAL fnorm,residual,dres;
  int icell;
  int nitmax;
  REAL factdens;
  REAL rloc;



  cudaEvent_t start,stop;
  cudaEvent_t startc,stopc;

  double t[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double tgat=0.,tcal=0.,tscat=0.,tall=0.;
  double tstart,tstop,tglob=0.,tup=0.;
  double tt;

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
  cudppPlan(cuparam.theCudpp,&(cuparam.scanplan), cuparam.config, stride, 1, 0);


  //======================== END CUDPP STUFF ========================/

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventCreate(&startc);
  cudaEventCreate(&stopc);


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
    tstart=MPI_Wtime();
    cudaEventRecord(start,0);
    // --------------- some inits for iterative solver
    if(iter==0){
      fnorm=0.;
      residual=0.;
    }
    else{
      residual=0.;
    }

    // --------------- setting the first oct of the level
    nextoct=firstoct[level-1];
    nreadtot=0;
    if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
      do {
	curoct=nextoct;
	nextoct=curoct->next; 
	
	// ------------ gathering the stencil value values


	// find a way to stick data on the device
	t[0]=MPI_Wtime();
	
	nextoct=gatherstencilgrav(curoct,stencil,stride,cpu, &nread);
	cudaMemcpy(cpu->dev_stencil,stencil,nread*sizeof(struct HGRID),cudaMemcpyHostToDevice);  
       

	t[3]=MPI_Wtime();
	// ------------ solving the hydro
	dim3 gridoct((nread/NTHREAD>1?nread/NTHREAD:1));
	dim3 blockoct(NTHREAD);

	cudaEventRecord(startc,0);

	dev_PoissonJacobi_single<<<gridoct,blockoct>>>(cpu->dev_stencil,level,cpu->rank,nread,stride,dxcur,(iter==0),factdens,resA);
	
	cudaEventRecord(stopc,0);
	cudaEventSynchronize(stopc);
	t[6]=MPI_Wtime();

	// ------------ computing the residuals
	
	if(iter>0){
	  cudppScan(cuparam.scanplan, resB, resA, nread);
	  cudaMemcpy(&rloc,resB,sizeof(REAL),cudaMemcpyDeviceToHost);
	  printf("rloc=%e\n",sqrt(rloc));
	  residual=(residual>rloc?residual:rloc);
	}

	// ------------ scatter back the data

	cudaMemcpy(stencil,cpu->dev_stencil,nread*sizeof(struct HGRID),cudaMemcpyDeviceToHost);  
	nextoct=scatterstencilgrav(curoct,stencil, nread, cpu);
	t[8]=MPI_Wtime();

	nreadtot+=nread;
	tgat+=(t[3]-t[0]);
	tcal+=(t[6]-t[3]);
	tscat+=(t[8]-t[6]);
	tall+=(t[8]-t[0]);

      }while(nextoct!=NULL);
    }

    /* cudaEventRecord(stop,0); */
    /* cudaEventSynchronize(stop); */
    /* cudaEventElapsedTime(&time,start,stop); */
    
    tt=MPI_Wtime();
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

    tstop=MPI_Wtime();
    tglob+=(tstop-tstart);
    tup+=(tstop-tt);

    if(iter>0){
      dres=sqrt(residual);
      if((dres)<param->poissonacc){
	if(level>=param->lcoarse) break;
      }
    }


  }
  printf("GPU | level=%d iter=%d res=%e tgat=%e tcal=%e tscat=%e tall=%e tup=%e tglob=%e\n",level,iter,dres,tgat/iter,tcal/iter,tscat/iter,tall/iter,tup/iter,tglob/iter);

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  cudppDestroyPlan(cuparam.scanplan);
  cudppDestroy(cuparam.theCudpp);
  cudaFree(resA);
  cudaFree(resB);

  return dres;
}
