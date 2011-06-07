
#include <cudpp.h>

#define NTHREAD 64


extern "C" void remove_valvec_GPU(float *vec, int nval, int stride, float avg, int level, int *vecl);
extern "C" float square_vec_GPU(float *vec, int nval, int stride, int level, int curcpu,int *vecl, int *veccpu, float *vec2, float *vecsum);
extern "C" float laplacian_vec2_GPU(float *vecden,float *vecpot,float *vecpotnew,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level,int curcpu, int nread,int stride,float dx,float factdens, float *vres, float *vecsum);
extern "C" int residual_vec2_GPU(float *vecden,float *vecpot,float *vres,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,float dx,float factdens);

//============================================================================
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

//============================================================================
__global__ void dev_remove_valvec_GPU(float *vec, int nval, int stride, float avg, int level, int *vecl)
{
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  
  int ioct;
  int icell;

  ioct=bx*blockDim.x+tx;

  for(icell=0;icell<8;icell++){
      vec[ioct+icell*stride]=(vec[ioct+icell*stride]-avg)*(level==vecl[ioct]);
  }
}


void remove_valvec_GPU(float *vec, int nval, int stride, float avg, int level, int *vecl){
  
  dim3 gridoct(stride/NTHREAD);
  dim3 blockoct(NTHREAD);

  dev_remove_valvec_GPU<<<gridoct,blockoct>>>(vec,nval,stride,1.,level,vecl);

}

//============================================================================
__global__ void dev_square_vec_GPU(float *vec, int nval, int stride, int level, int curcpu, int *vecl, int *veccpu, float *vec2)
{
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  
  int ioct;
  int icell;

  ioct=bx*blockDim.x+tx;

  for(icell=0;icell<8;icell++){
#ifdef WMPI
    vec2[ioct+icell*stride]=vec[ioct+icell*stride]*vec[ioct+icell*stride]*(level==vecl[ioct])*(curcpu==veccpu[ioct]);
#else
    vec2[ioct+icell*stride]=vec[ioct+icell*stride]*vec[ioct+icell*stride]*(level==vecl[ioct]);
#endif
  }
}


float square_vec_GPU(float *vec, int nval, int stride, int level, int curcpu, int *vecl, int *veccpu, float *vec2, float *vecsum)
{
  
  float res;

  // ==== some config for reduction
  CUDPPConfiguration config;
  CUDPPHandle planadd;
  config.algorithm = CUDPP_SCAN;
  config.op = CUDPP_ADD;
  config.datatype = CUDPP_FLOAT;
  CUDPPOption direction = CUDPP_OPTION_BACKWARD;
  CUDPPOption inclusivity = CUDPP_OPTION_INCLUSIVE;
  config.options = direction | inclusivity;
  cudppPlan(&planadd, config, stride*8, 1, 0);


  // ==== we compute the square of the array

  dim3 gridoct(stride/NTHREAD);
  dim3 blockoct(NTHREAD);

  cudaMemset(vecsum,0,stride*8*sizeof(float));
  dev_square_vec_GPU<<<gridoct,blockoct>>>(vec,nval,stride,level,curcpu,vecl,veccpu,vec2);
  cudppScan(planadd, vecsum, vec2, stride*8);
  cudaMemcpy(&res,vecsum,sizeof(float),cudaMemcpyDeviceToHost);

  //finish and return
  cudppDestroyPlan(planadd);

  return res;

}

// ================================================================================================================

__global__ void dev_laplacian_vec2_GPU(float *vecden,float *vecpot,float *vecpotnew,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int nread,int stride,float dx,float factdens, float *vres, int curcpu){

  int inei,icell,icellcoarse;
  float temp;
  float res;
  int vnei[6],vcell[6];
  int vneic[6],vcellc[6];
  float ominterp=0.2;
  
  int ioct;
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;

  ioct=bx*blockDim.x+tx;
  for(icell=0;icell<8;icell++){ // we scan the cells

    // we skip octs which do not belong to the current level
    if(vecl[ioct]!=level){
      vecpotnew[ioct+icell*stride]=vecpot[ioct+icell*stride]; // nevertheless we copy the (fixed) potential in the new potential
      vres[ioct+icell*stride]=0.;
      continue; 
    }

#ifdef WMPI
    // we skip octs which do not belong to the current cpu
    if(veccpu[ioct]!=curcpu){
      vecpotnew[ioct+icell*stride]=vecpot[ioct+icell*stride]; // nevertheless we copy the (fixed) potential in the new potential
      vres[ioct+icell*stride]=0.;
      continue; 
    }
#endif

    temp=0.;
    getcellnei_gpu(icell, vnei, vcell); // we get the neighbors

    // we compute the neighbor part of the laplacian
    for(inei=0;inei<6;inei++){ // we scan the neighbors
      if(vnei[inei]==6){ // the neighbor is in the current oct
	temp+=vecpot[ioct+vcell[inei]*stride];
      }
      else{
	if(vecl[vecnei[ioct+vnei[inei]*stride]]==vecl[ioct]){ // the octs share the same level
	  temp+=vecpot[vecnei[ioct+vnei[inei]*stride]+vcell[inei]*stride];
	}
	else{ // mixing values from two different levels
	  //
	  icellcoarse=vecicoarse[ioct];
	  getcellnei_gpu(icellcoarse,vneic,vcellc);
	  temp+=vecpot[vecnei[ioct+vnei[inei]*stride]+vcellc[inei]*stride]*(1.-ominterp)+vecpot[ioct+icell*stride]*ominterp;
	}
      }
    }

    //we setup the residual
    res=temp;
    
    // we finish the laplacian
    temp=temp/6.0f;
    temp=temp-dx*dx*vecden[ioct+icell*stride]/6.*factdens;

    // we finsih the residual
    res=res-6.0*vecpot[ioct+icell*stride];
    res=res/(dx*dx)-factdens*vecden[ioct+icell*stride];
    
    // we store the square of the residual
    vres[ioct+icell*stride]=res*res;
    
    // we store the new value
    vecpotnew[ioct+icell*stride]=temp;

    // ready for the next cell
  }

}


// ================================================================================================================

__global__ void dev_residual_vec2_GPU(float *vecden,float *vecpot,float *vres,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int nread,int stride,float dx,float factdens, int curcpu){

  int inei,icell,icellcoarse;
  float temp;
  float res;
  int vnei[6],vcell[6];
  int vneic[6],vcellc[6];
  float ominterp=0.2;
  
  int ioct;
  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;

  ioct=bx*blockDim.x+tx;
  for(icell=0;icell<8;icell++){ // we scan the cells

    // we skip octs which do not belong to the current level
    if(vecl[ioct]!=level){
      vres[ioct+icell*stride]=0.;
      continue; 
    }

#ifdef WMPI
    // we skip octs which do not belong to the current cpu
    if(veccpu[ioct]!=curcpu){
      vres[ioct+icell*stride]=0.;
      continue; 
    }
#endif


    temp=0.;
    getcellnei_gpu(icell, vnei, vcell); // we get the neighbors

    // we compute the neighbor part of the laplacian
    for(inei=0;inei<6;inei++){ // we scan the neighbors
      if(vnei[inei]==6){ // the neighbor is in the current oct
	temp+=vecpot[ioct+vcell[inei]*stride];
      }
      else{
	if(vecl[vecnei[ioct+vnei[inei]*stride]]==vecl[ioct]){ // the octs share the same level
	  temp+=vecpot[vecnei[ioct+vnei[inei]*stride]+vcell[inei]*stride];
	}
	else{ // mixing values from two different levels
	  //
	  icellcoarse=vecicoarse[ioct];
	  getcellnei_gpu(icellcoarse,vneic,vcellc);
	  temp+=vecpot[vecnei[ioct+vnei[inei]*stride]+vcellc[inei]*stride]*(1.-ominterp)+vecpot[ioct+icell*stride]*ominterp;
	}
      }
    }

    //we setup the residual
    res=temp;


    // we finsih the residual
    res=res-6.0*vecpot[ioct+icell*stride];
    res=res/(dx*dx)-factdens*vecden[ioct+icell*stride];
    
    // we store the residual
    vres[ioct+icell*stride]=res;
    
  }

}

//============================================================================

float laplacian_vec2_GPU(float *vecden,float *vecpot,float *vecpotnew,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,float dx,float factdens, float *vres, float *vecsum){
  
  float res;

  // ==== some config for reduction
  CUDPPConfiguration config;
  CUDPPHandle planadd;
  config.algorithm = CUDPP_SCAN;
  config.op = CUDPP_ADD;
  config.datatype = CUDPP_FLOAT;
  CUDPPOption direction = CUDPP_OPTION_BACKWARD;
  CUDPPOption inclusivity = CUDPP_OPTION_INCLUSIVE;
  config.options = direction | inclusivity;
  cudppPlan(&planadd, config, stride*8, 1, 0);


  // ==== we perform the laplacian

  dim3 gridoct(stride/NTHREAD);
  dim3 blockoct(NTHREAD);


  dev_laplacian_vec2_GPU<<<gridoct,blockoct>>>(vecden,vecpot,vecpotnew,vecnei,vecl,vecicoarse,veccpu,level,nread,stride,dx,factdens,vres,curcpu);

  cudppScan(planadd, vecsum, vres, stride*8);
  cudaMemcpy(&res,vecsum,sizeof(float),cudaMemcpyDeviceToHost);


  //finish and return
  cudppDestroyPlan(planadd);

  return res;
  
}



//============================================================================

int residual_vec2_GPU(float *vecden,float *vecpot,float *vres,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,float dx,float factdens){

  // ==== we perform the laplacian

  dim3 gridoct(stride/NTHREAD);
  dim3 blockoct(NTHREAD);


  dev_residual_vec2_GPU<<<gridoct,blockoct>>>(vecden,vecpot,vres,vecnei,vecl,vecicoarse,veccpu,level,nread,stride,dx,factdens,curcpu);
  cudaThreadSynchronize(); 

  return 0;
  
}
