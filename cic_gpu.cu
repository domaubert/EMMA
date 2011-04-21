#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <cutil.h>
#include <cudpp.h>

#define NPMAX 2097152
#define NOMAX 262144
#define NPBLOCK 64

extern "C" void call_cic_GPU(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu);

//#define CUERR() printf("\n %s \n",cudaGetErrorString(cudaGetLastError()))

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

// ====================================================================

__global__ void carte_mass(float dx, float *xp, float *yp, float *zp, float* md, float *mass, int offx, int offy, int offz)
{
	int thx=threadIdx.x;
	int bx=blockIdx.x;

	// Charge les positions depuis la mémoire globale vers la mémoire partagée
	int dataPos = bx*blockDim.x + thx; // position du thread dans data
	
	/* float x = xp[dataPos]/dx; */
	/* float y = yp[dataPos]/dx; */
	/* float z = zp[dataPos]/dx; */

	/* int xc=(int)x; */
	/* int yc=(int)y; */
	/* int zc=(int)z; */

	float tx,ty,tz;
	tx=(1.-offx)+(2.*offx-1.)*xp[dataPos]/dx;
	ty=(1.-offy)+(2.*offy-1.)*yp[dataPos]/dx;
	tz=(1.-offz)+(2.*offz-1.)*zp[dataPos]/dx;

	/* tx=(1.-offx)+(2.*offx-1.)*(x-xc); */
	/* ty=(1.-offy)+(2.*offy-1.)*(y-yc); */
	/* tz=(1.-offz)+(2.*offz-1.)*(z-zc); */
	
	mass[dataPos]=tx*ty*tz/(dx*dx*dx)*md[dataPos]; // 
}

// ====================================================================

__global__ void carte_flag_previous(int *d_data, unsigned int *d_flag  ) {

        unsigned int bx=blockIdx.x;
        unsigned int tx=threadIdx.x;
        unsigned int flag;

        int idglob=bx*blockDim.x+tx;
        uint current =d_data[idglob]; 
        uint previous=d_data[idglob-1]; 

        flag=((current!=previous)||(idglob==0)); 
        d_flag[idglob]=flag;
}


// ====================================================================

__global__ void carte_flag_next(int *d_data, unsigned int *d_flag, int np ) {

        unsigned int bx=blockIdx.x;
        unsigned int tx=threadIdx.x;
        unsigned int flag;

        int idglob=bx*blockDim.x+tx;
        uint current =d_data[idglob]; 
        uint next=d_data[idglob+1]; 

        flag=((current!=next)||(idglob==(np-1))); 
        d_flag[idglob]=flag;
}


// ====================================================================
// ====================================================================

void call_cic_GPU(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu){

  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  float dxcur;
  int ipart,ip;
  int i,j,k;
  
  //  ========================= CPU STUFF =====================

  int *keyp; //contains the particle keys
  unsigned int *flag; //contains the particle flags

  float *xp; //contains the particle x
  float *yp; //contains the particle y
  float *zp; //contains the particle z
  float *mp; //contains the particle mass
  float *massp2; //contains the particle cumulated mass

  //  int keyo[NOMAX]; // contains the oct key
  int keyoloc;
  struct OCT **keyodict; // dictionnary keyo -> pointer
  int ncart; // contains the max dimension along one direction


  //  ========================= GPU STUFF =====================
  // dimension configuration for Particle treatment
  dim3 dimGridPart(NPMAX/NPBLOCK); // np/64
  dim3 dimBlockPart(NPBLOCK); // 64

  unsigned int *flag_d;

  int *keyp_d;
  float *xp_d;
  float *yp_d;
  float *zp_d;
  float *m_d;
  float *mass_d;
  float *mass2_d;

  if(cpu->rank==0) printf("==> start CIC\n");

 // alloc 2

  keyp=(int*)calloc(NPMAX,sizeof(int));
  flag=(unsigned int*)calloc(NPMAX,sizeof(unsigned int));
  xp=(float*)calloc(NPMAX,sizeof(float));
  yp=(float*)calloc(NPMAX,sizeof(float));
  zp=(float*)calloc(NPMAX,sizeof(float));
  mp=(float*)calloc(NPMAX,sizeof(float));
  massp2=(float*)calloc(NPMAX,sizeof(float));
  keyodict=(struct OCT **)calloc(NOMAX,sizeof(struct OCT *));

  cudaMalloc((void **)&mass_d,sizeof(float)*NPMAX);
  cudaMalloc((void **)&mass2_d,sizeof(float)*NPMAX);
  cudaMalloc((void **)&flag_d,sizeof(unsigned int)*NPMAX);
  cudaMalloc((void **)&keyp_d,sizeof(int)*NPMAX);
  cudaMalloc((void **)&xp_d,sizeof(float)*NPMAX);
  cudaMalloc((void **)&yp_d,sizeof(float)*NPMAX);
  cudaMalloc((void **)&zp_d,sizeof(float)*NPMAX);
  cudaMalloc((void **)&m_d,sizeof(float)*NPMAX);

  //  ========================= GPU STUFF =====================




  // ========================== First we clean the density
  for(level=levelmax;level>=levelcoarse;level--)
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) curoct->cell[icell].density=0.;
	}while(nextoct!=NULL);
    }

  // =========================== Second start CIC
  

  for(ip=0;ip<NPMAX;ip++) keyp[ip]=-1; // initialize the particle keys
  for(ip=0;ip<NOMAX;ip++){
    keyodict[ip]=NULL;
  }

  ip=0;
  for(level=levelcoarse;level<=levelmax;level++)
    {
      dxcur=1./(1<<level); // size of a cell
      nextoct=firstoct[level-1];

      ncart=1<<(level-1); // number of octs along one dimension

      // ========================   sweeping the octs sort the oct keys
      if(nextoct==NULL) continue;
      do 
	{	  // =============== FIX LARGE OCT NUMBERS ===== !!!!

	  curoct=nextoct;
	  nextoct=curoct->next;

	  // we skip octs which do not belong to the current CPU (they will be considered through mpi)
	  if(curoct->cpu!=cpu->rank) continue;
	
	  // we compute the oct cartesian key
	  keyoloc=(int)(curoct->x/dxcur*0.5)+(int)(curoct->y/dxcur*0.5)*ncart+(int)(curoct->z/dxcur*0.5)*ncart*ncart; // the cartesian key of the oct
	  keyodict[keyoloc]=curoct; // building the dictionnary
	  
	}while(nextoct!=NULL);
	  
      
      // ==================== sweeping the dictionnary to key the particles in the right order
      
      ipart=0;
      for(ip=0;ip<NOMAX;ip++){
	if(keyodict[ip]==NULL) continue;
	curoct=keyodict[ip];
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
	    if(nexp!=NULL){ 
	      //printf("ip=%d\n",ip);
	      do{  
		curp=nexp; 
		nexp=curp->next;
		
		keyp[ipart]=ip*8+icell; // computing the key of each particle
		//printf("%d\n",keyp[ipart]);
		// we compute the relative position to ensure a consistency between cell and particle
		xp[ipart]=curp->x-(curoct->x+(icell%2)*dxcur);
		yp[ipart]=curp->y-(curoct->y+((icell/2)%2)*dxcur);
		zp[ipart]=curp->z-(curoct->z+(icell/4)*dxcur);
		mp[ipart]=curp->mass;
		if(curp->mass>0.5) printf("key=%d icell=%d %e %e %e %p\n",keyp[ipart],icell,curp->x,curp->y,curp->z,&curoct->cell);
		ipart++;
	      }while(nexp!=NULL); 
	    }
	  }
      }
      
      // at this stage each particle has been key-ed and are sorted


      // ==================================  launching the GPU horses =====================

      // --------- some config

      CUDPPConfiguration config;
      config.algorithm = CUDPP_SEGMENTED_SCAN;
      config.op = CUDPP_ADD;
      config.datatype = CUDPP_FLOAT;
      CUDPPOption direction = CUDPP_OPTION_FORWARD;
      CUDPPOption inclusivity = CUDPP_OPTION_INCLUSIVE;
      config.options = direction | inclusivity;
      
      CUDPPHandle plan;
      cudppPlan(&plan, config, ipart, 1, 0);

      // ---------- sending data

      CPU2GPU_INT(keyp_d,keyp,sizeof(int)*NPMAX);
      CPU2GPU(xp_d,xp,sizeof(float)*NPMAX);
      CPU2GPU(yp_d,yp,sizeof(float)*NPMAX);
      CPU2GPU(zp_d,zp,sizeof(float)*NPMAX);
      CPU2GPU( m_d,mp,sizeof(float)*NPMAX);

      // ---------- kernels start

      // flag segments
  
      carte_flag_next<<<dimGridPart,dimBlockPart>>>(keyp_d,flag_d,ipart); 
      cudaThreadSynchronize();   
      GPU2CPU_UINT(flag,flag_d,sizeof(unsigned int)*NPMAX);

      carte_flag_previous<<<dimGridPart,dimBlockPart>>>(keyp_d,flag_d); 
      cudaThreadSynchronize();   

      // scanning the 8 CIC calculations
      
      int cx,cy,cz;
      int ox,oy,oz;

      for(i=0;i<2;i++)
	{
	  for(j=0;j<2;j++)
	    {
	      for(k=0;k<2;k++)
		{
		  // segment scan
		  cudaMemset(mass_d,0,sizeof(float)*NPMAX);
		  carte_mass<<<dimGridPart,dimBlockPart>>>(dxcur,xp_d,yp_d,zp_d,m_d,mass_d,i,j,k);
		  cudaThreadSynchronize();   
		  //CUERR();
		  GPU2CPU(massp2,mass_d,sizeof(float)*NPMAX);
		  printf("mass=%e %e dx=%e i=%d j=%d k=%d %p\n",massp2[0],massp2[1],dxcur,i,j,k,mass_d);
		  cudppSegmentedScan(plan, mass2_d, mass_d, flag_d, ipart);
		  cudaThreadSynchronize();   
     
		  // ------------ getting the data back
		  
		  /* GPU2CPU(massp2,yp_d,sizeof(float)*NPMAX); */
		  /* printf("y=%e %e\n",massp2[0],massp2[1]); */
		  GPU2CPU(massp2,mass2_d,sizeof(float)*NPMAX);

		  // ------------ scatter in the tree
		  
		  int keyoct,icell;
		  
		  for(ip=0;ip<ipart;ip++){
		    if(flag[ip]==1){
		      keyoct=keyp[ip]>>3;
		      icell=keyp[ip]&7;
		      
		      cx=(icell&1)+i;
		      cy=((icell>>1)&1)+j;
		      cz=(icell>>2)+k;
		      
		      ox=keyoct%ncart;
		      oy=(keyoct/ncart)%ncart;
		      oz=keyoct/(ncart*ncart);
		      
		      if(cx==2){
			cx=0;
			ox=(ox+1)%ncart;
		      }

		      if(cy==2){
			cy=0;
			oy=(oy+1)%ncart;
		      }

		      if(cz==2){
			cz=0;
			oz=(oz+1)%ncart;
		      }
		      
		      keyoct=ox+oy*ncart+oz*ncart*ncart;
		      //		      if(mp[ip]>0.5) printf("keyoct=%d keyp=%d cx=%d cy=%d cz=%d ox=%d oy=%d oz=%d mp=%e\n",keyoct,keyp[ip],cx,cy,cz,ox,oy,oz,massp2[ip]);

		      icell=(cz<<2)+(cy<<1)+cx;
		      curoct=keyodict[keyoct];
		      //printf("curoct=%p keyoct=%d icell=%d mass=%e\n",curoct,keyoct,icell,massp2[ip]);
		      if(curoct!=NULL){
			curoct->cell[icell].density+=massp2[ip];
		      }
		    }
		  }
		  
		}
	    }
	}
      // ------------- Done

      cudppDestroyPlan(plan);

    }

  
  cudaFree(mass_d);
  cudaFree(mass2_d);
  cudaFree(xp_d);
  cudaFree(yp_d);
  cudaFree(zp_d);
  cudaFree(m_d);
  cudaFree(keyp_d);
  cudaFree(flag_d);

  free(flag);
  free(keyp);
  free(keyodict);
  free(xp);
  free(yp);
  free(zp);
  free(mp);
  free(massp2);
}

