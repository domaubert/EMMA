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


struct PLIGHT{
  int keyp;
  float xp;
  float yp;
  float zp;
  float mp;
};


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


struct OCT* cell2oct(struct CELL* cell)
{
  long int adress;
  struct OCT *poct;
  adress=(long int) cell;
  adress=adress-cell->idx*sizeof(struct CELL);
  poct=(struct OCT*) adress;

  return poct;
}

//------------------------------------------------------------------------

void getcellnei(int cindex, int *neip, int *cell)
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


// ================================================================

struct OCT * cic_nei_oct(struct OCT* curoct, int cx, int cy, int cz)
{
  struct OCT* woct;
  struct OCT* newoct;
  int ioct;
  
  // computing the oct "index" for neighbor cells
  ioct=(cx==2)+((cy==2)<<1)+((cz==2)<<2);
  
  woct=NULL;

  switch(ioct){
  case 0: //-------------------------------------
    woct=curoct;
    break;
  case 1://-------------------------------------
      woct=curoct->nei[1]->child;
    break;
  case 2://-------------------------------------
    woct=curoct->nei[3]->child;
    break;
  case 4://-------------------------------------
    woct=curoct->nei[5]->child;
    break;
  case 3://-------------------------------------
    if(curoct->nei[1]->child!=NULL){
      woct=curoct->nei[1]->child->nei[3]->child;
    }
    else if(curoct->nei[3]->child!=NULL){
      woct=curoct->nei[3]->child->nei[1]->child;
    }
    break;
  case 5://-------------------------------------
    if(curoct->nei[1]->child!=NULL){
      woct=curoct->nei[1]->child->nei[5]->child;
    }
    else if(curoct->nei[5]->child!=NULL){
      woct=curoct->nei[5]->child->nei[1]->child;
    }
    break;
  case 6://-------------------------------------
    if(curoct->nei[3]->child!=NULL){
      woct=curoct->nei[3]->child->nei[5]->child;
    }
    else if(curoct->nei[5]->child!=NULL){
      woct=curoct->nei[5]->child->nei[3]->child;
    }
    break;
  case 7://-------------------------------------
    if(curoct->nei[1]->child!=NULL){
      if(curoct->nei[1]->child->nei[3]->child!=NULL){
	woct=curoct->nei[1]->child->nei[3]->child->nei[5]->child;
      }
      else if(curoct->nei[1]->child->nei[5]->child!=NULL){
	woct=curoct->nei[1]->child->nei[5]->child->nei[3]->child;
      }
    }
    
    if((curoct->nei[3]->child!=NULL)&&(woct==NULL)){
      if(curoct->nei[3]->child->nei[1]->child!=NULL){
	woct=curoct->nei[3]->child->nei[1]->child->nei[5]->child;
      }
      else if(curoct->nei[3]->child->nei[5]->child!=NULL){
	woct=curoct->nei[3]->child->nei[5]->child->nei[1]->child;
      }
    }
    
    if((curoct->nei[5]->child!=NULL)&&(woct==NULL)){
      if(curoct->nei[5]->child->nei[1]->child!=NULL){
	woct=curoct->nei[5]->child->nei[1]->child->nei[3]->child;
      }
      else if(curoct->nei[5]->child->nei[3]->child!=NULL){
	woct=curoct->nei[5]->child->nei[3]->child->nei[1]->child;
      }
    }
    break;
  }

  return woct;
}
// ====================================================================

void getparticles(struct OCT *curoct, int *keyp, float *xp, float *yp, float *zp, float *mp, int *ipart, int root, float ix0, float iy0, float iz0, int idxoct,float idxcur, int iicell){ 
  
  int icell,icell2;
  struct PART *nexp;
  struct PART *curp;
  float x0,y0,z0;
  float xpc,ypc,zpc;
  int cx,cy,cz;
  int i,j,k;
  float dxcur;
  struct OCT * woct;
  int inx,iny,inz;

  for(icell=0;icell<8;icell++) // looping over cells in oct
    {
      if(root){ // we start the recursion
	dxcur=1./(1<<(curoct->level)); // size of a cell
	x0=(curoct->x+(icell%2)*dxcur); //corner coordinates
	y0=(curoct->y+((icell/2)%2)*dxcur);
	z0=(curoct->z+(icell/4)*dxcur);
	iicell=icell;
      }
      else{
	x0=ix0;
	y0=iy0;
	z0=iz0;
	dxcur=idxcur;
      }

      // ===============  looking for cic neighbors at levels >= current level
      for(i=0;i<2;i++){
	for(j=0;j<2;j++){
	  for(k=0;k<2;k++){
	    
	    if((!root)*((i+j+k)!=0)) break; // for higher levels particle the 8 cells should not be explored

	    cx=(icell&1)+i;
	    cy=((icell>>1)&1)+j;
	    cz=(icell>>2)+k;
      
	    // getting the neighbor oct
	    woct=cic_nei_oct(curoct,cx,cy,cz);

	    // at this stage we have the recipient oct
	    // we recompute the cell index
	    icell2=(cx&1)+((cy&1)<<1)+((cz&1)<<2);
	
	    // gathering particles
	    if(woct!=NULL){
	      if(woct->cell[icell2].child!=NULL){ // the cell is refined we go one level further
		getparticles(woct->cell[icell2].child, keyp,xp,yp,zp,mp,ipart,0,x0,y0,z0,idxoct,dxcur,iicell); 
	      }
	      else{
		nexp=woct->cell[icell2].phead; //sweeping the particles of the current cell */

		if(nexp!=NULL){ 
		  do{  
		    curp=nexp; 
		    nexp=curp->next;
		    
		    xpc=curp->x-dxcur*0.5;xpc=(xpc<0?1.+xpc:xpc);
		    ypc=curp->y-dxcur*0.5;ypc=(ypc<0?1.+ypc:ypc);
		    zpc=curp->z-dxcur*0.5;zpc=(zpc<0?1.+zpc:zpc);

		    // testing the particle (assuming particles are centered)
		    // is the lower left corner inside ?

		    inx=((xpc-x0)>=0.)*((xpc-x0)<dxcur);
		    iny=((ypc-y0)>=0.)*((ypc-y0)<dxcur);
		    inz=((zpc-z0)>=0.)*((zpc-z0)<dxcur);
		    //printf("m=%e inx=%d iny=%d inz=%d xp=%e x0=%e\n",curp->mass,inx,iny,inz,ypc,y0);
		    if((inx*iny*inz)!=1) continue;
	    
		    keyp[*ipart]=idxoct*8+iicell; // computing the key of each particle
		    
		    // we compute the relative position to ensure a consistency between cell and particle
		    xp[*ipart]=xpc-x0;
		    yp[*ipart]=ypc-y0;
		    zp[*ipart]=zpc-z0;
		    mp[*ipart]=curp->mass;
		    (*ipart)=(*ipart)+1;
		  }while(nexp!=NULL); 
		}
	      }
	    }
	  }
	}
      }
    }
}

// ====================================================================

void call_cic_GPU(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu){

  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *woct;
  int icell;
  /* struct PART *curp; */
  /* struct PART *nexp; */
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
  //int keyoloc;
  struct OCT **keyodict; // dictionnary keyo -> pointer
  //  int ncart; // contains the max dimension along one direction


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

  if(cpu->rank==0) printf("==> start CIC on GPU\n");

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

  int idxoct;
  for(level=levelcoarse;level<=levelmax;level++)
    {
      dxcur=1./(1<<level); // size of a cell
      nextoct=firstoct[level-1];

      //ncart=1<<(level-1); // number of octs along one dimension

      // ========================   sweeping the octs sort the oct keys
      if(nextoct==NULL) continue;
      idxoct=0;
      ipart=0;
      memset(keyp,0,sizeof(int)*NPMAX);
      memset(xp,0,sizeof(int)*NPMAX);
      memset(yp,0,sizeof(int)*NPMAX);
      memset(zp,0,sizeof(int)*NPMAX);
      memset(mp,0,sizeof(int)*NPMAX);
      memset(keyodict,0,sizeof(struct OCT* )*NOMAX);


      // gathering particles from levels>=current level
      do 
	{	  
	  // =============== FIX LARGE OCT NUMBERS ===== !!!!

	  curoct=nextoct;
	  nextoct=curoct->next;

	  // we skip octs which do not belong to the current CPU (they will be considered through mpi)
	  if(curoct->cpu!=cpu->rank) continue;
	  keyodict[idxoct]=curoct; // building the dictionnary

	  // gathering particles from all levels >= current level
	  getparticles(curoct, keyp,xp,yp,zp,mp,&ipart,1,0.,0.,0.,idxoct,0.,0); 
	  idxoct++;
	}while(nextoct!=NULL);
	  

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

      cudaMemset(mass_d,0,sizeof(float)*NPMAX);
      cudaMemset(mass2_d,0,sizeof(float)*NPMAX);
      cudaMemset(flag_d,0,sizeof(unsigned int)*NPMAX);
      // ---------- kernels start

      // flag segments
  
      carte_flag_next<<<dimGridPart,dimBlockPart>>>(keyp_d,flag_d,ipart); 
      cudaThreadSynchronize();   
      GPU2CPU_UINT(flag,flag_d,sizeof(unsigned int)*NPMAX);

      carte_flag_previous<<<dimGridPart,dimBlockPart>>>(keyp_d,flag_d); 
      cudaThreadSynchronize();   

      // scanning the 8 CIC calculations
      
      int cx,cy,cz;

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
		  cudppSegmentedScan(plan, mass2_d, mass_d, flag_d, ipart);
		  cudaThreadSynchronize();   
     
		  // ------------ getting the data back
		  
		  GPU2CPU(massp2,mass2_d,sizeof(float)*NPMAX);

		  // ------------ scatter in the tree
		  
		  int idxoct,icell;

		  for(ip=0;ip<ipart;ip++){
		    if(flag[ip]==1){
		      idxoct=keyp[ip]>>3; // oct index
		      icell=keyp[ip]&7; // cell index
		      curoct=keyodict[idxoct]; // the current oct

		      cx=(icell&1)+i;
		      cy=((icell>>1)&1)+j;
		      cz=(icell>>2)+k;

		      // getting the neighbor oct
		      woct=cic_nei_oct(curoct,cx,cy,cz);

		      // at this stage we have the recipitent oct
		      // we recompute the cell index
		      icell=(cx&1)+((cy&1)<<1)+((cz&1)<<2);
		      
		      if(woct!=NULL){
			woct->cell[icell].density+=massp2[ip];
		      }
		    }
		  }

		}
	    }
	}
      // ------------- Done

      cudppDestroyPlan(plan);

    } // going to next level

  
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

