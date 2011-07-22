
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hilbert.h"
#include "prototypes.h"
#include "vector.h"
#include "io.h"
#include "cic.h"
#include "oct.h"
#include "particle.h"
#include "amr.h"
#include "tools.h"
#include "segment.h"
#include "communication.h"
#include "mgrid.h"
#include "friedmann.h"
#include <time.h>
#include <mpi.h>

#ifdef WGPU
#include "interface.h"
#include "vector_gpu.h"
#include "cic_gpu.h"
#endif

#ifdef WHYDRO
#include "hydro_utils.h"
#endif

// ===============================================================================



//------------------------------------------------------------------------
 // the MAIN CODE
 //------------------------------------------------------------------------

#ifdef TESTCOSMO
float f_aexp(float aexp, float omegam, float omegav)
{
  return 1./sqrtf(omegam/aexp+omegav*aexp*aexp);
}
#endif

int main(int argc, char *argv[])
{
  struct OCT *grid;
  struct OCT **firstoct;
  struct OCT **lastoct;

  int level,levelcoarse,levelmax,levelmin;
  int nvcycles;
  int nrelax;
  int ngridmax,ngrid;
  int npartmax;
  int cur,curnext; // flat indexes with boundaries
  int i,il,ichild,icell,inext,ii,ip,j;
  int xp,yp,zp;
  int NBND=1,NBND2=2*NBND;
  float dx;
  int vnei[6],vcell[6]; // arrays to get neighbors
  int vnei2[6],vcell2[6]; // arrays to get neighbors
  int vnei3[6],vcell3[6]; // arrays to get neighbors
  int neip[7]; // contains the index of the six neighbors of the current parent cell +current parent
  int ci,cj,ck;
  int cinext,cjnext,cknext;
  float threshold;
  float tsim=0.;
  struct OCT oct;
  struct OCT* nextoct;
  struct OCT* curoct;
  struct OCT* desoct;
  struct CELL * parcell;
  struct CELL * newcell;
  struct CELL * newcell2;
  int tag;
  float dxcur;
  float *dens;
  int firstoct_currl;
  int nxoct;
  int lmap;
  int npart;
  struct PART *part;
  struct PART *nexploc, *curploc;

  struct OCT *freeoct; // the first free oct available for refinement
  struct OCT *newoct;
  int nref=0,ndes=0;

  float xc,yc,zc;
  int stride;
  float **vcomp;
  int ncomp;
  float acc;
  float dt;
  int ntot=0,nlev,noct;
  float ntotd=0.,nlevd=0.;

  float disp,mdisp;
  
  int dir;

  char filename[128]; 
  FILE *fd;
  struct PART *nexp;
  struct PART *nexp2;
  struct PART *curp;

  struct PART *lastpart;

  int curc;
  float dtnew=0.;
  int nbnd;

  float x,y,z;
  float vx,vy,vz;
  float mass,mtot;
  float idx;
  float faexp, faexp2;
  float aexp;
  unsigned key;

  struct CPUINFO cpu;

  struct PACKET **sendbuffer; 
  struct PACKET **recvbuffer; 

  struct PART_MPI **psendbuffer; 
  struct PART_MPI **precvbuffer; 

  struct RUNPARAMS param;

  size_t rstat;


  float omegam,omegav,Hubble,omegab=0.045;
  float avgdens; 
#ifdef PIC
  avgdens=1.;//we assume a unit mass in a unit length box
#else
  avgdens=0.;
#endif
  //========== TEST ZONE (IF REQUIRED)==========

/*   printf("size =%d\n",sizeof(struct CELL)); */
/*   printf("size =%d\n",sizeof(struct OCT)); */
/*   abort(); */


  //=========== some initial calls =============
  GetParameters(argv[1],&param); // reading the parameters file
    

#ifdef WMPI
  MPI_Status stat;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(cpu.nproc));
  MPI_Comm_rank(MPI_COMM_WORLD,&(cpu.rank));

  //========= creating a PACKET MPI type =======
  MPI_Datatype MPI_PACKET,oldtypes[3]; 
  int          blockcounts[3];
  
  /* MPI_Aint type used to be consistent with syntax of */
  /* MPI_Type_extent routine */
  MPI_Aint    offsets[3], extent;
  
  
  /* Setup description of the 8 MPI_FLOAT fields data */
  offsets[0] = 0;
  oldtypes[0] = MPI_FLOAT;
  blockcounts[0] = 8;
  
  /* Setup description of the 2 MPI_INT fields key, level */
  /* Need to first figure offset by getting size of MPI_FLOAT */
  MPI_Type_extent(MPI_FLOAT, &extent);
  offsets[1] = 8 * extent;
  oldtypes[1] = MPI_LONG;
  blockcounts[1] = 1;

  MPI_Type_extent(MPI_LONG, &extent);
  offsets[2] = offsets[1]+extent;
  oldtypes[2] = MPI_INT;
  blockcounts[2] = 1;

  /* Now define structured type and commit it */
  MPI_Type_struct(3, blockcounts, offsets, oldtypes, &MPI_PACKET);
  MPI_Type_commit(&MPI_PACKET);

#ifdef PIC
  //========= creating a PART MPI type =======
  MPI_Datatype MPI_PART;

  /* Setup description of the 7 MPI_FLOAT fields x,y,z,vx,vy,vz */
  offsets[0] = 0;
  oldtypes[0] = MPI_FLOAT;
  blockcounts[0] = 7;
  
  /* Setup description of the 4 MPI_INT fields idx key level icell*/
  /* Need to first figure offset by getting size of MPI_FLOAT */
  MPI_Type_extent(MPI_FLOAT, &extent);
  offsets[1] = 7 * extent;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 3;

  MPI_Type_extent(MPI_INT, &extent);
  offsets[2] = 3 * extent+offsets[1];
  oldtypes[2] = MPI_LONG;
  blockcounts[2] = 1;

  /* Now define structured type and commit it */
  MPI_Type_struct(3, blockcounts, offsets, oldtypes, &MPI_PART);
  MPI_Type_commit(&MPI_PART);
  
#endif
  //============================================

  cpu.MPI_PACKET=&MPI_PACKET;
#ifdef PIC
  cpu.MPI_PART=&MPI_PART;
#endif
  cpu.comm=MPI_COMM_WORLD;
#else
  cpu.rank=0;
  cpu.nproc=1;
#endif

  //=========== assigning values =============
  levelcoarse=param.lcoarse;
  levelmax=param.lmax;
  levelmin=param.mgridlmin;
  nvcycles=param.nvcycles;
  nrelax=param.nrelax;

  ngridmax=param.ngridmax;

#ifdef PIC
  npartmax=param.npartmax;
#ifdef PART2
  npart=2;
#else
  npart=128*128*128;
#endif

#ifdef PARTN
  npart=32768;
#endif
#endif

  threshold=param.amrthresh;
  lmap=param.levelmap;
  stride=fmax(8,param.stride);//pow(2,levelcoarse);
  ncomp=10;
  acc=param.poissonacc;
  dt=param.dt;
  cpu.maxhash=param.maxhash;
  cpu.levelcoarse=levelcoarse;

  //breakmpi();
  //========== allocations ===================
  
  //  if(cpu.rank==0) printf("Allocating %f GB cell=%f GB part=%f GB book=%f",(sizeof(struct OCT)*ngridmax+sizeof(struct PART)*npart+cpu.maxhash*sizeof(struct OCT*)+stride*ncomp*sizeof(float))/(1024*1024*1024.),sizeof(struct OCT)*ngridmax/(1024*1024*1024.),sizeof(struct PART)*npart/(1024*1024*1024.),(cpu.maxhash*sizeof(struct OCT*)+stride*ncomp*sizeof(float))/(1024.*1024.*1024.));

  int memsize=0.;
  grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT)); memsize+=ngridmax*sizeof(struct OCT);// the oct grid
#ifdef PIC
  part=(struct PART*)calloc(npartmax,sizeof(struct PART)); memsize+=npartmax*sizeof(struct PART);// the particle array
#endif
 
  if(cpu.rank==0){
    printf(" === alloc Memory ===\n");
    printf(" oct size=%f ngridmax=%d\n",sizeof(struct OCT)/1024./1024.,ngridmax);
    printf(" grid = %f MB\n",(ngridmax/(1024*1024.))*sizeof(struct OCT));
    printf(" part = %f MB\n",(npartmax/(1024*1024.))*sizeof(struct PART));
  }

  firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *)); memsize+=levelmax*sizeof(struct OCT *);// the firstoct of each level
  lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *)); memsize+=levelmax*sizeof(struct OCT *);// the last oct of each level
  cpu.htable=(struct OCT**) calloc(cpu.maxhash,sizeof(struct OCT *)); memsize+=cpu.maxhash*sizeof(struct OCT*);// the htable keys->oct address
  cpu.noct=(int *)calloc(levelmax,sizeof(int)); memsize+=levelmax*sizeof(int);// the number of octs per level

  lastpart=part-1; // the last particle points before the first at the very beginning

  vcomp=(float **)calloc(ncomp,sizeof(float*));
  for(i=0;i<ncomp;i++)
    {
      vcomp[i]=(float *)calloc(stride,sizeof(float));
    }

  memsize+=ncomp*stride*sizeof(float);


  //===================================================================================================

  // allocating the vectorized tree
  
  struct MULTIVECT vectors;
  
#ifdef NEWJACK
  vectors.vecpot=(float*)calloc(stride*8,sizeof(float));
  vectors.vecpotnew=(float*)calloc(stride*8,sizeof(float));
  vectors.vecden=(float*)calloc(stride*8,sizeof(float));

#ifdef WHYDRO
  vectors.vec_d=(float*)calloc(stride*8,sizeof(float));
  vectors.vec_u=(float*)calloc(stride*8,sizeof(float));
  vectors.vec_v=(float*)calloc(stride*8,sizeof(float));
  vectors.vec_w=(float*)calloc(stride*8,sizeof(float));
  vectors.vec_p=(float*)calloc(stride*8,sizeof(float));

  vectors.vec_dnew=(float*)calloc(stride*8,sizeof(float));
  vectors.vec_unew=(float*)calloc(stride*8,sizeof(float));
  vectors.vec_vnew=(float*)calloc(stride*8,sizeof(float));
  vectors.vec_wnew=(float*)calloc(stride*8,sizeof(float));
  vectors.vec_pnew=(float*)calloc(stride*8,sizeof(float));

#endif

  vectors.vecnei=(int *)calloc(stride*6,sizeof(int));
  vectors.vecl=(int *)calloc(stride,sizeof(int));
  vectors.veccpu=(int *)calloc(stride,sizeof(int));
  vectors.vecicoarse=(int *)calloc(stride,sizeof(int));
  memsize+= stride*32*4;
  if(cpu.rank==0) printf(" vect = %f MB\n",((stride*32)/(1024*1024.))*sizeof(float));
#endif 

#ifdef WGPU

  // assigning a GPU to each CPU
  int ngpu;
  int memgpu;
  initlocaldevice(cpu.rank,2);

  // main Tree data on GPU
  cudaMalloc((void **)&(vectors.vecl_d),sizeof(int)*stride);
  cudaMalloc((void **)&(vectors.veccpu_d),sizeof(int)*stride);
  cudaMalloc((void **)&(vectors.vecnei_d),sizeof(int)*stride*6);
  cudaMalloc((void **)&(vectors.vecicoarse_d),sizeof(int)*stride);
  cudaMalloc((void **)&(vectors.vecden_d),sizeof(float)*stride*8);
  cudaMalloc((void **)&(vectors.vecpot_d),sizeof(float)*stride*8);


  // temp GPU arrays
  cudaMalloc((void **)&(vectors.vecpotnew_d),sizeof(float)*stride*8);
  cudaMalloc((void **)&(vectors.vec2_d),sizeof(float)*stride*8);
  cudaMalloc((void **)&(vectors.vecsum_d),sizeof(float)*stride*8);

  memgpu=49*stride*4/(1024*1024);
  if(cpu.rank==0) printf("MEM GPU = %d MB allcoated\n",memgpu);


#endif
    
  if(cpu.rank==0) printf("Allocations %f GB done\n",memsize/(1024.*1024*1024));

  //========== setting up the parallel topology ===

  cpu.nsend=NULL;
  cpu.nrecv=NULL;

  // We segment the oct distributions at levelcoarse 
    cpu.bndoct=NULL;
    cpu.mpinei=NULL;
    cpu.dict=NULL;

    cpu.nbuff=param.nbuff;
    cpu.allkmin=(int*)calloc(cpu.nproc,sizeof(int));
    cpu.allkmax=(int*)calloc(cpu.nproc,sizeof(int));

    load_balance(levelcoarse,&cpu);

#ifdef WMPI
    MPI_Allgather(&cpu.kmin,1,MPI_INT,cpu.allkmin,1,MPI_INT,cpu.comm);
    MPI_Allgather(&cpu.kmax,1,MPI_INT,cpu.allkmax,1,MPI_INT,cpu.comm);
    MPI_Barrier(cpu.comm);
#else
    cpu.allkmin[0]=cpu.kmin;
    cpu.allkmax[0]=cpu.kmax;
#endif    
    

    
  //========== building the initial meshes ===

  if(cpu.rank==0) printf("building initial mesh\n");

  //breakmpi();
  // ZERO WE CREATE A ROOT CELL
  
  struct CELL root;
  root.child=grid;
  

  // FIRST WE POPULATE THE ROOT OCT
  grid->x=0.;
  grid->y=0.;
  grid->z=0.;

  grid->parent=NULL;
  grid->level=1;
  for(i=0;i<6;i++) grid->nei[i]=&root; //periodic boundary conditions
  grid->prev=NULL;
  grid->next=NULL;

  // setting the densities in the cells and the index
  for(icell=0;icell<8;icell++){ 
    grid->cell[icell].density=0.;
    grid->cell[icell].pot=0.;
    grid->cell[icell].temp=0.;
    grid->cell[icell].idx=icell;

#ifdef WHYDRO
    grid->cell[icell].d=0.;
    grid->cell[icell].u=0.;
    grid->cell[icell].v=0.;
    grid->cell[icell].w=0.;
    grid->cell[icell].p=0.;
    grid->cell[icell].a=0.;
#endif

  }

  grid->cpu=-1;
  grid->vecpos=-1;
  grid->border=0;

  // start the creation of the initial amr grid from level 1
  firstoct[0]=grid;
  lastoct[0]=grid;
  int noct2;
  int segok;

  newoct=grid+1;
  for(level=1;level<levelcoarse;level++){ // sweeping the levels from l=1 to l=levelcoarse
    dxcur=1./pow(2,level);
    nextoct=firstoct[level-1];
    noct2=0;
    if(nextoct==NULL) continue;
    do // sweeping level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++){ // sweeping the cells

	  segok=segment_cell(curoct,icell,&cpu,levelcoarse);// the current cell will be splitted according to a segmentation condition
	  if(segok==1){ 
	    //if(level==levelcoarse-1) printf(" segok=%d\n",segok);

	    noct2++;
	    // the newoct is connected to its mother cell
	    curoct->cell[icell].child=newoct;
	    
	    // a newoct is created
	    newoct->parent=&(curoct->cell[icell]);
	    newoct->level=curoct->level+1;
	    newoct->x=curoct->x+( icell   %2)*dxcur;
	    newoct->y=curoct->y+((icell/2)%2)*dxcur;
	    newoct->z=curoct->z+( icell   /4)*dxcur;

	    // filling the cells
	    for(ii=0;ii<8;ii++){
	      newoct->cell[ii].marked=0;
	      newoct->cell[ii].child=NULL;
	      newoct->cell[ii].density=0.;
	      newoct->cell[ii].idx=ii;
	      newoct->cell[ii].phead=NULL;
#ifdef WHYDRO    
	      newoct->cell[icell].d=0.;
	      newoct->cell[icell].u=0.;
	      newoct->cell[icell].v=0.;
	      newoct->cell[icell].w=0.;
	      newoct->cell[icell].p=0.;
	      newoct->cell[icell].a=0.;
#endif
	    }
	    
	    //the neighbours
	    getcellnei(icell, vnei, vcell);
	    for(ii=0;ii<6;ii++){
	      if((vnei[ii]!=6)){ 
		newoct->nei[ii]=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
	      }else{
		newoct->nei[ii]=&(curoct->cell[vcell[ii]]);
	      }
	    }

	    // vector data
	    newoct->vecpos=-1;
	    newoct->border=0;
	    
	    // preparing the next creations on level+1
	    newoct->next=NULL;
	    
	    if(firstoct[level]==NULL){
	      firstoct[level]=newoct;
	      newoct->prev=NULL;
	    }
	    else{
	      newoct->prev=lastoct[level];
	      lastoct[level]->next=newoct;
	    }
	    lastoct[level]=newoct;


	    // next oct ready
	    newoct++; 
	  }
 	}
      }while(nextoct!=NULL);
    if(cpu.rank==0) printf("level=%d noct=%d\n",level,noct2);
  }


 // ==================================== assigning CPU number to levelcoarse OCTS // filling the hash table // Setting up the MPI COMMS

  int newloadb=1;
  setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax,newloadb); // out of WMPI to compute the hash table
  newloadb=0;

#ifdef WMPI
  // allocating the communication buffers
  sendbuffer=(struct PACKET **)(calloc(cpu.nnei,sizeof(struct PACKET*)));
  recvbuffer=(struct PACKET **)(calloc(cpu.nnei,sizeof(struct PACKET*)));
  for(i=0;i<cpu.nnei;i++) {
    sendbuffer[i]=(struct PACKET *) (calloc(cpu.nbuff,sizeof(struct PACKET)));
    recvbuffer[i]=(struct PACKET *) (calloc(cpu.nbuff,sizeof(struct PACKET)));
  }

  psendbuffer=(struct PART_MPI **)(calloc(cpu.nnei,sizeof(struct PART_MPI*)));
  precvbuffer=(struct PART_MPI **)(calloc(cpu.nnei,sizeof(struct PART_MPI*)));
  for(i=0;i<cpu.nnei;i++) {
    psendbuffer[i]=(struct PART_MPI *) (calloc(cpu.nbuff,sizeof(struct PART_MPI)));
    precvbuffer[i]=(struct PART_MPI *) (calloc(cpu.nbuff,sizeof(struct PART_MPI)));
  }
#endif


  //===================================================================================================
  
  // ==== some initial dump

  sprintf(filename,"data/levstart.%05d.p%05d",0,cpu.rank);
  dumpcube(lmap,firstoct,0,filename,0.);
  /* sprintf(filename,"data/cpustart.%05d.p%05d",0,cpu.rank); */
  /* dumpcube(lmap,firstoct,3,filename,tsim); */

  // =====================  computing the memory location of the first freeoct and linking the freeocts

  freeoct=lastoct[levelcoarse-1]+1; //(at this stage the memory is perfectly aligned)
  freeoct->prev=NULL;
  freeoct->next=freeoct+1;
  for(curoct=freeoct+1;curoct<(grid+ngridmax);curoct++){ // sweeping all the freeocts
    curoct->prev=curoct-1;
    curoct->next=NULL;
    if(curoct!=(grid+ngridmax-1)) curoct->next=curoct+1;
  }
  //curoct->next=NULL; // for the last free oct
  //printf("freeoct=%p\n",freeoct);


#ifdef PIC
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================

  // ==================================== assigning particles to cells
  //breakmpi();
  if(cpu.rank==0) printf("==> starting part\n");
  firstoct_currl=0;
  for(il=1;il<levelcoarse;il++) firstoct_currl+=pow(pow(2,il-1),3); // the index of the first oct of current level
 
  // initialisation of particles
  

#ifdef PART2

  int ir,nr=2;
  ip=0;
  float dxcell=1./pow(2.,levelcoarse);
  float epsilon=0.;
  for(ir=0;ir<nr;ir++) {
    // first we read the position etc... (eventually from the file)
    if(ir==0){
      x=0.5;
      y=0.5;
      z=0.5;

      vx=0.;
      vy=0.;
      vz=0.;
      
      mass=1.0-epsilon;
    }
    else if(ir==1){

      x=0.5+0.1;
      y=0.5;
      z=0.5;

      vx=0.;
      vy=sqrt((1.-epsilon)/0.1)*.5;
      vz=0.;
      
      mass=epsilon;
    }
    
    // periodic boundary conditions
    
    x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
    y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
    z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 
    
    // it it belongs to the current cpu, we proceed and assign the particle to the particle array
    if(segment_part(x,y,z,&cpu,levelcoarse)){
      part[ip].x=x;
      part[ip].y=y;
      part[ip].z=z;
      
      part[ip].vx=vx;
      part[ip].vy=vy;
      part[ip].vz=vz;
      
      part[ip].mass=mass;
      lastpart=part+ip;
      part[ip].idx=ir;
      ip++;
    }
  }
  
  npart=ip; // we compute the localnumber of particle

#endif


#ifdef PARTN

  int ir,nr=32768;
  ip=0;
  float dxcell=1./pow(2.,levelcoarse);
  float th,ph,r;
  for(ir=0;ir<nr;ir++) {
    // first we read the position etc... (eventually from the file)
    if(ir==0){
      x=0.5;
      y=0.5;
      z=0.5;

      vx=0.;
      vy=0.;
      vz=0.;
      
      mass=1.;
    }
    else{


      th=acos(((float)(rand())/RAND_MAX*2-1.));
      ph=2*M_PI*(float)(rand())/RAND_MAX;
      r=(float)(rand())/RAND_MAX*0.3;

      x=r*sin(th)*cos(ph)+0.5;
      y=r*sin(th)*sin(ph)+0.5;
      z=r*cos(th)+0.5;

      vx=(float)(rand())/RAND_MAX*2.-1.;
      vy=(float)(rand())/RAND_MAX*2.-1.;
      vz=(float)(rand())/RAND_MAX*2.-1.;
      
      mass=0.;
    }
    
    // periodic boundary conditions
    
    x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
    y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
    z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 
    
    // it it belongs to the current cpu, we proceed and assign the particle to the particle array
    if(segment_part(x,y,z,&cpu,levelcoarse)){
      part[ip].x=x;
      part[ip].y=y;
      part[ip].z=z;
      
      part[ip].vx=vx;
      part[ip].vy=vy;
      part[ip].vz=vz;
      
      part[ip].mass=mass;
      lastpart=part+ip;
      part[ip].idx=ir;
      ip++;
    }
  }
  
  npart=ip; // we compute the localnumber of particle

#endif

#ifdef TESTPLUM
  int dummy;
  float dummyf;
  int npartf;

  //breakmpi();
  fd=fopen("utils/data.inp","r");
  if(fd==NULL) {
    printf("Error while reading particle file ABORT\n");
    abort();
  }
  fscanf(fd,"%d",&dummy);
  fscanf(fd,"%d",&npartf);
  fscanf(fd,"%f",&dummyf);

  ip=0.;
  for(i=0;i<npartf;i++)
    {
      //fscanf(fd,"%d %f %f %f %f %f %f %f",&part[i].idx,&part[i].mass,&(part[i].x),&(part[i].y),&(part[i].z),&(part[i].vx),&(part[i].vy),&(part[i].vz));
      fscanf(fd,"%d %f %f %f %f %f %f %f",&dummy,&mass,&x,&y,&z,&vx,&vy,&vz);
      
      x+=0.5;
      y+=0.5;
      z+=0.5;
      // periodic boundary conditions
    
      x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
      y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
      z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 

      // it it belongs to the current cpu, we proceed and assign the particle to the particle array
      if(segment_part(x,y,z,&cpu,levelcoarse)){
	part[ip].x=x;
	part[ip].y=y;
	part[ip].z=z;
	
	part[ip].vx=vx;
	part[ip].vy=vy;
	part[ip].vz=vz;
	
	part[ip].mass=mass;
	lastpart=part+ip;
	ip++;
      }
      
    }
  fclose(fd);
  npart=ip; // we compute the localnumber of particle

#endif  

#ifdef TESTCOSMO

#ifdef HPC
  long dummy;
#else
  int dummy;
#endif
  float dummyf;
  int npartf;

  int nploc;
  float munit;
  float ainit;
  float lbox;

#ifdef GRAFIC // ==================== read grafic file

  lastpart=read_grafic_part(part, &cpu, &munit, &ainit, &omegam, &omegav, &Hubble, &npart,omegab);

#ifdef SUPERCOMOV
  // at this stage we have to compute the conformal time
  tsim=-0.5*sqrt(omegam)*integ_da_dt_tilde(ainit,1.0,omegam,omegav,1e-8);
  
  // we compute the friedmann tables
  
  double tab_aexp[NCOSMOTAB];
  double tab_ttilde[NCOSMOTAB];
  double tab_t[NCOSMOTAB];
  
  compute_friedmann(ainit*0.95,NCOSMOTAB,omegam,omegav,tab_aexp,tab_ttilde,tab_t);
  
#endif

#else // ============================  explicit read here
#ifndef ZELDO
  fd=fopen("utils/IC.PM.0","rb");
#else
  fd=fopen("utils/ZEL.PM.0","rb");
#endif
  rstat=fread(&dummy,sizeof(dummy),1,fd);
  rstat=fread(&nploc,sizeof(int),1,fd);
  rstat=fread(&munit,sizeof(float),1,fd);
  rstat=fread(&ainit,sizeof(float),1,fd);
  rstat=fread(&lbox,sizeof(float),1,fd);
  rstat=fread(&omegam,sizeof(float),1,fd);
  rstat=fread(&omegav,sizeof(float),1,fd);
  rstat=fread(&Hubble,sizeof(float),1,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);

  if(cpu.rank==0) printf("%f %d %f %f\n",ainit,nploc,omegav,Hubble);
  
  mass=munit;
  tsim=ainit;
  aexp=ainit;

#ifdef SUPERCOMOV
  // at this stage we have to compute the conformal time
  tsim=-0.5*sqrt(omegam)*integ_da_dt_tilde(ainit,1.0,omegam,omegav,1e-8);

  // we compute the friedmann tables

  double tab_aexp[NCOSMOTAB];
  double tab_ttilde[NCOSMOTAB];
  double tab_t[NCOSMOTAB];

  compute_friedmann(ainit*0.95,NCOSMOTAB,omegam,omegav,tab_aexp,tab_ttilde,tab_t);

#endif

  // ------------------ dealing with particles

  float *pos;
  float *vel;
  int nread=(nploc<2097152?nploc:2097152); // we read by patch of 128^3
  int npatch=nploc/nread;
  int ipatch;

  pos=(float *)malloc(sizeof(float)*3*nread);
  vel=(float *)malloc(sizeof(float)*3*nread);
  int pstart=ftell(fd);

  ip=0.;
  for(ipatch=0;ipatch<npatch;ipatch++) {
    //    rstat=fread(&dummy,sizeof(dummy),1,fd); 
    //    fseek(fd,pstart,SEEK_SET);
    fseek(fd,pstart+(0*nploc+ipatch*nread)*sizeof(float)+1*sizeof(dummy),SEEK_SET);
    rstat=fread(pos,sizeof(float),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);

    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(1*nploc+ipatch*nread)*sizeof(float)+3*sizeof(dummy),SEEK_SET);
    rstat=fread(pos+nread,sizeof(float),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);

    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(2*nploc+ipatch*nread)*sizeof(float)+5*sizeof(dummy),SEEK_SET);
    rstat=fread(pos+2*nread,sizeof(float),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
  
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(3*nploc+ipatch*nread)*sizeof(float)+7*sizeof(dummy),SEEK_SET);
    rstat=fread(vel,sizeof(float),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);

    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(4*nploc+ipatch*nread)*sizeof(float)+9*sizeof(dummy),SEEK_SET);
    rstat=fread(vel+nread,sizeof(float),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);

    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(5*nploc+ipatch*nread)*sizeof(float)+11*sizeof(dummy),SEEK_SET);
    rstat=fread(vel+2*nread,sizeof(float),nread,fd);
    //rstat=fread(&dummy,sizeof(dummy),1,fd);
    
    for(i=0;i<nread;i++)
      {
	x=pos[i];
	y=pos[i+nread];
	z=pos[i+2*nread];
	vx=vel[i];
	vy=vel[i+nread];
	vz=vel[i+2*nread];
	// periodic boundary conditions
	
	x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
	y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
	z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 
	if(ip<4){
	  printf("x=%f y=%f z=%f\n",x,y,z);
	  printf("vx=%f\n",vx);
	}
	// it it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,&cpu,levelcoarse)){
	  if(ip<4) printf("passed\n");
	  part[ip].x=x;
	  part[ip].y=y;
	  part[ip].z=z;
	  
	  part[ip].vx=vx;
	  part[ip].vy=vy;
	  part[ip].vz=vz;
	  
	  part[ip].mass=mass;
	  part[ip].idx=i;
	  lastpart=part+ip;
	  ip++;
	}
      }
  }

  npart=ip;
  fclose(fd);
  free(pos);
  free(vel);
  if (cpu.rank==0) printf("cosmo readpart done with munit=%e\n",munit);
#endif
#endif

  // we set all the "remaining" particles mass to -1
  for(ii=npart;ii<npartmax;ii++) part[ii].mass=-1.0;


  // assigning particles to cells in coarse octs (assuming octs are aligned)

  if(cpu.rank==0) printf("start populating coarse grid with particles\n");
  struct PART* lastp[8]; // will contain the last particle of the 8 cells in each oct

  // FIRST WE CONSIDER THE LEVEL 1
  for(ii=0;ii<8;ii++) lastp[ii]=NULL; // we initialise the last part of each sub cell
  dxcur=0.5;
  for(i=0;i<npart;i++)
    {
      curc=(int)((part[i].x-grid[0].x)/dxcur)+(int)((part[i].y-grid[0].y)/dxcur)*2+(int)((part[i].z-grid[0].z)/dxcur)*4;
      
      if(grid[0].cell[curc].phead==NULL){
	grid[0].cell[curc].phead=&part[i];
	lastp[curc]=&part[i];
      }
      else{
	lastp[curc]->next=&part[i];
	part[i].prev=lastp[curc];
	lastp[curc]=&part[i];
      }
    }
  if(cpu.rank==0) printf("Part assigned root level ok\n");
  
  // WE POPULATE THE NEXT LEVELS BY SUBDIVISIONS
  int np=0;
  for(level=1;level<=levelcoarse-1;level++) // we stop at level coarse -1 because it will be assigned from levelcoarse-1
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  dxcur=1./pow(2,level+1); // size of a CELL at level +1
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if(curoct->cell[icell].child!=NULL){ // a child has been detected so we split the particle in 8 cells
		for(ii=0;ii<8;ii++) lastp[ii]=NULL; // we initialise the last part of each sub cell
		newoct=curoct->cell[icell].child; // temp oct for practical reasons
		nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
		if(nexp!=NULL){ 
		  do{  
		    curp=nexp; 
		    nexp=curp->next; 
		    
		    //curc is the index of the new cell at level+1
		    curc=(int)((curp->x-newoct->x)/dxcur)+(int)((curp->y-newoct->y)/dxcur)*2+(int)((curp->z-newoct->z)/dxcur)*4;
		    if(lastp[curc]==NULL){
		      // first particle in the current subcell
		      newoct->cell[curc].phead=curp;
		      curp->prev=NULL;
		      curp->next=NULL;
		      lastp[curc]=curp;
		    }
		    else{
		      // the current particle is linked to the last one in the current part
		      curp->prev=lastp[curc];
		      lastp[curc]->next=curp;
		      curp->next=NULL;
		      lastp[curc]=curp;
		    }
		  }while(nexp!=NULL); 
		  
		  // we empty the mother cell from particles
		  curoct->cell[icell].phead=NULL;
		  
		}
	      }
	    }
	}while(nextoct!=NULL);
    }


#if 1
  // ==================================== Check the number of particles and octs

  mtot=multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);

#endif	

  // ==================================== performing the CIC assignement
#ifndef WGPU
  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#else
  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#endif

#ifdef WMPI
    // ==================================== performing the CIC BOUNDARY CORRECTION 

  mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);

  // ======================================= Density boundary mpi update 

  mpi_exchange(&cpu,sendbuffer,recvbuffer,1,1);

#endif

  mtot=multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);

  /* sprintf(filename,"data/partstart.%05d.p%05d",0,cpu.rank); */
  /* dumppart(firstoct,filename,npart,levelcoarse,levelmax,tsim); */
  
  /* sprintf(filename,"data/denstart.%05d.p%05d",0,cpu.rank); */
  /* dumpcube(lmap,firstoct,1,filename,tsim); */
  /* abort(); */
#endif

  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================



#ifdef WHYDRO
  
#ifdef GRAFIC
  int ncellhydro;
  ncellhydro=read_grafic_hydro(&cpu,omegab);
  printf("%d hydro cell found in grafic file\n",ncellhydro);
#else
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================

  // initialisation of hydro quantities
  // Shock Tube

  struct Wtype WL, WR;
  float X0;
  if(cpu.rank==0) printf("Init Hydro\n");

   /* // TEST 1 */

  WL.d=1.;
  WL.u=0.;
  WL.v=0.;
  WL.w=0.75;
  WL.p =1.;

  WR.d=0.125;
  WR.u=0.;
  WR.v=0.;
  WR.w=0.;
  WR.p=0.1;
  X0=0.3;

   /* // SPHERICAL EXPLOSION */

  WL.d=10.;
  WL.u=0.;
  WL.v=0.;
  WL.w=0.;
  WL.p =1.;

  WR.d=0.125;
  WR.u=0.;
  WR.v=0.;
  WR.w=0.;
  WR.p=1.;
  X0=0.2;

  WL.a=sqrtf(GAMMA*WL.p/WL.d);
  WR.a=sqrtf(GAMMA*WR.p/WR.d);

  // ======================================================

  float dtot=0.;
  int nc=0;
  for(level=levelcoarse;level<=levelmax;level++) // (levelcoarse only for the moment)
    {
      dxcur=pow(0.5,level);
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      xc=curoct->x+( icell&1)*dxcur;
	      yc=curoct->y+((icell>>1)&1)*dxcur;
	      zc=curoct->z+((icell>>2))*dxcur;

	      float vrx=(((float)rand())/RAND_MAX)*0.02-0.01;
	      float vry=(((float)rand())/RAND_MAX)*0.02-0.01;
	      float vrz=(((float)rand())/RAND_MAX)*0.02-0.01;

	      /* if((zc>0.75)||(zc<0.25)){ */
	      /* 	curoct->cell[icell].d=1.0; */
	      /* 	curoct->cell[icell].u=-0.5+vrx; */
	      /* 	curoct->cell[icell].v=vry; */
	      /* 	curoct->cell[icell].w=vrz; */
	      /* 	curoct->cell[icell].p=2.5; */
	      /* 	curoct->cell[icell].a=sqrtf(GAMMA*2.5/1.); */
	      /* } */
	      /* else{ */
	      /* 	curoct->cell[icell].d=2.0; */
	      /* 	curoct->cell[icell].u=0.5+vrx; */
	      /* 	curoct->cell[icell].v=vry; */
	      /* 	curoct->cell[icell].w=vrz; */
	      /* 	curoct->cell[icell].p=2.5; */
	      /* 	curoct->cell[icell].a=sqrtf(GAMMA*2.5/1.); */
	      /* } */

	      
	      if((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5)<(X0*X0)){
	      	curoct->cell[icell].d=WL.d;
	      	curoct->cell[icell].u=WL.u;
	      	curoct->cell[icell].v=WL.v;
	      	curoct->cell[icell].w=WL.w;
	      	curoct->cell[icell].p=WL.p;
	      	curoct->cell[icell].a=WL.a;

	      }
	      else{
	      	curoct->cell[icell].d=WR.d;
	      	curoct->cell[icell].u=WR.u;
	      	curoct->cell[icell].v=WR.v;
	      	curoct->cell[icell].w=WR.w;
	      	curoct->cell[icell].p=WR.p;
	      	curoct->cell[icell].a=WR.a;
	      }
	      
	      /* if(zc<=X0){ */
	      /* 	curoct->cell[icell].d=WL.d; */
	      /* 	curoct->cell[icell].u=WL.u; */
	      /* 	curoct->cell[icell].v=WL.v; */
	      /* 	curoct->cell[icell].w=WL.w; */
	      /* 	curoct->cell[icell].p=WL.p; */
	      /* 	curoct->cell[icell].a=WL.a; */
	      /* } */
	      /* else{ */
	      /* 	curoct->cell[icell].d=WR.d; */
	      /* 	curoct->cell[icell].u=WR.u; */
	      /* 	curoct->cell[icell].v=WR.v; */
	      /* 	curoct->cell[icell].w=WR.w; */
	      /* 	curoct->cell[icell].p=WR.p; */
	      /* 	curoct->cell[icell].a=WR.a; */
	      /* } */
	      
	      if(level==levelcoarse) {
		dtot+=curoct->cell[icell].d;
		nc++;
	      }

	    }
	}while(nextoct!=NULL);
      
      //printf("level=%d avg=%e mind=%e maxd=%e\n",level,avg/ncell,mind,maxd);
    }
  
  avgdens+=dtot/nc;
  printf("avgdens=%e\n",avgdens);

#endif
  sprintf(filename,"data/denhydstart.%05d.p%05d",0,cpu.rank); 
  dumpcube(lmap,firstoct,101,filename,0.); 
  sprintf(filename,"data/prehydstart.%05d.p%05d",0,cpu.rank); 
  dumpcube(lmap,firstoct,105,filename,0.); 


  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================
  //===================================================================================================================================

#endif




  // =============================================== dumping information file


  FILE *fi;
  if(cpu.rank==0){
    fi=fopen("data/inforun.txt","w");
    fprintf(fi,"levelcoarse=%d\n",levelcoarse);
    fprintf(fi,"levelmax=%d\n",levelmax);
    fprintf(fi,"levelmap=%d\n",lmap);
    fprintf(fi,"nstepmax=%d\n",param.nsteps);
    fprintf(fi,"amr threshold=%e\n",threshold);
  }
  

  //================================================================================
  //================================================================================
  //================================================================================
  //
  //          AT THIS STAGE THE INITIAL SETUP HAS BEEN COMPLETED
  //
  //================================================================================
  //================================================================================
  //================================================================================


  int nsteps;
  int pass;
  int smark;
  int ismooth,nsmooth=2;
  int marker;

  float tmax;


#ifdef TESTCOSMO
#ifdef ZELDO
#ifdef SUPERCOMOV
  tmax=-0.5*sqrt(omegam)*integ_da_dt_tilde(0.5,1.0,omegam,omegav,1e-8);
#else
  tmax=0.5;
#endif
#else
#ifdef SUPERCOMOV
  tmax=0.; // in SC coordinates the maximal conformal time is 0.
#else
  tmax=1.;
#endif
#endif
#else
  tmax=1000.;
#endif

#ifdef WHYDRO
  tmax=0.21;
#endif

  FILE *fegy;

  //breakmpi();

#ifndef TESTCOSMO
    faexp=1.0;
    faexp2=1.0;
#endif


    //==================================== MAIN LOOP ================================================
    //===============================================================================================

    int ndumps=0;
    for(nsteps=0;(nsteps<=param.nsteps)*(tsim<=tmax);nsteps++){
    
#ifdef TESTCOSMO
#ifdef SUPERCOMOV
      aexp=interp_aexp(tsim,tab_aexp,tab_ttilde);
#else
      aexp=tsim;
#endif
      if(cpu.rank==0) printf("\n============== STEP %d aexp=%e z=%f ================\n",nsteps,aexp,1./aexp-1.);
#else
      if(cpu.rank==0) printf("\n============== STEP %d tsim=%e ================\n",nsteps,tsim);
#endif


#if 1
    // ==================================== marking the cells
    mark_cells(levelcoarse,levelmax,firstoct,nsmooth,threshold,&cpu,sendbuffer,recvbuffer);


    // ==================================== refining (and destroying) the octs

    curoct=refine_cells(levelcoarse,levelmax,firstoct,lastoct,freeoct,&cpu,grid+ngridmax);
    freeoct=curoct;
  
#endif

#ifdef WMPI
  // ==================================== after refinement we should remap the boundary cells
  setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax,newloadb);
#endif

  // ==================================== performing the CIC assignement
#ifndef WGPU
  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#else
  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#endif

#ifdef WMPI
  //------------- performing the CIC BOUNDARY CORRECTION 
  mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);
  //------------  Density boundary mpi update 
  mpi_exchange(&cpu,sendbuffer,recvbuffer,1,1);
#endif


  //======================================= cleaning the marks
  for(level=1;level<=levelmax;level++) // looping over levels
    {
      /* float maxd=0.,mind=1e30,avg=0.; */
      /* int ncell=0; */
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      curoct->cell[icell].marked=0.;
	    }
	}while(nextoct!=NULL);

      //printf("level=%d avg=%e mind=%e maxd=%e\n",level,avg/ncell,mind,maxd);
    }



  // ==================================== Check the number of particles and octs
    
    mtot=multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);
#ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&mtot,1,MPI_FLOAT,MPI_SUM,cpu.comm);
#endif

  int ltot,gtot=0,nomax,nomin;
  if(cpu.rank==0){
    printf("===================================================\n");
  }
  for(level=2;level<=levelmax;level++){
    ltot=cpu.noct[level-1];
    nomax=ltot;
    nomin=ltot;
    gtot+=ltot;
#ifdef WMPI
    MPI_Allreduce(&ltot,&nomax,1,MPI_INT,MPI_MAX,cpu.comm);
    MPI_Allreduce(&ltot,&nomin,1,MPI_INT,MPI_MIN,cpu.comm);
    MPI_Allreduce(MPI_IN_PLACE,&ltot,1,MPI_INT,MPI_SUM,cpu.comm);
#endif
    if(cpu.rank==0){
      printf("level=%2d noct=%9d min=%9d max=%9d\n",level,ltot,nomin,nomax);
    }
  }
#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&gtot,1,MPI_INT,MPI_MAX,cpu.comm);
#endif
  if(cpu.rank==0){
    printf("grid occupancy=%4.1f \n",(gtot/(1.0*ngridmax))*100.);
  }



#ifdef SELFGRAV
    // ==================================== POISSON Testing the jacobi iteration
      float res;
#ifdef TIME_JAC
  FILE *ft;
  //  sprintf(filename,"data/timejac.%05d.p%05d",nsteps+1,cpu.rank);
  //ft=fopen(filename,"w");
  
  double tg1,tg2;
  double tl1,tl2;
  double ta1,ta2;
  double ts1,ts2;
  double tm1,tm2;
  double t1,t2;
  double tc1;
  double tg,tl,ts,ta;
#endif

  if(cpu.rank==0){
    printf("=======================================\n");
    printf("==> Poisson Start with ncell coarse =%d \n",cpu.noct[levelcoarse-1]*8);
    printf("=======================================\n");
  }
  int icomp,iter,niter=param.niter;
  float norm_d;
  int nread;
  char compnei[8];
  char compvar[8];
  int ic;

  for(level=levelcoarse;level<=levelmax;level++)
    {
      norm_d=0.;
	
      /* ================================ 	// FINE LEVEL ONLY initial guess from parent cell */

      if(level>levelcoarse){
	nextoct=firstoct[level-1];
	if(nextoct==NULL){
	  //printf("Proc %d skipping level=%d\n",cpu.rank,level);
	}
	else{
	  do{
	    curoct=nextoct;
	    nextoct=curoct->next;
	    if(curoct->cpu!=cpu.rank) continue;
	    for(icell=0;icell<8;icell++){
	      curoct->cell[icell].pot=curoct->parent->pot;
	    }
	  }while(nextoct!=NULL);
	}

#ifdef WMPI
	mpi_exchange(&cpu,sendbuffer,recvbuffer,2,1);
#endif
      }


      // ====================== Here starts the new jacobian procedure ====================
      //===================================================================================
      //===================================================================================

#ifndef MULTIGRID
      // --- pure jacobi relaxation
      poisson_jacob(level,levelcoarse,levelmax,firstoct,&vectors,stride,&cpu,omegam,aexp,sendbuffer,recvbuffer,niter,acc,avgdens);
#else
      int igrid;
      if((level==levelcoarse)&&(levelcoarse!=levelmin)){
	// --- for coarse level two levels mgrid relaxation
	for(igrid=0;igrid<nvcycles;igrid++){
	  if(cpu.rank==0) printf("----------------------------------------\n");
	  clean_pot(levelcoarse-1,firstoct);	    // ------------- cleaning the vector positions
	  res=poisson_mgrid(level,levelcoarse,levelmax,levelmin,firstoct,&vectors,stride,&cpu,omegam,aexp,sendbuffer,recvbuffer,niter,acc,avgdens,nrelax);
	  if(res<acc) break;
	}
      }
      else{
	// --- pure jacobi relaxation
	poisson_jacob(level,levelcoarse,levelmax,firstoct,&vectors,stride,&cpu,omegam,aexp,sendbuffer,recvbuffer,niter,acc,avgdens);
      }

       if(cpu.rank==0) printf("----------------------------------------\n");
#endif

    }      
  
#ifdef TIME_JAC
  //fclose(ft);
#endif

#endif

  // ==================================== Force calculation and velocity update   // Corrector step
  if(cpu.rank==0) printf("Predictor\n");
#ifdef TESTCOSMO
#ifdef SUPERCOMOV
  faexp=1.0;
#else
  faexp=f_aexp(aexp,omegam,omegav);
#endif
#endif
  
#ifdef SUPERCOMOV
  forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dt*0.5*(nsteps!=0),&cpu,sendbuffer,recvbuffer);
#else
  forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dt*0.5*faexp*(nsteps!=0),&cpu,sendbuffer,recvbuffer);
#endif

#ifdef SELFGRAV
  // ==================================== hydro update   // Corrector step
  if(cpu.rank==0) printf("Predictor hydro\n");
  for(level=levelcoarse;level<=levelmax;level++)
    {
      nextoct=firstoct[level-1];
      correct_grav_hydro(nextoct, &cpu, dt*(nsteps!=0));
    }
#endif

  // ==================================== DUMP AFTER SYNCHRONIZATION

  if(nsteps%(param.ndumps)==0){
    if(cpu.rank==0) printf("Dumping .......\n");
    float tdump;

#ifdef TESTCOSMO
#ifdef SUPERCOMOV
    tdump=interp_aexp(tsim+0*(dt+dtnew)*0.5*(nsteps!=0),tab_aexp,tab_ttilde);
#else
    tdump=tsim+(dt+dtnew)*0.5*(nsteps!=0);
#endif
#else
    tdump=tsim+(dt+dtnew)*0.5*(nsteps!=0);
#endif
    // ===== Casting rays to fill a map
    sprintf(filename,"data/pot3d.%05d.p%05d",ndumps,cpu.rank);
    dumpcube(lmap,firstoct,2,filename,tdump);
    /* sprintf(filename,"data/lev3d.%05d.p%05d",ndumps,cpu.rank); */
    /* dumpcube(lmap,firstoct,0,filename,tdump); */
    sprintf(filename,"data/den3d.%05d.p%05d",ndumps,cpu.rank);
    dumpcube(lmap,firstoct,1,filename,tdump);
    sprintf(filename,"data/fx.%05d.p%05d",ndumps,cpu.rank);
    dumpcube(lmap,firstoct,6,filename,tdump);
    /* sprintf(filename,"data/fy.%05d.p%05d",ndumps,cpu.rank); */
    /* dumpcube(lmap,firstoct,7,filename,tdump); */
    /* sprintf(filename,"data/fz.%05d.p%05d",ndumps,cpu.rank); */
    /* dumpcube(lmap,firstoct,8,filename,tdump); */


#ifdef WMPI
    if(nsteps==0){
      sprintf(filename,"data/cpu3d.%05d.p%05d",ndumps,cpu.rank);
      dumpcube(lmap,firstoct,3,filename,tdump);
    }
#endif
  
  //==== Gathering particles for dump

#ifdef PIC
    sprintf(filename,"data/part.%05d.p%05d",ndumps,cpu.rank);
    dumppart(firstoct,filename,npart,levelcoarse,levelmax,tdump);
#endif

    // === Hydro dump
#ifdef WHYDRO
    sprintf(filename,"data/denhyd.%05d.p%05d",ndumps,cpu.rank); 
    dumpcube(lmap,firstoct,101,filename,tdump); 
    sprintf(filename,"data/velhyd.%05d.p%05d",ndumps,cpu.rank); 
    dumpcube(lmap,firstoct,102,filename,tdump); 
    sprintf(filename,"data/prehyd.%05d.p%05d",ndumps,cpu.rank); 
    dumpcube(lmap,firstoct,105,filename,tdump); 
#endif

    ndumps++;
  }


  // ==================================== New timestep
  dtnew=dt;
#ifdef TESTCOSMO
#ifdef SUPERCOMOV
  faexp2=1.0;
#else
  faexp2=f_aexp(aexp,omegam,omegav)/(aexp*aexp);
#endif
#endif

#ifdef PIC
  float dtpic;
  dtpic=comptstep(levelcoarse,levelmax,firstoct,faexp2,faexp,&cpu,param.dt);
  dtnew=(dtpic<dtnew?dtpic:dtnew);
#endif

#ifdef WHYDRO
  float dthydro;
  dthydro=comptstep_hydro(levelcoarse,levelmax,firstoct,faexp2,faexp,&cpu,param.dt);
  dtnew=(dthydro<dtnew?dthydro:dtnew);
  //if(cpu.rank==0) printf("dt=%e dthydro=%e\n",dtnew,dthydro);
#endif

#ifdef SELFGRAV
  float dtff;
  dtff=comptstep_ff(levelcoarse,levelmax,firstoct,aexp,&cpu,param.dt);
  dtnew=(dtff<dtnew?dtff:dtnew);
  if(cpu.rank==0) printf("dt=%e dthydro=%e dtpic=%e dtff=%e\n",dtnew,dthydro,dtpic,dtff);
#endif

  // ==================================== Force calculation and velocity update   // predictor step
#ifdef PIC
  if(cpu.rank==0) printf("Corrector\n");

#ifdef TESTCOSMO
#ifdef SUPERCOMOV
  faexp=1.0;
#else
  faexp=f_aexp(aexp,omegam,omegav);
#endif
#endif

#ifdef SUPERCOMOV
  forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dtnew*0.5,&cpu,sendbuffer,recvbuffer);
#else
  forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dtnew*0.5*faexp,&cpu,sendbuffer,recvbuffer);
#endif
#endif

  // ==================================== Moving Particles + Oct management
#ifdef PIC
  
#ifndef PARTN
  if(cpu.rank==0) printf("Moving particles\n");
    // Computing displacement (predictor)

#ifdef TESTCOSMO
#ifdef SUPERCOMOV
  faexp2=1.0;
#else
  faexp2=f_aexp(tsim+(dtnew)*0.5,omegam,omegav)/pow(tsim+0.5*(dtnew),2);
#endif
#endif
    
#ifdef SUPERCOMOV
  movepart(levelcoarse,levelmax,firstoct,dtnew,&cpu);
#else
  movepart(levelcoarse,levelmax,firstoct,dtnew*faexp2,&cpu);
#endif
    // Moving particles through cells (3 passes)
    
#ifdef WGPU
    partcellreorg(levelcoarse,levelmax,firstoct);
#else
    partcellreorg(levelcoarse,levelmax,firstoct);
#endif
    
#ifdef WMPI
    
    // Communication of particles
    int deltan;
    deltan=mpi_exchange_part(&cpu,psendbuffer,precvbuffer,&lastpart);
    
    // Recounting particles
    npart=npart+deltan;
#endif
#endif
#endif


#ifdef WHYDRO

    // ================================= performing hydro calculations
    
    double t0,t100,t20,t80;
    if(cpu.rank==0) printf("Start Hydro\n");
    for(level=levelcoarse;level<=levelmax;level++)
      {
	dxcur=pow(0.5,level);
	// ------------- cleaning the vector positions
	clean_vec(levelmax,firstoct);

	// --------------- setting the first oct of the level
	curoct=firstoct[level-1];
	

#ifdef WMPI
	int nolevel=(curoct==NULL);
	// --------------- do we need to calculate the level ?
	MPI_Allreduce(MPI_IN_PLACE,&nolevel,1,MPI_INT,MPI_SUM,cpu.comm);
	if(nolevel==cpu.nproc){
	  //we skip the calculation
	  epsilon=0.;
	  return epsilon;
	}
#endif
	
	if((curoct!=NULL)&&(cpu.noct[level-1]!=0)){
	  
	  // -------------  cleaning working arrays
	  
	  memset(vectors.vec_d,0,sizeof(float)*stride*8);
	  memset(vectors.vec_u,0,sizeof(float)*stride*8);
	  memset(vectors.vec_v,0,sizeof(float)*stride*8);
	  memset(vectors.vec_w,0,sizeof(float)*stride*8);
	  memset(vectors.vec_p,0,sizeof(float)*stride*8);
	  
	  memset(vectors.vec_dnew,0,sizeof(float)*stride*8);
	  memset(vectors.vec_unew,0,sizeof(float)*stride*8);
	  memset(vectors.vec_vnew,0,sizeof(float)*stride*8);
	  memset(vectors.vec_wnew,0,sizeof(float)*stride*8);
	  memset(vectors.vec_pnew,0,sizeof(float)*stride*8);
	  
	  memset(vectors.vecnei,0,sizeof(int)*stride*6);
	  memset(vectors.vecl,0,sizeof(int)*stride);
	  memset(vectors.veccpu,0,sizeof(int)*stride);
	  memset(vectors.vecicoarse,0,sizeof(int)*stride);
	  
      
	  // ------------------- setting neighbors to -1
	  for(i=0;i<stride*6;i++){vectors.vecnei[i]=-1;}
	  
	  // ---------------- gathering the neighbors
	  
	  int nocttotal,nread;
	  nocttotal=countvecocts(curoct,stride,&cpu,&nread); // we count the actual numbers of octs involved
	  stride=(nocttotal%64==0?nocttotal:(nocttotal/64+1)*64); // let us be wild ! we change stride for nread
	  //printf("new stride=%d\n",stride);

	  t0=MPI_Wtime();
	  nextoct=gathervecnei2(curoct,vectors.vecnei,stride,&cpu,&nread); // gathering the neighbours
	  
  
	  // ------------ gathering the PRIMITIVE values
	  /* nextoct=gathervec2(curoct,vectors.vec_d,10,vectors.vecl,vectors.vecicoarse,vectors.veccpu,stride,&cpu,&nread); // density  */
	  /* nextoct=gathervec2(curoct,vectors.vec_u,20,vectors.vecl,vectors.vecicoarse,vectors.veccpu,stride,&cpu,&nread); // u */
	  /* nextoct=gathervec2(curoct,vectors.vec_v,30,vectors.vecl,vectors.vecicoarse,vectors.veccpu,stride,&cpu,&nread); // v */
	  /* nextoct=gathervec2(curoct,vectors.vec_w,40,vectors.vecl,vectors.vecicoarse,vectors.veccpu,stride,&cpu,&nread); // w */
	  /* nextoct=gathervec2(curoct,vectors.vec_p,50,vectors.vecl,vectors.vecicoarse,vectors.veccpu,stride,&cpu,&nread); // p */
	  
	  nextoct=gathervechydro(curoct,&vectors,stride,&cpu, &nread);

	  // ------------ solving the hydro
	  
	  t20=MPI_Wtime();
	  hydrosolve(&vectors,level,cpu.rank,nread,stride,dxcur,dtnew);
	  t80=MPI_Wtime();
	  
	  
	  // ------------ scatter back the result
	  
	  nextoct=scattervechydro(curoct,&vectors, stride, &cpu);

	  t100=MPI_Wtime();

	  /*   nextoct=scattervec(curoct,vectors.vec_d,10,stride,&cpu,nread);  */
	  /* nextoct=scattervec(curoct,vectors.vec_u,20,stride,&cpu,nread);  */
	  /* nextoct=scattervec(curoct,vectors.vec_v,30,stride,&cpu,nread);  */
	  /* nextoct=scattervec(curoct,vectors.vec_w,40,stride,&cpu,nread);  */
	  /* nextoct=scattervec(curoct,vectors.vec_p,50,stride,&cpu,nread);  */

	}
      }

  if(cpu.rank==0) printf("Hydro done in %e (%e in hydro) sec\n",t100-t0,t80-t20);

#endif


#if 1
  // ==================================== Check the number of particles and octs
  multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);

#endif	 

#ifdef EGYCSV
  // ==================================== Energy Conservation Test
  
  float egy;
  egy=egypart(levelcoarse,levelmax,firstoct,&cpu);
  
#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&egy,1,MPI_FLOAT,MPI_SUM,cpu.comm);
#endif

  if(cpu.rank==0){
    printf("== Energy check etot=%e\n",egy);
    if(nsteps==0){
      fegy=fopen("egy.txt","w");
    }
    fprintf(fegy,"%e %e\n",tsim,egy);
  }

#endif  

  // =================================== dumping information

#ifdef PIC
  if(cpu.rank==0) fprintf(fi,"istep=%d t=%e dt=%e coarseres=%e\n",nsteps,tsim,(dt+dtnew)*0.5,res);
#endif

  //==================================== timestep completed, looping
  tsim+=dtnew;
  dt=dtnew;
  }

  fclose(fi);

#ifdef NEWJACK
      free(vectors.vecpot); //contains the potential in "stride" octs
      free(vectors.vecpotnew); //contains the potential in "stride" octs
      free(vectors.vecden); //contains the density   in "stride" octs

#ifdef WHYDRO
      free(vectors.vec_d); 
      free(vectors.vec_u); 
      free(vectors.vec_v); 
      free(vectors.vec_w); 
      free(vectors.vec_p); 

      free(vectors.vec_dnew); 
      free(vectors.vec_unew); 
      free(vectors.vec_vnew); 
      free(vectors.vec_wnew); 
      free(vectors.vec_pnew); 
#endif

      free(vectors.vecnei);//contains the cell neighbors of the octs
      free(vectors.vecl); // contains the level of the octs
      free(vectors.veccpu); // contains the level of the octs
      free(vectors.vecicoarse); // contains the level of the octs

#ifdef WGPU
      cudaFree(vectors.vecl_d);
      cudaFree(vectors.veccpu_d);
      cudaFree(vectors.vecden_d);
      cudaFree(vectors.vec2_d);
      cudaFree(vectors.vecsum_d);
      cudaFree(vectors.vecpot_d);
      cudaFree(vectors.vecpotnew_d);
      cudaFree(vectors.vecnei_d);
      cudaFree(vectors.vecicoarse_d);


#endif

#endif



  if(cpu.rank==0){
    printf("Done .....\n");
  }

#ifdef EGYCSV
  fclose(fegy);
#endif

#ifdef WMPI
  MPI_Finalize();
#endif

  return 0;
}
      
