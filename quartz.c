
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
#include "friedmann.h"
#include <time.h>
#include <mpi.h>

#ifdef WGPU
#include "interface.h"
#include "vector_gpu.h"
#include "cic_gpu.h"
#endif


#ifdef WHYDRO2
#include "hydro_utils.h"
#endif


#ifdef WGRAV
#include "poisson_utils.h"
#endif

void gdb_debug()
{
  int i = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for attach\n", getpid(), hostname);
  fflush(stdout);
  while (0 == i)
    sleep(5);
}

// ===============================================================================



//------------------------------------------------------------------------
 // the MAIN CODE
 //------------------------------------------------------------------------

#ifdef TESTCOSMO
REAL f_aexp(REAL aexp, REAL omegam, REAL omegav)
{
  return 1./sqrt(omegam/aexp+omegav*aexp*aexp);
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
  REAL dx;
  int vnei[6],vcell[6]; // arrays to get neighbors
  int vnei2[6],vcell2[6]; // arrays to get neighbors
  int vnei3[6],vcell3[6]; // arrays to get neighbors
  int neip[7]; // contains the index of the six neighbors of the current parent cell +current parent
  int ci,cj,ck;
  int cinext,cjnext,cknext;
  REAL threshold;
  REAL tsim=0.;
  struct OCT oct;
  struct OCT* nextoct;
  struct OCT* curoct;
  struct OCT* desoct;
  struct CELL * parcell;
  struct CELL * newcell;
  struct CELL * newcell2;
  int tag;
  REAL dxcur;
  REAL *dens;
  int firstoct_currl;
  int nxoct;
  int lmap;
  int npart;
  struct PART *part;
  struct PART *nexploc, *curploc;

  struct OCT *freeoct; // the first free oct available for refinement
  struct OCT *newoct;
  int nref=0,ndes=0;

  REAL xc,yc,zc;
  int stride;
  REAL **vcomp;
  int ncomp;
  REAL acc;
  REAL dt;
  int ntot=0,nlev,noct;
  REAL ntotd=0.,nlevd=0.;

  REAL disp,mdisp;
  
  int dir;

  char filename[128]; 
  FILE *fd;
  struct PART *nexp;
  struct PART *nexp2;
  struct PART *curp;

  struct PART *lastpart;

  int curc;
  REAL dtnew=0.;
  int nbnd;

  REAL x,y,z;
  REAL vx,vy,vz;
  REAL mass,mtot;
  REAL idx;
  REAL faexp, faexp2;
  REAL aexp;
  unsigned key;

  struct CPUINFO cpu;

  struct PACKET **sendbuffer; 
  struct PACKET **recvbuffer; 

  struct PART_MPI **psendbuffer; 
  struct PART_MPI **precvbuffer; 

  struct HYDRO_MPI **hsendbuffer;
  struct HYDRO_MPI **hrecvbuffer;

  struct FLUX_MPI **fsendbuffer;
  struct FLUX_MPI **frecvbuffer;

  struct RUNPARAMS param;

  size_t rstat;


  REAL omegam,omegav,Hubble,omegab=0.045;
  REAL avgdens; 
  REAL tmax;
#ifdef PIC
  avgdens=1.;//we assume a unit mass in a unit length box
#else
  avgdens=0.;
#endif


#ifdef TESTCOSMO
  double tab_aexp[NCOSMOTAB];
  double tab_ttilde[NCOSMOTAB];
  double tab_t[NCOSMOTAB];
  REAL ainit;
  REAL amax;

#endif

  //========== RIEMANN CHECK ====================/
#ifdef WHYDRO2
  int rtag=0;

#ifdef RIEMANN_EXACT
  rtag=1;
#endif

#ifdef RIEMANN_HLL
  rtag=2;
#endif

  if(rtag==0){
    printf("ERROR : RIEMANN SOLVER NOT SELECTED\n");
    abort();
  }
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
  
  
  /* Setup description of the 8 MPI_REAL fields data */
  offsets[0] = 0;
  oldtypes[0] = MPI_REAL;
  blockcounts[0] = 8;
  
  /* Setup description of the 2 MPI_INT fields key, level */
  /* Need to first figure offset by getting size of MPI_REAL */
  MPI_Type_extent(MPI_REAL, &extent);
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

  /* Setup description of the 7 MPI_REAL fields x,y,z,vx,vy,vz */
  offsets[0] = 0;
  oldtypes[0] = MPI_REAL;
  blockcounts[0] = 7;
  
  /* Setup description of the 4 MPI_INT fields idx key level icell*/
  /* Need to first figure offset by getting size of MPI_REAL */
  MPI_Type_extent(MPI_REAL, &extent);
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

#ifdef WHYDRO2
  //========= creating a WTYPE MPI type =======
  MPI_Datatype MPI_WTYPE;

  /* Setup description of the 6 MPI_REAL fields d,u,v,w,p,a */
  offsets[0] = 0;
  oldtypes[0] = MPI_REAL;
  blockcounts[0] = 6;

  /* Now define structured type and commit it */
  MPI_Type_struct(1, blockcounts, offsets, oldtypes, &MPI_WTYPE);
  MPI_Type_commit(&MPI_WTYPE);


  //========= creating a HYDRO MPI type =======
  MPI_Datatype MPI_HYDRO;

  /* Setup description of the 8 MPI_WTYPE fields one per oct*/
  offsets[0] = 0;
  oldtypes[0] = MPI_WTYPE;
  blockcounts[0] = 8;

  /* Setup description of the 2 MPI_INT fields key, level */
  /* Need to first figure offset by getting size of MPI_REAL */
  MPI_Type_extent(MPI_WTYPE, &extent);
  offsets[1] = 8 * extent;
  oldtypes[1] = MPI_LONG;
  blockcounts[1] = 1;

  MPI_Type_extent(MPI_LONG, &extent);
  offsets[2] = offsets[1]+extent;
  oldtypes[2] = MPI_INT;
  blockcounts[2] = 1;

  /* Now define structured type and commit it */
  MPI_Type_struct(3, blockcounts, offsets, oldtypes, &MPI_HYDRO);
  MPI_Type_commit(&MPI_HYDRO);



  //========= creating a FLUX MPI type =======
  MPI_Datatype MPI_FLUX;

  /* Setup description of the 8 MPI_WTYPE fields one per oct*/
  offsets[0] = 0;
  oldtypes[0] = MPI_REAL;
  blockcounts[0] = NFLUX*8;

  /* Setup description of the 2 MPI_INT fields key, level */
  /* Need to first figure offset by getting size of MPI_REAL */
  MPI_Type_extent(MPI_REAL, &extent);
  offsets[1] = NFLUX * 8 * extent;
  oldtypes[1] = MPI_LONG;
  blockcounts[1] = 1;

  MPI_Type_extent(MPI_LONG, &extent);
  offsets[2] = offsets[1]+extent;
  oldtypes[2] = MPI_INT;
  blockcounts[2] = 1;

  /* Now define structured type and commit it */
  MPI_Type_struct(3, blockcounts, offsets, oldtypes, &MPI_FLUX);
  MPI_Type_commit(&MPI_FLUX);


  
#endif


  //============================================

  cpu.MPI_PACKET=&MPI_PACKET;
#ifdef PIC
  cpu.MPI_PART=&MPI_PART;
#endif

#ifdef WHYDRO2
  cpu.MPI_HYDRO=&MPI_HYDRO;
  cpu.MPI_FLUX=&MPI_FLUX;
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
  
  //  if(cpu.rank==0) printf("Allocating %f GB cell=%f GB part=%f GB book=%f",(sizeof(struct OCT)*ngridmax+sizeof(struct PART)*npart+cpu.maxhash*sizeof(struct OCT*)+stride*ncomp*sizeof(REAL))/(1024*1024*1024.),sizeof(struct OCT)*ngridmax/(1024*1024*1024.),sizeof(struct PART)*npart/(1024*1024*1024.),(cpu.maxhash*sizeof(struct OCT*)+stride*ncomp*sizeof(REAL))/(1024.*1024.*1024.));

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

  vcomp=(REAL **)calloc(ncomp,sizeof(REAL*));
  for(i=0;i<ncomp;i++)
    {
      vcomp[i]=(REAL *)calloc(stride,sizeof(REAL));
    }

  memsize+=ncomp*stride*sizeof(REAL);


  //===================================================================================================




  // allocating the vectorized tree
  
  struct MULTIVECT vectors;
  
#ifdef NEWJACK
  vectors.vecpot=(REAL*)calloc(stride*8,sizeof(REAL));
  vectors.vecpotnew=(REAL*)calloc(stride*8,sizeof(REAL));
  vectors.vecden=(REAL*)calloc(stride*8,sizeof(REAL));



  vectors.vecnei=(int *)calloc(stride*6,sizeof(int));
  vectors.vecl=(int *)calloc(stride,sizeof(int));
  vectors.veccpu=(int *)calloc(stride,sizeof(int));
  vectors.vecicoarse=(int *)calloc(stride,sizeof(int));
  memsize+= stride*32*4;
  if(cpu.rank==0) printf(" vect = %f MB\n",((stride*32)/(1024*1024.))*sizeof(REAL));
#endif 
  
  // allocating the 6dim stencil
  struct HGRID *stencil;
  printf("stenci=%p\n",stencil);
  stencil=(struct HGRID*)calloc(stride,sizeof(struct HGRID));
  printf("stenci=%p mem=%f\n",stencil,stride*sizeof(struct HGRID)/(1024.*1024.));
  if(stencil==NULL) abort();
  

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
  cudaMalloc((void **)&(vectors.vecden_d),sizeof(REAL)*stride*8);
  cudaMalloc((void **)&(vectors.vecpot_d),sizeof(REAL)*stride*8);


  // temp GPU arrays
  cudaMalloc((void **)&(vectors.vecpotnew_d),sizeof(REAL)*stride*8);
  cudaMalloc((void **)&(vectors.vec2_d),sizeof(REAL)*stride*8);
  cudaMalloc((void **)&(vectors.vecsum_d),sizeof(REAL)*stride*8);

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


#ifdef WHYDRO2
    memset(&(grid->cell[icell].field),0,sizeof(struct Wtype));
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

#ifdef WHYDRO2
	      memset(&(newoct->cell[icell].field),0,sizeof(struct Wtype));
#endif

#ifdef WGRAV
	      memset(&(newoct->cell[icell].gdata),0,sizeof(struct Gtype));
	      memset(newoct->cell[icell].f,0,sizeof(REAL)*3);
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

#ifdef PIC
  psendbuffer=(struct PART_MPI **)(calloc(cpu.nnei,sizeof(struct PART_MPI*)));
  precvbuffer=(struct PART_MPI **)(calloc(cpu.nnei,sizeof(struct PART_MPI*)));
  for(i=0;i<cpu.nnei;i++) {
    psendbuffer[i]=(struct PART_MPI *) (calloc(cpu.nbuff,sizeof(struct PART_MPI)));
    precvbuffer[i]=(struct PART_MPI *) (calloc(cpu.nbuff,sizeof(struct PART_MPI)));
  }
#endif 

#ifdef WHYDRO2
  hsendbuffer=(struct HYDRO_MPI **)(calloc(cpu.nnei,sizeof(struct HYDRO_MPI*)));
  hrecvbuffer=(struct HYDRO_MPI **)(calloc(cpu.nnei,sizeof(struct HYDRO_MPI*)));
  for(i=0;i<cpu.nnei;i++) {
    hsendbuffer[i]=(struct HYDRO_MPI *) (calloc(cpu.nbuff,sizeof(struct HYDRO_MPI)));
    hrecvbuffer[i]=(struct HYDRO_MPI *) (calloc(cpu.nbuff,sizeof(struct HYDRO_MPI)));
  }

  fsendbuffer=(struct FLUX_MPI **)(calloc(cpu.nnei,sizeof(struct FLUX_MPI*)));
  frecvbuffer=(struct FLUX_MPI **)(calloc(cpu.nnei,sizeof(struct FLUX_MPI*)));
  for(i=0;i<cpu.nnei;i++) {
    fsendbuffer[i]=(struct FLUX_MPI *) (calloc(cpu.nbuff,sizeof(struct FLUX_MPI)));
    frecvbuffer[i]=(struct FLUX_MPI *) (calloc(cpu.nbuff,sizeof(struct FLUX_MPI)));
  }
#endif 



#endif


  //===================================================================================================
  
  // ==== some initial dump

  /* sprintf(filename,"data/levstart.%05d.p%05d",0,cpu.rank); */
  /* dumpcube(lmap,firstoct,0,filename,0.); */
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
  REAL dxcell=1./pow(2.,levelcoarse);
  REAL epsilon=0.;
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
      vy=sqrt((1.-epsilon)/0.1)*0.6;
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
  REAL dxcell=1./pow(2.,levelcoarse);
  REAL th,ph,r;
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


      th=acos(((REAL)(rand())/RAND_MAX*2-1.));
      ph=2*M_PI*(REAL)(rand())/RAND_MAX;
      r=(REAL)(rand())/RAND_MAX*0.3;

      x=r*sin(th)*cos(ph)+0.5;
      y=r*sin(th)*sin(ph)+0.5;
      z=r*cos(th)+0.5;

      vx=(REAL)(rand())/RAND_MAX*2.-1.;
      vy=(REAL)(rand())/RAND_MAX*2.-1.;
      vz=(REAL)(rand())/RAND_MAX*2.-1.;
      
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
  REAL dummyf;
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
  REAL dummyf;
  int npartf;

  int nploc;
  REAL munit;
  REAL lbox;

#ifdef GRAFIC // ==================== read grafic file

  lastpart=read_grafic_part(part, &cpu, &munit, &ainit, &omegam, &omegav, &Hubble, &npart,omegab);
  amax=1.;

#else // ============================  explicit read here
#ifndef ZELDO
  fd=fopen("utils/IC.PM.0","rb");
#else
  fd=fopen("utils/ZEL.PM.0","rb");
#endif
  rstat=fread(&dummy,sizeof(dummy),1,fd);
  rstat=fread(&nploc,sizeof(int),1,fd);
  rstat=fread(&munit,sizeof(REAL),1,fd);
  rstat=fread(&ainit,sizeof(REAL),1,fd);
  rstat=fread(&lbox,sizeof(REAL),1,fd);
  rstat=fread(&omegam,sizeof(REAL),1,fd);
  rstat=fread(&omegav,sizeof(REAL),1,fd);
  rstat=fread(&Hubble,sizeof(REAL),1,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);

  if(cpu.rank==0) printf("%f %d %f %f\n",ainit,nploc,omegav,Hubble);
  
  mass=munit;
  tsim=ainit;
  aexp=ainit;

  // at this stage we have to compute the conformal time
  tsim=-0.5*sqrt(omegam)*integ_da_dt_tilde(ainit,1.0,omegam,omegav,1e-8);

  // we compute the friedmann tables


  compute_friedmann(ainit*0.95,NCOSMOTAB,omegam,omegav,tab_aexp,tab_ttilde,tab_t);


  // ------------------ dealing with particles

  REAL *pos;
  REAL *vel;
  int nread=(nploc<2097152?nploc:2097152); // we read by patch of 128^3
  int npatch=nploc/nread;
  int ipatch;

  pos=(REAL *)malloc(sizeof(REAL)*3*nread);
  vel=(REAL *)malloc(sizeof(REAL)*3*nread);
  int pstart=ftell(fd);

  ip=0.;
  for(ipatch=0;ipatch<npatch;ipatch++) {
    //    rstat=fread(&dummy,sizeof(dummy),1,fd); 
    //    fseek(fd,pstart,SEEK_SET);
    fseek(fd,pstart+(0*nploc+ipatch*nread)*sizeof(REAL)+1*sizeof(dummy),SEEK_SET);
    rstat=fread(pos,sizeof(REAL),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);

    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(1*nploc+ipatch*nread)*sizeof(REAL)+3*sizeof(dummy),SEEK_SET);
    rstat=fread(pos+nread,sizeof(REAL),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);

    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(2*nploc+ipatch*nread)*sizeof(REAL)+5*sizeof(dummy),SEEK_SET);
    rstat=fread(pos+2*nread,sizeof(REAL),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
  
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(3*nploc+ipatch*nread)*sizeof(REAL)+7*sizeof(dummy),SEEK_SET);
    rstat=fread(vel,sizeof(REAL),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);

    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(4*nploc+ipatch*nread)*sizeof(REAL)+9*sizeof(dummy),SEEK_SET);
    rstat=fread(vel+nread,sizeof(REAL),nread,fd);
    //    rstat=fread(&dummy,sizeof(dummy),1,fd);

    //    rstat=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(5*nploc+ipatch*nread)*sizeof(REAL)+11*sizeof(dummy),SEEK_SET);
    rstat=fread(vel+2*nread,sizeof(REAL),nread,fd);
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


#ifdef WHYDRO2
  
#ifdef GRAFIC
  int ncellhydro;
  ncellhydro=read_grafic_hydro(&cpu,&ainit, &omegam, &omegav, &Hubble,omegab);

  printf("%d hydro cell found in grafic file with aexp=%e\n",ncellhydro,ainit);
  amax=1.0;
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
  REAL X0;
  if(cpu.rank==0) printf("Init Hydro\n");

   /* // TEST 1 */

  WL.d=1.;
  WL.u=0.;
  WL.v=0.;
  WL.w=0.;
  WL.p=1.0;

  WR.d=0.125;
  WR.u=0.;
  WR.v=0.;
  WR.w=0.;
  WR.p=0.1;

  X0=0.3;
  WL.a=sqrt(GAMMA*WL.p/WL.d);
  WR.a=sqrt(GAMMA*WR.p/WR.d);

  // ======================================================

  REAL dtot=0.;
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
	      xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
	      yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
	      zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

	      curoct->cell[icell].pot=GRAV*zc;

 	      /* RT INSTAB */

	      /* REAL amp=0.05; */
	      /* /\* REAL vrx=(((REAL)rand())/RAND_MAX)*2.*amp-amp; *\/ */
	      /* /\* REAL vry=(((REAL)rand())/RAND_MAX)*2.*amp-amp; *\/ */
	      /* /\* REAL vrz=(((REAL)rand())/RAND_MAX)*2.*amp-amp; *\/ */

	      /* REAL vrx=0.; */
	      /* REAL vry=0.; */
	      /* REAL vrz=-amp*(1.+cos(8.*M_PI*(xc-0.5)));//\*(1.+cos(8.*M_PI*(yc-0.5)))*(1.+cos(2.*M_PI*(zc-0.5)))/8.; */

	      /* curoct->cell[icell].field.u=vrx; */
	      /* curoct->cell[icell].field.v=vry; */
	      /* curoct->cell[icell].field.w=vrz; */

	     
	      
	      /* if(zc>0.75){ */
	      /* 	curoct->cell[icell].field.d=2.; */
	      /* } */
	      /* else{ */
	      /* 	curoct->cell[icell].field.d=1.; */
		
	      /* 	} */
	      	
	      /* curoct->cell[icell].field.p=2.5-curoct->cell[icell].field.d*GRAV*(zc-0.9); */
	      /* curoct->cell[icell].field.a=sqrt(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d); */


 	      /* /\* KH INSTAB *\/ */

	      /* REAL amp=0.05; */
	      /* /\* REAL vrx=(((REAL)rand())/RAND_MAX)*2.*amp-amp; *\/ */
	      /* /\* REAL vry=(((REAL)rand())/RAND_MAX)*2.*amp-amp; *\/ */
	      /* /\* REAL vrz=(((REAL)rand())/RAND_MAX)*2.*amp-amp; *\/ */

	      /* REAL vrx=amp*sin(2.*M_PI*xc); */
	      /* REAL vry=amp*sin(2.*M_PI*xc); */
	      /* REAL vrz=amp*sin(2.*M_PI*xc); */

	      /* if((zc>0.75)||(zc<0.25)){ */
	      /* 	curoct->cell[icell].field.d=1.0; */
	      /* 	curoct->cell[icell].field.u=0.5+vrx; */
	      /* 	curoct->cell[icell].field.v=vry; */
	      /* 	curoct->cell[icell].field.w=vrz; */
	      /* 	curoct->cell[icell].field.p=2.5; */
	      /* 	curoct->cell[icell].field.a=sqrt(GAMMA*2.5/1.); */
	      /* } */
	      /* else{ */
	      /* 	curoct->cell[icell].field.d=2.0; */
	      /* 	curoct->cell[icell].field.u=-0.5+vrx; */
	      /* 	curoct->cell[icell].field.v=vry; */
	      /* 	curoct->cell[icell].field.w=vrz; */
	      /* 	curoct->cell[icell].field.p=2.5; */
	      /* 	curoct->cell[icell].field.a=sqrt(GAMMA*2.5/2.); */
	      /* } */

	      /* SPHERICAL EXPLOSION */
	      
	      /* if((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5)<(X0*X0)){ */
	      /* 	memcpy(&(curoct->cell[icell].field),&WL,sizeof(struct Wtype)); */
	      /* } */
	      /* else{ */
	      /*  	memcpy(&(curoct->cell[icell].field),&WR,sizeof(struct Wtype)); */
	      /* } */
	      


	      /* ZELDOVICH PANCAKE */

	      REAL ZI=99;
	      REAL ZC=9;
	      omegam=1.0;
	      omegav=0.;
	      ainit=1./(1.+ZI);
	      amax=1./(1.+ZC);
	      curoct->cell[icell].field.d=1.+(1.+ZC)/(1.+ZI)*cos(2.*M_PI*(xc-0.5));
	      curoct->cell[icell].field.u=0.;//-(1.+ZC)/pow(1.+ZI,1.5)*sin(2.*M_PI*(xc-0.5))/(2.*M_PI);
	      curoct->cell[icell].field.v=0.;
	      curoct->cell[icell].field.w=0.;
	      curoct->cell[icell].field.p=1e-9;
	      curoct->cell[icell].field.a=sqrt(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d);



	      curoct->cell[icell].field.x=xc;
	      curoct->cell[icell].field.y=yc;
	      curoct->cell[icell].field.z=zc;



	      /* /\* SHOCK TUBE *\/ */
	      /* if(xc<=X0){ */

	      /* 	memcpy(&(curoct->cell[icell].field),&WL,sizeof(struct Wtype)); */
	      /* } */
	      /* else{ */
	      /* 	memcpy(&(curoct->cell[icell].field),&WR,sizeof(struct Wtype)); */
	      /* } */
	      
	      if(level==levelcoarse) {
	      	dtot+=curoct->cell[icell].field.d;
	      	nc++;
	      }

	    }
	}while(nextoct!=NULL);
      
      //printf("level=%d avg=%e mind=%e maxd=%e\n",level,avg/ncell,mind,maxd);
    }
  
  avgdens+=dtot/nc;
  printf("avgdens=%e\n",avgdens);

#endif

  /* sprintf(filename,"data/denhydstart.%05d.p%05d",0,cpu.rank);  */
  /* dumpcube(lmap,firstoct,101,filename,0.);  */
  /* sprintf(filename,"data/prehydstart.%05d.p%05d",0,cpu.rank);  */
  /* dumpcube(lmap,firstoct,105,filename,0.);  */


#ifdef TESTCOSMO
  tsim=ainit;
  aexp=ainit;

  // at this stage we have to compute the conformal time
  tsim=-0.5*sqrt(omegam)*integ_da_dt_tilde(ainit,1.0,omegam,omegav,1e-8);

  // we compute the friedmann tables

  compute_friedmann(ainit*0.95,NCOSMOTAB,omegam,omegav,tab_aexp,tab_ttilde,tab_t);

  tmax=-0.5*sqrt(omegam)*integ_da_dt_tilde(amax,1.0,omegam,omegav,1e-8);
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

#endif




#ifdef WGRAV
  /* REAL rc; */
  /* for(level=levelcoarse;level<=levelmax;level++) // (levelcoarse only for the moment) */
  /*   { */
  /*     dxcur=pow(0.5,level); */
  /*     nextoct=firstoct[level-1]; */
  /*     if(nextoct==NULL) continue; */
  /*     do // sweeping level */
  /* 	{ */
  /* 	  curoct=nextoct; */
  /* 	  nextoct=curoct->next; */
  /* 	  for(icell=0;icell<8;icell++) // looping over cells in oct */
  /* 	    { */
  /* 	      xc=curoct->x+( icell&1)*dxcur+dxcur*0.5; */
  /* 	      yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5; */
  /* 	      zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5; */
	    
  /* 	      rc=sqrt(pow(xc-0.5,2)+pow(yc-0.5,2)+pow(zc-0.5,2)); */
  
  /* 	      if(rc<=3.*dxcur){ */
  /* 		curoct->cell[icell].gdata.d=1.; */
  /* 		curoct->cell[icell].gdata.p=0.; */
  /* 	      } */

  /* 	    } */
  /* 	}while(nextoct!=NULL); */
      
  /*     //printf("level=%d avg=%e mind=%e maxd=%e\n",level,avg/ncell,mind,maxd); */
  /*   } */

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
  int ndumps=0;
  int pass;
  int smark;
  int ismooth,nsmooth=2;
  int marker;



#ifndef TESTCOSMO
#ifdef ZELDO
  tmax=-0.5*sqrt(omegam)*integ_da_dt_tilde(0.5,1.0,omegam,omegav,1e-8);
#else
  tmax=1000.;
#endif
#ifdef WHYDRO2
  tmax=0.2;
#endif
#endif

  FILE *fegy;

  //breakmpi();

#ifndef TESTCOSMO
    faexp=1.0;
    faexp2=1.0;
#endif

    


    //==================================== Restart =================================================

    if(param.nrestart!=0){
      FILE *fp;
      int ioct=0;
      struct OCT oct;
      struct OCT *ozero;

      sprintf(filename,"data/grid.%05d.p%05d",param.nrestart,cpu.rank); 

      printf("\n ====> Restarting from %s\n",filename);

      fp=fopen(filename,"rb");
      
      if(fp==NULL){
	printf("--ERROR -- Restart file not found\n");
	return 1;
      }
      
      // reading the time
      fread(&tsim,sizeof(REAL),1,fp);
      printf("tsim=%f\n",tsim);
      // reading the zero oct
      fread(&ozero,sizeof(struct OCT *),1,fp);

    
      fread(&oct,sizeof(struct OCT),1,fp);
      int curlev=-1;
      
      while(!feof(fp)){
	
	// copying raw data in the correct oct
	memcpy(grid+ioct,&oct,sizeof(struct OCT));
	
	//dealing with level per level data
	if(curlev!=oct.level){
	  // a new level has been found
	  firstoct[oct.level-1]=grid+ioct;
	  if(oct.level!=1){
	    lastoct[oct.level-2]=grid+ioct-1;
	  }
	  curlev=oct.level;
	  printf("curlev=%d\n",curlev);
	}


	//dealing with OCT pointers

	grid[ioct].next+=(grid-ozero);
	grid[ioct].prev+=(grid-ozero);
	grid[ioct].nexthash+=(grid-ozero);

	//dealing with CELL pointers
	long int adr;
	
	adr=(long int)grid[ioct].parent;
	grid[ioct].parent=(struct CELL *) (adr+(long int)(grid-ozero));
	
	for(i=0;i<6;i++){
	  adr=(long int)grid[ioct].nei[i];
	  grid[ioct].nei[i]=(struct CELL *) (adr+(long int)(grid-ozero));
	}


	// dealing with pointers in CELLS
	for(i=0;i<8;i++){
	  grid[ioct].cell[i].child+=(grid-ozero);
	}
	

#ifdef PIC
	//dealing with particles
	printf("ERROR RESTART NOT YET IMPLEMENTED WITH PARTICLES !\n");
	abort();
#endif

	ioct++;
	fread(&oct,sizeof(struct OCT),1,fp); //reading next oct
      }

      freeoct=grid+ioct;
      ndumps=param.nrestart+1;
    }




    //==================================== MAIN LOOP ================================================
    //===============================================================================================

    for(nsteps=0;(nsteps<=param.nsteps)*(tsim<=tmax);nsteps++){
    
#ifdef TESTCOSMO
      aexp=interp_aexp(tsim,tab_aexp,tab_ttilde);
      if(cpu.rank==0) printf("\n============== STEP %d aexp=%e z=%f t=%e tmax=%e================\n",nsteps,aexp,1./aexp-1.,tsim,tmax);
#else
      if(cpu.rank==0) printf("\n============== STEP %d tsim=%e ================\n",nsteps,tsim);
#endif

      
#if 1
      if(levelcoarse!=levelmax){
	
	// ==================================== marking the cells
	mark_cells(levelcoarse,levelmax,firstoct,nsmooth,threshold,&cpu,sendbuffer,recvbuffer);
	
	
	// ==================================== refining (and destroying) the octs
	
	curoct=refine_cells(levelcoarse,levelmax,firstoct,lastoct,freeoct,&cpu,grid+ngridmax);
	freeoct=curoct;
      }
    
#endif



#ifdef WMPI
      // ==================================== after refinement we should remap the boundary cells
      setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax,newloadb);
#endif
      
      
      // ==================================== performing the CIC assignement
#ifdef PIC
#ifndef WGPU
      call_cic(levelmax,levelcoarse,firstoct,&cpu);
#else
      
      call_cic(levelmax,levelcoarse,firstoct,&cpu);
#endif
#endif

#ifdef PIC
#ifdef WMPI
      //------------- performing the CIC BOUNDARY CORRECTION 
      mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);
      //------------  Density boundary mpi update 
      mpi_exchange(&cpu,sendbuffer,recvbuffer,1,1);
#endif
#endif
      
      //======================================= cleaning the marks
      
      clean_marks(levelmax,firstoct);
      
      // ==================================== Check the number of particles and octs
      mtot=multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);
      

      // =========== Grid Census : CLEAN THE NEXT FOLLOWING LINES ========================================
#ifdef WMPI
      MPI_Allreduce(MPI_IN_PLACE,&mtot,1,MPI_REAL,MPI_SUM,cpu.comm);
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
	  if(ltot!=0) printf("level=%2d noct=%9d min=%9d max=%9d\n",level,ltot,nomin,nomax);
	}
      }
#ifdef WMPI
      MPI_Allreduce(MPI_IN_PLACE,&gtot,1,MPI_INT,MPI_MAX,cpu.comm);
#endif
      if(cpu.rank==0){
	printf("grid occupancy=%4.1f \n",(gtot/(1.0*ngridmax))*100.);
      }

      // ==========================================================================================
      

#if 1
      // FOR RT INSTABILITY ONLY !
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
		  xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		  yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		  zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;
		  
		  curoct->cell[icell].pot=GRAV*zc;
		}
	    }while(nextoct!=NULL);
	}
#endif 
      

#ifdef WHYDRO2
#ifdef WMPI
    // ================================= exchange current state of hydro quantities 
    mpi_exchange_hydro(&cpu, hsendbuffer, hrecvbuffer,1);
    MPI_Barrier(cpu.comm);
#endif
#endif

#ifdef WGRAV 
    //==================================== Getting Density ====================================
    for(level=levelcoarse;level<=levelmax;level++)
      {
	FillDens(level,&param,firstoct,&cpu);
      }

    //===================================== JACOBI Poisson Solver ==============================
    if(cpu.rank==0) printf("Start Poisson Solver \n");
    REAL res;
    int igrid;
    struct Gtype Wi[8];
    struct CELL* curcell;
    int icell2;
    //breakmpi();
    for(level=levelcoarse;level<=levelmax;level++)
      {
	if((level==levelcoarse)&&(levelcoarse!=param.mgridlmin)){
	  for(igrid=0;igrid<param.nvcycles;igrid++){ // V-Cycles
	    printf("----------------------------------------\n");
	    res=PoissonMgrid(level,&param,firstoct,&cpu,stencil,stride,aexp);
	    if(res<param.poissonacc) break;
	  }
	}
	else{
	  PoissonJacobi(level,&param,firstoct,&cpu,stencil,stride,aexp);
	}
	
	//once done we propagate the solution to level+1

	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		curcell=&(curoct->cell[icell]);
		if(curcell->child!=NULL){
		  coarse2fine_grav(curcell,Wi);
		  for(icell2=0;icell2<8;icell2++){
		    memcpy(&(curcell->child->cell[icell2].gdata.p),&(Wi[icell2].p),sizeof(REAL));
		  }
		}
	      }
	  }while(nextoct!=NULL);
      }
    //==================================== End Poisson Solver ==========================
#endif



  // ==================================== Force calculation and velocity update   // Corrector step
#ifdef TESTCOSMO
  faexp=1.0;
#endif
  // -- force
#ifdef WGRAV  
  for(level=levelcoarse;level<=levelmax;level++)
    {
      PoissonForce(level,&param,firstoct,&cpu,stencil,stride);
    }
#endif

#ifndef NOCOUPLE
#ifdef WHYDRO2
#ifdef WGRAV
  // ================================== Gravitational source correction

  if(nsteps!=0){
  for(level=levelcoarse;level<=levelmax;level++)
    {
      struct Utype U0,U;
      struct Wtype W;
      struct Wtype Wnew;
      curoct=firstoct[level-1];
	if((curoct!=NULL)&&(cpu.noct[level-1]!=0)){
	  nextoct=curoct;
	  do{
	    curoct=nextoct;
	    nextoct=curoct->next;
	    if(curoct->cpu!=cpu.rank) continue; // we don't update the boundary cells
	    for(icell=0;icell<8;icell++){
	      int ref=0;

	      if(curoct->cell[icell].child==NULL){ // Leaf cell
		struct CELL *curcell;
		curcell=&(curoct->cell[icell]);

		memcpy(&Wnew,&(curcell->field),sizeof(struct Wtype));

#ifdef PRIMITIVE
		Wnew.u+=(-curcell->f[0]*dt*0.5);
		Wnew.v+=(-curcell->f[1]*dt*0.5);
		Wnew.w+=(-curcell->f[2]*dt*0.5);
#endif

#ifdef CONSERVATIVE
		W2U(&Wnew,&U);
		memcpy(&U0,&U,sizeof(struct Utype));
		// grav force correction

		U.du+=-(U0.d*curoct->cell[icell].f[0]*dt*0.5);
		U.dv+=-(U0.d*curoct->cell[icell].f[1]*dt*0.5);
		U.dw+=-(U0.d*curoct->cell[icell].f[2]*dt*0.5);
		U.E+=-(U.du*curoct->cell[icell].f[0]+U.dv*curoct->cell[icell].f[1]+U.dw*curoct->cell[icell].f[2])*dt*0.5;
		U2W(&U,&Wnew);
#endif	

		if(Wnew.p<0){
		  printf("oulah error in gravity coupling with pressure\n");
		  abort();
		}

		if(Wnew.d<0){
		  printf("oulah error in gravity coupling with pressure\n");
		  abort();
		}
		
		
		if((isnan(Wnew.p))||isnan(Wnew.d)){
		  printf("NAN\n");
		  abort();
		}
		

		memcpy(&(curcell->field),&Wnew,sizeof(struct Wtype));
	      }
	    }
	  }while(nextoct!=NULL);
	}
    }
  }
#endif
#endif
#endif

  // ==================================== DUMP AFTER SYNCHRONIZATION
  
  if(nsteps%(param.ndumps)==0){
    if(cpu.rank==0) printf("Dumping .......\n");
    REAL tdump=tsim;
    
    
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
    
#ifdef WHYDRO2
    printf("tdum=%f\n",tdump);
    sprintf(filename,"data/grid.%05d.p%05d",ndumps,cpu.rank); 
    dumpgrid(levelmax,firstoct,filename,tdump); 
#endif
    

    // === Gravity dump

#ifdef WGRAV
    printf("tdum=%f\n",tdump);
    sprintf(filename,"data/grid.%05d.p%05d",ndumps,cpu.rank); 
    dumpgrid(levelmax,firstoct,filename,tdump); 
    //abort();
#endif

    ndumps++;
  }


  // ============================================================================================
  // ==================================== New timestep ==========================================
  // ============================================================================================

  dtnew=dt;
#ifdef TESTCOSMO
  REAL dtcosmo;
  dtcosmo=-0.5*sqrt(omegam)*integ_da_dt_tilde(aexp*1.1,aexp,omegam,omegav,1e-8);
  dtnew=(dtcosmo<dtnew?dtcosmo:dtnew);
  printf("dtcosmo= %e ",dtcosmo);
#endif
  
#ifdef PIC
  REAL dtpic;
  dtpic=comptstep(levelcoarse,levelmax,firstoct,faexp2,faexp,&cpu,1e9);
  dtnew=(dtpic<dtnew?dtpic:dtnew);
  printf("dtpic= %e ",dtpic);
#endif
  

#ifdef WHYDRO2
  REAL dthydro;
  dthydro=comptstep_hydro(levelcoarse,levelmax,firstoct,faexp2,faexp,&cpu,1e9);
  //if(nsteps<5) dthydro/=5.;
  dtnew=(dthydro<dtnew?dthydro:dtnew);
  //dtnew=dthydro;
  //  if(cpu.rank==0) printf("dt=%e dthydro=%e\n",dtnew,dthydro);
  printf("dthydro= %e ",dthydro);

#ifdef WGRAV
  REAL dtff;
  dtff=comptstep_ff(levelcoarse,levelmax,firstoct,aexp,&cpu,1e9);
  dtnew=(dtff<dtnew?dtff:dtnew);
  printf("dtff= %e ",dtff);

  REAL dtforce;
  dtforce=comptstep_force(levelcoarse,levelmax,firstoct,aexp,&cpu,1e9);
  printf("dtforce= %e ",dtforce);
  dtnew=(dtforce<dtnew?dtforce:dtnew);


#endif
#endif



#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&dtnew,1,MPI_REAL,MPI_MIN,cpu.comm);
#endif

  printf("dtnew= %e \n",dtnew);


  


#ifdef WHYDRO2

  // ================================= performing hydro calculations
    
  double t0,t100,t20,t80;
  int nread,nreadtot;;
  REAL deltamax=0.;
  if(cpu.rank==0) printf("Start Hydro 2\n");
  //breakmpi();
  for(level=levelmax;level>=levelcoarse;level--)
    {
      
      //nreadtot=compute_hydro_fluxes(level,firstoct,stencil,stride,&cpu,aexp);

      dxcur=pow(0.5,level);

      // --------------- setting the first oct of the level
      nextoct=firstoct[level-1];
      //printf("lev=%d curoct=%p\n",level,nextoct);
      nreadtot=0;
      if((nextoct!=NULL)&&(cpu.noct[level-1]!=0)){
	do {
	  curoct=nextoct;
	  nextoct=curoct->next; 
	  //if(curoct->cpu!=cpu.rank) continue;
	  // -------------  cleaning working arrays
	    
	  memset(stencil,0,stride*sizeof(struct HGRID));
	    
	  // ------------ gathering the stencil value values
	  nextoct=gatherstencil(curoct,stencil,stride,&cpu, &nread);
	  //printf("nread=%d\n",nread);
	  // ------------ solving the hydro
	    
	  //t20=MPI_Wtime();
	  hydroM(stencil,level,cpu.rank,nread,stride,dxcur,dtnew);
	  //t80=MPI_Wtime();
	    	    
	  // ------------ scatter back the FLUXES
	    
	  nextoct=scatterstencil(curoct,stencil, nread, &cpu);
	  nreadtot+=nread;
	  //t100=MPI_Wtime();
	}while(nextoct!=NULL);
      }
      //printf("level=%d Nhydro=%d on proc %d\n",level,nreadtot,cpu.rank);

      // ---------------- at this stage we are ready to update the conservative variables
      int flx;
      REAL F[NFLUX];
      REAL Forg[NFLUX];
      REAL dtsurdx=dtnew/dxcur;
      REAL one;
      struct Utype U;
      struct Utype S;
      struct Utype U0;
      struct Wtype W;
      struct Wtype Wnew;
      struct CELL *neicell;
#ifdef DUAL_E
      REAL p,p0;
      REAL DE;
#endif

#ifdef WMPI
      // ================================= exchange current state of hydro quantities 
      MPI_Barrier(cpu.comm);
      mpi_exchange_flux(&cpu, fsendbuffer, frecvbuffer,1);
      MPI_Barrier(cpu.comm);
#endif

      //printf("dtsurdx=%e on proc %d at level=%d (dtnew=%e dxcur=%e)\n",dtsurdx,cpu.rank,level,dtnew,dxcur);

      if(nreadtot>0){

	curoct=firstoct[level-1];
	if((curoct!=NULL)&&(cpu.noct[level-1]!=0)){


	  nextoct=curoct;
	  do{
	    curoct=nextoct;
	    nextoct=curoct->next;
	    if(curoct->cpu!=cpu.rank) continue; // we don't update the boundary cells
	    for(icell=0;icell<8;icell++){
	      int ref=0;

	      if(curoct->cell[icell].child==NULL){ // Leaf cell
		struct CELL *curcell;
		curcell=&(curoct->cell[icell]);
		memcpy(&W,&(curcell->field),sizeof(struct Wtype));
		W2U(&W,&U);
		memcpy(&U0,&U,sizeof(struct Utype));
		

#ifdef DUAL_E
		DE=W.p/((GAMMA-1.)*U.E);
		p0=W.p;
		p=p0;
#endif
		
		memcpy(F,curcell->flux,sizeof(REAL)*NFLUX);// original fluxes
		
		// here we have to deal with coarse-fine boundaries

		if(level<levelmax){
		  int inei;
		  getcellnei(icell, vnei, vcell);

		  //loop over neighbours
		  for(inei=0;inei<6;inei++){
		    if(vnei[inei]!=6){
		    
		      if(curoct->nei[vnei[inei]]->child!=NULL){
			
#ifdef TRANSXP
			if(inei==1){
			  if((curoct->nei[inei]->child->x-curoct->x)<0.){
			    continue;
			  }
			}
#endif

#ifdef TRANSYP
			if(inei==3){
			  if((curoct->nei[inei]->child->y-curoct->y)<0.){
			    continue;
			  }
			}
#endif

#ifdef TRANSZP
			if(inei==5){
			  if((curoct->nei[inei]->child->z-curoct->z)<0.){
			    continue;
			  }
			}
#endif

#ifdef TRANSXM
			if(inei==0){
			  if((curoct->nei[inei]->child->x-curoct->x)>0.5){
			    continue;
			  }
			}
#endif

#ifdef TRANSYM
			if(inei==2){
			  if((curoct->nei[inei]->child->y-curoct->y)>0.5){
			    continue;
			  }
			}
#endif

#ifdef TRANSZM
			if(inei==4){
			  if((curoct->nei[inei]->child->z-curoct->z)>0.5){
			    continue;
			  }
			}
#endif
		      
			// the neighbor cell is at the same level or refined
			neicell=&(curoct->nei[vnei[inei]]->child->cell[vcell[inei]]);
			
			if(neicell->child!=NULL){
			  // the neighbor is split : fluxes must be averaged
			  int fcell[4];
			  getfcell(inei,fcell);
			  memset(F+NVAR*inei,0,NVAR*sizeof(REAL)); // reset the original flux
			
			  int iface;
			  REAL *Fnei;
			  int idxfnei[6]={1,0,3,2,5,4};
			  int j; 
			  // averaging the flux
			  for(iface=0;iface<4;iface++){
			    Fnei=neicell->child->cell[fcell[iface]].flux;
			    for(j=0;j<NVAR;j++) F[j+inei*NVAR]+=0.25*Fnei[j+idxfnei[inei]*NVAR];
			  }

			  /* if((neicell->child->level==6)&&(curoct->level==5)){ */
			  /*   printf("BOUND\n"); */
			  /*   for(flx=0;flx<6;flx++) printf("%e -> %e ** \n ",Forg[1+flx*5],F[1+flx*5]); */
			  /*   abort(); */
			  /* } */


			}
		      }
		    }
		  }
		}
		
		if((isnan(U.d))||isnan(U.E)){
		  printf("NAN0\n");
		  abort();
		}
		

		// ready to update
		one=1.;
		for(flx=0;flx<6;flx++){
		  U.d +=F[0+flx*NVAR]*dtsurdx*one;
		  U.du+=F[1+flx*NVAR]*dtsurdx*one;
		  U.dv+=F[2+flx*NVAR]*dtsurdx*one;
		  U.dw+=F[3+flx*NVAR]*dtsurdx*one;
		  U.E +=F[4+flx*NVAR]*dtsurdx*one;
#ifdef DUAL_E
		  p   +=F[5+flx*NVAR]*dtsurdx*one;
#endif
		  one*=-1.;
		}

		U2W(&U,&Wnew);

		if(Wnew.w!=0.) abort();
		/* if(Wnew.p<0){ */
		/*   printf("oulah error before gravity coupling\n"); */
		/*   abort(); */
		/* } */

#ifdef DUAL_E
		REAL divu=0.;
		REAL pavg=0.;
		REAL uloc;
		int inei;
		if(DE<1e-3){
		  // finish the pressure fix, needs the divergence of the velocity (1D for the moment)
		  // DOn't work with refinement !!!!!
		
		  getcellnei(icell, vnei, vcell);
		  
		  for(inei=0;inei<2;inei++){
		    // getting the neighbor
		    if(vnei[inei]!=6){
		      neicell=&(curoct->nei[vnei[inei]]->child->cell[vcell[inei]]);
		    }
		    else{
		      neicell=&(curoct->cell[vcell[inei]]);
		    }
		    uloc=neicell->field.u;
		    //pavg+=neicell->field.p;
		    divu+=pow(-1.,inei+1)*uloc;
		  }
		  divu/=(2.*dxcur);
		  //pavg/=2.;
		  pavg=p;
		  
		  p=p-(GAMMA-1.)*pavg*divu; // pressure is updated GAMMA OR GAMMA-1 ?
		}
#endif
		

#ifdef WGRAV
#ifndef NOCOUPLE
		/* // half gravitational force correction */

#ifdef CONSERVATIVE
		U.du+=-(U0.d*curoct->cell[icell].f[0]*dtnew*0.5);
		U.dv+=-(U0.d*curoct->cell[icell].f[1]*dtnew*0.5);
		U.dw+=-(U0.d*curoct->cell[icell].f[2]*dtnew*0.5);
		U.E+=-(U0.du*curoct->cell[icell].f[0]+U0.dv*curoct->cell[icell].f[1]+U0.dw*curoct->cell[icell].f[2])*dtnew*0.5;
		U2W(&U,&Wnew);
		

#endif

#ifdef PRIMITIVE
		Wnew.u+=(-curoct->cell[icell].f[0]*dtnew*0.5);
		Wnew.v+=(-curoct->cell[icell].f[1]*dtnew*0.5);
		Wnew.w+=(-curoct->cell[icell].f[2]*dtnew*0.5);
#ifdef DUAL_E
		if(DE<1e-3){
		  Wnew.p=p;
		  Wnew.a=sqrt(GAMMA*Wnew.p/Wnew.d);
		}

#endif
#endif

		


#endif	
#endif	
	
		/* if(Wnew.p<0){ */
		/*   printf("oulah error in gravity coupling with pressure\n"); */
		/*   abort(); */
		/* } */


		if(Wnew.p<0){
		  printf("oulah error in gravity coupling with pressure\n");
		  abort();
		}
		
		if((isnan(Wnew.p))||isnan(Wnew.d)){
		  printf("NAN\n");
		  abort();
		}

		Wnew.x=curcell->field.x;
		Wnew.y=curcell->field.y;
		Wnew.z=curcell->field.z;
		

		memcpy(&(curcell->field),&Wnew,sizeof(struct Wtype));


	      }
	      else{ // split cell : hydro quantities are averaged
		struct OCT *child;
		int i;
		child=curoct->cell[icell].child;
		memset(&W,0,sizeof(struct Wtype));
		for(i=0;i<8;i++){
		  W.d+=child->cell[i].field.d*0.125;
		  W.u+=child->cell[i].field.u*0.125;
		  W.v+=child->cell[i].field.v*0.125;
		  W.w+=child->cell[i].field.w*0.125;
		  W.p+=child->cell[i].field.p*0.125;
		}
		memcpy(&curoct->cell[icell].field,&W,sizeof(struct Wtype));

	      }
	    }
	  }while(nextoct!=NULL);
	}
      }
    }
  printf("deltamax=%e\n",deltamax);
  if(cpu.rank==0) printf("Hydro done in %e (%e in hydro) sec\n",t100-t0,t80-t20);

#endif

#ifdef WMPI
	MPI_Barrier(cpu.comm);
#endif

#ifdef PIC
  // ==================================== velocity update   // Predictor step
  if(cpu.rank==0) printf("Correcctor\n");
#ifdef TESTCOSMO
  faexp=1.0;
#endif

  // -- force
  
  for(level=levelcoarse;level<=levelmax;level++)
    {
      //PoissonForce(level,&param,firstoct,&cpu,stencil,stride);
      accelpart(level,firstoct,dtnew*0.5,&cpu,sendbuffer,recvbuffer);
    }
#endif

  // ==================================== Moving Particles + Oct management
#ifdef PIC
  
#ifndef PARTN
  if(cpu.rank==0) printf("Moving particles\n");
    // Computing displacement (predictor)

#ifdef TESTCOSMO
  faexp2=1.0;
#endif
    
  movepart(levelcoarse,levelmax,firstoct,dtnew,&cpu);
  
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


#if 1
  // ==================================== Check the number of particles and octs
  multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);

#endif	 

#ifdef EGYCSV
  // ==================================== Energy Conservation Test
  
  REAL egy;
  egy=egypart(levelcoarse,levelmax,firstoct,&cpu);
  
#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&egy,1,MPI_REAL,MPI_SUM,cpu.comm);
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

    if(cpu.rank==0) fclose(fi);




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
      
