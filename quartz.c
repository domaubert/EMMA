
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

#ifdef WRAD
#include "rad_utils.h"
#ifdef WCHEM
#include "chem_utils.h"
#include "atomic_data/Atomic.h"
#endif

#endif

#include "advanceamr.h"

#ifdef GPUAXL
#include "poisson_utils_gpu.h"
#include "hydro_utils_gpu.h"
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
  struct COSMOPARAM cosmo;
  param.cosmo=&cosmo;
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

#ifdef RIEMANN_HLLC
  rtag=3;
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

#ifdef TRANSXM
  printf("ERROR FIX GATHER FOR TRANSMISSIVE BNDARY\n");
  abort();
#endif

#ifdef TRANSXP
  printf("ERROR FIX GATHER FOR TRANSMISSIVE BNDARY\n");
  abort();
#endif

#ifdef TRANSYM
  printf("ERROR FIX GATHER FOR TRANSMISSIVE BNDARY\n");
  abort();
#endif

#ifdef TRANSYP
  printf("ERROR FIX GATHER FOR TRANSMISSIVE BNDARY\n");
  abort();
#endif

#ifdef TRANSZM
  printf("ERROR FIX GATHER FOR TRANSMISSIVE BNDARY\n");
  abort();
#endif

#ifdef TRANSZP
  printf("ERROR FIX GATHER FOR TRANSMISSIVE BNDARY\n");
  abort();
#endif


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
  npart=5;
#else
  npart=128*128*128;
#endif

#ifdef PARTN
  npart=32768;
#endif
#endif

  threshold=param.amrthresh;
  printf("DDD %d\n",param.stride);
  stride=fmax(8,param.stride);//pow(2,levelcoarse);
  ncomp=10;
  acc=param.poissonacc;
  dt=param.dt;
  cpu.maxhash=param.maxhash;
  cpu.levelcoarse=levelcoarse;

#ifdef GPUAXL
  cpu.nthread=param.nthread;
  cpu.nstream=param.nstream;
#endif
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
  cpu.noct =(int *)calloc(levelmax,sizeof(int)); memsize+=levelmax*sizeof(int);// the number of octs per level
  cpu.npart=(int *)calloc(levelmax,sizeof(int)); memsize+=levelmax*sizeof(int);// the number of particles per level

  lastpart=part-1; // the last particle points before the first at the very beginning



  //===================================================================================================




  // allocating the vectorized tree
  

  
  // allocating the 6dim stencil
  struct HGRID *stencil;
  struct STENGRAV gstencil;
  struct RGRID *rstencil;

#ifndef GPUAXL
  printf("stencil=%p with stride=%d\n",stencil,stride);
  stencil=(struct HGRID*)calloc(stride,sizeof(struct HGRID));
  printf("stenci=%p mem=%f\n",stencil,stride*sizeof(struct HGRID)/(1024.*1024.));

#ifdef WRAD
  rstencil=(struct RGRID*)calloc(stride,sizeof(struct RGRID));
#endif

#ifndef FASTGRAV
  struct GGRID *grav_stencil;
  grav_stencil=(struct GGRID*)calloc(stride,sizeof(struct GGRID));
  gstencil.stencil=grav_stencil;
  gstencil.res=(REAL *)calloc(stride*8,sizeof(REAL));
  gstencil.pnew=(REAL *)calloc(stride*8,sizeof(REAL));
  gstencil.resLR=(REAL *)calloc(stride,sizeof(REAL));
#else
  gstencil.d=(REAL *)calloc(stride*8,sizeof(REAL));
  gstencil.p=(REAL *)calloc(stride*8,sizeof(REAL));
  gstencil.pnew=(REAL *)calloc(stride*8,sizeof(REAL));

  gstencil.nei=(int *)calloc(stride*7,sizeof(int));
  gstencil.level=(int *)calloc(stride,sizeof(int));
  gstencil.cpu=(int *)calloc(stride,sizeof(int));
  gstencil.valid=(char *)calloc(stride,sizeof(char));
  
  gstencil.res=(REAL *)calloc(stride*8,sizeof(REAL));
  gstencil.resLR=(REAL *)calloc(stride,sizeof(REAL));
#endif
  
#endif

  

#ifdef GPUAXL
  // ================================== GPU ALLOCATIONS ===============
  countdevices(0);
  initlocaldevice(0,1);
  checkdevice(0);
#ifdef WGRAV
  create_pinned_gravstencil(&gstencil,stride);
#ifdef FASTGRAV
  struct STENGRAV dev_stencil;
  cpu.dev_stencil=&dev_stencil;
#endif
  create_gravstencil_GPU(&cpu,stride);
  printf("stencil created on device with adress %p\n",cpu.dev_stencil);
#endif

#ifdef WHYDRO2 

  printf("stencil=%p with stride=%d\n",stencil,stride);
  stencil=(struct HGRID*)calloc(stride,sizeof(struct HGRID));
  printf("stenci=%p mem=%f\n",stencil,stride*sizeof(struct HGRID)/(1024.*1024.));
  
  // UNCOMMENT BELOW FOR FASTHYDRO GPU
  create_pinned_stencil(&stencil,stride);  
  create_hydstencil_GPU(&cpu,stride); 
#endif 
  // ====================END GPU ALLOCATIONS ===============
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
    /* grid->cell[icell].density=0.; */
    /* grid->cell[icell].pot=0.; */
    /* grid->cell[icell].temp=0.; */
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
	      /* newoct->cell[ii].density=0.; */
	      newoct->cell[ii].idx=ii;
#ifdef PIC
	      newoct->cell[ii].phead=NULL;
#endif

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

  int ir,nr=5;
  ip=0;
  REAL dxcell=1./pow(2.,levelcoarse);
  REAL epsilon=0.;
  REAL r0=0.12;
  REAL vf=0.8;
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

      x=0.5+r0;
      y=0.5;
      z=0.5;

      vx=0.;
      vy=sqrt((1.-epsilon)/r0)*1.0; // this one is circular
      vz=0.;
      
      mass=epsilon;
    }
    else if(ir==2){

      x=0.5;
      y=0.5+r0*0.3;
      z=0.5;

      vy=0.;
      vx=-sqrt((1.-epsilon)/(r0*0.3))*1.0;//this one is circular
      vz=0.;
      
      mass=epsilon;
    }
    else if(ir==3){

      x=0.5-r0;
      y=0.5;
      z=0.5;

      vx=0.;
      vy=-sqrt((1.-epsilon)/r0)*vf;
      vz=0.;
      
      mass=epsilon;
    }
    else if(ir==4){

      x=0.5;
      y=0.5-r0;
      z=0.5;

      vy=0.;
      vx=sqrt((1.-epsilon)/r0)*vf;
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
      //r=(REAL)(rand())/RAND_MAX*0.3;
      r=0.12;

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

#ifdef TESTCOSMO // =================PARTICLE COSMOLOGICAL CASE

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

  lastpart=read_grafic_part(part, &cpu, &munit, &ainit, &npart, &param);
  amax=1.;
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
  ncellhydro=read_grafic_hydro(&cpu,&ainit, &param);

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
  struct Utype UL, UR;
  REAL X0;
  if(cpu.rank==0) printf("Init Hydro\n");

  /* /\*  /\\* // TEST 1 *\\/ *\/ */

  /* WL.d=1.; */
  /* WL.u=0.; */
  /* WL.v=0.; */
  /* WL.w=0.75; */
  /* WL.p=1.0; */

  /* WR.d=0.125; */
  /* WR.u=0.; */
  /* WR.v=0.; */
  /* WR.w=0.; */
  /* WR.p=0.1; */
 
  /* X0=0.3125; */
  /* tmax=0.15; */


  /*  /\* // TEST 123 *\/ */

  /* WL.d=1.; */
  /* WL.u=0.; */
  /* WL.v=0.; */
  /* WL.w=-2.0; */
  /* WL.p=0.4; */

  /* WR.d=1.; */
  /* WR.u=0.; */
  /* WR.v=0.; */
  /* WR.w=2.0; */
  /* WR.p=0.4; */
  /* X0=0.5; */
  /* tmax=0.15; */

  /*  /\* // TEST 3 *\/ */

  WL.d=1.;
  WL.u=0.;
  WL.v=0.;
  WL.w=0.;
  WL.p=1000.;

  WR.d=1.;
  WR.u=0.;
  WR.v=0.;
  WR.w=0.;
  WR.p=0.01;
  X0=0.5;
  tmax=0.012;


  /*  /\* // TEST 5 *\/ */

  /* WL.d=1.; */
  /* WL.u=-19.59745; */
  /* WL.v=0.; */
  /* WL.w=0.; */
  /* WL.p=1000.; */

  /* WR.d=1.; */
  /* WR.u=-19.59745; */
  /* WR.v=0.; */
  /* WR.w=0.; */
  /* WR.p=0.01; */
  /* X0=0.8; */
  /* tmax=0.012; */

  /*  /\* // TEST 4 *\/ */

  /* WR.d=5.99924; */
  /* WR.v=-19.5975; */
  /* WR.u=0.; */
  /* WR.w=0.; */
  /* WR.p=460.894; */

  /* WL.d=5.99242; */
  /* WL.v=6.19633; */
  /* WL.u=0.; */
  /* WL.w=0.; */
  /* WL.p=46.0950; */
  /* X0=0.6; */
  /* tmax=0.035; */

  /*  /\* // REVERSED TEST 1 *\/ */

  /* WL.d=0.125; */
  /* WL.u=0.; */
  /* WL.v=0.; */
  /* WL.w=0.; */
  /* WL.p=0.1; */

  /* WR.d=1.; */
  /* WR.u=-0.75; */
  /* WR.v=0.; */
  /* WR.w=0.; */
  /* WR.p=1.; */
  /* X0=0.7; */

  WL.a=sqrt(GAMMA*WL.p/WL.d);
  WR.a=sqrt(GAMMA*WR.p/WR.d);

  // ======================================================

  REAL dtot=0.;
  int nc=0;
  REAL dmax=0.;
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

	      /* curoct->cell[icell].pot=GRAV*zc; */

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
	      /* sZEL	curoct->cell[icell].field.d=2.; */
	      /* } */
	      /* else{ */
	      /* 	curoct->cell[icell].field.d=1.; */
		
	      /* 	} */
	      	
	      /* curoct->cell[icell].field.p=2.5-curoct->cell[icell].field.d*GRAV*(zc-0.9); */
	      /* curoct->cell[icell].field.a=sqrt(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d); */


	      /* /\* EVRARD ADIABATIC COLLAPSE *\/ */
	      
	      /* UR.d=1./sqrt((yc-0.5)*(yc-0.5)+(xc-0.5)*(xc-0.5)+(zc-0.5)*(zc-0.5)); */
	      /* UR.du=0.; */
	      /* UR.dv=0.; */
	      /* UR.dw=0.; */
	      /* UR.E=0.05; */
	      
	      /* U2W(&UR,&WR); */


	      /* /\* WR.d=1./sqrt((yc-0.5)*(yc-0.5)); *\/ */
	      /* /\* WR.u=0.; *\/ */
	      /* /\* WR.v=0.; *\/ */
	      /* /\* WR.w=0.; *\/ */
	      /* /\* WR.p=1e-1; *\/ */
	      /* /\* WR.a=sqrt(GAMMA*WR.p/WR.d); *\/ */


	      /* memcpy(&(curoct->cell[icell].field),&WR,sizeof(struct Wtype));  */
	      
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
	      


	      /* /\* /\\* ZELDOVICH PANCAKE *\\/ *\/ */

	      /* REAL ZI=100; */
	      /* REAL ZC=9; */
	      /* cosmo.om=1.0; */
	      /* cosmo.ov=0.; */
	      /* cosmo.ob=OMEGAB; */

	      /* ainit=1./(1.+ZI); */
	      /* amax=1.;//(1.+ZC); */
	      /* curoct->cell[icell].field.d=1.+(1.+ZC)/(1.+ZI)*cos(2.*M_PI*(xc-0.5)); */
	      /* curoct->cell[icell].field.u=0.-(1.+ZC)/pow(1.+ZI,1.5)*sin(2.*M_PI*(xc-0.5))/(2.*M_PI); */
	      /* curoct->cell[icell].field.v=0.;//-(1.+ZC)/pow(1.+ZI,1.5)*sin(2.*M_PI*(yc-0.5))/(2.*M_PI); */
	      /* curoct->cell[icell].field.w=0.;//-(1.+ZC)/pow(1.+ZI,1.5)*sin(2.*M_PI*(zc-0.5))/(2.*M_PI); */
	      /* curoct->cell[icell].field.p=1e-14; */
	      /* curoct->cell[icell].field.a=sqrt(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d); */
	      /* getE(&curoct->cell[icell].field); */



	      /* /\* SHOCK TUBE *\/ */
	      /* if(zc<=X0){ */

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
  printf("dmax=%e\n",dmax);
#endif

  /* sprintf(filename,"data/denhydstart.%05d.p%05d",0,cpu.rank);  */
  /* dumpcube(lmap,firstoct,101,filename,0.);  */
  /* sprintf(filename,"data/prehydstart.%05d.p%05d",0,cpu.rank);  */
  /* dumpcube(lmap,firstoct,105,filename,0.);  */




#endif

  //===================================================================================================================================
#ifdef WRAD

#ifdef WCHEM
  param.fudgecool=0.01;
  param.ncvgcool=0;
  if(NGRP!=NGRP_ATOMIC){
    printf("NGRP and NGRP_ATOMIC INCONSISTENT ! ERROR !\n");
    abort();
  }


#endif

#ifdef WRADTEST
  REAL X0=1./pow(2,levelcoarse);
  int igrp;

#ifndef TESTCOSMO
  param.unit.unit_l=13.2e3*PARSEC;
  param.unit.unit_v=LIGHT_SPEED_IN_M_PER_S;
  param.unit.unit_t=param.unit.unit_l/param.unit.unit_v;
  param.unit.unit_n=1.;
  ainit=1.;
#else
  ainit=1./(11.);;
  REAL om=1.;
  REAL ov=0.;
  REAL ob=0.045;
  REAL h0=0.73;
  REAL H0=h0*100*1e3/1e6/PARSEC;

  param.cosmo->om=om;
  param.cosmo->ov=ov;
  param.cosmo->ob=ob;
  param.cosmo->H0=h0*100.;

  REAL rstar= 13.2*1e6*PARSEC; // box size in m
  double rhoc=3.*H0*H0/(8.*M_PI*NEWTON_G); // comoving critical density (kg/m3)

  REAL rhostar=rhoc*om;
  REAL tstar=2./H0/sqrt(om); // sec
  REAL vstar=rstar/tstar; //m/s

  param.unit.unit_l=rstar;
  param.unit.unit_v=vstar;
  param.unit.unit_t=tstar;
  param.unit.unit_n=1.;
  param.unit.unit_mass=rhostar*pow(param.unit.unit_l,3);
  amax=1.;

#endif

  for(level=levelcoarse;level<=levelmax;level++) 
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
	      
	      for(igrp=0;igrp<NGRP;igrp++){
		curoct->cell[icell].rfield.e[igrp]=0.+EMIN; 
		if((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5)<(X0*X0)){ 
		  //curoct->cell[icell].rfield.src=5e8;
		  curoct->cell[icell].rfield.src=5e55/pow(X0,3)/8.*param.unit.unit_t/param.unit.unit_n*pow(ainit,2); 
		  printf("SRC=%e\n",curoct->cell[icell].rfield.src);
		}
		else{
		  curoct->cell[icell].rfield.src=0.; 
		}
		curoct->cell[icell].rfield.fx[igrp]=0.; 
		curoct->cell[icell].rfield.fy[igrp]=0.; 
		curoct->cell[icell].rfield.fz[igrp]=0.; 
#ifdef WCHEM
#ifndef COOLING
		curoct->cell[icell].rfield.temp=1e4; 
		curoct->cell[icell].rfield.xion=1.2e-3; 
#else
		curoct->cell[icell].rfield.temp=1e2; 
		curoct->cell[icell].rfield.xion=1e-5; 
#endif
		curoct->cell[icell].rfield.nh=0.2*pow(param.unit.unit_l,3)/param.unit.unit_n; 
#endif
	      }
	    }
	}while(nextoct!=NULL);
      
      
    }
#endif
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



#ifdef TESTCOSMO
  // ================== COMPUTATION OF FRIEDMANN TABLES

  tsim=ainit;
  aexp=ainit;

  // at this stage we have to compute the conformal time
  tsim=-0.5*sqrt(cosmo.om)*integ_da_dt_tilde(ainit,1.0,cosmo.om,cosmo.ov,1e-8);

  // we compute the friedmann tables

  compute_friedmann(ainit*0.95,NCOSMOTAB,cosmo.om,cosmo.ov,tab_aexp,tab_ttilde,tab_t);

  tmax=-0.5*sqrt(cosmo.om)*integ_da_dt_tilde(amax,1.0+1e-6,cosmo.om,cosmo.ov,1e-8);
  printf("tmax=%e\n",tmax);
  cosmo.tab_aexp=tab_aexp;
  cosmo.tab_ttilde=tab_ttilde;
#endif


  

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

#ifdef WHYDRO2
  tmax=5.;
#else
  tmax=1e15;
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
    
    // building the array of timesteps

    REAL *adt;
    adt=(REAL *)malloc(sizeof(REAL)*levelmax);
    for(level=1;level<=levelmax;level++) adt[level-1]=param.dt;

    int *ndt;
    ndt=(int *)malloc(sizeof(int)*levelmax);

    // preparing freeocts
    cpu.freeoct=freeoct;


    // Loop over time
    for(nsteps=0;(nsteps<=param.nsteps)*(tsim<=tmax);nsteps++){

      cpu.nsteps=nsteps;
      
#ifdef TESTCOSMO
      cosmo.aexp=interp_aexp(tsim,cosmo.tab_aexp,cosmo.tab_ttilde);
      cosmo.tsim=tsim;
      if(cpu.rank==0) printf("\n============== STEP %d aexp=%e z=%f t=%e tmax=%e================\n",nsteps,cosmo.aexp,1./cosmo.aexp-1.,tsim,tmax);
#else
#ifndef WRAD
      if(cpu.rank==0) printf("\n============== STEP %d tsim=%e ================\n",nsteps,tsim);
#else
      if(cpu.rank==0) printf("\n============== STEP %d tsim=%e [%e Myr] ================\n",nsteps,tsim,tsim*param.unit.unit_t/MYR);
#endif
#endif
      // Resetting the timesteps

      for(level=1;level<=levelmax;level++){
	//adt[level-1]=param.dt;
	ndt[level-1]=0;
      }

      //Recursive Calls over levels
      Advance_level(levelcoarse,adt,&cpu,&param,firstoct,lastoct,stencil,&gstencil,rstencil,stride,sendbuffer,recvbuffer,ndt,nsteps);
      

      // ==================================== dump
  
      if((nsteps%(param.ndumps)==0)||((tsim+dt)>=tmax)){

#ifndef TESTCOSMO
	REAL tdump=tsim+adt[levelcoarse-1];
#else
	REAL tdump=cosmo.aexp;
#endif
	// === Hydro dump
    
	//printf("tdum=%f\n",tdump);
	sprintf(filename,"data/grid.%05d.p%05d",ndumps,cpu.rank); 
	if(cpu.rank==0){
	  printf("Dumping .......");
	  printf("%s\n",filename);
	}
	dumpgrid(levelmax,firstoct,filename,tdump); 

#ifdef PIC
	sprintf(filename,"data/part.%05d.p%05d",ndumps,cpu.rank);
	if(cpu.rank==0){
	  printf("Dumping .......");
	  printf("%s\n",filename);
	}
	dumppart(firstoct,filename,npart,levelcoarse,levelmax,tdump);
#endif

	ndumps++;
      }

      //==================================== timestep completed, looping
      dt=adt[param.lcoarse-1];
      tsim+=dt;
    }
    

    // we are done let's free the ressources

    free(adt);
    free(ndt);


#ifdef GPUAXL

#ifdef WGRAV
    destroy_pinned_gravstencil(&gstencil,stride);
    destroy_gravstencil_GPU(&cpu,stride);
#endif

#ifdef WHYDRO2
    /* destroy_pinned_stencil(&stencil,stride); */
    /* destroy_hydstencil_GPU(&cpu,stride); */
#endif

#endif




    if(cpu.rank==0){
      printf("Done .....\n");
    }
    
    return 0;
}

