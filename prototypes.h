
#ifdef WMPI
#include <mpi.h>
#endif

#define GAMMA (1.4)
#define GRAV (0.25)

#ifdef SUPERCOMOV
#define NCOSMOTAB (262144)
#endif


//=======================================

struct RUNPARAMS{
  int npartmax; // the max particles number (per process)
  int ngridmax; // the max oct numbers (per process)
  int nbuff; // the mpi buffer size
  int ndumps; // the frequency of outputs
  int nsteps; // the maximal number of timesteps
  
  int lcoarse; // the coarse level
  int lmax; // the max level of refinement
  int levelmap; // the map projection level
  
  int niter; // the maximal number of iterations for the Poisson solver
  
  int stride; // the size of the stencil for vector based computations

  float dt; // the timsestep

  int maxhash; // the hash table size between hilbert keys and oct adress (should be typically = to (2^levelmax-1)^3
  
  float amrthresh; // the refinement criterion (refine if mcell>amrthresh)

  float poissonacc; // relaxation accuracy for Poisson equation
  int mgridlmin;    // coarsest level for multigrid relaxation
  int nvcycles; // number of vcycles for multigrid relaxation
  int nrelax; // number of smoothing cycles
};



//=======================================

// this structure exists for MPI communication protocol

struct PACKET{
  float data[8]; // the data to be transfered (8 since we transmit data per octs)
  long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};





//=======================================


struct CPUINFO{
  int rank;
  int nproc;

  int kmin;
  int kmax;
  int nkeys;

  float load;

  struct OCT **bndoct; // the list of external boundary octs

  int nebnd; // the number of external boundary octs
  int nnei; // the number of neighbors procs
  
  int *mpinei; // the ranks of the neighbors procs

  int *dict; // a hash table to converts ranks in local neighbour index

  struct OCT **htable; // the hashing table to recover the octs from hilbert keys

  int *allkmin;
  int *allkmax;

  int nbuff; // the max number of buffer cells

  int *nrecv;
  int *nsend;

#ifdef WMPI
  MPI_Datatype *MPI_PACKET; // the structured type for MPI messages (fields)
#ifdef PIC
  MPI_Datatype *MPI_PART; // the structured type for MPI messages (particles)
#endif
 
#ifdef WHYDRO2
  MPI_Datatype *MPI_HYDRO; // the structured type for MPI messages (particles)
  MPI_Datatype *MPI_FLUX; // the structured type for MPI messages (particles)
#endif

  MPI_Comm comm; // the communicator
#endif

  int maxhash; // the size of the hashtable between hilbert keys and oct adresses
  int *noct; // the number of octs per levels
  int levelcoarse; // the levelcoarse
};



// =================== Local types for hydro calculations
// W stands for primitive quantities
// U stands for conservative quantities
//


struct Wtype{
  float d;   // density
  float u;   // velocity
  float v;   // velocity
  float w;   // velocity
  float p;   // pressure
  float a;   // sound speed
};


struct Wtype_MPI{
  float d;   // density
  float u;   // velocity
  float v;   // velocity
  float w;   // velocity
  float p;   // pressure
  float a;   // sound speed
};


struct Utype{
  float d;    // density
  float du;   // momentum
  float dv;   // momentum
  float dw;   // momentum
  float E;    // Energy
};


struct Wtype1D{
  float d;   // density
  float u;   // velocity
  float p;   // pressure
  float a;   // sound speed
};

struct Utype1D{
  float d;    // density
  float du;   // momentum
  float E;    // Energy
};


//=======================================

struct PART
{
  float x;
  float y;
  float z;

  float vx;
  float vy;
  float vz;

#ifndef TESTCOSMO
  float fx;
  float fy;
  float fz;
#endif


  struct PART *next;
  struct PART *prev;

  float mass;

  int idx;

};

struct PART_MPI // For mpi communications
{
  float x;
  float y;
  float z;

  float vx;
  float vy;
  float vz;

  float mass;
  int idx;

  long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
  int icell; // the cell of destination
};

//=======================================



// this structure is for the communication of Hydro data
struct HYDRO_MPI{
  struct Wtype data[8]; // the data to be transfered (8 since we transmit data per octs)
  long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};


// this structure is for the communication of Flux data
struct FLUX_MPI{
  float data[8*30]; // the data to be transfered (8 since we transmit data per octs)
  long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};


//-------------------------------------
struct CELL
{
  struct OCT *child;
  float marked; // float for consistency with physical quantities during communications
  int idx; //index of the cell within the oct

  // the head particle
  struct PART * phead;

  // the physical quantities
  float density; // total density
  float pot;

  float temp;

#ifdef AXLFORCE
  float fx;
  float fy;
  float fz;
#endif

#ifdef WHYDRO
  float d; // gas density
  float u; // velx
  float v; // vely
  float w; // velz
  float p; // pressure
  float a; // sound speed
#endif

#ifdef WHYDRO2
  struct Wtype field;
  float flux[30]; // 6 fluxes of 5 variables each
#endif


};

// ----------------------------------------------------------------

struct CELLLIGHT
{
  struct OCT *child;
  float marked; // float for consistency with physical quantities during communications
  int idx; //index of the cell within the oct


#ifdef WHYDRO2
  struct Wtype field;
#endif


};



//-------------------------------------
//-------------------------------------------------------------------------
struct OCT
{
  // the cell properties
  struct CELL cell[8]; // MUSTN'T BE MOVED !!

  struct CELL *parent; // parent cell 
  struct CELL *nei[6];// neighbor cells at level - 1
  
  // the next two pointers are required for sweeps through a single level
  struct OCT *next; // next oct on the same level
  struct OCT *prev; // previous oct on the same level
  struct OCT *nexthash; // next oct in the hash list

  // the oct position (lowest left corner)
  float x;
  float y;
  float z;

  // parallel data
  int cpu; 
  int level;// level of the cells in the oct


  // vector info
  int vecpos;
  int border; // equal to zero if not a border

  
};


// ========================================
struct OCTLIGHT
{
  // the cell properties
  struct CELLLIGHT cell[8]; // MUSTN'T BE MOVED !!
};

// ========================================
struct OCTFLUX
{
  // the cell properties
  struct CELL cell[8]; // MUSTN'T BE MOVED !!
};

// ========================================
struct HGRID{
  struct OCTLIGHT oct[27];
  struct OCTFLUX new;
};

// ========================================
struct MULTIVECT{
  float *vecpot; //contains the potential in "stride" octs
  float *vecpotnew; //contains the potential in "stride" octs
  float *vecden; //contains the density   in "stride" octs

#ifdef WHYDRO
  float *vec_d; 
  float *vec_u; 
  float *vec_v; 
  float *vec_w;
  float *vec_p;

  float *vec_dnew; 
  float *vec_unew; 
  float *vec_vnew; 
  float *vec_wnew;
  float *vec_pnew;
#endif

  int *vecnei;//contains the cell neighbors of the octs
  int *vecl; // contains the level of the octs
  int *veccpu; // contains the level of the octs
  int *vecicoarse; // contains the level of the octs

#ifdef WGPU
  float *vecden_d;
  int *vecl_d;
  int *veccpu_d;
  float *vec2_d;
  float *vecsum_d;

  float *vecpot_d; 
  float *vecpotnew_d; 
  int *vecnei_d;
  int *vecicoarse_d; 
#endif
};

