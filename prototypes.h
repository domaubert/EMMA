

typedef double REAL;


#ifdef WMPI
#include <mpi.h>
#endif

#define GAMMA (5./3.)
#define CFL (0.4)
#define GRAV (0.25)
#define OMEGAB (0.045)

#define NCOSMOTAB (262144)


#ifdef DUAL_E
#define NVAR (6)
#else
#define NVAR (5)
#endif
#define NFLUX (6*NVAR)

//=======================================
#ifdef TESTCOSMO
struct COSMOPARAM{
  REAL aexp;
  REAL om;
  REAL ov;
  REAL ob;
  REAL *tab_aexp;
  REAL *tab_ttilde;
  REAL tsim;
};
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
  
  int niter; // the maximal number of iterations for the Poisson solver
  
  int stride; // the size of the stencil for vector based computations

  REAL dt; // the timsestep

  int maxhash; // the hash table size between hilbert keys and oct adress (should be typically = to (2^levelmax-1)^3
  
  REAL amrthresh; // the refinement criterion (refine if mcell>amrthresh)

  REAL poissonacc; // relaxation accuracy for Poisson equation
  int mgridlmin;    // coarsest level for multigrid relaxation
  int nvcycles; // number of vcycles for multigrid relaxation
  int nrelax; // number of smoothing cycles

  int nrestart; // the restart snapshot
  int nsubcycles; // number of subcyles in AMR advance procedure

#ifdef TESTCOSMO
  struct COSMOPARAM *cosmo;
#endif

  int nthread;
  int nstream;

};



//=======================================

// this structure exists for MPI communication protocol

struct PACKET{
  REAL data[8]; // the data to be transfered (8 since we transmit data per octs)
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

  REAL load;

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
  int *npart; // the number of particles per levels
  int levelcoarse; // the levelcoarse

  struct OCT *freeoct; // the location of the first free oct
  int nsteps; // the current coarse step index

#ifdef GPUAXL
  struct GGRID *dev_stencil;
  REAL *res;
  REAL *pnew;
  REAL *resLR;
  struct HGRID *hyd_stencil;

  int nstream;
  int nthread;
#endif

};





// =================== Local types for hydro calculations
// W stands for primitive quantities
// U stands for conservative quantities
//


struct Wtype{
  REAL d;   // density
  REAL u;   // velocity
  REAL v;   // velocity
  REAL w;   // velocity
  REAL p;   // pressure
  REAL a;   // sound speed
  REAL E;

};


struct Wtype_MPI{
  REAL d;   // density
  REAL u;   // velocity
  REAL v;   // velocity
  REAL w;   // velocity
  REAL p;   // pressure
  REAL a;   // sound speed
};


struct Utype{
  REAL d;    // density
  REAL du;   // momentum
  REAL dv;   // momentum
  REAL dw;   // momentum
  REAL E;    // Energy

#ifdef DUAL_E
  REAL eint; // internal energy
#endif

};


struct Wtype1D{
  REAL d;   // density
  REAL u;   // velocity
  REAL p;   // pressure
  REAL a;   // sound speed
};

struct Wtype1D_double{
  double d;   // density
  double u;   // velocity
  double p;   // pressure
  double a;   // sound speed
};


struct Utype1D{
  REAL d;    // density
  REAL du;   // momentum
  REAL E;    // Energy
};


//=======================================

struct PART
{
  REAL x;
  REAL y;
  REAL z;

  REAL vx;
  REAL vy;
  REAL vz;

#ifndef TESTCOSMO
  REAL fx;
  REAL fy;
  REAL fz;
#endif


  struct PART *next;
  struct PART *prev;

  REAL mass;

  int idx;
  int level;
  int is; // local timestep number per particle

#ifdef PART_EGY
  REAL epot;
  REAL ekin;
#endif

};

struct PART_MPI // For mpi communications
{
  REAL x;
  REAL y;
  REAL z;

  REAL vx;
  REAL vy;
  REAL vz;

  REAL mass;
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
  REAL data[8*NFLUX]; // the data to be transfered (8 since we transmit data per octs)
  long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};


//=========================================================

struct Gtype{
  REAL d;
  REAL p;
};

//-------------------------------------
struct CELL
{
  struct OCT *child;
  REAL marked; // REAL for consistency with physical quantities during communications
  int idx; //index of the cell within the oct

#ifdef PIC
  // the head particle
  struct PART * phead;

  // the physical quantities
  REAL density; // total density
  REAL temp;
#endif

#ifdef WGRAV
  struct Gtype gdata;
  REAL pnew; // new potential
  REAL res; // residual
  REAL f[3]; // the gravitational force component
#endif


#ifdef WHYDRO2
  struct Wtype field;
  struct Wtype fieldnew;
#endif
};


// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CELLFLUX
{

#ifdef WHYDRO2
  struct Wtype field;
  struct Utype deltaU;
  REAL flux[NFLUX]; // 6 fluxes of 5 variables each
  REAL divu;
#endif

};



// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CELLLIGHT
{
#ifdef WGRAV
  REAL f[3];
#endif

#ifdef WHYDRO2
  struct Wtype field; // hydrodynamical data
#endif
  char split;
};

// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CELLGRAV
{
  struct Gtype gdata; // gravitational data 
};

// ----------------------------------------------------------------
// ----------------------------------------------------------------



//-------------------------------------
//-------------------------------------------------------------------------
struct OCT
{
  // the cell properties
  struct CELL cell[8]; // MUSTN'T BE MOVED !!

  struct CELL *nei[6];// neighbor cells at level - 1
  struct CELL *parent; // parent cell 
 
  // the next two pointers are required for sweeps through a single level
  struct OCT *next; // next oct on the same level
  struct OCT *prev; // previous oct on the same level
  struct OCT *nexthash; // next oct in the hash list

  // the oct position (lowest left corner)
  REAL x;
  REAL y;
  REAL z;

  // parallel data
  int cpu; 
  int level;// level of the cells in the oct


  // ***************** CAN BE DELETED *****************//
  // vector info
  int vecpos;
  int border; // equal to zero if not a border
  // ***************** CAN BE DELETED *****************//

  int stenpos;

};


// ========================================
struct OCTLIGHT
{
  // the cell properties
  struct CELLLIGHT cell[8]; // MUSTN'T BE MOVED !!
};

// ========================================
struct OCTGRAV
{
  // the cell properties
  struct CELLGRAV cell[8]; 
};

// =======================================
#ifndef FASTGRAV
struct GGRID{
  struct OCTGRAV oct[27];
};

// =======================================
struct STENGRAV{
  struct GGRID *stencil;
  REAL *res;
  REAL *pnew;
  REAL *resLR;
};

#else
struct GGRID{
  int nei[27]; //pointers toward other neighbour octs in the stencil
  struct OCTGRAV oct; // the local one
};

// =======================================
struct STENGRAV{
  REAL *d; // density [8*stride]
  REAL *p; // potential [8*stride]
  REAL *pnew; // new potential [8*stride]

  REAL *res; //residual [8*stride]
  REAL *resLR; // low res residual for MG [stride]

  int *nei; // neighbour indexes [7*stride] (6 real neighbours plus the middle one)
  int *level; // oct levels [stride]
  int *cpu; // cpu of octs for MPI boundaries [stride]
  char *valid; // validity of oct (border=0 or inner=1);
};
#endif






// ========================================
struct OCTFLUX
{
  // the cell properties
  struct CELLFLUX cell[8]; // MUSTN'T BE MOVED !!
};

// ========================================
struct HGRID{
  struct OCTLIGHT oct[27];
  struct OCTFLUX New;
};

// ============================================================================
// ============================================================================




// ========================================
struct MULTIVECT{
  REAL *vecpot; //contains the potential in "stride" octs
  REAL *vecpotnew; //contains the potential in "stride" octs
  REAL *vecden; //contains the density   in "stride" octs


  int *vecnei;//contains the cell neighbors of the octs
  int *vecl; // contains the level of the octs
  int *veccpu; // contains the level of the octs
  int *vecicoarse; // contains the level of the octs

#ifdef WGPU
  REAL *vecden_d;
  int *vecl_d;
  int *veccpu_d;
  REAL *vec2_d;
  REAL *vecsum_d;

  REAL *vecpot_d; 
  REAL *vecpotnew_d; 
  int *vecnei_d;
  int *vecicoarse_d; 
#endif
};



//==============================================
