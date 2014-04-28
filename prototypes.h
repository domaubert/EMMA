typedef double REAL;

#ifdef WMPI
#include <mpi.h>
#endif


#define GAMMA (5./3.)
#define CFL (0.85)
#define GRAV (0.25)

#ifndef WHYDRO2 
#define OMEGAB (0.0)
#else
#define OMEGAB (0.049);
#define PMIN 1e-11
#endif

#define NCOSMOTAB (262144)


#ifdef DUAL_E

#ifndef WRADHYD
#define NVAR (6)
#else
#define NVAR (7)
#endif
#else

#define NVAR (5)
#endif



#define NFLUX (6*NVAR)


#ifdef WRAD
#define NVAR_R (5)
#define NGRP (1)
#define EMIN (0.)
#define NFLUX_R (6*NGRP*NVAR_R)
#endif


#define PI 3.14159265
// ================= PHYSICAL CONSTANTS ===================
#define LIGHT_SPEED_IN_M_PER_S (299792458.)
#define KBOLTZ (1.3806e-23) // J/K
#define PARSEC (3.085677e16) // in m
#define AVOGADRO (6.02214129e23) // mol-1
#define MYR (3.1536e13) // s
#define PROTON_MASS (1.67262158e-27) //kg
#define NEWTON_G (6.67384e-11) // SI
#define HELIUM_MASSFRACTION (0.24)
#ifdef WRADTEST
#define MOLECULAR_MU (1.0)
#else
#define MOLECULAR_MU (0.59)
#endif
//=======================================
#ifdef TESTCOSMO
struct COSMOPARAM{
  REAL aexp;
  REAL om;
  REAL ov;
  REAL ob;
  REAL H0;
  REAL *tab_aexp;
  REAL *tab_ttilde;
  REAL tsim;
  REAL unit_l;
  REAL tphy;
};
#endif

#ifdef STARS 
struct STARSPARAM{
  REAL overdensity_cond;// need overdensity_cond times the mean density to begin star formation
  REAL density_cond;	// Hydrogen density (m-3)
  REAL tcar;		// caracteristic time (yr)
  REAL tlife;		// life time of a radiative source (yr) 
  REAL feedback_eff;	// feedback efficiency
  REAL feedback_frac;	// fraction of thermal feedback over kinetic feedback
	
  REAL mstars;
  int  n;
  REAL thresh;

};
#endif


struct UNITS{
  REAL unit_l; // comoving length size of the box [meters]
  REAL unit_v; // unit velocity
  REAL unit_t; // unit time [seconds]
  REAL unit_n; // unit number [moles typically]
  REAL unit_mass; // unit mass [in kg, total mass is equal to one in unit codes]
};

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
  
  int gstride; // the size of the stencil for vector based computations (gravity)
  int hstride; // the size of the stencil for vector based computations (hydro)

  REAL dt; // the timsestep
  REAL tmax; // the simulation stops at tmax : corresponds to amax in cosmo
  REAL time_max; // for cosmo only : contains the time equivalent to amax (contained in tmax, yeah its obfuscated)

  int maxhash; // the hash table size between hilbert keys and oct adress (should be typically = to (2^levelmax-1)^3
  
  REAL amrthresh; // the refinement criterion (refine if mcell>amrthresh)
  int nsmooth; // the number of neighbour refinement steps

  REAL poissonacc; // relaxation accuracy for Poisson equation
  int mgridlmin;    // coarsest level for multigrid relaxation
  int nvcycles; // number of vcycles for multigrid relaxation
  int nrelax; // number of smoothing cycles

  int nrestart; // the restart snapshot
  int nsubcycles; // number of subcyles in AMR advance procedure

#ifdef TESTCOSMO
  struct COSMOPARAM *cosmo;
#endif

#ifdef STARS
  struct STARSPARAM *stars;
#endif

  int nthread;
  int nstream;

#ifdef WRAD
  REAL clight; // speed of light in units of the real one
  struct UNITS unit; // contains the units
  REAL fudgecool; // cooling fraction
  int ncvgcool; // cooling max iterations
  
  REAL denthresh; // density threshold to turn the sources on
  REAL tmpthresh; // temperature threshold to turn the sources on
  REAL srcint; // intensity of the sources
#endif

  REAL egy_rhs; // the right hand side of the energy conservation equation (0 in non cosmological case);
  REAL egy_0; // the initial energy
  REAL egy_last; // the last integrand for the energy equation (used for trapezoidal rule)
  REAL egy_timelast; // the last time for the integrand (used for trapezoidal rule)
  REAL egy_totlast; 
  FILE *fpegy; // the file with egy stats

};



//=======================================

// this structure exists for MPI communication protocol

struct PACKET{
  REAL data[8]; // the data to be transfered (8 since we transmit data per octs)
  unsigned long long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};





//=======================================


struct CPUINFO{
  int rank;
  int nproc;

  unsigned long long kmin;
  unsigned long long kmax;
  int nkeys;

  REAL load;

  struct PART *part; // the particle array
  struct OCT **bndoct; // the list of external boundary octs

  int nebnd; // the number of external boundary octs
  int nnei; // the number of neighbors procs
  
  int *mpinei; // the ranks of the neighbors procs

  int *dict; // a hash table to converts ranks in local neighbour index

  struct OCT **htable; // the hashing table to recover the octs from hilbert keys

  int *allkmin;
  int *allkmax;

  int nbuff; // the number of buffer cells = to the max of # of buffer cell from 1 neighbor 
  int nbufforg; // the max number of buffer cells (set from the parameter file)
  int nbuffpart;  // the number of particles to transmit

  int *nrecv; // the number of octs to be received by the local cpu, e.g cpu->nrecv[5] = nb of octs to be received from neigh # 5
  int *nsend; // the number of octs to be sent     by the local cpu, e.g cpu->nrecv[5] = nb of octs to be sent to       neigh # 5

  int *nrecv_coarse; // the number of l-1 octs to be received by the local cpu, e.g cpu->nrecv[5] = nb of l-1 octs to be received from neigh # 5
  int *nsend_coarse; // the number of l-1 octs to be sent     by the local cpu, e.g cpu->nrecv[5] = nb of l-1 octs to be sent to       neigh # 5

  struct PART *lastpart;
  struct PART *freepart;
  struct PART *firstpart;

#ifdef WMPI
  MPI_Datatype *MPI_PACKET; // the structured type for MPI messages (fields)
  struct PACKET **sendbuffer;
  struct PACKET **recvbuffer;
  struct PART_MPI **psendbuffer; 
  struct PART_MPI **precvbuffer;
  struct HYDRO_MPI **hsendbuffer; 
  struct HYDRO_MPI **hrecvbuffer;
  struct RAD_MPI **Rsendbuffer; 
  struct RAD_MPI **Rrecvbuffer;

#ifdef PIC
  MPI_Datatype *MPI_PART; // the structured type for MPI messages (particles)
#endif
 
#ifdef WHYDRO2
  MPI_Datatype *MPI_HYDRO; // the structured type for MPI messages (particles)
  MPI_Datatype *MPI_FLUX; // the structured type for MPI messages (particles)
#endif

#ifdef WRAD
  MPI_Datatype *MPI_RAD; // the structured type for MPI messages (particles)
#endif

  MPI_Comm comm; // the communicator
  
  
#endif

  int maxhash; // the size of the hashtable between hilbert keys and oct adresses
  int *noct; // the number of octs per levels
  int *npart; // the number of particles per levels
#ifdef STARS
  int *nstar;// the number of stars per levels
#endif

  int levelcoarse; // the levelcoarse

  struct OCT *freeoct; // the location of the first free oct
  struct OCT **firstoct; // the location of the first free oct
  int nsteps; // the current coarse step index
  REAL tinit; // the initial time
  int *ndumps; // the current dump number

#ifdef GPUAXL

#ifdef WGRAV
  
#ifndef FASTGRAV
  struct GGRID *dev_stencil;
  REAL *res;
  REAL *pnew;
  REAL *resLR;
#else
  struct STENGRAV *dev_stencil;
  struct STENGRAV *gpu_stencil;
#endif

#endif

#ifdef WHYDRO2
  struct HGRID *hyd_stencil;
#endif

#ifdef WRAD
  struct RGRID *rad_stencil;
#endif


  int nstream;
  int nthread;

  struct RUNPARAMS *dparam;


  REAL * gresA;
  REAL * gresB;
  unsigned long long int cuparam;
#endif
};





// =================== Local types for hydro calculations
// W stands for primitive quantities
// U stands for conservative quantities
//

#ifdef WRAD
struct Rtype{
  REAL e[NGRP];
  REAL fx[NGRP];
  REAL fy[NGRP];
  REAL fz[NGRP];
  REAL src;

#ifdef WCHEM
  REAL xion;
  REAL eint;
  REAL nh;
  REAL temp; // is a direct function of eint, nh and xion but stored for conveniency
#endif

};
#endif


struct Wtype{
  REAL d;   // density
  REAL u;   // velocity
  REAL v;   // velocity
  REAL w;   // velocity
  REAL p;   // pressure
  REAL a;   // sound speed
  REAL E;

#ifdef WRADHYD
  REAL X;
#endif
};


struct Wtype_MPI{
  REAL d;   // density
  REAL u;   // velocity
  REAL v;   // velocity
  REAL w;   // velocity
  REAL p;   // pressure
  REAL a;   // sound speed
  REAL E; 
#ifdef WRADHYD
  REAL X;
#endif
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

#ifdef WRADHYD
  REAL dX;
#endif

};


struct Wtype1D{
  REAL d;   // density
  REAL u;   // velocity
  REAL p;   // pressure
  REAL a;   // sound speed
#ifdef WRADHYD
  REAL X;
#endif
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
#ifdef WRADHYD
  REAL dX;
#endif
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

  REAL epot;
  REAL ekin;

#ifdef STARS
  int  isStar;
  REAL rhocell;
  REAL age;
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

#ifdef STARS
  REAL rhocell;
  REAL age;
#endif


  int idx;
  int level; // the level of the destination (to remove the key degeneracy)
  int icell; // the cell of destination
  int is;    // current step of particle


#ifdef STARS
  int isStar;
#endif


  unsigned long long key; // the destination hilbert key


};

//=======================================


// this structure is for the communication of Hydro data
struct HYDRO_MPI{
  struct Wtype data[8]; // the data to be transfered (8 since we transmit data per octs)
  unsigned long long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};

#ifdef WRAD
struct RAD_MPI{
  struct Rtype data[8]; // the data to be transfered (8 since we transmit data per octs)
  unsigned long long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};
#endif

// this structure is for the communication of Flux data
struct FLUX_MPI{
  REAL data[8*NFLUX]; // the data to be transfered (8 since we transmit data per octs)
  unsigned long long key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};


//=========================================================

struct Gtype{
  REAL d; //density
  REAL p; //pottential
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
  //REAL temp;
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


#ifdef WRAD
  struct Rtype rfield;
  struct Rtype rfieldnew;
#endif
};


// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CELLFLUX
{

#ifdef WHYDRO2
  //struct Wtype field;
  struct Utype deltaU;
  REAL flux[NFLUX]; // 6 fluxes of 5 variables each
  REAL divu;
#endif

#ifdef WRAD
  struct Rtype deltaR;
  struct Rtype rfield;
  struct Rtype rfieldnew;
  REAL rflux[NFLUX_R]; 
#endif

};

struct CELLFLUX_R
{

#ifdef WRAD
  struct Rtype deltaR;
  struct Rtype rfield;
  struct Rtype rfieldnew;
  REAL rflux[NFLUX_R]; 
#endif

};

struct CELLFLUX_H
{

#ifdef WHYDRO2
  //struct Wtype field;
  struct Utype deltaU;
  REAL flux[NFLUX]; // 6 fluxes of 5 variables each
  REAL divu;
#endif

};




// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CELLLIGHT
{

#ifdef WHYDRO2
  struct Wtype field; // hydrodynamical data
#endif

#ifdef WRAD
  struct Rtype rfield; // radiation data
#endif

#ifdef WGRAV
  REAL f[3];
#endif
  char split;
};

// =============================================================

struct CELLLIGHT_H
{

#ifdef WHYDRO2
  struct Wtype field; // hydrodynamical data
#endif

#ifdef WGRAV
  REAL f[3];
#endif
  char split;
};

// =========================================================

struct CELLLIGHT_R
{

#ifdef WRAD
  struct Rtype rfield; // radiation data
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
  /* int vecpos; */
  /* int border; // equal to zero if not a border */


};


// ========================================
struct OCTLIGHT_H
{
  // the cell properties
  struct CELLLIGHT_H cell[8]; // MUSTN'T BE MOVED !!
};


// ========================================
struct OCTLIGHT_R
{
  // the cell properties
  struct CELLLIGHT_R cell[8]; // MUSTN'T BE MOVED !!
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
  REAL *res2; 
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

struct OCTFLUX_H
{
  // the cell properties
  struct CELLFLUX_H cell[8]; // MUSTN'T BE MOVED !!
};

struct OCTFLUX_R
{
  // the cell properties
  struct CELLFLUX_R cell[8]; // MUSTN'T BE MOVED !!
};

// ========================================
struct HGRID{
  struct OCTLIGHT_H oct[27];
  struct OCTFLUX_H New;
};

// ============================================================================
// ============================================================================



// ========================================
struct RGRID{
  struct OCTLIGHT_R oct[7];
  struct OCTFLUX_R New;
};





// ========================================
struct MULTIVECT{
  REAL *vecpot; //contains the potential in "stride" octs
  REAL *vecpotnew; //contains the potential in "stride" octs
  REAL *vecden; //contains the density in "stride" octs


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

