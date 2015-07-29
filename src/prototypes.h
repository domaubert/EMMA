/**
  * \file prototypes.h
  *
  *
  */

#include <stdio.h>
#include <math.h>

#ifdef WMPI
#include <mpi.h>
#endif

#ifdef GSLRAND
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

#include "param.h"



#ifdef SINGLEPRECISION
// SINGLE PRECISION CASE

typedef float REAL;
#ifdef WMPI
#define MPI_REEL MPI_FLOAT
#endif

#define POW(A,B) powf(A,B)
#define SQRT(A) sqrtf(A)
#define EXP(A) expf(A)
#define FMIN(A,B) fminf(A,B)
#define FMAX(A,B) fmaxf(A,B)
#define FABS(A) fabsf(A)
#define CUDPP_REAL CUDPP_FLOAT

#else
// DOUBLE PRECISION CASE (BY DEFAULT)

typedef double REAL;
#ifdef WMPI
#define MPI_REEL MPI_DOUBLE
#endif

#define POW(A,B) pow(A,B)
#define SQRT(A) sqrt(A)
#define EXP(A) exp(A)
#define FMIN(A,B) fmin(A,B)
#define FMAX(A,B) fmax(A,B)
#define FABS(A) fabs(A)
#define CUDPP_REAL CUDPP_DOUBLE

#endif // SINGLEPRECISION

#ifdef DUAL_E
  #ifndef WRADHYD
    #define NVAR (6)
  #else
#ifdef HELIUM
#define NVAR (9)
#else
#define NVAR (7)
#endif
#endif
#else
#define NVAR (5)
#endif

#define NFLUX (6*NVAR)

#ifndef WHYDRO2
  #define OMEGAB (0.0)
#else
  #define OMEGAB (0.049); // 0.049 for PLANCK
//#define OMEGAB (0.31749); // 0.049 for PLANCK
#endif

#ifdef WRAD
  #define NFLUX_R (6*NGRP*NVAR_R)
#endif
//=======================================



#define LIFETIME_OF_STARS_IN_TEST (3e9)

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
  REAL overdensity_cond;///< need overdensity_cond times the mean density to begin star formation
  REAL density_cond;///< Hydrogen density (m-3)
  REAL efficiency;///< efficiency of star formation proccess
  REAL tlife;///< life time of a radiative source (yr)
  REAL mass_res;///< the mass resolution
  int  n;///< total number of stars
  REAL thresh;///< density threshold to allow star formation
#ifdef GSLRAND
  gsl_rng *rpoiss;
#endif

};
#endif

#ifdef SUPERNOVAE
struct SNPARAM{
  REAL feedback_eff;///< feedback efficiency
  REAL feedback_frac;///< fraction of kinetic feedback over thermal feedback
  REAL Esnfb;///<  total Energy of a SN
};
#endif // SUPERNOVAE


#ifdef MOVIE
struct MOVIEPARAM{
	char* folder;
	int lmap;
	REAL xmin;
	REAL xmax;
	REAL ymin;
	REAL ymax;
	REAL zmin;
	REAL zmax;

	float* map;
	float* map_reduce;
};
#endif

struct UNITS{

  REAL unit_l;///< comoving length size of the box [meters]
  REAL unit_v;///< unit velocity
  REAL unit_t;///< unit time [seconds]
  REAL unit_n;///< unit number [moles typically]
  REAL unit_mass;///< unit mass [in kg, total mass is equal to one in unit codes]
  REAL unit_d;///< density unit [typically Omegam*rhoc in kg/m3]
  REAL unit_N;///< number density unit [typically Omegam*rhoc/mp in 1./m3]
};

struct SCALE{
  REAL l;
  REAL v;
  REAL t;
  REAL d;
  REAL p;
  REAL E;
  REAL n;
  REAL mass;
  REAL N;
};

struct UVBACKGROUND{
  int N; ///< number of samples in the input background
  REAL *redshift; ///< size N
  REAL *Nphot; ///< size N
  REAL *value; ///< size NGRP
};

//=======================================
struct ATOMIC{
  char path[1024]; ///< path of the file containing the atomic data
  int n;
  int ngrp_space;
  int ngrp_time;

  REAL *space_bound;
  REAL *time_bound;

  REAL *hnu;
  REAL *alphae;
  REAL *alphai;
  REAL *factgrp;
};

struct SPECTRUM{
  char path[1024];
  int ntime;
  int nwave;
  REAL*time;
  REAL*wavelength;
  REAL**flux;
};


struct RUNPARAMS{
  int npartmax; ///< the max particles number (per process)
  int ngridmax; ///< the max oct numbers (per process)
  int nbuff; ///< the mpi buffer size
  int ndumps; ///< the frequency of outputs
  REAL dt_dump; ///< the physical time between 2 dumps in years
  int nsteps; ///< the maximal number of timesteps

  int lcoarse; ///< the coarse level
  int lmax; ///< the max level of refinement

  int niter; ///< the maximal number of iterations for the Poisson solver

  int gstride; ///< the size of the stencil for vector based computations (gravity)
  int hstride; ///< the size of the stencil for vector based computations (hydro)

  REAL dt; ///< the timsestep
  REAL tmax; ///< the simulation stops at tmax : corresponds to amax in cosmo
  REAL time_max; ///< for cosmo only : contains the time equivalent to amax (contained in tmax, yeah its obfuscated)

  int maxhash; ///< the hash table size between hilbert keys and oct adress (should be typically = to (2^levelmax-1)^3

  REAL amrthresh0;
  REAL amrthresh; ///< the refinement criterion (refine if mcell>amrthresh)

  int DM_res; ///< resolution of dark matter particle (equivalent level of lcoarse + DM_res)
  REAL dx_res; ///< maximum spatial resolution before blocking AMR in Parsec

  int nsmooth; ///< the number of neighbour refinement steps

  REAL poissonacc; ///< relaxation accuracy for Poisson equation
  int mgridlmin;    ///< coarsest level for multigrid relaxation
  int nvcycles; ///< number of vcycles for multigrid relaxation
  int nrelax; ///< number of smoothing cycles

  int nrestart; ///< the restart snapshot
  int nsubcycles; ///< number of subcyles in AMR advance procedure

  struct OUTPUTPARAM *out;

#ifdef TESTCOSMO
  struct COSMOPARAM *cosmo; ///< the cosmological parameters
#endif

#ifdef STARS
  struct STARSPARAM *stars; ///< the star formation parameters
#endif

#ifdef SUPERNOVAE
  struct SNPARAM *sn; ///< the supernovae parameters
#endif // SUPERNOVAE

  int nthread; ///< number of GPU threads
  int nstream; ///< number of GPU streams
  int ompthread; ///< numberf of OMP threads

  struct UNITS unit; ///< contains the units
  struct SCALE scale; ///< contains the scaling factor for units convertion (function of aexp)

#ifdef WRAD
  REAL clight; ///< speed of light in units of the real one
  REAL clightorg; ///< speed of light in units of the real one // saving the original value
  REAL fudgecool; ///< cooling fraction
  int ncvgcool; ///< cooling max iterations

  REAL denthresh; ///< density threshold to turn the sources on
  REAL tmpthresh; ///< temperature threshold to turn the sources on
  REAL srcint; ///< intensity of the sources
#ifdef HOMOSOURCE
  REAL bkg; ///< the uniform background intensity
#endif

#endif

  REAL egy_rhs; ///< the right hand side of the energy conservation equation (0 in non cosmological case);
  REAL egy_0; ///< the initial energy
  REAL egy_last; ///< the last integrand for the energy equation (used for trapezoidal rule)
  REAL egy_timelast; ///< the last time for the integrand (used for trapezoidal rule)
  REAL egy_totlast;
  FILE *fpegy; ///< the file with egy stats

  REAL rzoom; ///< the inner zoom radius
  REAL fzoom; ///< the scale factor for zoom radii (>1.)
  REAL lmaxzoom; ///< the maximal zoom level

#ifdef MOVIE
	struct MOVIEPARAM *movie; ///< the movie parameters
#endif

#ifdef UVBKG
  struct UVBACKGROUND uv; ///< the UV background
#endif

  struct ATOMIC atomic;
  struct SPECTRUM spectrum;
};



//=======================================

// this structure exists for MPI communication protocol

struct PACKET{
  REAL data[8]; ///< the data to be transfered (8 since we transmit data per octs)
  //unsigned long long key; // the destination hilbert key
  double key; ///< MODKEY
  int level; ///< the level of the destination (to remove the key degeneracy)
};





//=======================================


struct CPUINFO{
  int rank; ///< the local processor ID
  int nproc; ///< the total number of processors

  unsigned long long  kmin;
  unsigned long long  kmax;
  int nkeys;

  REAL load;

  struct PART *part; ///< the particle array
  struct OCT **bndoct; ///< the list of external boundary octs

  int nebnd; ///< the number of external boundary octs
  int nnei; ///< the number of neighbors procs

  int *mpinei; ///< the ranks of the neighbors procs

  int *dict; ///< a hash table to converts ranks in local neighbour index

  struct OCT **htable; ///< the hashing table to recover the octs from hilbert keys

  int *allkmin;
  int *allkmax;

  int nbuff; ///< the number of buffer cells = to the max of # of buffer cell from 1 neighbor
  int nbufforg; ///< the max number of buffer cells (set from the parameter file)
  int nbuffpart;  ///< the number of particles to transmit

  int *nrecv; ///< the number of octs to be received by the local cpu, e.g cpu->nrecv[5] = nb of octs to be received from neigh # 5
  int *nsend; ///< the number of octs to be sent     by the local cpu, e.g cpu->nrecv[5] = nb of octs to be sent to       neigh # 5

  int *nrecv_coarse; ///< the number of l-1 octs to be received by the local cpu, e.g cpu->nrecv[5] = nb of l-1 octs to be received from neigh # 5
  int *nsend_coarse; ///< the number of l-1 octs to be sent     by the local cpu, e.g cpu->nrecv[5] = nb of l-1 octs to be sent to       neigh # 5

  struct PART *lastpart;
  struct PART *freepart;
  struct PART *firstpart;

#ifdef WMPI
  MPI_Datatype *MPI_PACKET; ///< the structured type for MPI messages (fields)
  struct PACKET **sendbuffer;
  struct PACKET **recvbuffer;
  struct PART_MPI **psendbuffer;
  struct PART_MPI **precvbuffer;
  struct HYDRO_MPI **hsendbuffer;
  struct HYDRO_MPI **hrecvbuffer;
  struct RAD_MPI **Rsendbuffer;
  struct RAD_MPI **Rrecvbuffer;

  int mpiio_grid_offsets;
  int *mpiio_ncells;

  int mpiio_part_offsets;
  int *mpiio_nparts;
#ifdef STARS
  int mpiio_star_offsets;
  int *mpiio_nstars;
#endif // STARS

#ifdef PIC
  MPI_Datatype *MPI_PART; ///< the structured type for MPI messages (particles)
#endif

#ifdef WHYDRO2
  MPI_Datatype *MPI_HYDRO; ///< the structured type for MPI messages (particles)
  MPI_Datatype *MPI_FLUX; ///< the structured type for MPI messages (particles)
#endif

#ifdef WRAD
  MPI_Datatype *MPI_RAD; ///< the structured type for MPI messages (particles)
#endif

  MPI_Comm comm; ///< the communicator


#endif

  int maxhash; ///< the size of the hashtable between hilbert keys and oct adresses
  int *noct; ///< the number of octs per levels
  int *npart; ///< the number of particles per levels
#ifdef STARS
  int *nstar;///< the number of stars per levels
  int trigstar; ///< set to 1 after the first star has been formed
#endif

  int levelcoarse; ///< the levelcoarse

  struct OCT *freeoct; ///< the location of the first free oct
  struct OCT **firstoct; ///< the location of the first free oct
  int nsteps; ///< the current coarse step index
  REAL tinit; ///< the initial time
  int *ndumps; ///< the current dump number


  int *locNoct; ///< the local number of octs per levels
  struct OCT *** octList; ///< the dictionnary of all oct of the current cpu

#ifdef GPUAXL

#ifdef WGRAV

  struct GGRID *dev_stencil;
  REAL *res;
  REAL *pnew;
  REAL *resLR;

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
  unsigned long cuparam;
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
  REAL src[NGRP];

#ifdef SUPERNOVAE
  REAL snfb;
#endif // SUPERNOVAE


  REAL nhplus;
  REAL eint;
  REAL nh;
  REAL temp; // is a direct function of eint, nh and xion but stored for conveniency
  //  REAL deltaX; // the ionization variation (to track fronts)

#ifdef HELIUM
  REAL nheplus;
  REAL nhepplus;
#endif

};
#endif

/**
  * \stuct Wtype
  * \brief local primitive hydro quantities
  */

struct Wtype{
  REAL d;   ///< gas density in unit of average barionic density
  REAL u;   ///< velocity
  REAL v;   ///< velocity
  REAL w;   ///< velocity
  REAL p;   ///< pressure
  REAL a;   ///< sound speed
  REAL E;

#ifdef WRADHYD
  REAL dX;
#ifdef HELIUM
  REAL dXHE;
  REAL dXXHE;
#endif
#endif
};


struct Wtype_MPI{
  REAL d;   ///< gas density in unit of average barionic density
  REAL u;   ///< velocity
  REAL v;   ///< velocity
  REAL w;   ///< velocity
  REAL p;   ///< pressure
  REAL a;   ///< sound speed
  REAL E;
#ifdef WRADHYD
  REAL dX;
#ifdef HELIUM
  REAL dXHE;
  REAL dXXHE;
#endif

#endif
};

/**
  * \stuct Utype
  * \brief local conservative hydro quantities
  */

struct Utype{
  REAL d;    ///< gas density in unit of average barionic density
  REAL du;   ///< momentum
  REAL dv;   ///< momentum
  REAL dw;   ///< momentum
  REAL E;    ///< Energy

#ifdef DUAL_E
  REAL eint; ///< internal energy
#endif

#ifdef WRADHYD
  REAL dX;
#ifdef HELIUM
  REAL dXHE;
  REAL dXXHE;
#endif

#endif

};


struct Wtype1D{
  REAL d;   ///< gas density in unit of average barionic density
  REAL u;   ///< velocity
  REAL p;   ///< pressure
  REAL a;   ///< sound speed
#ifdef WRADHYD
  REAL dX;
#ifdef HELIUM
  REAL dXHE;
  REAL dXXHE;
#endif

#endif
};

struct Wtype1D_double{
  double d;   ///< gas density in unit of average barionic density
  double u;   ///< velocity
  double p;   ///< pressure
  double a;   ///< sound speed
};


struct Utype1D{
  REAL d;    ///< gas density in unit of average barionic density
  REAL du;   ///< momentum
  REAL E;    ///< Energy
#ifdef WRADHYD
  REAL dX;
#ifdef HELIUM
  REAL dXHE;
  REAL dXXHE;
#endif

#endif
};



//=======================================

// in fact I think there is a bug in the particle assignement scheme
// unitl now DFACT was set (incorrectly I think) to 2 but never encountered somthing massive
// in principle should be 1.
// lets try

#define DFACT (1.)

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
  int is; ///< local timestep number per particle

  REAL epot;
  REAL ekin;

#ifdef STARS
  int  isStar;
  int radiative_state;
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
  //unsigned long long key; // the destination hilbert key
  double key; ///< the destination hilbert key


  int idx;
  int level; ///< the level of the destination (to remove the key degeneracy)
  int icell; ///< the cell of destination
  int is;    ///< current step of particle


#ifdef STARS
  int isStar;
  int radiative_state;
#endif

};

//=======================================


// this structure is for the communication of Hydro data
struct HYDRO_MPI{
  struct Wtype data[8]; ///< the data to be transfered (8 since we transmit data per octs)
  //unsigned long long key; // the destination hilbert key
  double key;///< the destination hilbert key MODKEY
  int level;///< the level of the destination (to remove the key degeneracy)
};

#ifdef WRAD
struct RAD_MPI{
  struct Rtype data[8];///< the data to be transfered (8 since we transmit data per octs)
  //unsigned long long key; // the destination hilbert key MODKEY
  double key; ///< the destination hilbert key
  int level; ///< the level of the destination (to remove the key degeneracy)
};
#endif



//=========================================================

struct Gtype{
  REAL d; ///<contrast of matter density
  REAL p; ///<pottential
};

//-------------------------------------
struct CELL
{
  struct OCT *child;
  REAL marked; ///< REAL for consistency with physical quantities during communications
  int idx;///< index of the cell within the oct

#ifdef PIC
  struct PART * phead;///< the head particle
  REAL density;///< the physical quantities
  //REAL temp;
#endif

#ifdef WGRAV
  struct Gtype gdata;
  REAL pnew; ///< new potential
  REAL res;  ///< residual
  REAL f[3]; ///< the gravitational force component
#endif


#ifdef WHYDRO2
  struct Wtype field;
  struct Wtype fieldnew;
#endif


#ifdef WRAD
  struct Rtype rfield; ///< photons/s/m3
  struct Rtype rfieldnew; ///< photons/s/m3
#endif
};


struct LCELL
{
  REAL marked; ///< REAL for consistency with physical quantities during communications
  //int idx; //index of the cell within the oct
  int child;

#ifdef PIC
   // the physical quantities */
  float density; ///< total density
 #endif

#ifdef WGRAV
  //struct Gtype gdata;
  float den;
  float pot;
  float res;  ///< residual
  float f[3]; ///< the gravitational force component
#endif


#ifdef WHYDRO2
  //struct Wtype field;
  float d;
  float u;
  float v;
  float w;
  float p;
  float dX;
#ifdef HELIUM
  float dXHE;
  float dXXHE;
#endif

#endif


#ifdef WRAD
  //  struct Rtype rfield;
  double e[NGRP];
  double fx[NGRP];
  double fy[NGRP];
  float fz[NGRP];
  double src[NGRP];
  float snfb;
  double xion;
  double temp; ///< is a direct function of eint, nh and xion but stored for conveniency

#ifdef HELIUM
  double xHE;
  double xxHE;
#endif

#endif
//  REAL sfr;
};


// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CELLFLUX
{

#ifdef WHYDRO2
  //struct Wtype field;
  struct Utype deltaU;
  REAL flux[NFLUX]; ///< 6 fluxes of 5 variables each
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
  REAL flux[NFLUX]; ///< 6 fluxes of 5 variables each
  REAL divu;
#endif

};




// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CELLLIGHT
{

#ifdef WHYDRO2
  struct Wtype field; ///< hydrodynamical data
#endif

#ifdef WRAD
  struct Rtype rfield; ///< radiation data
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
  struct Wtype field; ///< hydrodynamical data
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
  struct Rtype rfield; ///< radiation data
#endif

  char split;
};

// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CELLGRAV
{
  struct Gtype gdata; ///< gravitational data
};

// ----------------------------------------------------------------
// ----------------------------------------------------------------



//-------------------------------------
//-------------------------------------------------------------------------
struct OCT
{
  /// the cell properties
  struct CELL cell[8]; // MUSTN'T BE MOVED !!

  struct CELL *nei[6];///< neighbor cells at level - 1
  struct CELL *parent; ///< parent cell

  ///< the next two pointers are required for sweeps through a single level
  struct OCT *next; ///< next oct on the same level
  struct OCT *prev; ///< previous oct on the same level
  struct OCT *nexthash; ///< next oct in the hash list

  ///< the oct position (lowest left corner)
  REAL x;
  REAL y;
  REAL z;

  // parallel data
  int cpu;
  int level;///< level of the cells in the oct

};

struct LOCT
{
/*   // the cell properties */
   struct LCELL cell[8]; // MUSTN'T BE MOVED !! */

   // the oct position (lowest left corner) */
   float x;
   float y;
   float z;

  // parallel data
  int cpu;
  int level;///< level of the cells in the oct

};

//struct OCT *SOCT;

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
  int nei[27]; ///< pointers toward other neighbour octs in the stencil
  struct OCTGRAV oct; ///< the local one
};

// =======================================
struct STENGRAV{
  REAL *d; ///< density [8*stride]
  REAL *p; ///< potential [8*stride]
  REAL *pnew; ///< new potential [8*stride]

  REAL *res; ///< residual [8*stride]
  REAL *res2;
  REAL *resLR; ///< low res residual for MG [stride]

  int *nei; ///< neighbour indexes [7*stride] (6 real neighbours plus the middle one)
  int *level; ///< oct levels [stride]
  int *cpu; ///< cpu of octs for MPI boundaries [stride]
  char *valid; ///< validity of oct (border=0 or inner=1);
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
  REAL *vecpot; ///< contains the potential in "stride" octs
  REAL *vecpotnew; ///< contains the potential in "stride" octs
  REAL *vecden; ///< contains the density in "stride" octs


  int *vecnei;///< contains the cell neighbors of the octs
  int *vecl; ///< contains the level of the octs
  int *veccpu; ///< contains the level of the octs
  int *vecicoarse; ///< contains the level of the octs

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


struct OUTPUTPARAM{
  int n_field;
  int n_field_tot;
  char *field_name[50];
  int field_id[50];
};
