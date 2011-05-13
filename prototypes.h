
#ifdef WMPI
#include <mpi.h>
#endif

//=======================================

struct RUNPARAMS{
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
  
};



//=======================================

// this structure exists for MPI communication protocol

struct PACKET{
  float data[8]; // the data to be transfered (8 since we transmit data per octs)
  int key; // the destination hilbert key
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

#ifdef WMPI
  MPI_Datatype *MPI_PACKET; // the structured type for MPI messages (fields)
  MPI_Datatype *MPI_PART; // the structured type for MPI messages (particles)
  MPI_Comm comm; // the communicator
#endif

  int maxhash; // the size of the hashtable between hilbert keys and oct adresses
  int *noct; // the number of octs per levels
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

  int key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
  int icell; // the cell of destination
};

//=======================================



//-------------------------------------
struct CELL
{
  struct OCT *child;
  float marked; // float for consistency with physical quantities during communications
  int idx; //index of the cell within the oct

  // the head particle
  struct PART * phead;

  // the physical quantities
  float density;
  float pot;
  float temp;

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
