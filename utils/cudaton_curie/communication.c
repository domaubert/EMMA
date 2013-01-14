#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WMPI
#include <mpi.h>
#include "Interface.h"

void get_elapsed(double *timein)
{
  *timein=MPI_Wtime();
  //printf("%lf %f\n",tloc,*timein);
}


int getlocalrank(void)
{
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  return mpi_rank;

}

void switchbuff(float *buff, int neighbor, int ndata)
{
  MPI_Status status;
  //  static int tag=0;
  int tag=0;
  int mpi_rank;
  

  FILE *ftest;
  char fnt[256];
  int dum=1;


  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Sendrecv_replace(buff,ndata,MPI_FLOAT,neighbor,tag,neighbor,tag,MPI_COMM_WORLD,&status);
  MPI_Barrier(MPI_COMM_WORLD);
   //tag++;
}

void mpisynch(void)
{
   MPI_Barrier(MPI_COMM_WORLD);
}

double mpireducemax(double *x)
{
  double xmax;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(x,&xmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  return xmax;
}


void topo_cartesian(int rank, int *dims, int *coords, int *neigh)
{
  int periods[3];
  
  MPI_Comm commcart;
  
  periods[0]=1;
  periods[1]=1;
  periods[2]=1;

  // creation of cartesian communicator
  MPI_Cart_create(MPI_COMM_WORLD,3,dims,periods,0,&commcart);

  // getting the cartesian position
  MPI_Cart_coords(commcart,rank,3,coords);
  
  // getting the neighbors
  MPI_Cart_shift(commcart,0,1,neigh+0,neigh+1); //X
  MPI_Cart_shift(commcart,1,1,neigh+2,neigh+3); //Y
  MPI_Cart_shift(commcart,2,1,neigh+4,neigh+5); //Z

  //printf(" proc #%d has coordinates %d %d %d and neighbors Xm=%d Xp=%d Ym=%d Yp=%d Zm=%d Zp=%d \n",rank,coords[0],coords[1],coords[2],neigh[0],neigh[1],neigh[2],neigh[3],neigh[4],neigh[5]);
  //  printf(" dims = %d %d %d\n",dims[0],dims[1],dims[2]);

}


#endif
