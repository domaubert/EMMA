#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  char nom[10];
  int mpi_size,mpi_rank;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

  gethostname(nom,10);
  printf("je suis %s avec le rank %d sur %d\n",nom,mpi_rank,mpi_size);
  MPI_Finalize();
  return 0;
}
