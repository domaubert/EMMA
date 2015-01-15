#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

// Exemple

//mpirun -np 3 mpi_parahop ./data 27 29 32


int main(int argc, char *argv[]){

  int i,i1,i2;
  int imax,imin,njobs,stride;
  int mpi_size,mpi_rank;
  int nx,ny,nz;
  int ncpu;

  char cparahop[255]="./parahop "; 
  char commande[512];
  char dir[255];
  char outname[255];
  char format[255];
  int full;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);


  sscanf(argv[2],"%d",&imin);
  sscanf(argv[3],"%d",&imax);
  sscanf(argv[4],"%d",&ncpu);

  printf("rank %d will process snaps #%d to #%d\n",mpi_rank,i1,i2);

  for(i=imin;i<=imax;i++){
    strcpy(format,cparahop);
    strcat(format,argv[1]);
    strcat(format," %d %d %d %d");
    sprintf(commande,format,i,ncpu,mpi_size,mpi_rank);
    printf("%s\n",commande);
    system(commande);
  }

  MPI_Finalize();
  return 0;
}
