#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

// Exemple
//  mpirun -np 3 ./a.out ./data_coarse_256_24MPC_alt_th 1 6 256


int main(int argc, char *argv[]){

  int i,i1,i2;
  int imax,imin,njobs,stride;
  int mpi_size,mpi_rank;
  int nx,ny,nz;
  int ncpu;

  char chop[255]="./hop -p 1. -in "; 
  char cregroup[255]="./regroup -f77 -douter 80. -dsaddle 200. -dpeak 240. -root "; 
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

  njobs=(imax-imin)+1;
  stride=njobs/mpi_size;

  i1=imin+mpi_rank*stride;
  i2=imin+(mpi_rank+1)*stride-1;
  if(mpi_rank==(mpi_size-1)) i2=imax;

  printf("rank %d will process snaps #%d to #%d\n",mpi_rank,i1,i2);

  for(i=i1;i<=i2;i++){
    strcpy(format,chop);
    strcat(format,argv[1]);
    strcat(format,"/part.%05d");
    strcat(format," -nf %d -o %s/hop.\%05d"); 
    sprintf(commande,format,i,ncpu,argv[1],i);
    printf("%s\n",commande);
    system(commande);


    strcpy(format,cregroup);
    strcat(format,"%s/hop.\%05d -o %s/reg.\%05d");
    sprintf(commande,format,argv[1],i,argv[1],i);
    printf("%s\n",commande);
    system(commande);


  }

  MPI_Finalize();
  return 0;
}
