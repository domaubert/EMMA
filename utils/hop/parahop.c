#include <stdio.h>
#include <stdlib.h>
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

  sscanf(argv[2],"%d",&i);
  sscanf(argv[3],"%d",&ncpu);
  sscanf(argv[4],"%d",&mpi_size);
  sscanf(argv[5],"%d",&mpi_rank);

  printf("rank %d will process snap #%d \n",mpi_rank,i);

  strcpy(format,chop);
  strcat(format,argv[1]);
  strcat(format,"/part.%05d");
  strcat(format," -nf %d -o %s/hop.\%05d.h\%05d -paranproc %d -pararank %d"); 
  sprintf(commande,format,i,ncpu,argv[1],i,mpi_rank,mpi_size,mpi_rank);
  printf("%s\n",commande);
  system(commande);
  
  
  strcpy(format,cregroup);
  strcat(format,"%s/hop.\%05d.h\%05d -o %s/reg.\%05d.h%05d");
  sprintf(commande,format,argv[1],i,mpi_rank,argv[1],i,mpi_rank);
  printf("%s\n",commande);
  system(commande);


  return 0;
}
