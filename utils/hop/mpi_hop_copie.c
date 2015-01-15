#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> // PAS OUBLIER
#include <string.h>

// Exemple
// compilation mpicc toto.c -o exec
// excution mpirun -np 16 exec 1 3
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

  MPI_Init(&argc,&argv); // envoi des parametres de la ligne de commande aux processeurs
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);


  sscanf(argv[2],"%d",&imin);
  sscanf(argv[3],"%d",&imax);

  // DECOUPAGE DES TACHES
  njobs=(imax-imin)+1; // NOMBRE DE TACHE GLOBAL A EXECUTER PAR TOUS LES PROCS
  stride=njobs/mpi_size; // NOMBRE DE TACHE PAR PROCESSEUR

  i1=imin+mpi_rank*stride; // DECOUPAGE DES TACHES
  i2=imin+(mpi_rank+1)*stride-1;
  if(mpi_rank==(mpi_size-1)) i2=imax;
  // FIN DECOUPAGE
  
  printf("rank %d will process snaps #%d to #%d\n",mpi_rank,i1,i2);

  for(i=i1;i<=i2;i++){

    // CONSTRUCTION DE LA LIGNE DE COMMANDE
    strcpy(format,chop);
    strcat(format,argv[1]);
    strcat(format,"/part.%05d");
    strcat(format," -nf %d -o %s/hop.\%05d"); 
    sprintf(commande,format,i,ncpu,argv[1],i);
    printf("%s\n",commande);
    /// === fin construction

    system(commande);


  }

  MPI_Finalize();
  return 0;
}
