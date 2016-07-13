/*
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

#include "prototypes.h"
#include "communication.h"

int main(int argc, char *argv[]){

//      setup

        MPI_Init(&argc,&argv);

        struct CPUINFO cpu;

        MPI_Datatype MPI_PACKET;
        MPI_Datatype MPI_PART;
        MPI_Datatype MPI_WTYPE;
        MPI_Datatype MPI_HYDRO;
        MPI_Datatype MPI_RTYPE;
        MPI_Datatype MPI_RAD;

        init_MPI(&cpu, &MPI_PACKET,&MPI_PART,&MPI_WTYPE,&MPI_HYDRO,&MPI_RTYPE,&MPI_RAD);

        printf("proc %d/%d\n", cpu.rank, cpu.nproc);

        printf("Test passed\n");
        return 0;
}
