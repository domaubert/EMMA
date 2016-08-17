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
        assert(cpu.rank < cpu.nproc);

        int loc_int = 1;
        MPI_Allreduce(MPI_IN_PLACE,&loc_int,1,MPI_INT,MPI_SUM,cpu.comm);
        assert(loc_int == cpu.nproc);

        // TODO add some message passing test to test MPI_Datatypes

        MPI_Barrier(cpu.comm);

        printf("Test passed\n");
        return 0;
}
