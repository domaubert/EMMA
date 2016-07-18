/*
 * Create a level 7 grid and test if the number of oct created is good on all level
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "prototypes.h"
#include "segment.h"
#include "amr.h"
#include "oct.h"
#include "communication.h"


int main(int argc, char *argv[]){

//      setup

        const int levelcoarse=6;
        const int levelmax=levelcoarse;
        const int ngridmax=POW(2,3*levelcoarse);

        struct RUNPARAMS param;
        param.lcoarse=levelcoarse;
        param.lmax=levelcoarse;
        param.nbuff=50000;

        struct CPUINFO cpu;

        MPI_Init(&argc,&argv);
        MPI_Datatype MPI_PACKET;
        MPI_Datatype MPI_PART;
        MPI_Datatype MPI_WTYPE;
        MPI_Datatype MPI_HYDRO;
        MPI_Datatype MPI_RTYPE;
        MPI_Datatype MPI_RAD;
        init_MPI(&cpu, &MPI_PACKET,&MPI_PART,&MPI_WTYPE,&MPI_HYDRO,&MPI_RTYPE,&MPI_RAD);

        load_balance(levelcoarse,&cpu);

        struct OCT *grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT));
        struct OCT **firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
        struct OCT **lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));

        struct CELL root;
        root = build_initial_grid(grid, firstoct, lastoct, &cpu, &param);

        const int val=(POW(2,param.lmax-1)<=512?POW(2,param.lmax-1):512); // limit to 2097152 octs in hash table i.e. 16e6 cells
        cpu.maxhash=POW(val,3);
        cpu.htable =	(struct OCT**) calloc(cpu.maxhash,sizeof(struct OCT *));

        cpu.nbuff=param.nbuff;
        cpu.nbufforg=param.nbuff;
        cpu.bndoct=(struct OCT**)calloc(cpu.nbufforg,sizeof(struct OCT*));

        cpu.mpinei=NULL;
        cpu.dict=NULL;

        cpu.sendbuffer=NULL;
        cpu.recvbuffer=NULL;
        cpu.psendbuffer=NULL;
        cpu.precvbuffer=NULL;
        cpu.hsendbuffer=NULL;
        cpu.hrecvbuffer=NULL;
        cpu.Rsendbuffer=NULL;
        cpu.Rrecvbuffer=NULL;
        setup_mpi(&cpu, firstoct, levelmax, levelcoarse, ngridmax, 1);

//      test

        int level;
        for(level=1;level<=levelmax;level++){

                // printf("level=%d\n",level);
                int noct=0;

                struct OCT  *nextoct = firstoct[level-1];
                if(nextoct==NULL) continue;
                do{
                        const struct OCT  *curoct=nextoct;
                        nextoct=curoct->next;
                        if(curoct->cpu!=cpu.rank) continue;
                        noct++;
                 }while(nextoct!=NULL);

                // printf("%d %d\n", noct, POW2(3*(level-1) ) );
                if (level>1){
                        assert( cpu.nproc==4);
                        assert( noct == POW2(3*(level-1))/4 );
                }

                MPI_Allreduce(MPI_IN_PLACE,&noct,1,MPI_INT,MPI_SUM,cpu.comm);
                // printf("%d %d\n", noct, POW2(3*(level-1) ) );
                assert( noct == POW2(3*(level-1)) );
        }

//      free

        free(grid);
        free(firstoct);
        free(lastoct);
        free(cpu.htable);
        free(cpu.bndoct);

        printf("Test passed\n");
        return 0;
}
