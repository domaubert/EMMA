/*
 * Create a grid, read graphic IC and assign part to grid
 * test if the number of read part is equal to 2^(3*Lcoarse)
 * test if all part are in (0,1,0,1,0,1) cube
 * test if the number of part on the grid is the same as the number of read part
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "prototypes.h"
#include "segment.h"
#include "amr.h"
#include "ic.h"
#include "communication.h"

int main(int argc, char *argv[]){

//      setup

        const int levelcoarse=6;
        const int levelmax=levelcoarse;
        const int ngridmax=POW(2,3*levelcoarse);
        const int npartmax=POW(2,3*levelcoarse)*2;

        struct RUNPARAMS param;
        param.lcoarse=levelcoarse;
        param.lmax=levelcoarse;

        struct COSMOPARAM cosmo;
        param.cosmo=&cosmo;

        struct CPUINFO cpu;
        cpu.rank=0;
        cpu.nproc=1;
        load_balance(levelcoarse,&cpu);
        cpu.levelcoarse=levelcoarse;


//        init

        struct OCT *grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT));
        struct OCT **firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
        struct OCT **lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));

        struct PART *part=(struct PART*)calloc(npartmax,sizeof(struct PART));

        struct CELL root;
        root = build_initial_grid(grid, firstoct, lastoct, &cpu, &param);

        REAL munit;
        REAL ainit;
        int npart;
        struct PART *lastpart=read_grafic_part(part, &cpu, &munit, &ainit, &npart, &param,levelcoarse);

        int val=(POW(2,param.lmax-1)<=512?POW(2,param.lmax-1):512); // limit to 2097152 octs in hash table i.e. 16e6 cells
        cpu.maxhash=POW(val,3);
        cpu.htable =	(struct OCT**) calloc(cpu.maxhash,sizeof(struct OCT *));
        setup_mpi(&cpu, firstoct, levelmax, levelcoarse, ngridmax, 1);

        part2grid(part,&cpu,npart);

//      test

        assert( npart == POW2(3*(levelcoarse)));

        int i;
        for(i=0;i<npart;i++){
                assert(part[i].x>=0);
                assert(part[i].x <1);
                assert(part[i].y>=0);
                assert(part[i].y <1);
                assert(part[i].z>=0);
                assert(part[i].z <1);
                assert(part[i].level == levelcoarse );
        }

        int count_part=0;

        int level;
        for(level=1;level<=levelmax;level++){

                struct OCT  *nextoct = firstoct[level-1];
                 do{
                        if(nextoct==NULL) 		continue;
                        struct OCT  *curoct=nextoct;
                        nextoct=curoct->next;

                        int icell;
                        for(icell=0;icell<8;icell++) {
                                struct CELL *curcell = &curoct->cell[icell];
                                struct PART *nexp=curcell->phead;
                                if(nexp!=NULL){
                                        do{
                                                struct PART *curp=nexp;
                                                nexp=curp->next;
                                                count_part++;

                                        }while(nexp!=NULL);
                                }
                        }
                }while(nextoct!=NULL);
        }

        // printf("count_part %d\n", count_part);
        assert( count_part == POW2(3*(levelcoarse)));

//      free

        free(grid);
        free(firstoct);
        free(lastoct);
        free(part);
        free(cpu.htable);

        printf("Test passed\n");
        return 0;
}
