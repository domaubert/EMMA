/*
 * create a grid, create a uniform particle distribution of particle
 * test if each cells get one and only one particle
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "prototypes.h"
#include "segment.h"
#include "amr.h"
#include "ic.h"
#include "communication.h"

int main(int argc, char *argv[]){

//      setup

        const int levelcoarse=6;
        const int levelmax=levelcoarse;
        const int ngridmax=POW2(3*levelcoarse);
        const int npartmax=POW2(3*levelcoarse)*2;
        const int nside = POW2(levelcoarse);

        struct RUNPARAMS param;
        param.lcoarse=levelcoarse;
        param.lmax=levelcoarse;

        struct CPUINFO cpu;
        cpu.rank=0;
        cpu.nproc=1;
        load_balance(levelcoarse,&cpu);
        cpu.levelcoarse=levelcoarse;

//        init

        struct OCT *grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT));
        struct OCT **firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
        struct OCT **lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
        struct CELL root;
        root = build_initial_grid(grid, firstoct, lastoct, &cpu, &param);


        const int val=(POW(2,param.lmax-1)<=512?POW(2,param.lmax-1):512); // limit to 2097152 octs in hash table i.e. 16e6 cells
        cpu.maxhash=POW(val,3);
        cpu.htable =	(struct OCT**) calloc(cpu.maxhash,sizeof(struct OCT *));
        setup_mpi(&cpu, firstoct, levelmax, levelcoarse, ngridmax, 1);


        struct PART *part=(struct PART*)calloc(npartmax,sizeof(struct PART));

        printf("Build uniform particle grid\n");
        int npart = POW2(3*levelcoarse);
        REAL dx = 1./nside;

        int k;
        for(k=0;k<nside;k++){
                int j;
                for(j=0;j<nside;j++){
                        int i;
                        for(i=0;i<nside;i++){
                                const int idx = i+j*nside+k*nside*nside;
                                part[idx].x = i*dx +dx/2;
                                part[idx].y = j*dx +dx/2;
                                part[idx].z = k*dx +dx/2;
                                part[idx].level = levelcoarse;
                        }
                }
        }

        printf("part2grid\n");
        part2grid(part,&cpu,npart);

//      test

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
        int part_in_cell_max=0;

        int level;
        for(level=1;level<=levelmax;level++){

                struct OCT  *nextoct = firstoct[level-1];
                 do{
                        if(nextoct==NULL) 		continue;
                        struct OCT  *curoct=nextoct;
                        nextoct=curoct->next;

                        int icell;
                        for(icell=0;icell<8;icell++) {
                                int count_part_cell=0;

                                struct CELL *curcell = &curoct->cell[icell];
                                struct PART *nexp=curcell->phead;
                                if(nexp!=NULL){
                                        do{
                                                struct PART *curp=nexp;
                                                nexp=curp->next;
                                                count_part++;
                                                count_part_cell++;
                                        }while(nexp!=NULL);
                                }
                                part_in_cell_max=fmax(part_in_cell_max,count_part_cell);
                        }
                }while(nextoct!=NULL);
        }



        // printf("count_part %d\n", count_part);
        assert( count_part == POW2(3*(levelcoarse)));

        // printf("part_in_cell_max %d\n", part_in_cell_max);
        assert( part_in_cell_max == 1 );

//      free

        free(grid);
        free(firstoct);
        free(lastoct);
        free(part);
        free(cpu.htable);

        printf("Test passed\n");
        return 0;
}
