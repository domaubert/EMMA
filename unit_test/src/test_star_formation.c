/*
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "prototypes.h"
#include "segment.h"
#include "amr.h"
#include "oct.h"
#include "ic.h"
#include "communication.h"
#include "stars.h"

int main(int argc, char *argv[]){

//      setup

        const int levelcoarse=2;
        const int levelmax=levelcoarse;
        const int ngridmax=POW2(3*levelcoarse);
        const int npartmax=POW2(3*(levelcoarse+1));
        const int nside = POW2(levelcoarse);

        srand(0000);

        struct RUNPARAMS param;
        param.lcoarse=levelcoarse;
        param.lmax=levelcoarse;
        param.nbuff=50000;
        param.npartmax=npartmax;
        param.unit.unit_mass = SOLAR_MASS;
        param.unit.unit_d=1.;
        param.unit.unit_t=1.;

        struct CPUINFO cpu;
        cpu.rank=0;
        cpu.nproc=1;
        load_balance(levelcoarse,&cpu);
        cpu.levelcoarse=levelcoarse;
        cpu.nbuff=param.nbuff;
        cpu.nbufforg=param.nbuff;
        cpu.nstar=	(int *)calloc(levelmax,sizeof(int));
        cpu.trigstar=0;

//        init

        struct OCT *grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT));
        struct OCT **firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
        struct OCT **lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
        struct CELL root;
        root = build_initial_grid(grid, firstoct, lastoct, &cpu, &param);

        printf("octlist\n");

        cpu.locNoct=(int *)calloc(levelmax,sizeof(int));
        cpu.octList=(struct OCT***)calloc(levelmax,sizeof(struct OCT**));

        int iLev;
        for(iLev = 0; iLev<levelcoarse; iLev++){
          cpu.locNoct[iLev] = (pow(2,3*(iLev+1))<ngridmax? pow(2,3*(iLev+1)):ngridmax) ;
          cpu.octList[iLev] = (struct OCT**)calloc(cpu.locNoct[iLev],sizeof(struct OCT*));
        }

        int level;
        for(level=1;level<=param.lmax;level++){
                setOctList(firstoct[level-1], &cpu, &param,level);
        }

        printf("hash\n");
        int val=(POW(2,param.lmax-1)<=512?POW(2,param.lmax-1):512); // limit to 2097152 octs in hash table i.e. 16e6 cells
        cpu.maxhash=POW(val,3);
        cpu.htable=(struct OCT**)calloc(cpu.maxhash,sizeof(struct OCT *));
        cpu.bndoct=(struct OCT**)calloc(cpu.nbufforg,sizeof(struct OCT*));

        cpu.mpinei=NULL;
        cpu.dict=NULL;
        printf("setup mpi\n");
        setup_mpi(&cpu, firstoct, levelmax, levelcoarse, ngridmax, 0);

        printf("part\n");
        struct PART *part = (struct PART*)calloc(npartmax,sizeof(struct PART));
        cpu.firstpart=part;

        struct PART * lastpart = part+1;
        cpu.lastpart=lastpart;

        struct PART * freepart = part+2;
        cpu.freepart=freepart;
        freepart->prev=NULL;
        freepart->next=freepart+1;

        struct PART *curp;
        for(curp=freepart+1;curp<(part+npartmax);curp++){
          curp->prev=curp-1;
          curp->next=NULL;
          if(curp!=(part+npartmax-1)) curp->next=curp+1;
        }


        printf("setting uniform density\n");

        REAL mean_dens= SQRT(3.*M_PI/(32.*NEWTON_G)) * 128;

        int iOct;
        for(iOct=0; iOct<cpu.locNoct[levelcoarse-1]; iOct++){
                struct OCT *curoct=cpu.octList[levelcoarse-1][iOct];

                int icell;
                for(icell=0;icell<8;icell++) {
                        struct CELL *curcell = &curoct->cell[icell];
                        curcell->field.d=mean_dens;
                }
        }

        struct STARSPARAM stars;
        stars.n=0;
        param.stars=&stars;
        param.stars->thresh = 0;
        param.stars->mass_res = 150;
        param.stars->efficiency=1.;

        REAL dt = 1.;
        REAL aexp = 1.;
        int is = 0;
        Stars(&param, &cpu, dt, aexp, levelcoarse, is);

//      test

        int count_part=0;
        int part_in_cell_max=0;
        REAL tot_mass = 0;

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

                                REAL dv = 1./POW2(3*level);
                                REAL cell_mass=curcell->field.d * dv;

                                if(nexp!=NULL){
                                        do{
                                                struct PART *curp=nexp;
                                                nexp=curp->next;
                                                count_part++;
                                                count_part_cell++;

                                                assert( curp->isStar == 6 );

                                                assert(curp->x>=0);
                                                assert(curp->x <1);
                                                assert(curp->y>=0);
                                                assert(curp->y <1);
                                                assert(curp->z>=0);
                                                assert(curp->z <1);
                                                assert(curp->level == levelcoarse );

                                                assert( curp->mass == 150 );
                                                cell_mass+=curp->mass;

                                        }while(nexp!=NULL);
                                }
                                part_in_cell_max=fmax(part_in_cell_max,count_part_cell);

                                if (level==levelcoarse){
                                // printf("cell=%f mean=%f\n", cell_mass,  mean_dens*dv);
                                assert( cell_mass ==  mean_dens*dv );
                                }
                        }
                }while(nextoct!=NULL);
        }

        printf("count_part %d\n", count_part);
        assert( count_part == POW2(3*(levelcoarse)));

        printf("part_in_cell_max %d\n", part_in_cell_max);
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
