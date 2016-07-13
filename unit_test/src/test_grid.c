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


int main(int argc, char *argv[]){

//      setup

        const int levelcoarse=6;
        const int levelmax=levelcoarse;
        const int ngridmax=POW(2,3*levelcoarse);

        struct RUNPARAMS param;
        param.lcoarse=levelcoarse;

        struct CPUINFO cpu;
        cpu.rank=0;
        cpu.nproc=1;
        load_balance(levelcoarse,&cpu);

        struct OCT *grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT));
        struct OCT **firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
        struct OCT **lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));

        cpu.locNoct=(int *)calloc(levelmax,sizeof(int));
        cpu.octList=(struct OCT***)calloc(levelmax,sizeof(struct OCT**));

        int iLev;
        for(iLev = 0; iLev<levelcoarse; iLev++){
          cpu.locNoct[iLev] = (pow(2,3*(iLev+1))<ngridmax? pow(2,3*(iLev+1)):ngridmax) ;
          cpu.octList[iLev] = (struct OCT**)calloc(cpu.locNoct[iLev],sizeof(struct OCT*));
        }

        struct CELL root;
        root = build_initial_grid(grid, firstoct, lastoct, &cpu, &param);

//      test

        int level;
        for(level=1;level<=levelmax;level++){
                int noct=0;

                struct OCT  *nextoct = firstoct[level-1];
                 do{
                        if(nextoct==NULL) 		continue;
                        struct OCT  *curoct=nextoct;
                        nextoct=curoct->next;
                        noct++;
                 }while(nextoct!=NULL);

                assert( noct == POW2(3*(level-1)));
        }

        for(level=1;level<=param.lmax;level++){
                setOctList(firstoct[level-1], &cpu, &param,level);
        }

        for(level=1;level<levelmax;level++){
                FILE *fp;
                fp = fopen ("filename.txt","r");
                if (fp!=NULL)
                {
                        fscanf(fp,"Some String\n", &var);
                        fclose (fp);
                }printf("%d\t%d\t%d\n",level,cpu.locNoct[level-1],POW2(3*level) );
                assert(cpu.locNoct[level-1] == POW2(3*level) );
        }

//      free

        free(grid);
        free(firstoct);
        free(lastoct);

        printf("Test grid creation passed\n");
        return 0;
}
