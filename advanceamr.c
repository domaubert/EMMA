#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"
#include "amr.h"
#include "hydro_utils.h"
#include "friedmann.h"


void Advance_level(int level,REAL *adt, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct, struct HGRID *stencil, int stride, REAL aexp,struct PACKET **sendbuffer, struct PACKET **recvbuffer){
  
  struct OCT *curoct;
  REAL dtnew;
  REAL dt;
  REAL dtmax;
  int npart=0;
  int mtot;

  dt=0.;

  do{
    // ================= I we refine the current level
    if((param->lmax!=param->lcoarse)&&(level<param->lmax)){
    
      // marking the cells of the current level
      L_mark_cells(level,param,firstoct,param->nrelax,param->amrthresh,cpu,sendbuffer,recvbuffer);
    
      // refining (and destroying) octs
      curoct=L_refine_cells(level,param,firstoct,lastoct,cpu->freeoct,cpu,firstoct[0]+param->ngridmax);
      cpu->freeoct=curoct;
    }
    
    //======================================= cleaning the marks
    
    clean_marks(param->lmax,firstoct);
    
    // ==================================== Check the number of particles and octs
    mtot=multicheck(firstoct,npart,param->lcoarse,param->lmax,cpu->rank,cpu->noct);
        
    // =========== Grid Census ========================================
      
    grid_census(param,cpu);
  
    // == Ready to advance

    // ================= II We compute the timestep of the current level
    dtnew=adt[level-1];
#ifdef TESTCOSMO
    REAL dtcosmo;
    dtcosmo=-0.5*sqrt(omegam)*integ_da_dt_tilde(aexp*1.1,aexp,omegam,omegav,1e-8);
    dtnew=(dtcosmo<dtnew?dtcosmo:dtnew);
    printf("dtcosmo= %e ",dtcosmo);
#endif
  
#ifdef WHYDRO2
    REAL dthydro;
    dthydro=L_comptstep_hydro(level,param,firstoct,1.0,1.0,cpu,1e9);
    dtnew=(dthydro<dtnew?dthydro:dtnew);
    printf("dthydro= %e ",dthydro);
#endif
    
    adt[level-1]=dtnew;
    if(level==param->lcoarse) adt[level-2]=dtnew; // we synchronize coarser levels with the coarse one
    
    // ================= III Recursive call to finer level
    
    if(level<param->lmax){
      if(cpu->noct[level]>0) Advance_level(level+1,adt,cpu,param,firstoct,lastoct,stencil,stride,aexp,sendbuffer,recvbuffer);
    }
    
    // finer level must be synchronized now
    
    // ================= IV advance solution at the current level
    
    hydro(level,param,firstoct,cpu,stencil,stride,adt[level-1]);

    dt+=fmin(adt[level-1],adt[level-2]-dt);
  }while(dt<adt[level-2]);
  

}
