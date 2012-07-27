#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"
#include "amr.h"
#include "hydro_utils.h"
#include "friedmann.h"

// ===============================================================
// ===============================================================

void dispndt(struct RUNPARAMS *param, struct CPUINFO *cpu, int *ndt){
  
  int level;
  int ndtmax=pow(param->nsubcycles,param->lmax-param->lcoarse+1);
  int i,j,na;

  printf("\n");
  for(j=0;j<ndtmax;j++)printf("#");
  printf("\n");

  for(level=param->lcoarse;level<=param->lmax;level++){
    na=pow(2,param->lmax-level+1);

    for(j=0;j<ndt[level-1];j++){
      //one arrow
      for(i=0;i<(na-1);i++) printf("~");
      printf(">");
    }
    for(j=ndt[level-1];j<ndtmax;j++) printf(" ");
    printf("\n");

  }
  

  for(j=0;j<ndtmax;j++)printf("#");
  printf("\n \n");
  
}

// ===============================================================
// ===============================================================

REAL Advance_level(int level,REAL *adt, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct, struct HGRID *stencil, int stride, REAL aexp,struct PACKET **sendbuffer, struct PACKET **recvbuffer,int *ndt){
 

  if(cpu->rank==0){
    printf("\n === entering level =%d with stride=%d sten=%p aexp=%e\n",level,stride,stencil,aexp);
  }
  struct OCT *curoct;
  REAL dtnew;
  REAL dt;
  REAL dtfine;
  int npart=0;
  int mtot;

  dt=0.;
  int is=0;
  
  if(level==param->lcoarse){
    // ==================================== Check the number of particles and octs
    mtot=multicheck(firstoct,npart,param->lcoarse,param->lmax,cpu->rank,cpu->noct);
        
    // =========== Grid Census ========================================
      
    grid_census(param,cpu);
  }


  do{
    printf("----\n");
    printf("subscyle #%d subt=%e\n",is,dt);
    

    // ================= I we refine the current level
    if((param->lmax!=param->lcoarse)&&(level<param->lmax)){
      // refining (and destroying) octs
      curoct=L_refine_cells(level,param,firstoct,lastoct,cpu->freeoct,cpu,firstoct[0]+param->ngridmax);
      cpu->freeoct=curoct;
    }
    
    
    // ==================================== Check the number of particles and octs
    mtot=multicheck(firstoct,npart,param->lcoarse,param->lmax,cpu->rank,cpu->noct);


#ifdef WHYDRO2
#ifdef NOFLUX
    // =============================== cleaning the updated values of hydro quantities
    clean_new_hydro(level,param,firstoct,cpu);
#endif
#endif

    // == Ready to advance

    // ================= II We compute the timestep of the current level
    dtnew=param->dt*(cpu->nsteps>2);
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
    printf("sum=%e\n",adt[level-1]+dtnew);
    adt[level-1]=dtnew;

    if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one
    adt[level-1]=fmin(adt[level-1],adt[level-2]-dt);// we force synchronization

    // ================= III Recursive call to finer level
    
    if(level<param->lmax){
      if(cpu->noct[level]>0){
	dtfine=Advance_level(level+1,adt,cpu,param,firstoct,lastoct,stencil,stride,aexp,sendbuffer,recvbuffer,ndt);
	// coarse and finer level must be synchronized now
	adt[level-1]=dtfine;
	if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one
      }
    }
    


    
    // ================= IV advance solution at the current level
    
#ifndef NOFLUX
    hydro(level,param,firstoct,cpu,stencil,stride,adt[level-1]);
#else
    advancehydro(level,param,firstoct,cpu,stencil,stride,adt[level-1]);
#endif

    if((param->lmax!=param->lcoarse)&&(level<param->lmax)){
      
      // cleaning the marks
      L_clean_marks(level,firstoct);
      // marking the cells of the current level
      L_mark_cells(level,param,firstoct,param->nrelax,param->amrthresh,cpu,sendbuffer,recvbuffer);
    }

    dt+=adt[level-1]; // advance local time
    is++;
    ndt[level-1]++;

    // Some Eye candy for timesteps display

    if(cpu->rank==0) dispndt(param,cpu,ndt);
    

  }while((dt<adt[level-2])&&(is<param->nsubcycles));
  
  // ================= V Computing the new refinement map
  
  if(cpu->rank==0){
    printf("--\n");
    printf("exiting level =%d\n",level);
  }

  return dt;

}
