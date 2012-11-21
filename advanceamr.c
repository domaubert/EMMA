#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"
#include "amr.h"
#include "hydro_utils.h"
#include "friedmann.h"
#include "cic.h"
#include "particle.h"

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
#ifdef WHYDRO2
REAL L_comptstep_hydro(int level, struct RUNPARAMS *param,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax){
  
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icell;
  REAL aa;
  REAL va,vx,vy,vz;
  REAL dt;
  REAL Smax=0.,S1;

  //Smax=fmax(Smax,sqrt(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);
  // Computing new timestep
  dt=tmax;
  // setting the first oct
  
  nextoct=firstoct[level-1];
  
  if(nextoct!=NULL){
    dxcur=pow(0.5,level); // +1 to protect level change
    do // sweeping through the octs of level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    vx=curoct->cell[icell].field.u; 
	    vy=curoct->cell[icell].field.v; 
	    vz=curoct->cell[icell].field.w; 
	    va=sqrt(vx*vx+vy*vy+vz*vz); 
	    aa=sqrt(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d); 
	    Smax=fmax(Smax,va+aa); 

	  }
      }while(nextoct!=NULL);
  }

  dt=fmin(dxcur*CFL/(Smax*3.),dt);

  /* #ifdef WMPI */
  /*   // reducing by taking the smallest time step */
  /*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_REAL,MPI_MIN,cpu->comm); */
  /* #endif   */

  return dt;
}
#endif

// ===============================================================
#ifdef WGRAV
REAL L_comptstep_ff(int level,struct RUNPARAMS *param,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax){
  
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icell;
  REAL dtloc;
  REAL dt;

  dt=tmax;
  // setting the first oct
      
  nextoct=firstoct[level-1];
      
  if(nextoct!=NULL){
    dxcur=pow(0.5,level); // +1 to protect level change
    do // sweeping through the octs of level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    dtloc=0.1*sqrt(2.*M_PI/(3.*curoct->cell[icell].gdata.d*aexp));
	    dt=fmin(dt,dtloc);
	  }
      }while(nextoct!=NULL);
  }
  return dt;
}
#endif

// ===============================================================
// ===============================================================

REAL Advance_level(int level,REAL *adt, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct, struct HGRID *stencil, struct STENGRAV *gstencil, int stride,struct PACKET **sendbuffer, struct PACKET **recvbuffer,int *ndt, int nsteps){
 
#ifdef TESTCOSMO
  struct COSMOPARAM *cosmo;
  cosmo=param->cosmo;
#endif
  struct OCT *curoct;
  REAL dtnew;
  REAL dt=0.;
  REAL dtold;
  REAL dtfine;
  int npart=0;
  int mtot;
  int is=0;
  REAL tloc;
  REAL aexp;

#ifdef TESTCOSMO
  tloc =cosmo->tsim;
  aexp=cosmo->aexp;
#else
  tloc=0.;
  aexp=1.0;
#endif

  if(cpu->rank==0){
    printf("\n === entering level =%d with stride=%d sten=%p aexp=%e\n",level,stride,stencil,aexp);
  }

  // ==================================== Check the number of particles and octs
  mtot=multicheck(firstoct,npart,param->lcoarse,param->lmax,cpu->rank,cpu);

  if(level==param->lcoarse){
        
    // =========== Grid Census ========================================
    grid_census(param,cpu);
  }


  do{
    printf("----\n");
    printf("subscyle #%d subt=%e\n",is,dt);
    
#ifdef TESTCOSMO
    cosmo->aexp=interp_aexp(tloc,cosmo->tab_aexp,cosmo->tab_ttilde);
    aexp=cosmo->aexp;
#endif


    // ================= I we refine the current level
    if((param->lmax!=param->lcoarse)&&(level<param->lmax)){
      // refining (and destroying) octs
      curoct=L_refine_cells(level,param,firstoct,lastoct,cpu->freeoct,cpu,firstoct[0]+param->ngridmax);
      cpu->freeoct=curoct;
    }
    



    // =============================== cleaning 
#ifdef WHYDRO2
    clean_new_hydro(level,param,firstoct,cpu);
#endif

#ifdef PIC
    L_clean_dens(level,param,firstoct,cpu);
#endif


    // == Ready to advance


    // ================= IV advance solution at the current level


    printf("ndt=%d nsteps=%d\n",ndt[param->lcoarse-1],nsteps);


#ifdef PIC
    // ==================================== performing the CIC assignement
    L_cic(level,firstoct,param,cpu);
#endif
 
    /* //==================================== Getting Density ==================================== */
    FillDens(level,param,firstoct,cpu); 

#ifdef WGRAV

    /* //====================================  Poisson Solver ========================== */
    PoissonSolver(level,param,firstoct,cpu,gstencil,stride,aexp); 

    /* //====================================  Force Field ========================== */
    PoissonForce(level,param,firstoct,cpu,gstencil,stride,aexp);

#endif




    // ================= II We compute the timestep of the current level
    dtnew=param->dt;//*(cpu->nsteps>2);
    
#ifdef TESTCOSMO
    REAL dtcosmo;
    dtcosmo=-0.5*sqrt(cosmo->om)*integ_da_dt_tilde(aexp*1.1,aexp,cosmo->om,cosmo->ov,1e-8);
    dtnew=(dtcosmo<dtnew?dtcosmo:dtnew);
    printf("dtcosmo= %e ",dtcosmo);

#ifdef WGRAV
    REAL dtff;
    dtff=L_comptstep_ff(level,param,firstoct,aexp,cpu,1e9);

    dtnew=(dtff<dtnew?dtff:dtnew);
    printf("dtff= %e ",dtff);
#endif

#endif
  
#ifdef WHYDRO2
    REAL dthydro;
    dthydro=L_comptstep_hydro(level,param,firstoct,1.0,1.0,cpu,1e9);
    dtnew=(dthydro<dtnew?dthydro:dtnew);
    printf("dthydro= %e ",dthydro);
#endif

#ifdef PIC
    REAL dtpic;
    dtpic=L_comptstep(level,param,firstoct,1.0,1.0,cpu,1e9);
    printf("dtpic= %e ",dtpic);
    dtnew=(dtpic<dtnew?dtpic:dtnew);
#endif



    printf("\n");
    dtold=adt[level-1];
    adt[level-1]=dtnew;

    if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one
    adt[level-1]=fmin(adt[level-1],adt[level-2]-dt);// we force synchronization




   // ================= III Recursive call to finer level
    if(level<param->lmax){
      if(cpu->noct[level]>0){
	dtfine=Advance_level(level+1,adt,cpu,param,firstoct,lastoct,stencil,gstencil,stride,sendbuffer,recvbuffer,ndt,nsteps);
	// coarse and finer level must be synchronized now
	adt[level-1]=dtfine;
	if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one
      }
    }


#ifdef PIC
    //================ Part. Update ===========================

    if(cpu->rank==0) printf("Start PIC on %d part with dt=%e on level %d\n",cpu->npart[level-1],adt[level-1],level);
    

    // =============================== Leapfrog Euler Step # 0
    if((is==0)&&(cpu->nsteps==0)){
      L_accelpart(level,firstoct,adt[level-1]*0.5,cpu); // computing the particle acceleration and velocity
    }
    else{
      L_accelpart(level,firstoct,(adt[level-1]+dtold)*0.5,cpu); // computing the particle acceleration and velocity
    }

    L_movepart(level,firstoct,adt[level-1],cpu); // moving the particles

    L_partcellreorg(level,firstoct); // reorganizing the particles of the level throughout the mesh


#endif


#ifdef WHYDRO2
      //=============== Hydro Update ======================
    HydroSolver(level,param,firstoct,cpu,stencil,stride,adt[level-1]);

      // ================================= gravitational correction for Hydro

    grav_correction(level,param,firstoct,cpu,adt[level-1]);
#endif


    //}
  // ================= V Computing the new refinement map
    if((param->lmax!=param->lcoarse)&&(level<param->lmax)){
      
      // cleaning the marks
      L_clean_marks(level,firstoct);
      // marking the cells of the current level
      L_mark_cells(level,param,firstoct,param->nrelax,param->amrthresh,cpu,sendbuffer,recvbuffer);
    }


    // ====================== VI Some bookkeeping ==========
    dt+=adt[level-1]; // advance local time
    tloc+=adt[level-1]; // advance local time
    is++;
    ndt[level-1]++;

    // Some Eye candy for timesteps display

    if(cpu->rank==0) dispndt(param,cpu,ndt);
    
    // === Loop

  }while((dt<adt[level-2])&&(is<param->nsubcycles));
  
  
  if(cpu->rank==0){
    printf("--\n");
    printf("exiting level =%d\n",level);
  }

  return dt;

}
