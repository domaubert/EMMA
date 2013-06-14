#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"
#include "amr.h"
#include "hydro_utils.h"
#ifdef WRAD
#include "rad_utils.h"
#include "src_utils.h"
#endif
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

#ifdef WRAD
REAL L_comptstep_rad(int level, struct RUNPARAMS *param,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax){
  
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
    dt=CFL*dxcur/(aexp*param->clight*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v)/3.0; // UNITS OK
  }

  return dt;
}
#endif





// ===============================================================
// ===============================================================

REAL Advance_level(int level,REAL *adt, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct, struct HGRID *stencil, struct STENGRAV *gstencil, struct RGRID *rstencil,struct PACKET **sendbuffer, struct PACKET **recvbuffer,int *ndt, int nsteps,REAL tloc){
 
#ifdef TESTCOSMO
  struct COSMOPARAM *cosmo;
  cosmo=param->cosmo;
#endif
  struct OCT *curoct;
  REAL dtnew;
  REAL dt=0.;
  REAL dtold;
  REAL dtvel;
  REAL dtfine;
  int npart=0;
  int mtot;
  int is;
  REAL aexp;
  int nsource;
  int hstride=param->hstride;
  int gstride=param->gstride;

#ifdef TESTCOSMO
  aexp=cosmo->aexp;
#else
  aexp=1.0;
#endif

  if(cpu->rank==0){
    printf("\n === entering level =%d with gstride=%d hstride=%d sten=%p aexp=%e\n",level,gstride,hstride,stencil,aexp);
  }


  // ==================================== Check the number of particles and octs
  mtot=multicheck(firstoct,npart,param->lcoarse,param->lmax,cpu->rank,cpu);

  if(level==param->lcoarse){
    // =========== Grid Census ========================================
    grid_census(param,cpu);
  }
  
  is=0;
#ifdef PIC
  //reset substep index of the particles
  L_reset_is_part(level,firstoct);
#endif


  do{
    printf("----\n");
    printf("subscyle #%d subt=%e\n",is,dt);
    
#ifdef TESTCOSMO
    aexp=interp_aexp(tloc,cosmo->tab_aexp,cosmo->tab_ttilde);
    //aexp=cosmo->aexp;
#endif


    
    // =============================== cleaning 
#ifdef WHYDRO2
    clean_new_hydro(level,param,firstoct,cpu);
#endif


#ifdef PIC
    L_clean_dens(level,param,firstoct,cpu);
#endif

#ifdef WRAD
    //sanity_rad(level,param,firstoct,cpu,aexp);
    clean_new_rad(level,param,firstoct,cpu,aexp);
#endif

    // == Ready to advance


  // ================= I we refine the current level
  if((param->lmax!=param->lcoarse)&&(level<param->lmax)){
    // refining (and destroying) octs
#ifdef WRAD
    //sanity_rad(level,param,firstoct,cpu,aexp);
#endif
    curoct=L_refine_cells(level,param,firstoct,lastoct,cpu->freeoct,cpu,firstoct[0]+param->ngridmax,aexp);
    cpu->freeoct=curoct;
  }

  // ==================================== Check the number of particles and octs
  mtot=multicheck(firstoct,npart,param->lcoarse,param->lmax,cpu->rank,cpu);

  // ================= IV advance solution at the current level


    printf("ndt=%d nsteps=%d\n",ndt[param->lcoarse-1],nsteps);


#ifdef PIC
    // ==================================== performing the CIC assignement
    L_cic(level,firstoct,param,cpu);
#ifdef WMPI
    mpi_cic_correct(cpu, cpu->sendbuffer, cpu->recvbuffer, 0);

    mpi_exchange(cpu,cpu->sendbuffer, cpu->recvbuffer,1,1);
    
#endif
#endif
 
    /* //==================================== Getting Density ==================================== */
#ifdef WGRAV
    FillDens(level,param,firstoct,cpu);  // Here Hydro and Gravity are coupled

    /* //====================================  Poisson Solver ========================== */
    PoissonSolver(level,param,firstoct,cpu,gstencil,gstride,aexp); 

    /* //====================================  Force Field ========================== */
    PoissonForce(level,param,firstoct,cpu,gstencil,gstride,aexp);
#endif

#ifdef PIC
    // Mid point rule correction step
    if((is>0)||(cpu->nsteps>0)){
      L_accelpart(level,firstoct,adt,-1,cpu); // computing the particle acceleration and velocity
      // here -1 forces the correction : all particles must be corrected
    }
    L_levpart(level,firstoct,is); // assigning all the particles to the current level
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

#ifdef WRAD
    REAL dtrad;
    dtrad=L_comptstep_rad(level,param,firstoct,aexp,cpu,1e9);
    printf("dtrad= %e ",dtrad);
    dtnew=(dtrad<dtnew?dtrad:dtnew);
#endif



    printf("\n");

#ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&dtnew,1,MPI_DOUBLE,MPI_MIN,cpu->comm);
#endif

    /// ================= Assigning a new timestep for the current level
    dtold=adt[level-1];
    adt[level-1]=dtnew;
    printf("inital dtnew=%e\n",dtnew);

    if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one

    // the coarser level may be shorter than the finer one
    if(adt[level-2]<adt[level-1]){
      adt[level-1]=adt[level-2];
    }
    else{ // otherwise standard case
      adt[level-1]=fmin(adt[level-1],adt[level-2]-dt);// we force synchronization
    }

   // ================= III Recursive call to finer level
    if(level<param->lmax){
      if(cpu->noct[level]>0){
	dtfine=Advance_level(level+1,adt,cpu,param,firstoct,lastoct,stencil,gstencil,rstencil,sendbuffer,recvbuffer,ndt,nsteps,tloc);
	// coarse and finer level must be synchronized now
	adt[level-1]=dtfine;
	if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one
      }
    }

    // ==================================== Check the number of particles and octs
    mtot=multicheck(firstoct,npart,param->lcoarse,param->lmax,cpu->rank,cpu);

#ifdef PIC
    //================ Part. Update ===========================

    if(cpu->rank==0) printf("Start PIC on %d part with dt=%e on level %d\n",cpu->npart[level-1],adt[level-1],level);

#ifdef PART_EGY
    //computing energy
    L_egypart(level,firstoct); // computing the particle acceleration and velocity
#endif

    // predictor step
    L_accelpart(level,firstoct,adt,is,cpu); // computing the particle acceleration and velocity

#ifndef PARTN
    L_movepart(level,firstoct,adt,is,cpu); // moving the particles
#endif

    L_partcellreorg(level,firstoct); // reorganizing the particles of the level throughout the mesh

#ifdef WMPI
    int deltan;
    deltan=mpi_exchange_part(cpu, cpu->psendbuffer, cpu->precvbuffer, &(cpu->lastpart));

    //printf("proc %d receives %d particles\n",cpu->rank,deltan);
    //update the particle number within this process
    npart=npart+deltan;
    
    mtot=multicheck(firstoct,npart,param->lcoarse,param->lmax,cpu->rank,cpu);
#endif


#endif




#ifdef WHYDRO2

#ifdef WMPI
    mpi_exchange_hydro(cpu,cpu->hsendbuffer,cpu->hrecvbuffer,1);
#endif

      //=============== Hydro Update ======================
    HydroSolver(level,param,firstoct,cpu,stencil,hstride,adt[level-1]);

      // ================================= gravitational correction for Hydro

#ifdef WGRAV
    grav_correction(level,param,firstoct,cpu,adt[level-1]); // Here Hydro and Gravity are coupled
#endif
#endif


#ifdef WRAD
    //=============== Building Sources and counting them ======================
    nsource=FillRad(level,param,firstoct,cpu,(level==param->lcoarse)&&(nsteps==0),aexp);  // Computing source distribution and filling the radiation fields
 
#ifdef WMPI
    mpi_exchange_rad(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,1);
#endif
    //=============== Radiation Update ======================
    RadSolver(level,param,firstoct,cpu,rstencil,hstride,adt[level-1],aexp);
    //sanity_rad(level,param,firstoct,cpu,aexp);
#endif


    //}
  // ================= V Computing the new refinement map
    if((param->lmax!=param->lcoarse)&&(level<param->lmax)){
      
      // cleaning the marks
      L_clean_marks(level,firstoct);
      // marking the cells of the current level
      L_mark_cells(level,param,firstoct,2,param->amrthresh,cpu,sendbuffer,recvbuffer);
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
