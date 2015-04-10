#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "prototypes.h"
#include "amr.h"
#include "oct.h"
#include "io.h"
#include "hydro_utils.h"
#include "poisson_utils.h"
#include "tools.h"
#ifdef WRAD
#include "rad_utils.h"
#include "src_utils.h"
#endif
#include "friedmann.h"
#include "cic.h"
#include "particle.h"
#include "communication.h"

#ifdef STARS
#include "stars.h"
#endif

#ifdef SUPERNOVAE
#include "supernovae.h"
#endif

#ifdef MOVIE
#include "movie.h"
#endif // MOVIE


 #ifdef CURIE
 #include "ccc_user.h"
 #endif // CURIE

// ===============================================================
// ===============================================================

void dispndt(struct RUNPARAMS *param, struct CPUINFO *cpu, int *ndt){

  int level;
  int ndtmax=POW(param->nsubcycles,param->lmax-param->lcoarse+1);
  int i,j,na;
  int nl;

  if(cpu->rank==RANK_DISP){
    printf("\n");
    for(j=0;j<ndtmax;j++)printf("#");
    printf("\n");
  }
  for(level=param->lcoarse;level<=param->lmax;level++){
    na=POW(2,param->lmax-level+1);
    nl=ndt[level-1];

#ifdef WMPI
    int nlmax;
    MPI_Allreduce(&nl,&nlmax,1,MPI_INT,MPI_MAX,cpu->comm);
    //printf("levef=%d rank=%d nl=%d nmax=%d\n",level,cpu->rank,nl,nlmax);
    nl=nlmax;
#endif

    if(cpu->rank==RANK_DISP){
      for(j=0;j<nl;j++){
      //one arrow
      for(i=0;i<(na-1);i++) printf("~");
      printf(">");
      }
      for(j=nl;j<ndtmax;j++) printf(" ");
      printf("\n");
    }
  }


  if(cpu->rank==RANK_DISP){
  for(j=0;j<ndtmax;j++)printf("#");
  printf("\n \n");
  }
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
  //  REAL amax=0.;

  //Smax=FMAX(Smax,SQRT(Warray[i].u*Warray[i].u+Warray[i].v*Warray[i].v+Warray[i].w*Warray[i].w)+Warray[i].a);
  // Computing new timestep
  dt=tmax;
  // setting the first oct

  nextoct=firstoct[level-1];

  if(nextoct!=NULL){
    dxcur=POW(0.5,level); // +1 to protect level change
    do // sweeping through the octs of level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	if(curoct->cpu!=cpu->rank) continue;

	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    vx=curoct->cell[icell].field.u;
	    vy=curoct->cell[icell].field.v;
	    vz=curoct->cell[icell].field.w;
	    va=SQRT(vx*vx+vy*vy+vz*vz);
	    aa=SQRT(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d);
	    //amax=FMAX(aa,amax);
	    //if(level==7) printf("v=%e\n",va+aa);
	    Smax=FMAX(Smax,va+aa);
	  }
      }while(nextoct!=NULL);
  }

  if(Smax>0.) dt=FMIN(dxcur*CFL/(Smax*3.),dt);
  //printf("Smax=%e amax=%e dt=%e\n",Smax,amax,dt);

  //if(level==7) printf("DT HY=%e\n",dt);
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
    dxcur=POW(0.5,level); // +1 to protect level change
    do // sweeping through the octs of level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	if(curoct->cpu!=cpu->rank) continue;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    dtloc=0.1*SQRT(2.*M_PI/(3.*(curoct->cell[icell].gdata.d+1.)*aexp));
 	    /* if(curoct->cell[icell].gdata.d<0.){ */
	    /*   printf("ouhla %e\n",dtloc); */
	    /*   abort(); */
	    /* } */
	    dt=FMIN(dt,dtloc);
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
#ifdef STARS
#ifdef ACCEL_RAD_STAR
  if((cpu->trigstar==0)){
    param->clight=1e-4;
  }
  else{
    param->clight=param->clightorg;
    //    printf("SWITCH VEL %e %e\n",param->clight,param->clightorg);
  }
#endif
#endif

  //nextoct=firstoct[level-1];

  //if(nextoct!=NULL){
    //dxcur=POW(0.5,level); // +1 to protect level change
    //dt=CFL*dxcur/(aexp*cfact*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v)/3.0; // UNITS OK
    dt=CFL/(aexp*param->clight*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v)/3.0/(REAL)(1<<level); // UNITS OK
    //}

    //    printf("dxcur=%e aexp=%e CFL=%e\n",1./(REAL)(1<<level),aexp,CFL);
    //abort();
  return dt;
}

#if 0
REAL L_comptstep_radfront(int level, struct RUNPARAMS *param,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax){

  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icell;
  REAL aa;
  REAL vf;
  REAL dt;
  REAL Smax;



  if(level==param->lcoarse){
    Smax=1e-3*aexp*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v; // initial velocity is c=1e-3
    // Computing new timestep
    //dt=tmax;

    dxcur=POW(0.5,level);

    // setting the first oct
    nextoct=firstoct[level-1];
    REAL clight=(aexp*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v);

    if(nextoct!=NULL){
      do // sweeping through the octs of level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu!=cpu->rank) continue;

	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {

	      REAL deltaX;
	      int ip=icell+(icell%2?-1:1);
	      REAL XP=curoct->cell[ip].rfield.nhplus/curoct->cell[ip].rfield.nh;
	      REAL X0=curoct->cell[icell].rfield.nhplus/curoct->cell[icell].rfield.nh;
	      if(XP<=X0) {
		deltaX=curoct->cell[ip].rfield.deltaX;
	      }
	      else{
		deltaX=curoct->cell[icell].rfield.deltaX;
	      }

	      REAL alpha= FABS(XP-X0)/dxcur;
	      REAL contrast= FABS(XP-X0)/(XP+X0)*2.;
	      if(contrast>0.1){
	      if((FABS(XP-0.5)<0.4)||(FABS(X0-0.5)<0.4))
		{
		  vf=deltaX/alpha*1.;

		  if(vf>clight) printf("vf/c=%3.1f deltaX=%e XP=%e X0=%e alpha=%e c=%e contrast=%e\n",vf/clight,deltaX,XP,X0,alpha,clight,contrast);
		  Smax=FMAX(Smax,FMIN(clight,vf));
		}
	    }

	  }
	}while(nextoct!=NULL);
    }

    dt=CFL*dxcur/(Smax*3.);
    Smax=Smax/(aexp*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v); // switch back to natural units for speed of light

#ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&Smax,1,MPI_REEL,MPI_MAX,cpu->comm);
    if(cpu->rank==RANK_DISP) printf("Smax=%f\n",Smax);
#endif
    param->clightorg=Smax;
  }


  if(NGRP>1) printf("ERROR NGRP>1\n");

  Smax=param->clightorg*(aexp*LIGHT_SPEED_IN_M_PER_S/param->unit.unit_v); // switch back to natural units for speed of light;
  // Check
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
    do // sweeping through the octs of level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	if(curoct->cpu!=cpu->rank) continue;

	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    struct Rtype *R;
	    R=&curoct->cell[icell].rfield;

	    REAL F=SQRT(POW(R->fx[0],2)+POW(R->fy[0],2)+POW(R->fz[0],2));
	    REAL E=R->e[0];

	    if(F>(Smax*E)) {
	      REAL corg=F/E;
	      R->fx[0]*=Smax/corg;
	      R->fy[0]*=Smax/corg;
	      R->fz[0]*=Smax/corg;
	      //printf("SOL VIOLATION ABORT F=%e E=%e C=%e CE=%e\n",F,E,Smax,Smax*E);
	      //abort();
	    }
	  }
      }while(nextoct!=NULL);
  }



  return dt;
}
#endif
#endif





// ===============================================================
// ===============================================================

REAL Advance_level(int level,REAL *adt, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct, struct HGRID *stencil, struct STENGRAV *gstencil, struct RGRID *rstencil,int *ndt, int nsteps, REAL tloc){

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
  int nsub=param->nsubcycles;
  int ip;
  //int *ptot = (int*)calloc(2,sizeof(int));
  int ptot[2];
  int deltan[2];
  ptot[0]=0;
  ptot[1]=0;

  REAL ekp,epp,ein;
  REAL RHS;
  REAL delta_e;
  REAL drift_e;
  REAL htilde;


#ifdef TESTCOSMO
  aexp=cosmo->aexp;
	param->cosmo->tphy	= a2t(param, aexp);
#else
  aexp=1.0;
#endif

  if(cpu->rank==RANK_DISP){
    printf("\n === entering level =%d with gstride=%d hstride=%d sten=%p aexp=%e adt=%e\n",level,gstride,hstride,stencil,aexp,adt[level-1]);
  }

  // ==================================== Check the number of particles and octs
  ptot[0]=0; for(ip=1;ip<=param->lmax;ip++){
    ptot[0]+=cpu->npart[ip-1]; // total of local particles

 //    if((level==12)&&(cpu->rank==217)){
 //     printf("l=%ip n=%ip\n",ip,cpu->npart[ip-1]);
 //    }
  }
#ifdef STARS
  ptot[1]=0; for(ip=1;ip<=param->lmax;ip++) ptot[1]+=cpu->nstar[ip-1];
#endif


  mtot=multicheck(firstoct,ptot,param->lcoarse,param->lmax,cpu->rank,cpu,param,0);

  if(level==param->lcoarse){
    // =========== Grid Census ========================================
    grid_census(param,cpu);
  }

  is=0;
#ifdef PIC
  //reset substep index of the particles
  L_reset_is_part(level,firstoct);
#endif



  // ========================== subcycling starts here ==================

  do{

    if(cpu->rank==RANK_DISP){
      printf("----\n");
      printf("subscyle #%d subt=%e nsub=%d\n",is,dt,nsub);
    }

#ifdef TESTCOSMO
    aexp=interp_aexp(tloc,(double *)cosmo->tab_aexp,(double *)cosmo->tab_ttilde);
    //aexp=cosmo->aexp;
#endif


    // == Ready to advance

  // ================= I we refine the current level
  if((param->lmax!=param->lcoarse)&&(level<param->lmax)){

    if(ndt[level-1]%2==0){
      // enforcing the 2 levels rule
      L_check_rule(level,param,firstoct,cpu);

#ifdef WMPI
      mpi_cic_correct(cpu, cpu->sendbuffer, cpu->recvbuffer, 3);
      mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,3,1,level); // propagate the rule check
#ifdef WHYDRO2
      mpi_exchange_hydro(cpu,cpu->hsendbuffer,cpu->hrecvbuffer,0); // propagate hydro for refinement
#endif
#ifdef WRAD
      mpi_exchange_rad_level(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,0,level); // propagate rad for refinement
#endif
#endif
      // refining (and destroying) octs
      curoct=L_refine_cells(level,param,firstoct,lastoct,cpu->freeoct,cpu,firstoct[0]+param->ngridmax,aexp);
      cpu->freeoct=curoct;
    }
    // ==================================== Check the number of particles and octs
    ptot[0]=0; for(ip=1;ip<=param->lmax;ip++){
      ptot[0]+=cpu->npart[ip-1]; // total of local particles
    }

#ifdef STARS
  ptot[1]=0; for(ip=1;ip<=param->lmax;ip++) ptot[1]+=cpu->nstar[ip-1];
#endif


    mtot=multicheck(firstoct,ptot,param->lcoarse,param->lmax,cpu->rank,cpu,param,1);

    ptot[0]=0; for(ip=1;ip<=param->lmax;ip++){
      ptot[0]+=cpu->npart[ip-1]; // total of local particles
      /* if((level==11)&&(cpu->rank==217)){ */
      /* 	printf("AP l=%ip n=%ip\n",ip,cpu->npart[ip-1]); */
      /* } */

    }

#ifdef WMPI
    //reset the setup in case of refinement
    setup_mpi(cpu,firstoct,param->lmax,param->lcoarse,param->ngridmax,1); // out of WMPI to compute the hash table
    MPI_Barrier(cpu->comm);
#endif
  }

    // =============================== cleaning
#ifdef WHYDRO2
    clean_new_hydro(level,param,firstoct,cpu);
#endif


#ifdef PIC
    L_clean_dens(level,param,firstoct,cpu);
#endif

#ifdef WRAD
    clean_new_rad(level,param,firstoct,cpu,aexp);
#endif




 // ================= IV advance solution at the current level


    if(cpu->rank==RANK_DISP) printf("ndt=%d nsteps=%d\n",ndt[param->lcoarse-1],nsteps);


#ifdef PIC
    // ==================================== performing the CIC assignement
    L_cic(level,firstoct,param,cpu);
    MPI_Barrier(cpu->comm);
    if(cpu->rank==RANK_DISP) printf("Local CIC done\n");
#ifdef WMPI
    mpi_cic_correct(cpu, cpu->sendbuffer, cpu->recvbuffer, 0);
    mpi_exchange(cpu,cpu->sendbuffer, cpu->recvbuffer,1,1);
#endif
    MPI_Barrier(cpu->comm);
    if(cpu->rank==RANK_DISP) printf("CIC done\n");
#endif

    /* //==================================== Getting Density ==================================== */
#ifdef WGRAV

    FillDens(level,param,firstoct,cpu);  // Here Hydro and Gravity are coupled

    /* //====================================  Poisson Solver ========================== */
    setOctList(firstoct[level-1], cpu, param,level);
    PoissonSolver(level,param,firstoct,cpu,gstencil,gstride,aexp);




    /* //====================================  Force Field ========================== */
#ifdef WMPI
    mpi_exchange(cpu,cpu->sendbuffer, cpu->recvbuffer,2,1);
#endif
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

    /* //====================================  I/O======= ========================== */
    // at this stage particles are synchronized at aexp
    // ready to dump particles
    // (note fields are dumped in quartz.c

#ifndef EDBERT

      int cond1 = nsteps%param->ndumps==0;
      int cond2 = 0;
      int cond3 = tloc>=param->time_max;

#ifdef TESTCOSMO

      if (param->dt_dump){
        cond1=0;

        int offset=0;
        if (nsteps==0) offset = (int)(param->cosmo->tphy/param->dt_dump);

        REAL a=param->cosmo->tphy;
        REAL b=(int)(*(cpu->ndumps)+offset)*param->dt_dump;
        cond2=a>b;
        if(cpu->rank==RANK_DISP)printf("t=%.2e yrs next dump at %.2e yrs\n",a,b+cond2*param->dt_dump);
      }
#endif // TESTCOSMO
    if(level==param->lcoarse){
if(cond1||cond2||cond3){
	if(cpu->rank==RANK_DISP) printf(" tsim=%e adt=%e\n",tloc,adt[level-1]);
	dumpIO(tloc,param,cpu,firstoct,adt,1);
	//dumpIO(tloc,param,cpu,firstoct,adt,0);
	//abort();
      }
    }
#endif



#if 0
    // =============== Computing Energy diagnostics

    ekp=0.;
    ein=0.;
#ifdef PIC
    //    if(level==param->lcoarse){
      egypart(cpu,&ekp,&epp,param,aexp);
#endif

      /* /\* /\\* // lets try to compute the potential from the grid *\\/ *\/ */
      epp=0.;
      REAL potloc=0., einloc=0., ekploc=0.;
      struct OCT *nextoct;
      int icell;
      REAL dx;
      int levelin;
      REAL u,v,w;
      //      if(cpu->rank==RANK_DISP) printf("get pop\n");

      for(levelin=param->lcoarse;levelin<=param->lmax;levelin++){
      	nextoct=firstoct[levelin-1];
      	dx=1./POW(2,levelin);
      	if(nextoct!=NULL){
      	  do // sweeping level
      	    {
      	      curoct=nextoct;
      	      nextoct=curoct->next;
      	      if(curoct->cpu!=cpu->rank) continue;
      	      for(icell=0;icell<8;icell++) // looping over cells in oct
      		{
      		  if(curoct->cell[icell].child==NULL){
      		    potloc+=aexp*dx*dx*dx*(curoct->cell[icell].gdata.d)*(curoct->cell[icell].gdata.p)*0.5;
#ifdef WHYDRO2
		    einloc+=dx*dx*dx*(curoct->cell[icell].field.p)/(GAMMA-1.);

		    u=curoct->cell[icell].field.u;
		    v=curoct->cell[icell].field.v;
		    w=curoct->cell[icell].field.w;

		    ekploc+=dx*dx*dx*(curoct->cell[icell].field.d)*(u*u+v*v+w*w)*0.5;
#endif
      		  }
      		}
      	    }while(nextoct!=NULL);
      	}
      }

      epp=potloc;
      ekp=ekp+ekploc;
      ein=ein+einloc;
#ifdef WMPI
      REAL sum_ekp,sum_epp,sum_ein;
      MPI_Allreduce(&ekp,&sum_ekp,1,MPI_REEL,MPI_SUM,cpu->comm);
      MPI_Allreduce(&epp,&sum_epp,1,MPI_REEL,MPI_SUM,cpu->comm);
      MPI_Allreduce(&ein,&sum_ein,1,MPI_REEL,MPI_SUM,cpu->comm);
      ekp=sum_ekp;
      epp=sum_epp;
      ein=sum_ein;
#endif


      if(nsteps==1){
#ifdef TESTCOSMO
	htilde=2./SQRT(param->cosmo->om)/faexp_tilde(aexp,param->cosmo->om,param->cosmo->ov)/aexp;
	param->egy_last=epp*htilde;
#else
	param->egy_last=0.;
	htilde=0.;
#endif
	param->egy_rhs=0.;
  	param->egy_0=ekp+epp+ein;
	param->egy_timelast=aexp;
	param->egy_totlast=ekp+epp+ein;
      }


    if((level>=param->lcoarse)&&(nsteps>1)) {
#ifdef TESTCOSMO
      htilde=2./SQRT(param->cosmo->om)/faexp_tilde(aexp,param->cosmo->om,param->cosmo->ov)/aexp;
      RHS=param->egy_rhs;

      //RHS=RHS+0.5*(epp*htilde+param->egy_last)*(tloc-param->egy_timelast); // trapezoidal rule
      RHS=RHS+0.5*(epp/aexp+param->egy_last)*(aexp-param->egy_timelast); // trapezoidal rule

      //delta_e=(((ekp+epp)-param->egy_totlast)/(0.5*(epp*htilde+param->egy_last)*(tloc-param->egy_timelast))-1.);
      delta_e=(ekp+epp+ein-param->egy_totlast-0.5*(epp/aexp+param->egy_last)*(aexp-param->egy_timelast));
      drift_e=(((ekp+epp+ein)-param->egy_0)/RHS-1.);
      //      param->egy_last=epp*htilde;
      param->egy_last=epp/aexp;
#else
      RHS=0.;
      delta_e=((ekp+epp+ein)-param->egy_totlast)/param->egy_totlast;
      drift_e=((ekp+epp+ein)-param->egy_0)/param->egy_0;
#endif
      param->egy_rhs=RHS;
      param->egy_timelast=aexp;
      param->egy_totlast=ekp+epp+ein;

      if(cpu->rank==RANK_DISP){
	FILE *fpe;
	fpe=fopen("energystat.txt","a");
	fprintf(fpe,"%e %e %e %e %e %e %e %e %e %e %d\n",aexp,delta_e,drift_e,ekp,epp,RHS,ein,adt[level-1],param->egy_0,htilde,level);
	fclose(fpe);
	printf("Egystat rel. err= %e drift=%e\n",delta_e,drift_e);
      }
    }

#endif

    // ================= II We compute the timestep of the current level

    REAL adtold=adt[level-1]; // for energy conservation
    dtnew=param->dt;//*(cpu->nsteps>2);


    // Overshoot tmax
    dtnew=((param->time_max-tloc)<dtnew?(param->time_max-tloc):dtnew);
    if(cpu->rank==RANK_DISP){    printf("dtnew=%e %e %e",dtnew,param->time_max,tloc);
      printf("aexp=%e tloc=%e param->tmax=%e dtnew=%e ",aexp,tloc,param->time_max,dtnew);
    }

    // Free Fall
#ifdef WGRAV
    REAL dtff;
#ifdef TESTCOSMO
    dtff=L_comptstep_ff(level,param,firstoct,aexp,cpu,1e9);
#else
    dtff=L_comptstep_ff(level,param,firstoct,1.0,cpu,1e9);
#endif
    dtnew=(dtff<dtnew?dtff:dtnew);
    if(cpu->rank==RANK_DISP) printf("dtff= %e ",dtff);
#endif


    // Cosmo Expansion
#ifdef TESTCOSMO
    REAL dtcosmo;
    dtcosmo=-0.5*SQRT(cosmo->om)*integ_da_dt_tilde(aexp*1.02,aexp,cosmo->om,cosmo->ov,1e-8);
    dtnew=(dtcosmo<dtnew?dtcosmo:dtnew);
    if(cpu->rank==RANK_DISP) printf("dtcosmo= %e ",dtcosmo);
#endif




    // Courant Condition Hydro
#ifdef WHYDRO2
    REAL dthydro;
    dthydro=L_comptstep_hydro(level,param,firstoct,1.0,1.0,cpu,1e9);
    if(cpu->rank==RANK_DISP) printf("dthydro= %e ",dthydro);
    dtnew=(dthydro<dtnew?dthydro:dtnew);
#endif

    // Courant Condition Particle
#ifdef PIC
    REAL dtpic;
    dtpic=L_comptstep(level,param,firstoct,1.0,1.0,cpu,1e9);
    if(cpu->rank==RANK_DISP) printf("dtpic= %e ",dtpic);
    dtnew=(dtpic<dtnew?dtpic:dtnew);
#endif

    // Courant Condition Radiation
#ifdef WRAD
    REAL dtrad;
    dtrad=L_comptstep_rad(level,param,firstoct,aexp,cpu,1e9);
    if(cpu->rank==RANK_DISP) printf("dtnew=%e dtrad= %e\n",dtnew,dtrad);
#ifdef RADSTEP
    dtnew=(dtrad<dtnew?dtrad:dtnew);
#endif

#ifdef WRADTEST
#ifndef TESTCLUMP
#ifndef WRADHYD
    dtnew=(dtrad<dtnew?dtrad:dtnew);
#endif
#endif
#endif

#endif

#ifdef FLOORDT
    // REALLY WEIRD ===
    // Apparently there are some truncation errors in REDUCTION operation on double
    // that makes multi-processor dt different than the mono processor ones !
    // the few lines below removes this...

    REAL dtf=floor(dtnew*1e10)/1e10;
    dtnew=dtf;
#endif

    // SOME hack to force refinement before integration
    /* if(nsteps==0){ */
    /*   printf("skip udpate\n"); */
    /*   dtnew=0.; */
    /* } */


    //printf("dtnew %e before\n",dtnew);

    /// ================= Assigning a new timestep for the current level
    dtold=adt[level-1];
    adt[level-1]=dtnew;
    if(dtnew==0.){
      printf("WARNING NULL dt rank %d n %d dt %e\n",cpu->rank,cpu->noct[level-1],dtnew);
      //abort();
    }

    if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one

    // the coarser level may be shorter than the finer one
    if(adt[level-2]<adt[level-1]){
      adt[level-1]=adt[level-2];
    }
    else{ // otherwise standard case
      adt[level-1]=FMIN(adt[level-1],adt[level-2]-dt);// we force synchronization
    }

#ifdef WMPI
    REAL tdum=0.;
    REAL tdum2=0.;
    MPI_Allreduce(adt+level-1,&tdum,1,MPI_REEL,MPI_MIN,cpu->comm);
    MPI_Allreduce(adt+level-2,&tdum2,1,MPI_REEL,MPI_MIN,cpu->comm);
    adt[level-1]=tdum;
    adt[level-2]=tdum2;
#endif

    //printf("dtnew %e\n",tdum);


   // ================= III Recursive call to finer level

  double tt2,tt1;
    MPI_Barrier(cpu->comm);
    tt1=MPI_Wtime();

    int nlevel=cpu->noct[level]; // number at next level
#ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&nlevel,1,MPI_INT,MPI_SUM,cpu->comm);
#endif

#if 1
    if(level<param->lmax){
      if(nlevel>0){
	dtfine=Advance_level(level+1,adt,cpu,param,firstoct,lastoct,stencil,gstencil,rstencil,ndt,nsteps,tloc);
	// coarse and finer level must be synchronized now
	adt[level-1]=dtfine;
	if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one
      }
    }
#endif
#ifdef TESTCOSMO
	param->cosmo->tphy	= a2t(param, aexp);
#endif

#ifdef WMPI
    tdum=0.;
    tdum2=0.;
    MPI_Allreduce(adt+level-1,&tdum,1,MPI_REEL,MPI_MIN,cpu->comm);
    MPI_Allreduce(adt+level-2,&tdum2,1,MPI_REEL,MPI_MIN,cpu->comm);
    adt[level-1]=tdum;
    adt[level-2]=tdum2;
#endif

    MPI_Barrier(cpu->comm);
    tt2=MPI_Wtime();

    /* if(param->lcoarse==level){ */
    /*   if(cpu->rank==RANK_DISP){ */
    /* 	int lev; */
    /* 	for(lev=7;lev<=10;lev++){ */
    /* 	  printf("lev %d dt=%e\n",lev,adt[lev-1]); */
    /* 	} */
    /*   } */
    /* } */

    if(cpu->rank==RANK_DISP) printf("SUBLEVEL TIME lev=%d %e\n",level,tt2-tt1);
    // ==================================== Check the number of particles and octs
    ptot[0]=0; for(ip=1;ip<=param->lmax;ip++) ptot[0]+=cpu->npart[ip-1]; // total of local particles
#ifdef STARS
    ptot[1]=0; for(ip=1;ip<=param->lmax;ip++) ptot[1]+=cpu->nstar[ip-1];
#endif
    mtot=multicheck(firstoct,ptot,param->lcoarse,param->lmax,cpu->rank,cpu,param,2);



#ifdef PIC
    //================ Part. Update ===========================
#ifdef WMPI
    int maxnpart;
    MPI_Allreduce(cpu->npart+level-1,&maxnpart,1,MPI_INT,MPI_SUM,cpu->comm);
#endif
    if(cpu->rank==RANK_DISP) printf("Start PIC on %d part with dt=%e on level %d\n",maxnpart,adt[level-1],level);

    //printf("1 next=%p on proc=%d\n",firstoct[0]->next,cpu->rank);

    // moving particles around

    // predictor step
    L_accelpart(level,firstoct,adt,is,cpu); // computing the particle acceleration and velocity

#ifndef PARTN
    L_movepart(level,firstoct,adt,is,cpu); // moving the particles
#endif

    L_partcellreorg(level,firstoct); // reorganizing the particles of the level throughout the mesh

#ifdef WMPI
    //reset the setup in case of refinement
    //printf("2 next=%p on proc=%d\n",firstoct[0]->next,cpu->rank);
    setup_mpi(cpu,firstoct,param->lmax,param->lcoarse,param->ngridmax,1); // out of WMPI to compute the hash table
    MPI_Barrier(cpu->comm);
#endif

#ifdef WMPI

    //printf("3 next=%p on proc=%d\n",firstoct[0]->next,cpu->rank);
    cpu->firstoct = firstoct;
    int deltan[2];
    mpi_exchange_part(cpu, cpu->psendbuffer, cpu->precvbuffer,deltan,level);

    //printf("4 next=%p on proc=%d\n",firstoct[0]->next,cpu->rank);
    //printf("proc %d receives %d particles %d stars freepart=%p\n",cpu->rank,deltan[0],deltan[1],cpu->freepart);
    //update the particle number within this process
    //npart=npart+deltan;

    ptot[0]=deltan[0]; for(ip=1;ip<=param->lmax;ip++) ptot[0]+=cpu->npart[ip-1]; // total of local particles
#ifdef STARS
    ptot[1]=deltan[1]; for(ip=1;ip<=param->lmax;ip++) ptot[1]+=cpu->nstar[ip-1]; // total of local stars
#endif


    mtot=multicheck(firstoct,ptot,param->lcoarse,param->lmax,cpu->rank,cpu,param,3);
#endif


#endif




    //=============== Hydro Update ======================
#ifdef WHYDRO2
    double th[10];
    MPI_Barrier(cpu->comm);
    th[0]=MPI_Wtime();
#ifdef WMPI
    mpi_exchange_hydro_level(cpu,cpu->hsendbuffer,cpu->hrecvbuffer,1,level);
#endif

    MPI_Barrier(cpu->comm);
    th[1]=MPI_Wtime();
    HydroSolver(level,param,firstoct,cpu,stencil,hstride,adt[level-1]);

    MPI_Barrier(cpu->comm);
    th[2]=MPI_Wtime();



#ifdef WGRAV
    // ================================= gravitational correction for Hydro
    grav_correction(level,param,firstoct,cpu,adt[level-1]); // Here Hydro and Gravity are coupled
#endif

    MPI_Barrier(cpu->comm);
    th[3]=MPI_Wtime();

#ifdef WMPI
    if(level>param->lcoarse){
      mpi_hydro_correct(cpu,cpu->hsendbuffer,cpu->hrecvbuffer,level);

#ifdef WMPI
    MPI_Barrier(cpu->comm);
    mpi_exchange_hydro_level(cpu,cpu->hsendbuffer,cpu->hrecvbuffer,1,level);
#endif

    }
    MPI_Barrier(cpu->comm);
#endif

    MPI_Barrier(cpu->comm);
    th[4]=MPI_Wtime();

    if(cpu->rank==RANK_DISP) printf("HYD -- Ex=%e HS=%e GCorr=%e HCorr=%e Tot=%e\n",th[1]-th[0],th[2]-th[1],th[3]-th[2],th[4]-th[3],th[4]-th[0]);

#endif


    // ===================================== RADIATION
#ifdef WRAD
    double tcomp[10];
    int chemonly;
    MPI_Barrier(cpu->comm);
    tcomp[0]=MPI_Wtime();
    //=============== Building Sources and counting them ======================
    nsource=FillRad(level,param,firstoct,cpu,(level==param->lcoarse)&&(nsteps==0),aexp, tloc);  // Computing source distribution and filling the radiation fields

#ifndef COARSERAD
    chemonly=0;
#else
    chemonly=1;
    if(cpu->rank==RANK_DISP) printf("Dynamical Cooling\n");
#endif

#ifdef WMPI
    //printf("cpu #%d ready 1\n",cpu->rank);
    MPI_Barrier(cpu->comm);
    tcomp[1]=MPI_Wtime();
    mpi_exchange_rad_level(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,1,level);
    //printf("cpu #%d ready 2\n",cpu->rank);
    MPI_Barrier(cpu->comm);
    tcomp[2]=MPI_Wtime();
#endif

#ifdef WMPI
    MPI_Barrier(cpu->comm);
    //printf("cpu #%d ready 3\n",cpu->rank);
#endif

    //=============== Radiation Update ======================
#ifndef COARSERAD
    REAL time_rad;
    time_rad=RadSolver(level,param,firstoct,cpu,rstencil,hstride,adt[level-1],aexp,chemonly);
    if(cpu->rank==RANK_DISP){
#ifndef GPUAXL
      printf("==== CPU RAD TOTAL TIME =%e\n",time_rad);
#else
      printf(" === GPU RAD TOTAL TIME =%e\n",time_rad);
#endif
    }
#endif

    //printf("cpu #%d ready 4\n",cpu->rank);

#ifdef WMPI
    MPI_Barrier(cpu->comm);
#endif

#ifdef WMPI
    //printf("proc %d waiting\n",cpu->rank);
    MPI_Barrier(cpu->comm);
    tcomp[3]=MPI_Wtime();
    if(level>param->lcoarse){
      mpi_rad_correct(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,level);

      MPI_Barrier(cpu->comm);
      tcomp[5]=MPI_Wtime();
      mpi_exchange_rad_level(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,1,level);

    }
    MPI_Barrier(cpu->comm);
    tcomp[4]=MPI_Wtime();
#endif

    if(cpu->rank==RANK_DISP) printf("RAD -- Fill=%e Ex=%e RS=%e Corr=%e Tot=%e\n",tcomp[1]-tcomp[0],tcomp[2]-tcomp[1],tcomp[3]-tcomp[2],tcomp[4]-tcomp[3],tcomp[5]-tcomp[0]);


#ifdef WRADHYD
    if(cpu->rank==RANK_DISP) printf("TRADHYD l=%d Total=%e Noct=%d\n",level,tcomp[5]-th[0],cpu->noct[level-1]);
#endif

#endif

    /* //===================================creating new stars=================================// */

#ifdef STARS
#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI
#ifdef ZOOM
    if(level>=param->lmaxzoom)
#endif //ZOOM
      {
	Stars(param,cpu, adt[level-1], aexp, level, is);
      }
#endif // STARS

    /* //===================================Supernovae=========================================// */

#ifdef SUPERNOVAE
#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI
#ifndef SNTEST
    supernovae(param,cpu, adt[level-1], aexp, level, is);
#else //ifdef SNTEST
    supernovae(param,cpu, adt[level-1], tloc, level, is);
#endif // SNTEST
#endif // SUPERNOVAE

  // ================= V Computing the new refinement map
    REAL dxkpc=0.;
    REAL dxnext=POW(0.5,level+1)*aexp;
#ifdef KPCLIMIT
    dxkpc=1e2*PARSEC/param->unit.unit_l;
#endif

    if(dxnext>dxkpc){ // ENFORCE Kennicut scale
      if((param->lmax!=param->lcoarse)&&(level<param->lmax)){

#ifndef ZOOM
	if((ndt[level-1]%2==1)||(level==param->lcoarse)){
#else
	  if((ndt[level-1]%2==1)||(level>=param->lmaxzoom)){
#endif
	    L_clean_marks(level,firstoct);
	    // marking the cells of the current level
	    L_mark_cells(level,param,firstoct,param->nsmooth,param->amrthresh,cpu,cpu->sendbuffer,cpu->recvbuffer);
	  }
	}
    }
    else{
      //if(cpu->rank==RANK_DISP)
	printf("Blocking refinement to level %d : dx[%d]=%e dxkpc=%e\n",level+1,level+1,dxnext,dxkpc);
    }

    // ====================== VI Some bookkeeping ==========
    dt+=adt[level-1]; // advance local time
    tloc+=adt[level-1]; // advance local time
    is++;
    ndt[level-1]++;

    if(cpu->rank==RANK_DISP) printf("is=%d ndt=%d dt=%le tloc=%le %le\n",is,ndt[level-1],dt,tloc,adt[level-2]);

    // Some Eye candy for timesteps display

    //dispndt(param,cpu,ndt);

    // === Loop


#ifdef MOVIE
	if (level==param->lcoarse)	dumpMovie(firstoct, param, cpu, level, (float)aexp);
#endif

  }while((dt<adt[level-2])&&(is<nsub));


  if(cpu->rank==RANK_DISP){
    printf("--\n");
    printf("exiting level =%d\n",level);
  }


  return dt;

}

  // ============================================================================================================
  // ============================================================================================================

#ifdef WRAD
#ifdef COARSERAD

  REAL Advance_level_RAD(int level,REAL dtmax, REAL *adt, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct, struct HGRID *stencil, struct STENGRAV *gstencil, struct RGRID *rstencil, int nsteps, REAL tloc, int nrad){

#ifdef TESTCOSMO
  struct COSMOPARAM *cosmo;
  cosmo=param->cosmo;
#endif
   struct OCT *curoct;
  REAL dtnew;
  REAL dt=0.;
  int npart=0;
  int mtot;
  int is;
  REAL aexp;
  int nsource;
  int hstride=param->hstride;
  int gstride=param->gstride;
  int nsub=param->nsubcycles;
  int ip;
  int ptot[2];
  int deltan[2];
  ptot[0]=0;
  ptot[1]=0;

  REAL ekp,epp,ein;
  REAL RHS;
  REAL delta_e;
  REAL drift_e;
  REAL htilde;

  is=0;

  do{

    /* if(cpu->rank==RANK_DISP){ */
    /*   printf("----\n"); */
    /*   printf("subscyle #%d subt=%e nsub=%d\n",is,dt,nsub); */
    /* } */

#ifdef TESTCOSMO
    aexp=interp_aexp(tloc,(double *)cosmo->tab_aexp,(double *)cosmo->tab_ttilde);
    //aexp=cosmo->aexp;
#else
    aexp=1.0;
#endif



    // == Ready to advance

#ifdef WRAD
    clean_new_rad(level,param,firstoct,cpu,aexp);
#endif

    // ================= II We compute the timestep of the current level

    dtnew=param->dt;//*(cpu->nsteps>2);


    // Overshoot tmax
    dtnew=((param->time_max-tloc)<dtnew?(param->time_max-tloc):dtnew);


    // Courant Condition Radiation
    REAL dtrad;
    dtrad=L_comptstep_rad(param->lcoarse,param,firstoct,aexp,cpu,1e9);
    //if(cpu->rank==RANK_DISP) printf("dtmax=%e dtrad= %e\n",dtmax,dtrad);

    dtnew=(dtrad<dtmax?dtrad:dtmax); // we cannot go further than the outer timestep


    /// ================= Assigning a new timestep for the current level
    adt[level-1]=dtnew;

    if(dtnew==0.){
      printf("WARNING NULL dt rank %d n %d dt %e\n",cpu->rank,cpu->noct[level-1],dtnew);
      //abort();
    }

    if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one

    // the coarser level may be shorter than the finer one
    if(adt[level-2]<adt[level-1]){
      adt[level-1]=adt[level-2];
    }
    else{ // otherwise standard case
      adt[level-1]=FMIN(adt[level-1],adt[level-2]-dt);// we force synchronization
    }

#ifdef WMPI
    REAL tdum=0.;
    REAL tdum2=0.;
    MPI_Allreduce(adt+level-1,&tdum,1,MPI_REEL,MPI_MIN,cpu->comm);
    MPI_Allreduce(adt+level-2,&tdum2,1,MPI_REEL,MPI_MIN,cpu->comm);
    adt[level-1]=tdum;
    adt[level-2]=tdum2;
#endif

    // ================= III Recursive call to finer level

    double tt2,tt1;
    MPI_Barrier(cpu->comm);
    tt1=MPI_Wtime();

    int nlevel=cpu->noct[level]; // number at next level
#ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&nlevel,1,MPI_INT,MPI_SUM,cpu->comm);
#endif


    if(level<param->lmax){
      if(nlevel>0){
	REAL dtfine;
	dtfine=Advance_level_RAD(level+1,dtmax,adt,cpu,param,firstoct,lastoct,stencil,gstencil,rstencil,nsteps,tloc,nrad);
	// coarse and finer level must be synchronized now
	adt[level-1]=dtfine;
	if(level==param->lcoarse) adt[level-2]=adt[level-1]; // we synchronize coarser levels with the coarse one
      }
    }


#ifdef WMPI
    tdum=0.;
    tdum2=0.;
    MPI_Allreduce(adt+level-1,&tdum,1,MPI_REEL,MPI_MIN,cpu->comm);
    MPI_Allreduce(adt+level-2,&tdum2,1,MPI_REEL,MPI_MIN,cpu->comm);
    adt[level-1]=tdum;
    adt[level-2]=tdum2;
#endif

    MPI_Barrier(cpu->comm);
    tt2=MPI_Wtime();

    /* if(param->lcoarse==level){ */
    /*   if(cpu->rank==RANK_DISP){ */
    /* 	int lev; */
    /* 	for(lev=7;lev<=10;lev++){ */
    /* 	  printf("lev %d dt=%e\n",lev,adt[lev-1]); */
    /* 	} */
    /*   } */
    /* } */

    // ===================================== RADIATION

    double tcomp[10];
    MPI_Barrier(cpu->comm);
    tcomp[0]=MPI_Wtime();

    //=============== Building Sources and counting them ======================



    nsource=FillRad(level,param,firstoct,cpu,(level==param->lcoarse)&&(nsteps==0),aexp, tloc);  // Computing source distribution and filling the radiation fields

#ifdef HOMOSOURCE
    homosource(param,firstoct,cpu,level); // FOR HOMOGENOUS SRC
#endif

#ifdef WMPI
      MPI_Barrier(cpu->comm);
      tcomp[1]=MPI_Wtime();
      mpi_exchange_rad_level(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,1,level);
      MPI_Barrier(cpu->comm);
      tcomp[2]=MPI_Wtime();
#endif

#ifdef WMPI
      MPI_Barrier(cpu->comm);
#endif

      //=============== Radiation Update ======================
      int chemonly=0;
      REAL time_rad;
      time_rad=RadSolver(level,param,firstoct,cpu,rstencil,hstride,adt[level-1],aexp,chemonly);

      if(nrad%10==0){
	if(cpu->rank==RANK_DISP){
#ifndef GPUAXL
	  printf("==== CPU RAD TOTAL TIME =%e\n",time_rad);
#else
	  printf(" === GPU RAD TOTAL TIME =%e\n",time_rad);
#endif
	}
      }

#ifdef CURIE
REAL time_remain;
int error;
error = ccc_tremain(&time_remain)
 if (!error) printf("Time remaining :%lf\n", time_remain);
#endif // CURIE


#ifdef WMPI
      MPI_Barrier(cpu->comm);
#endif

#ifdef WMPI
      //printf("proc %d waiting\n",cpu->rank);
      MPI_Barrier(cpu->comm);
      tcomp[3]=MPI_Wtime();
      if(level>param->lcoarse){
	mpi_rad_correct(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,level);

	MPI_Barrier(cpu->comm);
	tcomp[5]=MPI_Wtime();
	mpi_exchange_rad_level(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,1,level);

      }
      MPI_Barrier(cpu->comm);
      tcomp[4]=MPI_Wtime();

#endif


    // ====================== VI Some bookkeeping ==========
    dt+=adt[level-1]; // advance local time
    tloc+=adt[level-1]; // advance local time
    is++;

    //if(cpu->rank==RANK_DISP) printf("is=%d dt=%le tloc=%le %le\n",is,dt,tloc,adt[level-2]);

    // Some Eye candy for timesteps display

    //dispndt(param,cpu,ndt);

    // === Loop


  }while((dt<adt[level-2])&&(is<nsub));


  /* if(cpu->rank==RANK_DISP){ */
  /*   printf("--\n"); */
  /*   printf("RAD exiting level =%d\n",level); */
  /* } */

  return dt;

}
#endif
#endif
