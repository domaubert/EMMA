// ----------------------------------------------------------
// ----------------------------------------------------------
/** \file stars.c
  * \brief contain the stars functions
  * \author Nicolas Deparis
  *
  * Need the STARS preprocessor flag
  *
  * Documentation on the implementation can be found at :
  * https://github.com/domaubert/EMMA/wiki/Stars
  */
// ----------------------------------------------------------
// ----------------------------------------------------------

#ifdef STARS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prototypes.h"
#include "particle.h" //findlastpart
#include "hydro_utils.h" // W2U and U2W
#include "tools.h" // rdm and gpoiss

#ifdef WMPI
#include <mpi.h>
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initStar(struct Wtype *field, struct PART *star, struct RUNPARAMS *param, int level, REAL xc, REAL yc, REAL zc,int idx, REAL aexp, int is, REAL dt,REAL dx, REAL mstar) {
// ----------------------------------------------------------//
/// compute the initial state of star particle
// ----------------------------------------------------------//

  // some parameters
	star->next = NULL;
	star->idx = -1000;
	star->level = level;
	star->is = is;
	star->isStar = 6;
	star->rhocell = field->d;

  // set star position to cell center
	star->x = xc ;
	star->y = yc ;
	star->z = zc ;

  // set star velocity to fluid velocity
 	star->vx = field->u;
	star->vy = field->v;
	star->vz = field->w;

#define RDM_STARS
#ifdef RDM_STARS
  // random position
	star->x += rdm(-0.5,0.5)*dx;
	star->y += rdm(-0.5,0.5)*dx;
	star->z += rdm(-0.5,0.5)*dx;

  // compute random velocity component
	REAL r = rdm(0,1) * field->a;
	REAL theta  = acos(rdm(-1,1));
	REAL phi = rdm(0,2*M_PI);

  // add random velocity component
	star->vx += r * sin(theta) * cos(phi);
	star->vy += r * sin(theta) * sin(phi);
	star->vz += r * cos(theta) ;
#endif // RDM_STARS

  //mass
	star->mass = mstar;

  //energy
	star->epot = 0.0;
	star->ekin = 0.5 * star->mass * (POW(star->vx,2) + POW(star->vy,2) + POW(star->vz,2));

  // age
#ifdef TESTCOSMO
	star->age = param->cosmo->tphy;
#else
/// TODO fix age of star ifndef TESTCOSMO
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int testCond(struct CELL *cell, struct RUNPARAMS *param, REAL aexp, int level){
// ----------------------------------------------------------//
/// test if a cell is elligible to form star
// ----------------------------------------------------------//

	if (cell->child != NULL) return 0;
	int A=1,B=1,C=1;

	// test if density is over the threshold
	A = 	cell->field.d > param->stars->thresh;

#ifdef WGRAV
#ifdef JEANSCRIT

	REAL dx = POW(0.5,level);

	REAL rho_m = (cell->gdata.d+1.) / param->stars->thresh;

	REAL fact_rho = POW(aexp,3)/param->unit.unit_d;
	REAL fact_t = POW(aexp,2) * param->unit.unit_t;

	// local free fall time in seconde in code unit
	REAL t_ff = SQRT(3.*M_PI/(32.*NEWTON_G * rho_m/ fact_rho));
	t_ff /= fact_t;

	// local Jeans time in second in code unit
	REAL t_j = dx/cell->field.a;

	B = t_j > t_ff;
#endif // JEANSCRIT
#endif // WGRAV

#ifdef WRAD
#ifdef TEMPCRIT
  C= (cell->rfield.temp < 2e4);
#endif // TEMPCRIT
#endif // WRAD


	return A && B && C;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void conserveField(struct Wtype *field, struct RUNPARAMS *param, struct PART *star, REAL dx, REAL aexp, REAL mstar){
// ----------------------------------------------------------//
/// compute the conservations equations after a star is created
// ----------------------------------------------------------//

	REAL drho = mstar / POW(dx,3.);
	struct Utype U;
	struct Wtype W;

	memcpy(&W,field,sizeof(struct Wtype));

	W2U(&W, &U);

//	density
	U.d    -= drho;

#ifdef WRADHYD
	REAL xion=W.dX/W.d;
	U.dX = U.d*xion;
#endif

//	momentum
	U.du -= star->vx * drho;
	U.dv -= star->vy * drho;
	U.dw -= star->vz * drho;

//	internal energy
#ifdef DUAL_E
	U.eint=U.eint*(1.-drho/W.d); // assuming T and x remain constant
#endif

	U2W(&U, &W);

//	total energy
	getE(&W);
	W.a=SQRT(GAMMA*W.p/W.d);
	W.p=FMAX(W.p,PMIN);
	memcpy(field,&W,sizeof(struct Wtype));

	if(isnan(U.du)){
	  printf("drho=%e vx=%e\n",drho,star->vx);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL getSFR(struct CELL *cell, struct RUNPARAMS *param, REAL aexp, int level){
/**
  * Compute the local Star Formation Rate
  * typically return the local density over a caracteristic time
  *
  */

#ifdef SCHAYE
	/* REAL A=1.515e-4; // Ms/yr/kpc2 */
	/* A=A*(2e30)/(3600.*24.*365)/POW(1e3*PARSEC,2.); // kg/sec/m2 */
	REAL A=1.0028e-20; // Kg/sec/m2
	REAL E=1.; // Ms/pc2
	E=E*SOLAR_MASS/POW(PARSEC,2.);

	REAL P=cell->field.p/POW(aexp,5)*param->unit.unit_n*param->unit.unit_d*POW(param->unit.unit_v,2); // Physical Pressure (S.I)
	REAL geff=5./3.;

	REAL tstars 	= 1./(A*POW(E,-1.4)*POW(geff/NEWTON_G*P,(1.4-1.)/2.));
	//printf("A=%e E=%e P=%e p=%e c=%e tstars=%e\n",A,E,P,cell->field.p,param->unit.unit_d,tstars/(3600.*24*265*1e9));
	//	abort();

	REAL tstartilde = tstars / POW(aexp,2)/param->unit.unit_t; //tstars in code unit

	REAL SFR=cell->field.d / tstartilde*param->stars->efficiency; // NOTE the efficiency is neutralized

#else
#ifndef SIMPLESTAR

//	REAL rho_m = cell->gdata.d+1.;
	REAL rho_m = cell->field.d;

	REAL fact_rho = POW(aexp,3)/param->unit.unit_d;
	REAL fact_t = POW(aexp,2) * param->unit.unit_t;

	// local free fall time in seconde in code unit
	REAL t_ff = SQRT(3.*M_PI/(32.*NEWTON_G * rho_m/ fact_rho)); /// TODO find the expression in the case of a cosmological Poisson equation
	t_ff /= fact_t;

	// star formation rate in kg/s/m3 in code unit
	REAL SFR = param->stars->efficiency * cell->field.d  / t_ff;

#else

	// DOM HACK STAR FORMATION
	REAL tcarac=param->stars->efficiency; // Gyrs
	REAL fact_t=POW(aexp,2) * param->unit.unit_t;

	REAL tstars=tcarac*(3600.*24*365*1e9)/SQRT(cell->field.d/param->stars->thresh);
	tstars /= fact_t;

  REAL SFR = cell->field.d/tstars;
#endif // SIMPLESTAR
#endif // SCHAYE

  //printf("SFR=%e\n",SFR);
	return SFR; // kg.yr-1.m-3 in code unit
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL getMstars2create(struct CELL *cell, struct RUNPARAMS *param, REAL dt, REAL aexp, int level){

  REAL dx = POW(0.5,level);
  REAL dv = POW(dx,3);
  REAL SFR=getSFR(cell,param,aexp,level);

//average mass to create in cell
	REAL m_sfr =  SFR  *dt *dv;

	/// TODO drawn in a statistical law

	return m_sfr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int getNstars2create(struct CELL *cell, struct RUNPARAMS *param, REAL dt, REAL aexp, int level, REAL mstar){
// ----------------------------------------------------------//
/// Compute the number of stars to create in a given cell\n
/// And do a random draw in a Poisson law
// ----------------------------------------------------------//

  REAL dv = POW(0.5,3*level);
  const REAL SFR=getSFR(cell,param,aexp,level);
  // Average number of stars created
	const REAL lambda =  SFR  / mstar * dt * dv;

//#define NOPOISS
#ifndef NOPOISS
#ifdef GSLRAND
	unsigned int N = gsl_ran_poisson (param->stars->rpoiss, (double)lambda);
#else
	unsigned int N = gpoiss(lambda); //Poisson drawing
#endif // GSLRAND
#else
  unsigned int N = round(lambda);
#endif // NOPOISS

	if(N){
	  //printf("tstar=%e lambda=%e\n",t0/(3600.*24.*365.*1e9),lambda);
	  //printf("AVG star creation =%e /eff %d  SFR=%e\n",lambda,N,SFR);
	}

	REAL M_in_cell = cell->field.d * dv; // mass of the curent cell in code unit
	if(N * mstar >= M_in_cell) N = 0.9*M_in_cell / mstar ; // 0.9 to prevent void cells

#ifdef CONTINUOUSSTARS
  N=1;
#endif // CONTINUOUSSTARS

  return N;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addStar(struct CELL *cell, struct Wtype *init_field, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dttilde, struct RUNPARAMS *param, REAL aexp, int is,  int nstars, REAL mstar){
// ----------------------------------------------------------//
/// Add a stellar particle in the double linked list\n
/// Call the initialisation and the conservation functions
// ----------------------------------------------------------//

#ifdef PIC
	struct PART *star = cpu->freepart;

	if (star==NULL){
    printf("\n");
    printf("----------------------------\n");
    printf("No more memory for particles\n");
    printf("----------------------------\n");
	//	exit(0);
	}else{

		cpu->freepart = cpu->freepart->next;
		cpu->freepart->prev = NULL;

		if (cell->phead==NULL){
			cell->phead = star;
			star->prev = NULL;
			star->next =NULL;
		}else{
			struct PART *lasp = findlastpart(cell->phead);
			lasp->next = star;
			star->prev = lasp;
			star->next = NULL;
		}

		cpu->npart[level-1]++;
		cpu->nstar[level-1]++;

		REAL dx = POW(2.0,-level);

		initStar(init_field, star, param, level, xc, yc, zc, -1, aexp, is, dttilde, dx,mstar);

		conserveField(&cell->field   , param, star,  dx, aexp,mstar);
		conserveField(&cell->fieldnew, param, star,  dx, aexp,mstar);

		//printPart(star);
	}
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initThresh(struct RUNPARAMS *param,  REAL aexp){
// ----------------------------------------------------------//
/// Compute the density threshold of star formation
// ----------------------------------------------------------//

#ifdef TESTCOSMO
#ifdef SCHAYE
	// std value for overdensity = 55.7
	param->stars->thresh = FMAX(1e6*POW(aexp,3.) *PROTON_MASS/param->unit.unit_d, param->stars->overdensity_cond* (param->cosmo->ob/param->cosmo->om));
#else

  REAL   k=-1;                                      // Comoving density case
  if (param->stars->density_cond>0) k=POW(aexp,3);  // Physical density case

  // Hydrogen atom per cubic meter in code unit
   REAL thresh_1 = k * param->stars->density_cond * PROTON_MASS / param->unit.unit_d*param->unit.unit_N;

  // density in kg.m-3 in code unit
  //REAL thresh_1 = k * param->stars->density_cond / param->unit.unit_d*param->unit.unit_N;

  // overdensity
   REAL thresh_2 = param->stars->overdensity_cond * (param->cosmo->ob/param->cosmo->om);

  // final threshold
  param->stars->thresh = FMAX(thresh_1,thresh_2);
#endif

#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int setStarsState(struct RUNPARAMS *param, struct CPUINFO *cpu, int level){
// ----------------------------------------------------------//
/**
  Define the state of a particle in function of his age\n
  State are defined as follow:\n

  stars radiation  supernovae

  000 -> 0 -> DM
  001 -> 1 -> Not allowed
  010 -> 2 -> Not allowed
  011 -> 3 -> Not allowed
  100 -> 4 -> Dead Star
  101 -> 5 -> SN
  110 -> 6 -> Rad
  111 -> 7 -> SN + rad
*/
// ----------------------------------------------------------//

  if (param->stars->n){

    //if(cpu->rank == RANK_DISP) printf("setting states");

    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *curoct=cpu->octList[level-1][iOct];

      int icell;
      for(icell=0;icell<8;icell++) {
        struct CELL *curcell = &curoct->cell[icell];
        struct PART *nexp=curcell->phead;
        do{
          if(nexp==NULL) continue;
          struct PART *curp=nexp;
          nexp=curp->next;

#ifdef TESTCOSMO
          const REAL age =  param->cosmo->tphy - curp->age;
#else
          // TODO considere non socmo case
          const REAL age =0;
#endif // TESTCOSMO


          if( curp->isStar && !(curp->isStar==4) ){ // isStar and not dead
            if(age>=0){ // for inter-level communications

//------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------//
#ifdef DECREASE_EMMISIVITY_AFTER_TLIFE
              const REAL tlife_rad = 100*param->stars->tlife;
#else
              const REAL tlife_rad = param->stars->tlife;
#endif

#ifdef SUPERNOVAE
              const REAL tlife_sn  = param->sn->tlife;
#else
              const REAL tlife_sn  = param->stars->tlife;
#endif // SUPERNOVAE
              const int SNdead=(curp->isStar==8);
              const int isSN =age>=tlife_sn ?1:0;
              const int isRAD=age>=tlife_rad?0:1;

//	      printf("age=%e  tsn=%e\n",age,tlife_sn);

	      int isAGN=0;
#ifdef AGN
	      if(curp->isStar==100) isAGN=1;
#endif

	      if(!isAGN){
		if(SNdead){
		  if(isRAD){
		    curp->isStar=8;
		  }
		  else{
		    curp->isStar=4;
		  }
		}
		else{
		  if( (isRAD==0) && (isSN==0)) curp->isStar=4;
		  if( (isRAD==0) && (isSN==1)) curp->isStar=5;
		  if( (isRAD==1) && (isSN==0)) curp->isStar=6;
		  if( (isRAD==1) && (isSN==1)) curp->isStar=7;
		}
	      }

//------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------//

            }
          }

#ifdef WRAD
          int igrp_time;
          for(igrp_time=0;igrp_time<NGRP_TIME;igrp_time++){
            if (age>param->atomic.time_bound[igrp_time]){
              curp->radiative_state=igrp_time;
              break;
            }
          }
#endif // WRAD

 //------------------------------------------------//

        }while(nexp!=NULL);
      }
    }
    //if(cpu->rank == RANK_DISP) printf(" done\n");
  }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL setmStar(struct RUNPARAMS *param,int level){
// ----------------------------------------------------------//
/// Compute the stellar particle mass
// ----------------------------------------------------------//
  REAL mstars_level=0; // mass of stars at level

	//mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*level)*param->stars->overdensity_cond; // variable mass
	//mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*param->lcoarse)*param->stars->overdensity_cond; // coarse mass+ overdensity
	//mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse)); // coarse mass

  REAL mlevel=0;
  REAL res=param->stars->mass_res;

  if(res>100){
    mstars_level = param->stars->mass_res*SOLAR_MASS/param->unit.unit_mass;
  }else{

    if(res>=0){
      mlevel=param->lcoarse; // Fix mass
    }else{
      mlevel=level-1; //Adaptative mass
      res*=-1;
    }
#ifdef TESTCOSMO
    mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(mlevel+res));
#endif // TESTCOSMO
// TODO considere ifndef TESTCOSMO
  }

  return mstars_level;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void setID(struct RUNPARAMS *param, struct CPUINFO *cpu,int level){
/**
  * Compute the Id of newly formed stars
  * We need to know the total numbers of new stars for all process
  */


  const int debug=0;
  if (debug) printf("setting ID \n");

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI

  //count old and new stars
  int locNOld=0;
  int locNNew=0;

  int iOct;
  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    struct OCT *curoct=cpu->octList[level-1][iOct];

    int icell;
    for(icell=0;icell<8;icell++) {
      struct CELL *curcell = &curoct->cell[icell];

      struct PART * nexp=curoct->cell[icell].phead;
      if(nexp!=NULL){
        do{
          struct PART *curp=nexp;
          nexp=curp->next;

            if(curp->idx==-1000) locNNew++;
            else                 locNOld++;

        }while(nexp!=NULL);
      }
    }
  }

  if (debug) printf("locNNew %d\n", locNNew);
  if (debug) printf("locNOld %d\n", locNOld);


#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI

  //reduce
  int noldtot;
#ifdef WMPI
  MPI_Allreduce(&locNOld,&noldtot,1,MPI_INT,MPI_SUM,cpu->comm);
#else
  noldtot = locNOld;
#endif // WMPI


  int* NNew = (int*)calloc(cpu->nproc,sizeof(int));
#ifdef WMPI
  MPI_Allgather(&locNNew,1,MPI_INT,NNew ,1,MPI_INT, cpu->comm);
#else
  NNew[0]=locNNew;
#endif // WMPI

  int cpuoffset = 0;
  int i;
  for (i=0;i<cpu->rank;i++){
    cpuoffset += NNew[i];
  }

  if (debug) printf("cpuoffset %d\n", cpuoffset);

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI

  // setting ID
  int curID=0;

  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    struct OCT *curoct=cpu->octList[level-1][iOct];

    int icell;
    for(icell=0;icell<8;icell++){
      struct CELL *curcell = &curoct->cell[icell];

      struct PART *nexp=curcell->phead;
      if(nexp!=NULL){
        do{
          struct PART *curp=nexp;
          nexp=curp->next;

          if(curp->idx==-1000){
            //if (debug) printf("%d %d %d\n",noldtot, cpuoffset, curID);
            curp->idx = noldtot+ cpuoffset+ curID++;
            //printf("%d \n",curp->idx);
          }

        }while(nexp!=NULL);
      }
    }
  }

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI

  free(NNew);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Stars(struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){
// ----------------------------------------------------------//
/// The stars creation function.\n
/// Scan if cell is allowed to form star\n
/// If true, compute how many star are needed\n
/// and add them to the linked list\n
// ----------------------------------------------------------//

const int debug=0;

#ifdef GSLRAND
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  param->stars->rpoiss=r;
#endif

  if(cpu->rank == RANK_DISP) printf("STARS\n");

  if(debug) printf("setting star states\n");
  setStarsState(param, cpu, level);
  if(debug) printf("setting star states ok\n");

  REAL mmax = 0;
  REAL percentvol =0;

  int n_unit_stars = 0;
  int n_part_stars = 0;
  int nstarsmax_in_one_cell=0;
  int nstarsmin_in_one_cell=1e3;

  if(debug) printf("init threshold\n");
  initThresh(param, aexp);

  if(debug) printf("init mstar\n");
  REAL mstars_level = setmStar(param,level);
  if(debug) printf("mstars_level = %f\n",mstars_level );

  if(debug) printf("main loop\n");
  int iOct;
  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    struct OCT *curoct=cpu->octList[level-1][iOct];

    int icell;
    for(icell=0;icell<8;icell++) {
      struct CELL *curcell = &curoct->cell[icell];

      if( testCond(curcell, param, aexp, level) ) {
	REAL dx = POW(2.0,-level);
	REAL xc=curoct->x+( icell    & 1)*dx+dx*0.5;
	REAL yc=curoct->y+((icell>>1)& 1)*dx+dx*0.5;
	REAL zc=curoct->z+( icell>>2    )*dx+dx*0.5;

	percentvol += POW(dx,3);

	int N = getNstars2create(curcell, param, dt, aexp, level,mstars_level);
  if(N && debug) printf("N stars = %d\n",N);

	nstarsmax_in_one_cell = FMAX(nstarsmax_in_one_cell,N);
	nstarsmin_in_one_cell = FMIN(nstarsmin_in_one_cell,N);

	if(N>0.2*param->npartmax){
    printf("You are about to create more stars than 20 percent of npartmax N=%d nmax=%d on proc %d\n",(int)N,param->npartmax,cpu->rank);
    if(N>param->npartmax){
      printf("You are about to create more stars than npartmax N=%d nmax=%d on proc %d--> ABORT\n",(int)N,param->npartmax,cpu->rank);
      abort();
    }
  }

  //if(N) printf("N_Rho_Temp_Seuil_z\t%d\t%e\t%e\t%e\t%e\n", N, curcell->field.d, curcell->rfield.temp, param->stars->thresh,1.0/aexp - 1.0  );


#ifdef MULTIPLESTARS
  struct Wtype init_state=curcell->field;
	int ipart;
	for (ipart=0;ipart< N; ipart++){
	  addStar(curcell, &init_state, level, xc, yc, zc, cpu, dt, param, aexp, is,n_part_stars++, mstars_level);
	  n_unit_stars++;
  }
#else
#ifdef CONTINUOUSSTARS

  REAL m = getMstars2create(curcell, param, dt, aexp, level);
  //if(m>mstars_level)
  {
    addStar(curcell, &curcell->field, level, xc, yc, zc, cpu, dt, param, aexp, is,n_part_stars++, m);
    n_unit_stars++;
  }
#else
  if(N){
    addStar(curcell, &curcell->field, level, xc, yc, zc, cpu, dt, param, aexp, is,n_part_stars++, N*mstars_level);
    n_unit_stars++;
  }
#endif // CONTINUOUSSTARS
#endif // MULTIPLESTARS

/*
if (N){
printWtype(&init_state);
printWtype(&curcell->field);
}
*/
      }
      mmax = FMAX(curcell->field.d, mmax);
    }
  }

//    printf("Adding %d stars on cpu %d\n",n_part_stars, cpu->rank);


#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&n_unit_stars,1,MPI_INT,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&n_part_stars,1,MPI_INT,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&mmax,1,MPI_REEL,MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&percentvol,1,MPI_REEL,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&nstarsmax_in_one_cell,1,MPI_INT,MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&nstarsmin_in_one_cell,1,MPI_INT,MPI_MIN,cpu->comm);
#endif

  param->stars->n += n_part_stars;
  if(cpu->rank==RANK_DISP) {
    printf("Mmax=%e\tthresh=%e\tvol=%e\n", mmax, param->stars->thresh, percentvol );
    if (n_unit_stars){
      printf("%d stars added in %d particles on level %d --> min=%d max=%d \n", n_unit_stars, n_part_stars, level,nstarsmin_in_one_cell,nstarsmax_in_one_cell);
      printf("%d stars particles in total\n",param->stars->n);
    }
    int nsta;
#ifndef AGN
    nsta=param->stars->n;
#else
    nsta=param->stars->n-param->stars->nagn;
#endif
    if(cpu->trigstar==0 && nsta>0) printf("FIRST_STARS at z=%e\n",1./aexp-1.);
    if(param->stars->n>0) cpu->trigstar=1;
  }

#ifdef GSLRAND
  gsl_rng_free (r);
#endif

  setID(param, cpu,level);

}
#endif//STARS
