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

void initStar(struct CELL * cell, struct PART *star, struct RUNPARAMS *param, int level, REAL xc, REAL yc, REAL zc,int idx, REAL aexp, int is, REAL dt,REAL dx, REAL mlevel) {
// ----------------------------------------------------------//
/// compute the initial state of star particle
// ----------------------------------------------------------//

  // some parameters
	star->next = NULL;
	star->idx = idx;
	star->level = level;
	star->is = is;
	star->isStar = 1;
  star->rhocell = cell->field.d;

  // random position
	star->x = xc + rdm(-0.5,0.5) * dx;
	star->y = yc + rdm(-0.5,0.5) * dx;
	star->z = zc + rdm(-0.5,0.5) * dx;

  // set star velocity to fluid velocity
 	star->vx = cell->field.u;
	star->vy = cell->field.v;
	star->vz = cell->field.w;

  // random component
	REAL r = rdm(0,1) * cell->field.a ;
	REAL theta  = acos(rdm(-1,1));
	REAL phi = rdm(0,2*M_PI);

  // add random component
	star->vx += r * sin(theta) * cos(phi);
	star->vy += r * sin(theta) * sin(phi);
	star->vz += r * cos(theta) ;

  //mass
	star->mass = mlevel;

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
	int A,B;

	// test if density is over the threshold
	A = 	cell->field.d > param->stars->thresh;
#ifdef WGRAV
  // test the Jeans criterion
	//B = cell->field.a/POW(2.,-level) > SQRT(6.*aexp * cell->gdata.d +1.) ;
	REAL dx = POW(2.,-level);

  REAL jeans = M_PI/(16.*NEWTON_G) * POW(cell->field.a/dx,2);
  jeans *= POW(aexp,3) /(param->unit.unit_d);

  B = cell->field.d > jeans;
  B = 1;
#else
	B = 1;
#endif

	return A && B;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void conserveField(struct Wtype *field, struct RUNPARAMS *param, struct PART *star, REAL dx, REAL aexp, REAL mlevel){
// ----------------------------------------------------------//
/// compute the conservations equations after a star is created
// ----------------------------------------------------------//

	REAL drho = mlevel / POW(dx,3.);
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
	U.eint=U.eint*(1.-drho/W.d); // assuming T and x remain constant

	U2W(&U, &W);

//	total energy
	getE(&(W));
	W.a=SQRT(GAMMA*W.p/W.d);
	W.p=FMAX(W.p,PMIN);
	memcpy(field,&W,sizeof(struct Wtype));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int getNstars2create(struct CELL *cell, struct RUNPARAMS *param, REAL dt, REAL aexp, int level, REAL mlevel){
// ----------------------------------------------------------//
/// Compute the number of stars to create in a given cell\n
/// And do a random draw in a Poisson law
// ----------------------------------------------------------//

#ifdef SCHAYE
	/* REAL A=1.515e-4; // Ms/yr/kpc2 */
	/* A=A*(2e30)/(3600.*24.*365)/POW(1e3*PARSEC,2.); // kg/sec/m2 */
	REAL A=1.0028e-20; // Kag/sec/m2
	REAL E=1.; // Ms/pc2
	E=E*SOLAR_MASS/POW(PARSEC,2.);

	REAL P=cell->field.p/POW(aexp,5)*param->unit.unit_n*param->unit.unit_d*POW(param->unit.unit_v,2); // Physical Pressure (S.I)
	REAL geff=5./3.;

	REAL tstars 	= 1./(A*POW(E,-1.4)*POW(geff/NEWTON_G*P,(1.4-1.)/2.));
	//printf("A=%e E=%e P=%e p=%e c=%e tstars=%e\n",A,E,P,cell->field.p,param->unit.unit_d,tstars/(3600.*24*265*1e9));
	//	abort();

	REAL tstartilde = tstars / POW(aexp,2)/param->unit.unit_t; //tstars in code unit

	REAL M_in_cell 	= cell->field.d * POW(2.0,-3.0*level); // mass of the curent cell in code unit

	REAL lambda =  param->stars->efficiency * M_in_cell / mlevel * dt/ tstartilde; // Average number of stars created
#else
	// local free fall time in seconde in code unit
	// REAL t_ff = 1. / SQRT(6*aexp*cell->gdata.d);
	//printf("Local SFR=%e M0/yr/Mpc3\n", SFR/SOLAR_MASS*31556926*POW(PARSEC,3));
	
	REAL dx = POW(0.5,level);
	REAL dv = POW(0.5,3*level);
	
	REAL rho_m = (cell->gdata.d+1.) / param->stars->thresh;
	REAL rho_b =  cell->field.d     / param->stars->thresh;
	
	REAL fact_rho = POW(aexp,3)/param->unit.unit_d;
	REAL fact_t = POW(aexp,2) * param->unit.unit_t;
	
	// local free fall time in seconde in code unit
	REAL t_ff = SQRT(3.*M_PI/(32.*NEWTON_G * rho_m/ fact_rho)); /// TODO find the expression in the case of a cosmological Poisson equation
	t_ff /= fact_t;
  
  // local Jeans time in seconde in code unit
	REAL t_j = dx/cell->field.a;
	
	// star formation rate in kg/s/m3 in code unit
	REAL SFR = param->stars->efficiency * cell->field.d  / t_ff;

	// Jeans efficiency
	//SFR *= t_j/t_ff;
	
  // Average number of stars created
	REAL lambda =  SFR  / mlevel * dt * dv;
	//printf("rho=%e tff=%e tj=%e SFR=%e tstar=%e\n",cell->field.d, t_ff, t_j, SFR,t_ff*t_ff/t_j/param->stars->efficiency*fact_t/(3600.*24.*365.*1e9));

#endif //SCHAYE

#ifdef GSLRAND
	int N = gsl_ran_poisson (param->stars->rpoiss, lambda);
#else
	int N = gpoiss(lambda); //Poisson drawing
#endif
	//printf("AVG star creation =%e /eff %d\n",lambda,N);
	REAL M_in_cell = cell->field.d * POW(2.0,-3.0*level); // mass of the curent cell in code unit
	if(N * mlevel >= M_in_cell) N = 0.9*M_in_cell / mlevel ; // 0.9 to prevent void cells

	return N;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addStar(struct CELL *cell, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dttilde, struct RUNPARAMS *param, REAL aexp, int is,  int nstars, REAL mlevel){
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

		initStar(cell, star, param, level, xc, yc, zc, param->stars->n+nstars, aexp, is, dttilde, dx,mlevel);

		conserveField(&(cell->field   ), param, star,  dx, aexp,mlevel);
		conserveField(&(cell->fieldnew), param, star,  dx, aexp,mlevel);
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
/// Define the state of a particle in function of his age\n
/// State are defined as follow:\n
///
///  state = 0 -> Dark Matter\n
///  state = 1 -> Radiative Star\n
///  state = 2 -> Supernovae\n
///  state = 3 -> Supernovae + decreasing luminosity\n
///  state = 4 -> Decreasing luminosity\n
///  state = 5 -> Dead star
// ----------------------------------------------------------//

  if (param->stars->n){

    if(cpu->rank == RANK_DISP) printf("setting states\n");

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

          //------------------------------------------------//
          if(curp->isStar && curp->isStar < 5){ // Star not dead

            REAL t0 =  param->cosmo->tphy - curp->age;
            if(t0>=0){ // for inter-level communications
              REAL tlife = param->stars->tlife;

              if( (curp->isStar==4) && (t0>=100*tlife) ){
                curp->isStar=5; // decreasing luminosity -> dead star
              }

              if( curp->isStar==3){
                curp->isStar=4; //Supernovae + decreasing luminosity -> decreasing luminosity
                //curently supernovae are instantaneous
                /// TODO implement slow feedback
              }

              if(curp->isStar==2){
                curp->isStar=5; // supernovae -> dead star
                //curently supernovae are instantaneous
              }

              if( (curp->isStar==1) && (t0>=tlife) ){
#ifdef SUPERNOVAE
  #ifdef DECREASE_EMMISIVITY_AFTER_TLIFE
                curp->isStar=3; // radiative -> supernovae + decreasing luminosity
  #else
                curp->isStar=2; // radiative -> supernovae
  #endif // DECREASE_EMMISIVITY_AFTER_TLIFE
#else
                curp->isStar=4; // radiative -> decreasing luminosity
#endif // SUPERNOVAE
              }
            }
          }
          //------------------------------------------------//

        }while(nexp!=NULL);
      }
    }
    if(cpu->rank == RANK_DISP) printf("setting states done\n");
  }
  return 0;
}

REAL setmStar(struct RUNPARAMS *param,int level){
// ----------------------------------------------------------//
/// Compute the stellar particle mass
// ----------------------------------------------------------//
  REAL mstars_level=0; // mass of stars at level
#ifdef TESTCOSMO
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
    mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(mlevel+res));
  }
#endif // TESTCOSMO

// TODO considere ifndef TESTCOSMO


  return mstars_level;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Stars(struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){
// ----------------------------------------------------------//
/// The stars creation function.\n
/// Scan if cell is allowed to form star\n
/// If true, compute how many star are needed\n
/// and add them to the linked list\n
// ----------------------------------------------------------//

#ifdef GSLRAND
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  param->stars->rpoiss=r;
#endif

  if(cpu->rank == RANK_DISP) printf("STARS\n");
  setStarsState(param, cpu, level);

  REAL mmax = 0;
  REAL percentvol =0;
  int nstars = 0;

  initThresh(param, aexp);
  REAL mstars_level = setmStar(param,level);

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
	if(N>0.2*param->npartmax){
	  printf("You are about to create more stars than 20 percent of npartmax N=%d nmax=%d on proc %d\n",N,param->npartmax,cpu->rank);
	}

	if(N>param->npartmax){
	  printf("You are about to create more stars than npartmax N=%d nmax=%d on proc %d--> ABORT\n",N,param->npartmax,cpu->rank);
	  abort();
	}

	//	if(N) printf("N_Rho_Temp_Seuil_z\t%d\t%e\t%e\t%e\t%e\n", N, curcell->field.d, curcell->rfield.temp, param->stars->thresh,1.0/aexp - 1.0  );
	int ipart;
	for (ipart=0;ipart< N; ipart++){
          addStar(curcell, level, xc, yc, zc, cpu, dt, param, aexp, is,nstars++, mstars_level);
        }
      }
      mmax = FMAX(curcell->field.d, mmax);
    }
  }

#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&nstars,1,MPI_INT,   MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&mmax,  1,MPI_REEL,	MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&percentvol,  1,MPI_REEL,	MPI_SUM,cpu->comm);
#endif

  param->stars->n += nstars ;
  if(cpu->rank==RANK_DISP) {
    printf("Mmax=%e\tthresh=%e\tvol=%e\n", mmax, param->stars->thresh, percentvol );
    if (nstars){
      printf("%d stars added on level %d \n", nstars, level);
      printf("%d stars in total\n",param->stars->n);
    }
    if(cpu->trigstar==0 && param->stars->n>0) printf("FIRST_STARS at z=%e\n",1./aexp-1.);
    if(param->stars->n>0) cpu->trigstar=1;
  }

#ifdef GSLRAND
  gsl_rng_free (r);
#endif

}
#endif//STARS
