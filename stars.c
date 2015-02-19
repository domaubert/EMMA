#ifdef STARS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "prototypes.h"
#include "particle.h"
#include "oct.h"
#include "hydro_utils.h"

#include "cic.h"
#include "segment.h"

#ifdef WMPI
#include <mpi.h>
#endif


/*
known t_car :
 2.1 Gyr (Kennicutt 1998)
 2.4 Gyr (Rownd & Young 1999)
 3.0 Gyr (Rasera & Teyssier 2006)
 3.2 Gyr (Springel & Hernquist 2003)
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 		SOME FONCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef TESTCOSMO
REAL a2t(struct RUNPARAMS *param, REAL az ){
//  	from Ned Wright http://www.astro.ucla.edu/~wright/CC.python

	REAL age = 0.;
	REAL a, adot;
	REAL h = param->cosmo->H0/100;

	REAL or = 4.165e-5/(h*h);
	REAL ok = 1. - param->cosmo->om - or - param->cosmo->ov;

	int i, n=1000;
	for (i=0; i<n;i++){
		a = az*(i+0.5)/n;
		adot = SQRT(ok+(param->cosmo->om / a)+(or/ (a*a) )+ (param->cosmo->ov*a*a) );
		age = age + 1./adot;
	}
	REAL zage = az*age/n;
	REAL zage_Gyr = (977.8 /param->cosmo->H0)*zage;

	return zage_Gyr*1e9;
}
#endif


REAL rdm(REAL a, REAL b){
	return 	(rand()/(REAL)RAND_MAX ) * (b-a) + a ;
}

int gpoiss(REAL lambda){
	int k=1;
	REAL p = rdm(0,1);
	REAL P = exp(-lambda);
	REAL sum=P;
	if (sum>=p){
	  k=0;
	}
	else{
	  do { 	P*=lambda/(REAL)k;
	    sum+=P;
	    if (sum>=p) break;
	    k++;
	  }while(k<1e6);
	}
	//printf("k=%d lambda=%e sum=%e p=%e\n",k,lambda,sum,p);
	return k;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		STARS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int countstar(struct PART* phead)
{
  struct PART* curp;
  struct PART* nexp;
  int npart=0;

  curp=NULL;
  nexp=phead; //sweeping the particles of the current cell */
  if(nexp!=NULL){
    do{
      curp=nexp;
      nexp=curp->next;
      npart+=(curp->isStar==1);
    }while(nexp!=NULL);
  }

  return npart;
}


void initStar(struct CELL * cell, struct PART *star, struct RUNPARAMS *param, int level, REAL xc, REAL yc, REAL zc,int idx, REAL aexp, int is, REAL dt,REAL dx, REAL mlevel) {

	star->next  = NULL;
	star->idx   = idx;
	star->level = level;
	star->is    = is;

	star->x = xc + rdm(-0.5,0.5) * dx;
	star->y = yc + rdm(-0.5,0.5) * dx;
	star->z = zc + rdm(-0.5,0.5) * dx;

 	star->vx = cell->field.u;
	star->vy = cell->field.v;
	star->vz = cell->field.w;

	REAL r 		= rdm(0,1) * cell->field.a ;
	REAL theta   	= acos(rdm(-1,1));
	REAL phi 	= rdm(0,2*PI);

	star->vx += r * sin(theta) * cos(phi);
	star->vy += r * sin(theta) * sin(phi);
	star->vz += r * cos(theta) ;

	star->mass  = mlevel;

	star->epot = 0.0;
	star->ekin = 0.5 * star->mass * (POW(star->vx,2) + POW(star->vy,2) + POW(star->vz,2));

	star->isStar = 1;
#ifdef TESTCOSMO
	star->age = param->cosmo->tphy;
#endif
	star->rhocell = cell->field.d;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int testCond(struct CELL *cell, REAL dttilde, REAL dxtilde, struct RUNPARAMS *param, REAL aexp, int level){

	if (cell->child != NULL) return 0;
	int A,B;
	A = 	cell->field.d > param->stars->thresh;
#ifdef WGRAV
	B = 0?	cell->field.a/POW(2.,-level) > SQRT(6.*aexp * cell->gdata.d +1.) 				: 1;
#else
	B = 1;
#endif

	return 		A && B;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void conserveField(struct Wtype *field, struct RUNPARAMS *param, struct PART *star, REAL dx, REAL aexp, REAL mlevel){

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
//	???????????????

	U.eint=U.eint*(1.-drho/W.d);


	U2W(&U, &W);

//	total energy
	getE(&(W));
	W.a=SQRT(GAMMA*W.p/W.d);

	memcpy(field,&W,sizeof(struct Wtype));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int getNstars2create(struct CELL *cell, struct RUNPARAMS *param, REAL dttilde, REAL aexp, int level, REAL mlevel){

	REAL gas_efficiency = 1e0;	// maybe need to be passed in param??


#ifdef SCHAYE
	/* REAL A=1.515e-4; // Ms/yr/kpc2 */
	/* A=A*(2e30)/(3600.*24.*365)/POW(1e3*PARSEC,2.); // kg/sec/m2 */
	REAL A=1.0028e-20; // Kag/sec/m2
	REAL E=1.; // Ms/pc2
	E=E*2e30/POW(PARSEC,2.);

	REAL P=cell->field.p/POW(aexp,5)*param->unit.unit_n*param->unit.unit_d*POW(param->unit.unit_v,2); // Physical Pressure (S.I)
	REAL geff=5./3.;

	REAL tstars 	= 1./(A*POW(E,-1.4)*POW(geff/NEWTON_G*P,(1.4-1.)/2.));
	//printf("A=%e E=%e P=%e p=%e c=%e tstars=%e\n",A,E,P,cell->field.p,param->unit.unit_d,tstars/(3600.*24*265*1e9));
	//	abort();
#else
	REAL tstars 	= param->stars->tcar * 31556926 / SQRT(cell->field.d / param->stars->thresh );
#endif //SCHAYE

#ifdef WRADHYD
	REAL tstartilde = tstars / POW(aexp,2)/param->unit.unit_t;
#else
	REAL tstartilde = 1;
#endif

	REAL M_in_cell 	= cell->field.d * POW(2.0,-3.0*level);

	REAL lambda =  gas_efficiency * M_in_cell / mlevel * dttilde/ tstartilde; // Average number of stars created

	int N 		= gpoiss(lambda);

	//printf("AVG star creation =%e /eff %d\n",lambda,N);

	if(N * mlevel >= M_in_cell ) N = 0.9*M_in_cell / mlevel ; // 0.9 to prevent void cells

	return N;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void addStar(struct CELL *cell, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dttilde, struct RUNPARAMS *param, REAL aexp, int is,  int nstars, REAL mlevel){
#ifdef PIC
	struct PART *star = cpu->freepart;

	if (star==NULL){
    printf("\n");
    printf("----------------------------\n");
    printf("No more memory for particles\n");
    printf("----------------------------\n");
		exit(0);
	}else{

		cpu->freepart 		= cpu->freepart->next;
		cpu->freepart->prev 	= NULL;

		if (cell->phead==NULL){
			cell->phead 		= star;
			star->prev 		= NULL;
			star->next =NULL;
		}else{
			struct PART *lasp 	= findlastpart(cell->phead);
			lasp->next 		= star;
			star->prev 		= lasp;
			star->next=NULL;
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

#ifdef TESTCOSMO
	REAL k =(param->stars->density_cond >0)? 1 : -POW(aexp,3.0);

	REAL rhocrittilde 	= param->stars->density_cond * PROTON_MASS;

#ifdef SCHAYE
	// std value for overdensity = 55.7
	param->stars->thresh = FMAX(1e6*POW(aexp,3.) *PROTON_MASS/param->unit.unit_d, param->stars->overdensity_cond* (param->cosmo->ob/param->cosmo->om));
#else
	param->stars->thresh    = FMAX( k * rhocrittilde / param->unit.unit_d, param->stars->overdensity_cond * (param->cosmo->ob/param->cosmo->om));
#endif

#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int setStarsState(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, int level){

  /*
  state = 0 -> Dark Matter
  state = 1 -> Radiative Star
  state = 2 -> Supernovae
  state = 3 -> Supernovae + decreasing luminosity
  state = 4 -> decreasing luminosity
  state = 5 -> Dead star
  */

  struct OCT  *curoct;
	struct OCT  *nextoct=firstoct[level-1];

	do{
    if(nextoct==NULL) 		continue;
	  curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu != cpu->rank) 	continue;
    int icell;
	  for(icell=0;icell<8;icell++) {
	    struct CELL *curcell = &curoct->cell[icell];
      struct PART *curp;
      struct PART *nexp=curcell->phead;
      do{
        if(nexp==NULL) continue;
        curp=nexp;
        nexp=curp->next;

        //------------------------------------------------//
        if(curp->isStar){
          REAL t0 =  param->cosmo->tphy - curp->age;
          if(t0>=0){ // for inter-level communications
            REAL tlife = param->stars->tlife;

            if( (curp->isStar==4) && (t0>=50*tlife) ){
              curp->isStar=5; /// decreasing luminosity -> dead star
            }

            if( curp->isStar==3){
              curp->isStar=4; ///Supernovae + decreasing luminosity -> decreasing luminosity
              //curently supernovae are instantaneous
            }

            if(curp->isStar==2){
              curp->isStar=5; /// supernovae -> dead star
              //curently supernovae are instantaneous
            }

            if( (curp->isStar==1) && (t0>=tlife) ){
  #ifdef DECREASE_EMMISIVITY_AFTER_TLIFE
              curp->isStar=3; /// radiative -> supernovae + decreasing luminosity
  #else
              curp->isStar=2; /// radiative -> supernovae
  #endif

            }
          }
        }
        //------------------------------------------------//

      }while(nexp!=NULL);
    }
  }while(nextoct!=NULL);
  return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void createStars(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){

	if(cpu->rank == RANK_DISP) printf("STARS\n");

	struct OCT  *nextoct=firstoct[level-1];

	REAL dx = POW(2.0,-level);
	REAL mmax = 0;
  REAL mstars_level=0; // mass of stars at level
	int nstars = 0;
	int nsl = 0;

	initThresh(param, aexp);

#ifdef TESTCOSMO
	//mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*level)*param->stars->overdensity_cond; // variable mass
	//mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*param->lcoarse)*param->stars->overdensity_cond; // coarse mass+ overdensity
	//mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse)); // coarse mass

	REAL res = param->stars->mass_res;
	REAL mlevel=0;
	if(res>=0){mlevel=param->lcoarse;}
	else{mlevel=level+1;res*=-1.;}
  mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(mlevel+res));
#endif

	do {	if(nextoct==NULL) 		continue;
	  struct OCT  *curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu != cpu->rank) 	continue;

    int icell;
	  for(icell=0;icell<8;icell++) {
	    struct CELL *curcell = &curoct->cell[icell];

#ifndef SNTEST
	    if( testCond(curcell, dt, dx, param, aexp, level) ) {
	      REAL xc=curoct->x+( icell    & 1)*dx+dx*0.5;
	      REAL yc=curoct->y+((icell>>1)& 1)*dx+dx*0.5;
	      REAL zc=curoct->z+( icell>>2    )*dx+dx*0.5;

	      int N = getNstars2create(curcell, param, dt, aexp, level,mstars_level);

	      //	if(N) printf("N_Rho_Temp_Seuil_z\t%d\t%e\t%e\t%e\t%e\n", N, curcell->field.d, curcell->rfield.temp, param->stars->thresh,1.0/aexp - 1.0  );
        int ipart;
	      for (ipart=0;ipart< N; ipart++){
#ifdef PIC
		addStar(curcell, level, xc, yc, zc, cpu, dt, param, aexp, is,nstars++, mstars_level);
#endif //PIC
	      }
	    }
#endif //SNTEST
	    mmax = FMAX(curcell->field.d, mmax);
	  }
	}while(nextoct!=NULL);


#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&nstars,1,MPI_INT,   MPI_SUM,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&mmax,  1,MPI_REEL,	MPI_MAX,cpu->comm);
#endif

	param->stars->n += nstars ;
	if(cpu->rank==RANK_DISP) {
		printf("Mmax=%e\tthresh=%e\n", mmax, param->stars->thresh );
		if (nstars){
      printf("%d stars added on level %d \n", nstars, level);
      printf("%d stars in total\n",param->stars->n);
    }
		if(cpu->trigstar==0 && param->stars->n>0) printf("FIRST_STARS\t%e\n",1./aexp-1.);
		if(param->stars->n>0) cpu->trigstar=1;
//		printf("\n");
	}

  int l;
	for (l = param->lcoarse; l<=param->lmax;l++){
#ifdef WMPI
	  MPI_Allreduce(&cpu->nstar[l-1],&nsl,1,MPI_INT,   MPI_SUM,cpu->comm);
	  MPI_Barrier(cpu->comm);
#endif
	  //if(cpu->rank==RANK_DISP && nsl) {	printf("%d stars on level %d \n", nsl, l);	}
	}


//	if(cpu->rank==RANK_DISP) {	printf("\n");}
}
#endif
