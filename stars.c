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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void feedback(struct CELL *cell, struct RUNPARAMS *param, REAL aexp){

	REAL rhoE 	= 3.7e-15 * 1.0e-7 * 1e3; // erg.g-1  en J.kg-1   Kay 2002

	REAL H0 	= param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL rhostar 	= 3.0 * pow(H0,2.0) * param->cosmo->om /(8.0*PI*NEWTON_G);

	REAL Mtot 	= rhostar * pow(param->cosmo->unit_l,3.0);
	REAL M 		= param->stars->mstars * Mtot;

	cell->field.E  += M * rhoE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL rdm(){ return rand()/pow(2.0,32.0); }

int fac(int x){
	return 		(x>1)?(fac(x-1)*x):1;
}

int gpoiss(REAL m){
	int  a = -1;
	REAL s = -1;
	REAL v =  0;
	REAL q = exp(-m);
	REAL u = rdm();

	while( s <= u){ 
		a++;
		v +=  pow(m,a)/fac(a); 
		s = v*q;
	}
	return a ;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initStar(struct CELL * cell, struct PART *star, struct RUNPARAMS *param, int level, REAL xc, REAL yc, REAL zc,int idx, REAL aexp) {

	star->next  = NULL;
	star->idx   = idx;
	star->level = level;
	star->is    = 0;

	star->x = xc;
	star->y = yc;
	star->z = zc;

	star->vx = cell->field.u + ( rdm()*2.0-1.0 ) *cell->field.a;
	star->vy = cell->field.v + ( rdm()*2.0-1.0 ) *cell->field.a;
	star->vz = cell->field.w + ( rdm()*2.0-1.0 ) *cell->field.a;

	star->mass  = param->stars->mstars;
	star->rhocell = cell->field.d;

	star->epot = 0.0;
	star->ekin = 0.0;

	star->isStar = 1;
	star->age = 1.0/aexp-1.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int testCond(struct CELL *cell,REAL dttilde, REAL dxtilde, struct RUNPARAMS *param, REAL aexp, int level){

//	printf("%e\t%e\n", cell->field.d, n*param->stars->mstars);

	return cell->field.d > param->amrthresh/2.0 * (param->cosmo->ob/param->cosmo->om) / pow(2.0,-3.0*level);

	/*
	REAL DVcoarse 	=  pow(2.0,-3.0*param->lcoarse);
	REAL DV		=  pow(2.0,-3.0*param->level);

	REAL rhobarre 	= (param->cosmo->ob/param->cosmo->om);

	REAL N 		= param->amrthresh / pow(2.0,-3.0*param->lcoarse);

	REAL rhocrit 	= N * rhobarre * DVcoarse/DV ;
	*/


//	return cell->field.d > param->stars->overdensity_cond * (param->cosmo->ob/param->cosmo->om) * pow(aexp,-alpha0) ;

/*	REAL rhocrit 	= param->stars->density_cond * PROTON_MASS;

	REAL H0 	= param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL a3rhostar 	= pow(aexp,3.0) *  3.0 * pow(H0,2.0) * param->cosmo->om /(8.0*PI*NEWTON_G);

//	printf("H0\t%e\t\t a3rs\t%e\t\t a^-3\t%e\n",H0 , a3rhostar, pow(aexp,-alpha0) );
//	printf("rho\t%e\t\t rcrit\t%e\t\t seuil\t%e\n",cell->field.d , rhocrit,rhocrit/a3rhostar*pow(aexp,-alpha0) );
	return cell->field.d > rhocrit/a3rhostar*pow(aexp,-alpha0);
*/
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void removeMfromgas(struct CELL *cell, struct PART *star, REAL drho){

	cell->field.d    -= drho;		cell->fieldnew.d -= drho;

	struct Utype U;				struct Utype Unew;
	W2U( &(cell->field), &U);		W2U( &(cell->fieldnew), &Unew);

	U.du -= star->vx * star->mass;		Unew.du -= star->vx * star->mass;
	U.dv -= star->vy * star->mass;		Unew.dv -= star->vy * star->mass;
	U.dw -= star->vz * star->mass;		Unew.dw -= star->vz * star->mass;

	U2W(&U, &(cell->field));		U2W(&Unew, &(cell->fieldnew));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int getNstars2create(struct CELL *cell, struct RUNPARAMS *param, REAL dttilde, REAL aexp, int level){
	
	REAL H0 = param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL a2ts = pow(aexp,2.0) * 2.0/(H0 *sqrt(param->cosmo->om));

	REAL rho0 = param->stars->overdensity_cond * (param->cosmo->ob / param->cosmo->om);
	REAL tstars = param->stars->tcar * 31556926 / sqrt(cell->field.d / rho0 );

	REAL tstartilde = tstars / a2ts;

	REAL lambda = cell->field.d * pow(2.0,-3.0*level) / param->stars->mstars * dttilde/ tstartilde;

	return gpoiss(lambda);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addStar(struct CELL *cell, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dttilde, struct RUNPARAMS *param, REAL aexp){

	struct PART *star = cpu->freepart;

	if (star==NULL && cpu->rank==0){		
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
		}else{
			struct PART *lasp 	= findlastpart(cell->phead);
			lasp->next 		= star;
			star->prev 		= lasp;
		}

		removeMfromgas(cell, star, param->stars->mstars/pow(2.0,-3.0*level) );
		initStar(cell, star, param, level, xc, yc, zc, 0, aexp );
		
		cpu->npart[level-1]++;
		cpu->nstar[level-1]++;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void createStars(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level){

	param->stars->mstars	= (param->cosmo->ob/param->cosmo->om) * pow(2.0,-3.0*param->lcoarse);

	REAL rhocrittilde 	= param->amrthresh/pow(2.0,-3.0*level) * (param->cosmo->ob/param->cosmo->om);
	REAL H0 		= param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL a3rhostar 		= pow(aexp,3.0) *  3.0 * pow(H0,2.0) * param->cosmo->om /(8.0*PI*NEWTON_G);
	REAL rhocrit 		= rhocrittilde * a3rhostar ; 

	struct OCT  *curoct;
	struct OCT  *nextoct=firstoct[level-1];
	struct CELL *curcell;

	REAL xc, yc, zc;
	REAL dx = pow(2.0,-level);

	int icell, ipart;
	int nstars = 0;

	srand(time(NULL));
	REAL mmax=0;

	if(cpu->rank==0){
		printf("\n");
		printf("================================\n");
		printf("   Starting Add Stars routine   \n");
		printf("================================\n");	
	}

	do {	if(nextoct==NULL) continue;
		if(nextoct->cpu != cpu->rank) continue;

		curoct=nextoct;
		nextoct=curoct->next;

	      	for(icell=0;icell<8;icell++) {
			curcell = &curoct->cell[icell];

			if(curcell->child==NULL && testCond(curcell, dt, dx, param, aexp, level) ) {
				xc=curoct->x+( icell    & 1)*dx+dx*0.5; 
				yc=curoct->y+((icell>>1)& 1)*dx+dx*0.5;
				zc=curoct->z+( icell>>2    )*dx+dx*0.5; 										

				for (ipart=0;ipart< getNstars2create(curcell, param, dt, aexp, level); ipart++){
					addStar(curcell, level, xc, yc, zc, cpu, dt, param, aexp);
					nstars++;
				}
			}
			if (curcell->field.d > mmax) mmax = curcell->field.d;
		}
	}while(nextoct!=NULL);

#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&nstars,1,MPI_INT,MPI_SUM,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&mmax,1,MPI_FLOAT,MPI_MAX,cpu->comm);
#endif

	param->stars->n += nstars ;
	if(cpu->rank==0) {
		printf("%d stars added on level %d \n", nstars, level);
		printf("%d stars in total\n",param->stars->n);
		printf("rho max\t%e\nthreshold\t%e\n",mmax,  param->amrthresh * (param->cosmo->ob/param->cosmo->om)/pow(2.0,-3.0*level) );
		printf("rho crit [m-3] \t %e\n", rhocrit/PROTON_MASS);
		printf("\n");
	}
}
#endif
