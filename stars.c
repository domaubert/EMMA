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

REAL rdm(){
	return rand()/pow(2.0,32.0);
}

int fac(int x)
{
	int res ;
	if (x>1)	res = (fac(x-1)*x) ;
	else      	res = 1 ;

	return res ;
}

int gpoiss(REAL m){

	int  a = -1 ;
	REAL s = -1 ;
	REAL v = 0 ;
	REAL q = exp(-m) ; 
	REAL u = rdm() ;   

	while( s <= u){ 
		a++;
		v +=  pow(m,a)/fac(a) ; 
		s = v*q ;
	}
	return a ;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initStar(struct CELL * cell, struct PART *star, int level, REAL m ,REAL xc, REAL yc, REAL zc,int idx, REAL aexp) {

	star->next = NULL;

	star->idx   = idx;

	star->level = level;
	star->is    = 0;

	star->x = xc;
	star->y = yc;
	star->z = zc;


	star->vx = cell->field.u + ( rdm()*2.0-1.0 ) *cell->field.a;
	star->vy = cell->field.v + ( rdm()*2.0-1.0 ) *cell->field.a;
	star->vz = cell->field.w + ( rdm()*2.0-1.0 ) *cell->field.a;

	star->mass  = m;
//	star->rhocell = cell->field.d;

	star->epot = 0.0;
	star->ekin = 0.0;

	star->isStar = 1;
	star->age = 1.0/aexp-1.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int testCond(struct CELL *cell,REAL dttilde, REAL dxtilde, struct RUNPARAMS *param, REAL aexp){

	REAL H0 = param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL om = param->cosmo->om;
	REAL L  = param->cosmo->unit_l;

	int cond_overdensity = 0;
	int cond_density     = 0;
	int cond_jeans       = 0;

	REAL tdyn = 2.0*PI / sqrt(6.0*aexp*(cell->gdata.d+1.0));

//------OVERDENSITY--------------------------------------------------------------------------------------------//

	REAL rhotilde 	 = cell->field.d * (param->cosmo->om / param->cosmo->ob);
	cond_overdensity = rhotilde > param->stars->overdensity_cond;

	return (cond_overdensity)?1:0;

/*
//------DENSITY------------------------------------------------------------------------------------------------//

	if (cond_overdensity){
		REAL N = param->stars->density_cond;	// minimum hydrogene atome per cm^3	

		REAL mH = 1.6710e-27  ; // (kg)
		REAL rhocrit 	= N * mH * 1e6 ; //

		REAL rhostar 	= 3.0 * pow(H0,2) * om /(8.0*PI*NEWTON_G);
		REAL rho 	= rhotilde / pow(a,3) * rhostar;

		cond_density = (rho > rhocrit); 

	}else{return 0;}
	////////////	A CHECKER semble problematique
*/

//------JEANS--------------------------------------------------------------------------------------------------//
/*
	if (cond_overdensity){
		cond_jeans = (dxtilde/cell->field.a > tdyn);
	}else{return 0;}
*/
//------RANDOM-------------------------------------------------------------------------------------------------//	

/*	if(cond_overdensity) {
		REAL p = 1.0 - exp(- param->stars->eff * dttilde/tdyn);
		return (rdm() < p)?1:0;
	}else{return 0;}
*/


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
REAL getdrho(struct CELL *cell, REAL dttilde, struct RUNPARAMS *param, REAL aexp){

	REAL tdyn = 2.0*PI / sqrt( 6.0*aexp *(cell->gdata.d+1.0) );
	return cell->field.d *param->stars->eff *dttilde /tdyn;
}
*/

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

int getNstars2create(struct CELL *cell, struct RUNPARAMS *param, REAL dttilde, REAL mstars, REAL aexp, int level){

	REAL tstars = param->stars->tcar * 31556926 / sqrt(cell->field.d / param->cosmo->ob/param->cosmo->om);

	REAL H0 = param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL ts = 2.0/(H0 *sqrt(param->cosmo->om));
	REAL tstartilde = tstars /( pow(aexp,2.0) * ts);

	REAL lambda = cell->field.d * pow(2.0,-3.0*level) / mstars * dttilde/ tstartilde;

	return gpoiss(lambda);
}



int addStar(struct CELL *cell, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dttilde, struct RUNPARAMS *param, REAL aexp, REAL mstars){

	struct PART *star = cpu->freepart;

	if (star==NULL && cpu->rank==0){		
		printf("\n");
		printf("----------------------------\n");
		printf("No more memory for particles\n");
		printf("----------------------------\n");
		return 1;
	}else{

		cpu->freepart = cpu->freepart->next;
		cpu->freepart->prev = NULL;

		if (cell->phead==NULL){
			cell->phead = star;
			star->prev = NULL;				
		} else{
			struct PART *lasp = findlastpart(cell->phead);
			lasp->next = star;
			star->prev = lasp;
		}


	//	printf("N : %d\tl : %e\t m : %e\tt : %e\tms : %e\n",N,lambda, cell->field.d * pow(2.0,-3.0*level) / mstars, dttilde/ tstartilde, mstars);

		initStar(cell, star, level, mstars, xc, yc, zc, cpu->npart[level-1]++, aexp );
		removeMfromgas(cell, star,  mstars/pow(2.0,-3.0*level) );
		return 0;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void createStars(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int curlevel){

	struct OCT  *nextoct;
	struct OCT  *curoct;
	struct CELL *curcell;

	REAL xc, yc, zc;
	REAL dx;
	REAL mstars = param->stars->mstars;

	int level, icell, i;
	int nstars;

	if(cpu->rank==0){
	
		srand(time(NULL));

		printf("\n");
		printf("================================\n");
		printf("   Starting Add Stars routine   \n");
		printf("================================\n");
	}



	for(level=curlevel;level<=param->lmax;level++) {

		nextoct=firstoct[level-1];
 		dx=pow(2.0,-level);
		nstars = 0;

		do {
			if(nextoct==NULL) continue;
			curoct=nextoct;
			nextoct=curoct->next;
			if(curoct->cpu!=cpu->rank) continue;

		      	for(icell=0;icell<8;icell++) {
				curcell = &curoct->cell[icell];

				if(curcell->child==NULL && testCond(curcell, dt, dx, param, aexp) ) {

					xc=curoct->x+( icell    & 1)*dx+dx*0.5; 
					yc=curoct->y+((icell>>1)& 1)*dx+dx*0.5;
					zc=curoct->z+( icell>>2    )*dx+dx*0.5; 										
					
					for (i=0;i< getNstars2create(curcell, param, dt, mstars, aexp, level) ;i++){
						addStar(curcell, level, xc, yc, zc, cpu, dt, param, aexp, mstars);
						nstars++;
					}
				}
			}
		}while(nextoct!=NULL);

#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&nstars,1,MPI_INT,MPI_SUM,cpu->comm);
#endif
	if(cpu->rank==0) printf("%d stars added on level %d \n", nstars, level);
	}
param->stars->n += nstars ;
if(cpu->rank==0) printf("%d stars in total\n",param->stars->n);
}
#endif
