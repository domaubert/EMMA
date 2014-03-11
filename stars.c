#ifdef STARS 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 

#include "prototypes.h"
#include "particle.h"
#include "oct.h"


#include "cic.h"
#include "segment.h"


#ifdef WMPI
#include <mpi.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initStar(struct CELL * cell, struct PART *star, struct PART *prev , int level, REAL m ,REAL xc, REAL yc, REAL zc,int idx, REAL aexp) {

	star->next = NULL;
	star->prev = prev;
	
	star->idx   = idx;

	star->level = level;
	star->is    = 0;

	star->x = xc;
	star->y = yc;
	star->z = zc;

	REAL n = 100000;
	star->vx = cell->field.u + ( ((float)(rand()%n)*2.0-1.0 )/n *cell->field.a;
	star->vy = cell->field.v + ( ((float)(rand()%n)*2.0-1.0 )/n *cell->field.a;
	star->vz = cell->field.w + ( ((float)(rand()%n)*2.0-1.0 )/n *cell->field.a;

	star->mass = m;

	star->epot = 0.0;
	star->ekin = 0.0 ;

	star->isStar = 1;
	star->age = 1.0/aexp-1.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int countStars(struct CELL * cell){
	int nstars =0;

	struct PART * pcur;
	struct PART * pnext;

	pnext = cell->phead;
	
	do {
		pcur=pnext;
		pnext = pcur->next;
		if (pcur->isStar) nstars++;
	}while(pnext != NULL);

	return nstars;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int testCond(struct CELL *cell,REAL dttilde, REAL dxtilde,struct RUNPARAMS *param,  REAL aexp){

	REAL H0 = param->cosmo->H0 *1000/1e6/PARSEC;
	REAL om = param->cosmo->om;
	REAL L  = param->cosmo->unit_l;

	int cond_overdensity = 0;
	int cond_density = 0;
	int cond_jeans = 0;
	REAL jeans = 1.0/(2.0*PI*sqrt(6.0*aexp*(cell->gdata.d+1.0) ));;

//------OVERDENSITY--------------------------------------------------------------------------------------------//

	REAL cond = param->stars->overdensity_cond;	// stars formation begin at minimum cond times the mean density 
	REAL rhotilde 	= cell->field.d * (param->cosmo->om/ param->cosmo->ob) ;	// 
//	printf("%lf\n",rhotilde);
	cond_overdensity = (rhotilde > cond);

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


//------JEANS--------------------------------------------------------------------------------------------------//


	if (cond_density){
		cond_jeans = (dxtilde/cell->field.a > jeans);
//			if (!cond_jeans) abort();  					////////////	A CHECKER
	}else{return 0;}

//------RANDOM-------------------------------------------------------------------------------------------------//	
*/
	if(cond_overdensity) {
		REAL p = 1.0 - exp(- param->stars->eff * dttilde/jeans);
		REAL random = ((float)(rand()%100000))/100000;
		return (random < p)?1:0;	
		return 1;
	}else{return 0;}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL getdrho(struct CELL *cell, REAL dttilde, struct RUNPARAMS *param, REAL aexp){

	REAL Ho = param->cosmo->H0 *1000/1e6/PARSEC; // (s)
	REAL tstar = 2.0 /( Ho * sqrt(param->cosmo->om) ); 
	REAL dt = dttilde * pow(aexp,2) * tstar ; 

	REAL tcar = param->stars->t_car * 31556926.0; // (s) 
	REAL eff = param->stars->eff;
	
	return (cell->field.d) * 2.0 * PI / sqrt(6.0*aexp*(cell->gdata.d+1.0)) * eff * dt/tcar ; 
//	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void removeMfromgas(struct CELL * cell, REAL drho){
	cell->fieldnew.d -= drho;
	cell->field.d    -= drho;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addStar(struct CELL *cell, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dt, struct RUNPARAMS *param, REAL aexp, REAL drho){

	struct PART * lasp;
	struct PART * prev;
	struct PART * star = cpu->freepart;

	if (cpu->freepart->next){	
		cpu->freepart->next->prev = NULL;
		cpu->freepart = cpu->freepart->next;
	}else{
		printf("no more memory for particles\n");
		abort();
	}

	if (cell->phead!=NULL){
		lasp = findlastpart(cell->phead);
		prev = lasp;
		lasp->next  = star;
	}else{
		prev = NULL;
		cell->phead = star;
	}	

	initStar(cell, star,  prev, level, drho*pow(2.0,-3*level), xc, yc, zc, cpu->npart[level-1]++, aexp );
	removeMfromgas(cell, drho);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void accretion(struct CELL * cell, struct RUNPARAMS *param, REAL dt, int nstarsincell, REAL dx, REAL aexp){

	REAL drho = getdrho(cell, dt, param, aexp );		
	REAL m = drho * pow(dx,3) / nstarsincell;
	REAL fact ;

	struct PART * pcur;
	struct PART * pnext = cell->phead;

	do {
		pcur=pnext;
		pnext = pcur->next;
		if (pcur->isStar) {
			fact = sqrt(m / pcur->mass );

			pcur->vx += fact * cell->field.u;
			pcur->vy += fact * cell->field.v;
			pcur->vz += fact * cell->field.w;
	
			pcur->mass += drho * pow(dx,3) / nstarsincell;	
		}
	}while(pnext != NULL);

	removeMfromgas(cell, drho);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void createStars(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp){

	struct OCT  *nextoct;
	struct OCT  *curoct;
	struct CELL *curcell;

	REAL xc, yc, zc;
	REAL dx;
	REAL drho;
	REAL drhomin = (param->cosmo->ob/ param->cosmo->om) / pow(2,3*param->lmax);

	int level, icell;
	int nstars = 0;
	int nstarsincell;

	if(cpu->rank==0){
		printf("\n");
		printf("================================\n");
		printf("   Starting Add Stars routine   \n");
		printf("================================\n");
	}

	if(cpu->rank==0) srand(time(NULL));

	for(level=param->lcoarse;level<=param->lmax;level++) {

		nextoct=firstoct[level-1];
		 

		dx=1./pow(2.0,level);		
		nstars = 0;

		do {
			if(nextoct==NULL) continue;
			curoct=nextoct;
			nextoct=curoct->next;
			if(curoct->cpu!=cpu->rank) continue;

		      	for(icell=0;icell<8;icell++) {
				curcell = &curoct->cell[icell];

				if(curcell->child == NULL && testCond(curcell,dt, dx, param, aexp) ) {

					drho = getdrho(curcell, dt, param, aexp);

					if (drho > drhomin){

						xc=curoct->x+( icell    & 1)*dx+dx*0.5; 
						yc=curoct->y+((icell>>1)& 1)*dx+dx*0.5;
						zc=curoct->z+( icell>>2    )*dx+dx*0.5; 
					
	/*					nstarsincell = countStars(curcell);
						if(nstarsincell){
							accretion(curcell, param, dt, nstarsincell, dx);					
						}else{
							addStar(curcell, level, xc, yc, zc, cpu, dt, param);
							nstars++;
						}
	*/

					
						
						addStar(curcell, level, xc, yc, zc, cpu, dt, param, aexp, drho);
						nstars++;
					}else{
						printf("blocked stars\n");
						abort();
					}
				}			

			}

		}while(nextoct!=NULL);
#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&nstars,1,MPI_INT,MPI_SUM,cpu->comm);
#endif
	if(cpu->rank==0) printf("%d stars added on level %d \n", nstars, level);
	

	}
if(cpu->rank==0) printf("\n");
}

#endif
