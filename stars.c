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
		adot = sqrt(ok+(param->cosmo->om / a)+(or/ (a*a) )+ (param->cosmo->ov*a*a) );
		age = age + 1./adot;
	}
	REAL zage = az*age/n;
	REAL zage_Gyr = (977.8 /param->cosmo->H0)*zage;

	return zage_Gyr*1e9;
}

REAL rdm(REAL a, REAL b){
	return 	(rand()/(REAL)RAND_MAX ) * (b-a) + a ;
}

int fac(int x){
	return 	(x>1)?(fac(x-1)*x):1;
}

int gpoiss(REAL lambda, REAL max_k){

	int k=0;                
	REAL p = rdm(0,1); 	
	REAL P = exp(-lambda);  
	REAL sum=P;               
	if (sum>=p) return 0;     
	for (k=1; k<max_k; ++k) { 
		P*=lambda/(REAL)k;
		sum+=P;           
		if (sum>=p) break;
  	}

	return k;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL getFeedbackEgy(struct RUNPARAMS *param, REAL aexp){

	REAL rhoE 	= 3.7e-15; 	// erg.g-1 Kay 2002 // 4e48 erg.Mo-1 springel hernquist 2003 -> OK	
	rhoE 		*= 1.0e-7 * 1e3; 	// erg.g-1  en J.kg-1


	
	REAL H0 	= param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL rho0 	= 3.0 * pow(H0,2.0) * param->cosmo->om /(8.0*PI*NEWTON_G) ;

	REAL Mtot 	= rho0 * pow(param->cosmo->unit_l,3.0);
	REAL M 		= param->stars->mstars * Mtot;	

	REAL E 		= M * rhoE;


	REAL ts 	= 2.0/(H0 *sqrt(param->cosmo->om));
	REAL rs		= param->cosmo->unit_l;
	REAL es		= pow(rs/ts,2);

	E 		*=  pow(aexp,2.0) / es;

	return E;
}

void kineticFeedback(struct CELL *cell, struct RUNPARAMS *param, REAL E){

	int vnei[6],vcell[6];
	getcellnei(cell->idx, vnei, vcell); 

	struct CELL *curcell;
	struct OCT  *parent;
	E *= (1.-param->stars->feedback_frac) ;


	int icell;
	for(icell=0; icell<6;icell++){
		parent = cell2oct(cell);

		curcell = parent->nei[vnei[icell]];
	//	if (curcell->child!=NULL)

	//	->cell[vcell[icell]];		
	

	}



}

int feedback(struct CELL *cell, struct RUNPARAMS *param, REAL aexp, REAL t, REAL dt){

	REAL s8 	= 2e7 ;		// life time of a 8 M0 star
	s8 		*= 31556926; 	// years en s 

	REAL E		= getFeedbackEgy(param,aexp)  * param->stars->feedback_eff ;

	REAL H0 	= param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL ts 	= 2.0/(H0 *sqrt(param->cosmo->om));
	dt 		*= pow(aexp,2.0) * ts; 


	REAL t0;
	struct PART *nexp;
	struct PART *curp;

	nexp=cell->phead;
	if(nexp==NULL) return 0;
 	do{ 	curp=nexp;
		nexp=curp->next;

		if (curp->isStar){
			t0 = t - curp->age;
			if (  t0>0 && t0<2e8 ) {
				cell->field.E  += param->stars->feedback_frac * E * exp( -t0/s8 ) *dt/s8 ;    //thermal Feedback
				kineticFeedback(cell, param, E);
			}
		}
	}while(nexp!=NULL);

	return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initStar(struct CELL * cell, struct PART *star, struct RUNPARAMS *param, int level, REAL xc, REAL yc, REAL zc,int idx, REAL aexp, int is, REAL dt,REAL dx ) {

	star->next  = NULL;
	star->idx   = idx;
	star->level = level;
	star->is    = is;

	star->x = xc + rdm(-0.5,0.5) * dx;
	star->y = yc + rdm(-0.5,0.5) * dx;
	star->z = zc + rdm(-0.5,0.5) * dx;

	REAL vx = cell->field.u + rdm(-1.,1.) * cell->field.a;
	REAL vy = cell->field.v + rdm(-1.,1.) * cell->field.a;
	REAL vz = cell->field.w + rdm(-1.,1.) * cell->field.a;

	REAL vmax = pow(2.0,-level)/dt;

	star->vx = (fabs(vx)<vmax)?vx:vmax;
	star->vy = (fabs(vy)<vmax)?vy:vmax;
	star->vz = (fabs(vz)<vmax)?vz:vmax;

	star->mass  = param->stars->mstars;

	star->epot = 0.0;
	star->ekin = 0.5 * star->mass * (pow(vx,2) + pow(vy,2) + pow(vz,2));

	star->isStar = 1;
	star->age = a2t(param,aexp);
	star->rhocell = cell->field.d;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int testCond(struct CELL *cell, REAL dttilde, REAL dxtilde, struct RUNPARAMS *param, REAL aexp, int level){

	if (cell->child != NULL) return 0;

	int A = param->stars->overdensity_cond? 		(cell->field.d > param->stars->overdensity_cond * (param->cosmo->ob/param->cosmo->om) )  	: 1;
	int B = param->stars->density_cond?			(cell->field.d > param->stars->thresh)   							: 1;
	int C = 0?						cell->field.a/pow(2.,-level) > sqrt(6.*aexp * cell->gdata.d +1.) 				: 1;

	return 		A && B && C;
;

	

//	printf("%e\t%e\n", cell->field.d, n*param->stars->mstars);

//	return cell->field.d > param->amrthresh * (param->cosmo->ob/param->cosmo->om) / pow(2.0,-3.0*level);

//	printf("%e\t%e\n", (cell->gdata.d + 1.0)*pow(2.0,-3.0*level) , param->amrthresh);

//	return 	(cell->gdata.d+1.0)*pow(2.0,-3.0*level) > param->amrthresh;


/*	REAL DVcoarse 	=  pow(2.0,-3.0*param->lcoarse);
	REAL DV		=  pow(2.0,-3.0*param->level);

	REAL rhobarre 	= (param->cosmo->ob/param->cosmo->om);

	REAL N 		= param->amrthresh / pow(2.0,-3.0*param->lcoarse);

	REAL rhocrit 	= N * rhobarre * DVcoarse/DV ;
	*/


//	return cell->field.d > param->stars->overdensity_cond * (param->cosmo->ob/param->cosmo->om) * pow(aexp,-alpha0) ;

/*	REAL rhocrit 	= param->stars->density_cond * PROTON_MASS;

	REAL H0 	= param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL a3rhostar 	= pow(aexp,	3.0) *  3.0 * pow(H0,2.0) * param->cosmo->om /(8.0*PI*NEWTON_G);

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

	REAL tstars 	= param->stars->tcar * 31556926 / sqrt(cell->field.d / param->stars->thresh );
	REAL tstartilde = tstars / a2ts;

	REAL lambda = cell->field.d * pow(2.0,-3.0*level) / param->stars->mstars * dttilde/ tstartilde;

	return gpoiss(lambda, 0.75* cell->field.d * pow(2.0,-3.0*level) / param->stars->mstars );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addStar(struct CELL *cell, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dttilde, struct RUNPARAMS *param, REAL aexp, int is,  int nstars){

	struct PART *star = cpu->freepart;

	if (star==NULL){
		if( cpu->rank==0){
			printf("\n");
			printf("----------------------------\n");
			printf("No more memory for particles\n");
			printf("----------------------------\n");
		}
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

		cpu->npart[level-1]++;
		cpu->nstar[level-1]++;
	
		removeMfromgas(cell, star, param->stars->mstars/pow(2.0,-3.0*level) );
		initStar(cell, star, param, level, xc, yc, zc, param->stars->n+nstars, aexp, is, dttilde, pow(2.0,-level));	
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initThresh(struct RUNPARAMS *param,  REAL aexp){
	REAL k	= 0;
	if ( param->stars->density_cond > 0 ){
		k = 1;
	}else{
		k = -pow(aexp,3.0);
	}	

	REAL H0 		= param->cosmo->H0 *1000.0/1e6/PARSEC;
	REAL rhostar 		= 3.0 * pow(H0,2.0) * param->cosmo->om /(8.0*PI*NEWTON_G);

	REAL rhocrittilde 	= param->stars->density_cond * PROTON_MASS;

	param->stars->thresh    = k * rhocrittilde / rhostar;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void createStars(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){

	struct OCT  *curoct;
	struct OCT  *nextoct=firstoct[level-1];
	struct CELL *curcell;

	REAL xc, yc, zc;
	REAL dx = pow(2.0,-level);
	REAL mmax = 0;
	REAL mmin = 1e9;

	REAL t = a2t(param,aexp) ;

	int l,icell, ipart;
	int nstars = 0;
	int nsl = 0;
	int N;

	srand(time(NULL));
	initThresh(param, aexp);
	param->stars->mstars	= (param->cosmo->ob/param->cosmo->om) * pow(2.0,-3.0*param->lcoarse);
	param->cosmo->tphy	= a2t(param, aexp);


	if(cpu->rank==0){
		printf("\n");
		printf("================================\n");
		printf("   Starting Add Stars routine   \n");
		printf("================================\n");	
	}

	do {	if(nextoct==NULL) 		continue;
		curoct=nextoct;
		nextoct=curoct->next;
		if(curoct->cpu != cpu->rank) 	continue;

	      	for(icell=0;icell<8;icell++) {
			curcell = &curoct->cell[icell];

			if( testCond(curcell, dt, dx, param, aexp, level) ) {
				xc=curoct->x+( icell    & 1)*dx+dx*0.5; 
				yc=curoct->y+((icell>>1)& 1)*dx+dx*0.5;
				zc=curoct->z+( icell>>2    )*dx+dx*0.5; 										

				N = getNstars2create(curcell, param, dt, aexp, level);

			//	if(N) printf("N_Rho_Temp_Seuil_z\t%d\t%e\t%e\t%e\t%e\n", N, curcell->field.d, curcell->rfield.temp, param->stars->thresh,1.0/aexp - 1.0  );

				for (ipart=0;ipart< N; ipart++){
						addStar(curcell, level, xc, yc, zc, cpu, dt, param, aexp, is,nstars++);						
				}
			}

			feedback(curcell, param, aexp, t, dt);

			if (curcell->field.d > mmax) mmax = curcell->field.d;
			if (curcell->field.d < mmin) mmin = curcell->field.d;
		}
	}while(nextoct!=NULL);

#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&nstars,1,MPI_INT,   MPI_SUM,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&mmax,  1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&mmin,  1,MPI_DOUBLE,MPI_MIN,cpu->comm);
#endif

	param->stars->n += nstars ;
	if(cpu->rank==0) {
		printf("Mmax_thresh_overthresh\t%e\t%e\t%e\n", mmax, param->stars->thresh, param->stars->overdensity_cond * (param->cosmo->ob/param->cosmo->om) );
		printf("%d stars added on level %d \n", nstars, level);
		printf("%d stars in total\n",param->stars->n);
		printf("\n");
	}

	for (l = param->lcoarse; l<=param->lmax;l++){
#ifdef WMPI
		MPI_Allreduce(&cpu->nstar[l-1],&nsl,1,MPI_INT,   MPI_SUM,cpu->comm);
		MPI_Barrier(cpu->comm);
#endif
		if(cpu->rank==0 && nsl) {	printf("%d stars on level %d \n", nsl, l);	}
	}
	if(cpu->rank==0) {	printf("\n");}
}
#endif
