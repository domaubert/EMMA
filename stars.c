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
//		FEEDBACK
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedback(struct CELL *cell, struct RUNPARAMS *param, REAL t0, REAL aexp, int level, REAL dt ){

	REAL s8 	 = 2e7;		// life time of a 8 M0 star
	     s8 	*= 31556926; 	// years en s

	REAL dv 	= pow( pow(2.,-level) * aexp * param->cosmo->unit_l, 3.); 

	REAL E 		= param->stars->Esnfb/dv  * param->stars->feedback_frac;
	     E	       *= exp( -t0/s8 )/s8;

	cell->rfield.snfb += E;	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int feedback(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt){

	cell->rfield.snfb = 0;

	if(param->stars->feedback_eff){

		int Nsn = 0;
		REAL t0;

		struct PART *nexp; //
		struct PART *curp;

		nexp=cell->phead;
		if(nexp==NULL) return 0;
	 	do{ 	curp=nexp;
			nexp=curp->next;

			if (curp->isStar && curp->isStar < 3){
				t0 =  param->cosmo->tphy - curp->age;
				if(t0>=0){ // for inter-level communications
					if ( t0 >= param->stars->tlife 	)  	curp->isStar = 2; 
					if ( t0 >= 2e8 			)  	curp->isStar = 3;
	
					thermalFeedback(cell, param, t0*31556926, aexp, level, dt );
					Nsn++;
				}
			}
		}while(nexp!=NULL);
		return Nsn;
	}
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		STARS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initStar(struct CELL * cell, struct PART *star, struct RUNPARAMS *param, int level, REAL xc, REAL yc, REAL zc,int idx, REAL aexp, int is, REAL dt,REAL dx ) {

	star->next  = NULL;
	star->idx   = -1;
	star->level = level;
	star->is    = is;

	star->x = xc + rdm(-0.5,0.5) * dx;
	star->y = yc + rdm(-0.5,0.5) * dx;
	star->z = zc + rdm(-0.5,0.5) * dx;


	star->vx = cell->field.u;
	star->vy = cell->field.v;
	star->vz = cell->field.w;

// Random component
/*	REAL r 		= rdm(0,1) * cell->field.a ;
	REAL theta   	= acos(rdm(-1,1));
	REAL phi 	= rdm(0,2*PI);

	star->vx += r * sin(theta) * cos(phi);
	star->vy += r * sin(theta) * sin(phi);
	star->vz += r * cos(theta) ;
*/
	star->mass  = param->stars->mstars;

	star->epot = 0.0;
	star->ekin = 0.5 * star->mass * (pow(star->vx,2) + pow(star->vy,2) + pow(star->vz,2));

	star->isStar = 1;
	star->age = param->cosmo->tphy;
	star->rhocell = cell->field.d;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int testCond(struct CELL *cell, REAL dttilde, REAL dxtilde, struct RUNPARAMS *param, REAL aexp, int level){

	if (cell->child != NULL) return 0;

	int A = 	cell->field.d > param->stars->thresh;
	int B = 0?	cell->field.a/pow(2.,-level) > sqrt(6.*aexp * cell->gdata.d +1.) 				: 1;

	return 		A && B;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void conserveField(struct Wtype *field, struct RUNPARAMS *param, struct PART *star, REAL dx, REAL aexp){

	REAL drho = param->stars->mstars / pow(dx,3.);
	struct Utype U;
	struct Wtype W;

	memcpy(&W,field,sizeof(struct Wtype));

	W2U(&W, &U);

//	density
	U.d    -= drho;

//	momentum
	U.du -= star->vx * drho;		
	U.dv -= star->vy * drho;		
	U.dw -= star->vz * drho;		

//	internal energy
//	???????????????

	U2W(&U, &W);

//	total energy
	getE(&(W));
	W.a=sqrt(GAMMA*W.p/W.d);

	memcpy(field,&W,sizeof(struct Wtype));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int getNstars2create(struct CELL *cell, struct RUNPARAMS *param, REAL dttilde, REAL aexp, int level){
	
	REAL gas_efficiency = 1.;	// maybe need to be passed in param??

	REAL tstars 	= param->stars->tcar * 31556926 / sqrt(cell->field.d / param->stars->thresh );
	REAL tstartilde = tstars / pow(aexp,2)/param->unit.unit_t;

	REAL M_in_cell 	= cell->field.d * pow(2.0,-3.0*level);

	REAL lambda 	= gas_efficiency * M_in_cell / param->stars->mstars * dttilde/ tstartilde; // Average number of stars created

	int N 		= gpoiss(lambda);
	
	if(N * param->stars->mstars >= M_in_cell ) N = M_in_cell / param->stars->mstars ;

	//if (N) printf("N %d \t M in cell %e \t Mstars %e \t dt %e\t dtt %e\t E %e \t l %e  \n",N,M_in_cell, param->stars->mstars, dttilde, tstartilde, cell->rfieldnew.eint, lambda);
	//if (N) printf("rfield.eint %e \t rfieldnew.eint %e \n ", cell->rfield.eint, cell->rfieldnew.eint);
	//printf("M in cell %e \t Mstars %e \t dt %e\t dtt %e\t d %e \t l %e  \n N=%e",M_in_cell, param->stars->mstars, dttilde, tstartilde, cell->field.d, lambda,N);
	return N;
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
			star->next =NULL;
		}else{
			struct PART *lasp 	= findlastpart(cell->phead);
			lasp->next 		= star;
			star->prev 		= lasp;
			star->next=NULL;
		}

		cpu->npart[level-1]++;
		cpu->nstar[level-1]++;
	
		REAL dx = pow(2.0,-level);

		initStar(cell, star, param, level, xc, yc, zc, param->stars->n+nstars, aexp, is, dttilde, dx);	

		conserveField(&(cell->field   ), param, star,  dx, aexp);
		conserveField(&(cell->fieldnew), param, star,  dx, aexp);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initThresh(struct RUNPARAMS *param,  REAL aexp){

	REAL k =(param->stars->density_cond >0)? 1 : -pow(aexp,3.0);
	
	REAL rhostar		= param->unit.unit_mass/pow(param->unit.unit_l,3) ;
	REAL rhocrittilde 	= param->stars->density_cond * PROTON_MASS;

	param->stars->thresh    = fmax( k * rhocrittilde / rhostar, param->stars->overdensity_cond * (param->cosmo->ob/param->cosmo->om));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void createStars(struct OCT **octList, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){

	struct OCT  *curoct;
	struct CELL *curcell;
	struct PART *nexp;
	struct PART *curp;

	REAL xc, yc, zc;
	REAL dx = pow(2.0,-level);
	REAL mmax = 0;
	REAL mmin = 1e9;

	int l,icell, ipart, icpu;
	int nstars = 0;
	int nsl = 0;
	int N;
	int Nsn=0;


	initThresh(param, aexp);
	param->cosmo->tphy	= a2t(param, aexp);

	int iOct;
	for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
		struct OCT *curoct=octList[iOct]; 
	  for(icell=0;icell<8;icell++) {
			curcell = &curoct->cell[icell];

			if( testCond(curcell, dt, dx, param, aexp, level) ) {
				xc=curoct->x+( icell    & 1)*dx+dx*0.5; 
				yc=curoct->y+((icell>>1)& 1)*dx+dx*0.5;
				zc=curoct->z+( icell>>2    )*dx+dx*0.5; 										

				N = getNstars2create(curcell, param, dt, aexp, level);

			//	if(N) printf("N_Rho_Temp_Seuil_z\t%d\t%e\t%e\t%e\t%e\n", N, curcell->field.d, curcell->rfield.temp, param->stars->thresh,1.0/aexp - 1.0  );

				for (ipart=0;ipart< N; ipart++){
						addStar(curcell, level, xc, yc, zc, cpu, dt, param, aexp, is, ++nstars );	
				}
			}

			Nsn += feedback(curcell, param, cpu, aexp,level,dt);
			mmax = fmax(curcell->field.d, mmax);
		}
	  }


	// set IDX

	int *nStarCpu = (int *)calloc(cpu->nproc, sizeof(int));
#ifdef WMPI
	MPI_Allgather(&nstars, 1, MPI_INT, nStarCpu, 1, MPI_INT,cpu->comm);
#endif
	int nbefore = 0;	for(icpu=0; icpu<cpu->rank; icpu++) nbefore += nStarCpu[icpu];
	int istar   = 0;


	for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
		struct OCT *curoct=octList[iOct]; 
  	for(icell=0;icell<8;icell++) {
			curcell = &curoct->cell[icell];
			nexp=curoct->cell[icell].phead;
			if(nexp!=NULL) 
			do{	
				curp=nexp;
				nexp=curp->next;
				if(curp->isStar && curp->idx==-1){
					curp->idx = param->stars->n + nbefore +  istar;
				//	printf("\t param->stars->n = %d\tnbefore = %d\tistar = %d\tIDX = %d\n",param->stars->n,nbefore,istar,curp->idx);
					istar++;
				}
			}while(nexp!=NULL);
		}
	}





#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&nstars,1,MPI_INT,   MPI_SUM,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&mmax,  1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&Nsn,   1,MPI_INT,   MPI_SUM,cpu->comm);

	for (l = param->lcoarse; l<=param->lmax;l++){
	  MPI_Allreduce(&cpu->nstar[l-1],&nsl,1,MPI_INT,   MPI_SUM,cpu->comm);
	  MPI_Barrier(cpu->comm);
	  //if(cpu->rank==0 && nsl) {	printf("%d stars on level %d \n", nsl, l);	}
	}
#endif

	param->stars->n += nstars ;
	if(cpu->rank==0) {
		printf("STARS.......................................\n");
		printf("Mmax_thresh\t%e\t%e\n", mmax, param->stars->thresh );
		printf("\t%d\tstars added on level %d \n", nstars, level);
		printf("\t%d\tstars in total\n",param->stars->n);
		printf("\t%d\tActive SN\n",Nsn);
	}



}
#endif
