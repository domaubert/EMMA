#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prototypes.h"
#include "oct.h"

#ifdef WMPI
#include <mpi.h>
#endif

void cleanSNFBfield(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, int level){
#ifdef WRAD
	struct OCT  *nextoct=firstoct[level-1];

	do {	if(nextoct==NULL) 		continue;
	  struct OCT  *curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu != cpu->rank) 	continue;

		int icell;
	  for(icell=0;icell<8;icell++) {
	    struct CELL *curcell = &curoct->cell[icell];

			curcell->rfield.snfb = 0;
//			curcell->rfield.src   =0.;
//			curcell->rfieldnew.src=0.;

		}
	}while(nextoct!=NULL);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedback(struct CELL *cell,  REAL E){
#ifdef WRAD
 	struct OCT* oct = cell2oct(cell);
	int i;
	for(i=0;i<8;i++){
		struct CELL* curcell = &oct->cell[i];
		cell->rfield.snfb += E/8.;
	}
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kineticFeedback(struct CELL *cell, REAL E){

	float dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};
	float dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
	float dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

	REAL e = E/8.;

	struct OCT* oct = cell2oct(cell);

	int i;
	for(i=0;i<8;i++){
		struct CELL* curcell = &oct->cell[i];

		REAL V = SQRT(2.*e/curcell->field.d);

		printf("--> %e \t ",  curcell->field.u);
		curcell->field.u += dir_x[i] * 0.52532198881 * V;
		curcell->field.v += dir_y[i] * 0.52532198881 * V;
		curcell->field.w += dir_z[i] * 0.52532198881 * V;
		printf(" %e \t %e \n ",  curcell->field.u , V);
	
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL computeFeedbackEnergy(struct RUNPARAMS *param, REAL t0, REAL aexp, int level, REAL mstar){

	REAL s8 	 = param->stars->tlife;		// life time of a massive star (~20 Myr for 8 M0 star)
	s8 	*= 31556926; 	// years en s

	REAL dv = POW( POW(2.,-level) * aexp * param->unit.unit_l, 3.);   
	REAL E  = mstar*param->unit.unit_mass*SN_EGY*param->stars->feedback_eff/dv;

	E	     *= exp( -t0*31556926/s8 )/s8;

	return E;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int SN_TMP_PARAM = 1;

int feedback(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt){

#ifndef SNTEST
	cell->rfield.snfb = 0;

	if(param->stars->feedback_eff){

		int Nsn = 0;
		REAL t0;

		struct PART *nexp;
		struct PART *curp;

		nexp=cell->phead;
		if(nexp==NULL) return 0;
	 	do{ 	curp=nexp;
			nexp=curp->next;

			if (curp->isStar && curp->isStar < 3){
				t0 =  param->cosmo->tphy - curp->age;
				if(t0>=0){ // for inter-level communications
				  if ( t0 >= (param->stars->tlife*11)){
				    curp->isStar = 3;
				  }
				  else if ( t0 >= param->stars->tlife){
				    curp->isStar = 2; 

						REAL E = computeFeedbackEnergy(param, t0, aexp, level, curp->mass) ;

				    thermalFeedback(cell, E*(1.-param->stars->feedback_frac));
				    kineticFeedback(cell, E*(   param->stars->feedback_frac));

				    Nsn++;
				  }
				}
			}
		}while(nexp!=NULL);
		return Nsn;
	}
	return 0;

#else

	struct OCT* oct = cell2oct(cell);

	if (oct->x == 0.5 && oct->y == 0.5 && oct->z == 0.5 && cell->idx == 0){

		REAL in_yrs = param->unit.unit_t/MYR *1e6;

		REAL t = aexp * in_yrs;

//		if ( t >= param->stars->tlife && t < param->stars->tlife*11){
		if ( t >= param->stars->tlife && SN_TMP_PARAM){
				SN_TMP_PARAM = 0;

				printf("SN active at t = %e\n", t);

				REAL s8 	 = param->stars->tlife; //yrs

				REAL dv = POW( POW(2.,-level) * param->unit.unit_l, 3.); //m-3

			
				REAL msn = 8. *2e30 ; //kg

				REAL E  = msn *SN_EGY*param->stars->feedback_eff/dv; //J.m-3 

				

	//			E	     *= exp( -(t-param->stars->tlife)/s8 ) * dt*in_yrs/s8; //J.m-3

				printf("total feedback EGY %e\n",E);

		    kineticFeedback(cell, E*(   param->stars->feedback_frac));
		    thermalFeedback(cell, E*(1.-param->stars->feedback_frac));
		}
	}
	return 1;

#endif

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		FEEDBACK
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void supernovae(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){

	int Nsn = 0;
	struct OCT  *nextoct=firstoct[level-1];

	cleanSNFBfield(firstoct, param, cpu, level);

	do {	if(nextoct==NULL) 		continue;
	  struct OCT* curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu != cpu->rank) 	continue;

		int icell;
	  for(icell=0;icell<8;icell++) {
	    struct CELL *curcell = &curoct->cell[icell];
		
			Nsn += feedback(curcell, param, cpu, aexp, level, dt);

		}
	}while(nextoct!=NULL);

#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&Nsn,   1,MPI_INT,   MPI_SUM,cpu->comm);
#endif

	if(cpu->rank==RANK_DISP) {	printf("%d\tActive SN\n",Nsn);}

}


