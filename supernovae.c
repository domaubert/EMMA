#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prototypes.h"
#include "oct.h"

#ifdef WMPI
#include <mpi.h>
#endif

#ifdef SUPERNOVAE
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

//			curcell->rfield.snfb = 0;
//			curcell->rfield.src   =0.;
//			curcell->rfieldnew.src=0.;

		}
	}while(nextoct!=NULL);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedback(struct CELL *cell,  REAL E){
#ifdef WRAD
    cell->field.E += E;
    cell->field.p += E*(GAMMA-1.);
    //cell->rfield.snfb =1;
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kineticFeedback(struct CELL *cell, REAL E){

    int type = 0;   //oct
//    int type = 1;   //cardinal
//    int type = 2;   //cardinal

    switch(type){

        case 0:{
            struct OCT* oct = cell2oct(cell);

            int i;
            for(i=0;i<8;i++){
                struct CELL* curcell = &oct->cell[i];

                REAL e = E/8.;
                REAL V = SQRT(2.*e/curcell->field.d);

                float dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};
                float dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
                float dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

                curcell->field.u += dir_x[i] * V /2.; // 1/2 = cos45 * cos45
                curcell->field.v += dir_y[i] * V /2.;
                curcell->field.w += dir_z[i] * V /2.;

            }
        }
        break;

        case 1:{

            struct CELL* neicell[6];
            getneicell(cell, neicell);

            int i;
            for(i=0;i<6;i++){
                struct CELL* curcell = neicell[i];
                printf("neicell = %p\n",neicell[i]);

                if (curcell!=NULL) {
                    REAL dir_x[]={-1., 1., 0., 0., 0., 0.};
                    REAL dir_y[]={ 0., 0.,-1., 1., 0., 0.};
                    REAL dir_z[]={ 0., 0., 0., 0.,-1., 1.};

                    REAL e = E/6.;
                    REAL dv = SQRT(2.*e/curcell->field.d);

                    curcell->field.u += dir_x[i]*dv;
                    curcell->field.v += dir_y[i]*dv;
                    curcell->field.w += dir_z[i]*dv;

//                    curcell->field.u += 1;
//                    curcell->field.v += 1;
//                    curcell->field.w += 1;



//                REAL u = curcell->field.u;
//                REAL v = curcell->field.v;
//                REAL w = curcell->field.w;
//
//                REAL x = cell2oct(neicell[i])->x *64;
//                REAL y = cell2oct(neicell[i])->y *64;
//                REAL z = cell2oct(neicell[i])->z *64;
//
//                REAL dx = dir_x[i];
//                REAL dy = dir_y[i];
//                REAL dz = dir_z[i];

                REAL stop = 0;
                }
            }
        }
        break;

        case 2:{

            struct OCT* curoct = cell2oct(cell);
            struct CELL* curcell = NULL;

            REAL e = E/6.;

            //xm
            curcell =  &(curoct->nei[0]->child->cell[1]);
            REAL dv = SQRT(2.*e/curcell->field.d);
            curcell->field.u += -dv;

            //xp
            curcell =  &(curoct->cell[1]);
            curcell->field.u += dv;


            curcell =  &(curoct->nei[2]->child->cell[2]);
            curcell->field.v += -dv;

            curcell =  &(curoct->cell[2]);
            curcell->field.v += dv;

            curcell =  &(curoct->nei[4]->child->cell[4]);
            curcell->field.w += -dv;

            curcell =  &(curoct->cell[4]);
            curcell->field.w += dv;

        }
        break;
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL computeFeedbackEnergy(struct RUNPARAMS *param, REAL t0, REAL aexp, int level, REAL mstar){
#ifdef STARS
	REAL s8 	 = param->stars->tlife;		// life time of a massive star (~20 Myr for 8 M0 star)
	s8 	*= 31556926; 	// years en s

	REAL dv = POW( POW(2.,-level) * aexp * param->unit.unit_l, 3.);
	REAL E  = mstar*param->unit.unit_mass*SN_EGY*param->stars->feedback_eff/dv;

	E	     *= exp( -t0*31556926/s8 )/s8;

	return E;
#endif // STARS
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int SN_TMP_PARAM = 1;
int feedback(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt){
#ifndef SNTEST
#ifdef PIC
	cell->rfield.snfb = 0;

	if(param->stars->feedback_eff){

		int Nsn = 0;
		REAL t0;

		struct PART *nexp;
		struct PART *curp;

		nexp=cell->phead;
		if(nexp==NULL) return 0;
	 	do{ curp=nexp;
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

				    thermalFeedback(cell, E*(1.-param->sn->feedback_frac));
				    kineticFeedback(cell, E*(   param->sn->feedback_frac));

				    Nsn++;
				  }
				}
			}
		}while(nexp!=NULL);
		return Nsn;
	}
	return 0;
#endif // PIC
#else // ifdef SNTEST

	struct OCT* oct = cell2oct(cell);

	if (oct->x == 0 && oct->y == 0 && oct->z == 0 && cell->idx == 0){

		REAL in_yrs = param->unit.unit_t/MYR *1e6;
		REAL t = aexp * in_yrs;

		if(!SN_TMP_PARAM){
            printf("======================================\n");
            printf("===SN EXPLODE=========================\n");
            printf("======================================\n");
        }

//		if ( t >= param->stars->tlife && t < param->stars->tlife*11){
		if ( t >= LIFETIME_OF_STARS_IN_TEST && SN_TMP_PARAM ){

// src http://cdsads.u-strasbg.fr/abs/2009A%26A...495..389B
// S20
// mass = 7.6e6
// L = 4.6e43 erg/s
// L = 4.6e43 / 1e7 J/s
// L = 4.6e43 / 1e7 / (29.609895722*1.6022e-19)
//
//S100
//mass = 9.5e8 M0
//L = 5.75e45 erg/s

			SN_TMP_PARAM = 0;
			printf("SN active at t = %e\n", t);


			REAL dx = POW(2.,-level) * param->unit.unit_l; //m
			REAL dv = POW( dx, 3.); //m3

			REAL msn = 7.6e6 * 2e30 ; //kg

			REAL E  = msn *SN_EGY/dv; //J.m-3

//				REAL s8 	 = param->stars->tlife; //yrs
//				E	     *= exp( -(t-param->stars->tlife)/s8 ) * dt*in_yrs/s8; //J.m-3

			printf("total feedback EGY %e erg\n",E*dv*1e7);

			dv = POW( 2., -3.*level);
			printf("eblast= %e erg   \n", 		E*dv);

		    kineticFeedback(cell, E);
		    //thermalFeedback(cell, E);
		}
	}
	return 1;

#endif

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		FEEDBACK
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void checksupernovae(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is, int *N){

	struct OCT  *nextoct=firstoct[level-1];

	do {	if(nextoct==NULL) 		continue;
	  struct OCT* curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu != cpu->rank) 	continue;

      int icell;
	  for(icell=0;icell<8;icell++) {
            struct CELL *curcell = &curoct->cell[icell];

            if (curcell->field.u || curcell->field.v || curcell->field.w ){

                REAL dx = POW(2.0,-curoct->level);

                REAL xc=curoct->x+( curcell->idx    & 1)*dx;
                REAL yc=curoct->y+((curcell->idx>>1)& 1)*dx;
                REAL zc=curoct->z+( curcell->idx>>2    )*dx;

                printf("%f %f %f\t",xc,yc,zc);

                REAL vx = curcell->field.u;
                REAL vy = curcell->field.v;
                REAL vz = curcell->field.w;

                printf("%f %f %f\n",vx,vy,vz);

                (*N)++;
                REAL stop = 0;

            }
		}
	}while(nextoct!=NULL);
}


void supernovae(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){
	if(cpu->rank==RANK_DISP) printf("Supernovae\n");

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


#endif //SUPERNOVAE
