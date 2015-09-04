// ----------------------------------------------------------
// ----------------------------------------------------------
/** \file supernovae.c
  * \brief contain functions for supernovae and stellar wind
  * \author Nicolas Deparis
  *
  * Need the SUPERNOVAE preprocessor flag
  *
  * More documentation can be found on the wiki :
  * https://github.com/domaubert/EMMA/wiki/Supernovae
  */
// ----------------------------------------------------------
// ----------------------------------------------------------

#ifdef SUPERNOVAE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef WMPI
#include <mpi.h>
#endif

#include "prototypes.h"
#include "oct.h" //cell2oct
#include "hydro_utils.h" // W2U and U2W
#include "convert.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedbackCell(struct CELL *cell,  REAL E){
// ----------------------------------------------------------
/// Inject an energy "E" in the cell "cell" on thermal form.
//----------------------------------------------------------
#ifdef SNTEST
  printf("injecting Energy in thermal form within a cell\n");
#endif // SNTEST

  cell->field.E += E;
  cell->field.p += E*(GAMMA-1.);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedbackOct(struct CELL *cell,  REAL E){
// ----------------------------------------------------------
/// Inject an energy "E" in all cells of the oct contening
/// the cell "cell" uniformly on thermal form.
// ----------------------------------------------------------
#ifdef SNTEST
  printf("injecting Energy in thermal form within a oct\n");
#endif // SNTEST

    struct OCT* oct = cell2oct(cell);
    int i;
    for(i=0;i<8;i++){
        struct CELL* curcell = &oct->cell[i];
        REAL e = E/8.;
        curcell->field.E += e;
        curcell->field.p += e*(GAMMA-1.);
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kineticFeedback(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, REAL aexp, int level, REAL E){
// ----------------------------------------------------------//
/// Inject an energy "E" in all cells of the oct contening
/// the cell "cell" on kinetic form, radially to the center
/// of the oct and uniformly in all cells
// ----------------------------------------------------------//

#ifdef SNTEST
  //REAL msn = param->unitary_stars_test->mass * SOLAR_MASS /param->unit.unit_mass;
  //REAL mtot_feedback = msn * param->sn->ejecta_proportion;
  REAL mtot_feedback =0;
#else
  REAL mtot_feedback = curp->mass* param->sn->ejecta_proportion;
#endif // SNTEST

  struct OCT* oct = cell2oct(cell);

  int i;
  for(i=0;i<8;i++){
    struct CELL* curcell = &oct->cell[i];
    REAL dv = POW(2.,-3*level);// curent volume

    REAL e = E/8.; //uniform energy distribution

    REAL me = mtot_feedback/8.;

    if (mtot_feedback==0){
      // if there's no ejecta, we injecte the energy directly to the gas
      me = curcell->field.d * dv;
    }

    REAL rho_i = curcell->field.d; //initial density
    REAL rho_e = me/dv; // density ejecta

    REAL vxi = curcell->field.u; // initial velocity
    REAL vyi = curcell->field.v;
    REAL vzi = curcell->field.w;

    REAL ve = SQRT(2.*e/rho_e);// ejecta velocity

    REAL dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};// diagonal projection
    REAL dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
    REAL dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

#ifdef SNTEST
  REAL  x_src=param->unitary_stars_test->src_pos_x,
        y_src=param->unitary_stars_test->src_pos_y,
        z_src=param->unitary_stars_test->src_pos_z;

  if ( (x_src==0) && (y_src==0) && (z_src==0)) {
      REAL x[]={1., 0., 0., 0., 0., 0., 0., 0.};// diagonal projection
      REAL y[]={1., 0., 0., 0., 0., 0., 0., 0.};
      REAL z[]={1., 0., 0., 0., 0., 0., 0., 0.};
      int i;
      for (i=0;i<8; i++){
        dir_x[i]=x[i];
        dir_y[i]=y[i];
        dir_z[i]=z[i];
      }
    }
#endif // SNTEST


#ifdef PIC
    REAL vx0 = curp->vx;
    REAL vy0 = curp->vy;
    REAL vz0 = curp->vz;
#else
		REAL vx0 = curcell->field.u;
	  REAL vy0 = curcell->field.v;
	  REAL vz0 = curcell->field.w;
#endif // PIC

    REAL vxe = vx0 + ve*dir_x[i]/SQRT(3.); // projection on axis in the particle framework
    REAL vye = vy0 + ve*dir_y[i]/SQRT(3.);
    REAL vze = vz0 + ve*dir_z[i]/SQRT(3.);

#ifdef PIC
    if(mtot_feedback!=0)  curp->mass -= me; // new particle mass
#else
    curcell->field.d -= rho_e;
#endif // PIC

    if(mtot_feedback==0) {
      rho_i -= rho_e;
      //TODO verif if needed;
     // curcell->field.d -= rho_e;
    }

    curcell->field.u = (vxi*rho_i + vxe*rho_e)/(rho_i+rho_e); //new velocity
    curcell->field.v = (vyi*rho_i + vye*rho_e)/(rho_i+rho_e);
    curcell->field.w = (vzi*rho_i + vze*rho_e)/(rho_i+rho_e);

    curcell->field.d += rho_e; //new density

    //Energy conservation
#ifdef DUAL_E
    struct Utype U; // conservative field structure
    W2U(&curcell->field, &U); // primitive to conservative
#ifdef DUAL_E
    U.eint*=1.+rho_e/curcell->field.d; // compute new internal energy
#endif // DUAL_E

    U2W(&U, &curcell->field); // back to primitive
#else
    curcell->field.p*=1.+rho_e/curcell->field.d; // compute new internal energy
#endif
    getE(&curcell->field); //compute new total energy
    curcell->field.p=FMAX(curcell->field.p,PMIN);
    curcell->field.a=SQRT(GAMMA*curcell->field.p/curcell->field.d); // compute new sound speed
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kineticFeedback_test(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, REAL aexp, int level, REAL E){

  struct OCT* oct = cell2oct(cell);
  int i;
  for(i=0;i<8;i++){
    struct CELL* curcell = &oct->cell[i];

    REAL e = E/8.;
    REAL v = SQRT(2.*e/curcell->field.d);

    REAL dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};// diagonal projection
    REAL dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
    REAL dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

    curcell->field.u = v*dir_x[i]/SQRT(3.);
    curcell->field.v = v*dir_y[i]/SQRT(3.);
    curcell->field.w = v*dir_z[i]/SQRT(3.);

    //Energy conservation
    curcell->field.p=FMAX(curcell->field.p,PMIN);
    getE(&curcell->field); //compute new total energy
    curcell->field.a=SQRT(GAMMA*curcell->field.p/curcell->field.d); // compute new sound speed
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL computeFeedbackEnergy(struct RUNPARAMS *param, REAL aexp, int level, REAL mstar){
// ----------------------------------------------------------//
/// Compute the total feedback energy
// ----------------------------------------------------------//
  REAL egy = param->sn->sn_egy; // j/kg
  egy *= mstar/POW(0.5,3*level) * POW(aexp,2.)/POW(param->unit.unit_v,2.); // j/m3 in code unit
  return egy;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SNTEST
#ifdef PIC
int feedback(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt){
// ----------------------------------------------------------//
/// Scan all particles of cell "cell" and look for supernovae
/// If SN found, compute energy and inject it
// ----------------------------------------------------------//

  int Nsn = 0;

  struct PART *nexp=cell->phead;

  if(nexp==NULL) return 0;
  do{
    struct PART *curp=nexp;
    nexp=curp->next;

    if (curp->isStar==2 || curp->isStar==3){ // if curp is in SN state

      REAL E = computeFeedbackEnergy(param, aexp, level, curp->mass);

      //printf("Energy injected =%e mass=%e aexp=%e cellE=%e\n",E,curp->mass,aexp, cell->field.E);

      if(param->sn->feedback_frac){
        // oct feedback
        // if there is kinetic feedback the energy is injected in a oct
        thermalFeedbackOct(cell, E*(1.-param->sn->feedback_frac));
        kineticFeedback(param, cell,curp,aexp,level, E*param->sn->feedback_frac);
      }else{
        // cell feedback
        // if there is only thermal feedback the energy can be injected in just a cell
        thermalFeedbackCell(cell, E*(1.-param->sn->feedback_frac));
      }
      Nsn++;
    }
  }while(nexp!=NULL);
  return Nsn;
}
#endif // PIC

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#else // ifdef SNTEST
int SN_TMP_PARAM = 1;
int feedback(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt){
// ----------------------------------------------------------//
/// For sedov test
/// A supernovae explode in a uniform medium
// ----------------------------------------------------------//

	struct OCT* oct = cell2oct(cell);
  REAL dx = POW(2.0,-level);

  REAL  x_src=param->unitary_stars_test->src_pos_x,
        y_src=param->unitary_stars_test->src_pos_y,
        z_src=param->unitary_stars_test->src_pos_z;


  REAL x=(oct->x+( cell->idx    & 1)*dx);
  REAL y=(oct->y+((cell->idx>>1)& 1)*dx);
  REAL z=(oct->z+( cell->idx>>2    )*dx);

	/*
	REAL x=(oct->x+( cell->idx    & 1)*dx+dx*0.5);
	REAL y=(oct->y+((cell->idx>>1)& 1)*dx+dx*0.5);
	REAL z=(oct->z+( cell->idx>>2    )*dx+dx*0.5);
  */

  REAL dX = x-x_src;
  REAL dy = y-y_src;
  REAL dz = z-z_src;

	REAL R=SQRT(dX*dX+dy*dy+dz*dz);
  REAL rmax = 0.6 * POW(0.5,param->lcoarse);

	if (R<=rmax){



//    if ( (x==x_src) && (y==y_src) && (z==z_src) && cell->child==NULL){


		REAL in_yrs = param->unit.unit_t/MYR *1e6;
		REAL t = aexp * in_yrs;


		if ( t >= param->unitary_stars_test->lifetime && SN_TMP_PARAM ){

      printf("======================================\n");
      printf("===SN EXPLODE=========================\n");
      printf("======================================\n");

      //printf("x_src=%e y_src=%e, z_src=%e\n", param->unitary_stars_test->src_pos_x, param->unitary_stars_test->src_pos_y, param->unitary_stars_test->src_pos_z);

      //REAL msn = param->unitary_stars_test->mass * SOLAR_MASS;
      //REAL E = computeFeedbackEnergy(param, 1, level, msn/param->unit.unit_mass);

      REAL e=1. ;
      REAL E=e /POW( 2.,-3.*level);

    //  if ( (x_src==0) && (y_src==0) && (z_src==0)) E/=8.;

    //  printf("msn=%e\n",  param->unitary_stars_test->mass );

      printf("cell egy t0=%e\n",  cell->field.E);

      //thermalFeedbackCell(cell, E);
      //thermalFeedbackOct(cell, E);
      kineticFeedback_test(param, cell,NULL,aexp,level, E);

      printf("cell egy t1=%e\n",  cell->field.E);

      printf("SN active at t = %e\n", t);
      printf("SN pos =>  x=%e y=%e z=%e \n", x,y,z);

      printf("eblast= %e \n", E*POW( 2.,-3.*level));

      return 1;
    }
    return 0;
  }
  return 0;
}
#endif // SNTEST

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void supernovae(struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){
// ----------------------------------------------------------//
/// Call the feedback function for all cells of the grid
// ----------------------------------------------------------//
  if(param->sn->feedback_eff){
    if(cpu->rank==RANK_DISP) printf("SUPERNOVAE on level %d\n", level);

    int Nsn = 0;

    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *curoct=cpu->octList[level-1][iOct];

      int icell;
      for(icell=0;icell<8;icell++) {
        struct CELL *curcell = &curoct->cell[icell];
        Nsn += feedback(curcell, param, cpu, aexp, level, dt);
      }
    }


  #ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&Nsn,   1,MPI_INT,   MPI_SUM,cpu->comm);
  #endif

    if(cpu->rank==RANK_DISP && Nsn) {printf("%d\tActive SN\n",Nsn);}
  }
}


#endif // SUPERNOVAE
