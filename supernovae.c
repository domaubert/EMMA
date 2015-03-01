#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef WMPI
#include <mpi.h>
#endif

#include "prototypes.h"
#include "oct.h" //cell2oct
#include "hydro_utils.h" // W2U and U2W

/// ----------------------------------------------------------//
/// ----------------------------------------------------------//

/// This file contain functions for supernovae and stellar wind
/// More documentation can be found on the wiki :
/// https://github.com/domaubert/EMMA/wiki/Supernovae

/// ----------------------------------------------------------//
/// ----------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedbackCell(struct CELL *cell,  REAL E){
/// ----------------------------------------------------------//
/// Inject an energy "E" in the cell "cell" on thermal form.
/// ----------------------------------------------------------//

#ifdef WRAD
        cell->field.E += E;
        cell->field.p += E*(GAMMA-1.);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedbackOct(struct CELL *cell,  REAL E){
/// ----------------------------------------------------------//
/// Inject an energy "E" in all cells of the oct contening
/// the cell "cell" uniformly on thermal form.
/// ----------------------------------------------------------//

#ifdef WRAD
    struct OCT* oct = cell2oct(cell);
    int i;
    for(i=0;i<8;i++){
        struct CELL* curcell = &oct->cell[i];
        REAL e = E/8.;
        curcell->field.E += e;
        curcell->field.p += e*(GAMMA-1.);
    }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kineticFeedback(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, REAL aexp, int level, REAL E){
/// ----------------------------------------------------------//
/// Inject an energy "E" in all cells of the oct contening
/// the cell "cell" on kinetic form, radially to the center
/// of the oct and uniformly in all cells
/// The proportion of ejecta follow fig 109 of Starburst99 model:
/// http://www.stsci.edu/science/starburst99/figs/mass_inst_e.html
/// ----------------------------------------------------------//

  REAL ejecta_proportion = 0.5260172663907063;
  REAL mtot_feedback = curp->mass* ejecta_proportion;

  int i;
  for(i=0;i<8;i++){
    struct OCT* oct = cell2oct(cell);
    struct CELL* curcell = &oct->cell[i];

    REAL e = E/8.; //uniform energy distribution
    REAL me = mtot_feedback/8.; //uniform ejecta distribution

    REAL dv = POW(2.,-3*level);// curent volume

    REAL rho_i = curcell->field.d; //initial density
    REAL rho_e = me/dv; // density ejecta

    REAL vxi = curcell->field.u; // initial velocity
    REAL vyi = curcell->field.v;
    REAL vzi = curcell->field.w;

    REAL ve = SQRT(2.*e/rho_e);//velocity ejecta

    REAL dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};// diagonal projection
    REAL dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
    REAL dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

    REAL vxe = curp->vx + ve*dir_x[i]/2.; // projection on axis in the particle framework
    REAL vye = curp->vy + ve*dir_y[i]/2.; // cos45*cos45 = 1/2
    REAL vze = curp->vz + ve*dir_z[i]/2.;

    curcell->field.d += rho_e; //new density
    curp->mass -= me; // new particle mass

    curcell->field.u = (vxi*rho_i + vxe*rho_e)/(rho_i+rho_e); //new velocity
    curcell->field.v = (vyi*rho_i + vye*rho_e)/(rho_i+rho_e);
    curcell->field.w = (vzi*rho_i + vze*rho_e)/(rho_i+rho_e);

    //Energy conservation

    struct Utype U; // conservative field structure
    W2U(&curcell->field, &U); // primitive to conservative
    U.eint*=1.+rho_e/curcell->field.d; // compute new internal energy
    U2W(&U, &curcell->field); // back to primitive

    getE(&curcell->field); //compute new total energy
    curcell->field.a=SQRT(GAMMA*curcell->field.p/curcell->field.d); // compute new sound speed
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL computeFeedbackEnergy(struct RUNPARAMS *param, REAL aexp, int level, REAL mstar){
/// ----------------------------------------------------------//
/// Compute the total feedback energy following fig 115 of the
/// Starburst99 model :
/// http://www.stsci.edu/science/starburst99/figs/energy_inst_e.html
/// ----------------------------------------------------------//

  //REAL egy = 1.936421963946603e+55 *1e-7/(1e6*SOLAR_MASS); // erg/1e6M0 -> j/kg

/*
  REAL dx = POW(0.5,level) * aexp * param->unit.unit_l;
  REAL dv = POW(dx,3);

  REAL egy = 9.73565592733335e11; // j/kg
  egy *= mstar * param->unit.unit_mass; // j
  egy /= dv; // j/m3
  egy *= POW(aexp,5.)/(param->unit.unit_d *POW(param->unit.unit_v,2.)) ; // j/m3 in code unit
*/



  //REAL dv = POW(0.5,3*level) *  POW(aexp,3) *  POW(param->unit.unit_l,3);

  REAL egy = 9.73565592733335e11; // j/kg
  egy *= mstar* POW(aexp,2.) / ( POW(0.5,3*level) *  POW(param->unit.unit_v,2.) ) ; // j/m3 in code unit
  return egy;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SNTEST
#ifdef PIC
int feedback(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt){
/// ----------------------------------------------------------//
/// Scan all particles of cell "cell" and look for supernovae
/// If SN found, compute energy and inject it
/// ----------------------------------------------------------//

  int Nsn = 0;

  struct PART *nexp=cell->phead;

  if(nexp==NULL) return 0;
  do{
    struct PART *curp=nexp;
    nexp=curp->next;

    if (curp->isStar==2 || curp->isStar==3){ // if curp is in SN state

      REAL E = computeFeedbackEnergy(param, aexp, level, curp->mass);

      printf("Energy injected =%e mass=%e aexp=%e\n",E,curp->mass,aexp);
      //abort();

      thermalFeedbackOct(cell, E*(1.-param->sn->feedback_frac));
      kineticFeedback(param, cell,curp,aexp,level, E*param->sn->feedback_frac);

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
/// ----------------------------------------------------------//
/// For sedov test
/// A supernovae explode in a uniform medium
/// ----------------------------------------------------------//

	struct OCT* oct = cell2oct(cell);

  if (!SN_TMP_PARAM && cpu->rank==RANK_DISP){
    printf("======================================\n");
    printf("===SN EXPLODE=========================\n");
    printf("======================================\n");
  }

	if (oct->x == 0.5 && oct->y == 0.5 && oct->z == 0.5 && cell->idx == 0){

		REAL in_yrs = param->unit.unit_t/MYR *1e6;
		REAL t = aexp * in_yrs;

		if ( t >= LIFETIME_OF_STARS_IN_TEST && SN_TMP_PARAM ){
			SN_TMP_PARAM = 0;

// src http://cdsads.u-strasbg.fr/abs/2009A%26A...495..389B
// L(M) = 9.315613314066386e+16 photon/s/kg

//    REAL msn = 7.6e6 * 2e30 ; //kg L=1.42e61
      REAL msn = 26.84 * 2e30 ; //kg L=5e48

      REAL E = computeFeedbackEnergy(param, 0, 1, level, msn/param->unit.unit_mass );

      REAL  fKIN = compute_fkin(param,cell,E,level,1.);
      thermalFeedback(cell, E*(1.-fKIN));
      kineticFeedback(cell, E*(   fKIN));

      printf("SN active at t = %e\n", t);
      printf("eblast= %e \n", E*POW( 2.,-3.*level));
    }
    return 1;
  }
}
#endif // SNTEST

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void supernovae(struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){
/// ----------------------------------------------------------//
/// Call the feedback function for all cells of the grid
/// ----------------------------------------------------------//

  if(param->sn->feedback_eff){
    if(cpu->rank==RANK_DISP) printf("SUPERNOVAE\n");

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
