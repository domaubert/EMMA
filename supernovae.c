#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prototypes.h"
#include "oct.h"
#include "hydro_utils.h"

#ifdef WMPI
#include <mpi.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedbackCell(struct CELL *cell,  REAL E){
// ----------------------------------------------------------//
// Inject an energy "E" in the cell "cell" on thermal form.
// ----------------------------------------------------------//

#ifdef WRAD
        cell->field.E += E;
        cell->field.p += E*(GAMMA-1.);
#endif
}

void thermalFeedbackOct(struct CELL *cell,  REAL E){
// ----------------------------------------------------------//
// Inject an energy "E" in all cells of the oct contening the cell "cell" on thermal form.
// ----------------------------------------------------------//

#ifdef WRAD
    struct OCT* oct = cell2oct(cell);
    int i;
    for(i=0;i<8;i++){
        struct CELL* curcell = &oct->cell[i];

        REAL e = E/8.;
        curcell->field.E += E;
        curcell->field.p += E*(GAMMA-1.);
    }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kineticFeedback(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, REAL aexp, int level, REAL E){
// ----------------------------------------------------------//
// Inject an energy "E" in all cells of the oct contening the cell "cell" on kinetic form, radially to the center of the oct and uniformly in all cells
// ----------------------------------------------------------//

  REAL mtot_feedback = curp->mass*0.5260172663907063;  // http://www.stsci.edu/science/starburst99/figs/mass_inst_e.html
  curp->mass -= mtot_feedback;

  int i;
  for(i=0;i<8;i++){
    struct OCT* oct = cell2oct(cell);
    struct CELL* curcell = &oct->cell[i];

    REAL e = E/8.; //uniform energy distribution
    REAL me = mtot_feedback/8.; //uniform ejecta distribution

    REAL dv = POW(2.,-3*level);

    REAL rho_i = curcell->field.d; //initial density
    REAL rho_e = me/dv; // density ejecta

    REAL vxi = curcell->field.u; // initial velocity
    REAL vyi = curcell->field.v;
    REAL vzi = curcell->field.w;

    REAL ve = SQRT(2.*e/rho_e);//velocity ejecta

    REAL dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};
    REAL dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
    REAL dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

    REAL vxe = curp->vx + ve*dir_x[i]/2.; // projection on axis in the particle framework
    REAL vye = curp->vy + ve*dir_y[i]/2.; // cos45*cos45 = 1/2
    REAL vze = curp->vz + ve*dir_z[i]/2.;

    curcell->field.d += rho_e; //new density

    curcell->field.u = (vxi*rho_i + vxe*rho_e)/(rho_i+rho_e); //new velocity
    curcell->field.v = (vyi*rho_i + vye*rho_e)/(rho_i+rho_e);
    curcell->field.w = (vzi*rho_i + vze*rho_e)/(rho_i+rho_e);

    //Energy conservation

    struct Utype U; // conservative field structure
    W2U(&curcell->field, &U); // primitive to conservative
    U.eint*=1.+me/curcell->field.d; // compute new internal energy
    U2W(&U, &curcell->field); // back to primitive

    getE(&curcell->field); //compute new total energy
    curcell->field.a=SQRT(GAMMA*curcell->field.p/curcell->field.d); // compute new sound speed
  //  curcell->field.p=curcell->field.E*(GAMMA-1.); // compute new pressure
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL computeFeedbackEnergy(struct RUNPARAMS *param, REAL aexp, int level, REAL mstar){
// ----------------------------------------------------------//
// Compute the total Enegy of feedback
// ----------------------------------------------------------//

  REAL egy = 1.936421963946603e+55; // erg/1e6M0 from http://www.stsci.edu/science/starburst99/figs/energy_inst_e.html
  egy *= 1e-7/(1e6*1.989e30); // j/kg

  // egy *= a5/(d*v2)*mass ;  // j/kg in code unit
  // mass/d = l3
  // egy *= a5/v2*l3;
  // l3/v2 = l3/l2*t2 = l*t2
  // egy *= a5*l*t2

  egy *= POW(aexp,5)*param->unit.unit_l*POW(param->unit.unit_l,2);  // j/kg in code unit
  egy *= mstar; // j in code unit
  egy /= POW(2.,-3*level); // j/m3 in code unit

  return egy;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int SN_TMP_PARAM = 1;
#ifndef SNTEST
#ifdef PIC
int feedback(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt){
// ----------------------------------------------------------//
// Parcour all particles and look for supernovae
// If found compute energy and inject it
// ----------------------------------------------------------//

  int Nsn = 0;
  REAL t0;

  struct PART *nexp=cell->phead;
  struct PART  *curp;

  if(nexp==NULL) return 0;
  do{ curp=nexp;
    nexp=curp->next;

    if (curp->isStar==2 || curp->isStar==3){

      REAL E = computeFeedbackEnergy(param, aexp, level, curp->mass);
      fKIN = param->sn->feedback_frac;
      //printf("Energy injected =%e mass=%e aexp=%e\n",E,curp->mass,aexp);
      thermalFeedbackOct(cell, E*(1.-fKIN));
      kineticFeedback(param, cell,curp,aexp,level, E*fKIN);

      Nsn++;
    }
  }while(nexp!=NULL);
  return Nsn;
}
#endif // PIC
#else // ifdef SNTEST
int feedback(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt){
// ----------------------------------------------------------//
// For sedov test
// A supernovae explode in a uniform medium
// ----------------------------------------------------------//


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

void supernovae(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level, int is){
// ----------------------------------------------------------//
// Call the feedback function for all cells of the grid
// ----------------------------------------------------------//


  if(param->sn->feedback_eff){
    if(cpu->rank==RANK_DISP) printf("SUPERNOVAE\n");

    int Nsn = 0;
    struct OCT  *nextoct=firstoct[level-1];

    do {
      if(nextoct==NULL) continue;
      struct OCT* curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu != cpu->rank) continue;

      int icell;
      for(icell=0;icell<8;icell++) {
        struct CELL *curcell = &curoct->cell[icell];

        Nsn += feedback(curcell, param, cpu, aexp, level, dt);
      }
    }while(nextoct!=NULL);

  #ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&Nsn,   1,MPI_INT,   MPI_SUM,cpu->comm);
  #endif

    if(cpu->rank==RANK_DISP && Nsn) {printf("%d\tActive SN\n",Nsn);}
  }
}
