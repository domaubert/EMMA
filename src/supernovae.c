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
#include <float.h>

#ifdef WMPI
#include <mpi.h>
#endif

#include "prototypes.h"
#include "oct.h" //cell2oct
#include "hydro_utils.h" // W2U and U2W
#include "convert.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REAL interp_egy(struct RUNPARAMS *param,REAL t){
  int i;
  for(i=0;i<10000;i++){
    if (param->sn->egy_loss_t[i]>t){
      break;
    }
  }
  REAL dx = param->sn->egy_loss_t[i]  -param->sn->egy_loss_t[i-1];
  REAL dy = param->sn->egy_loss_egy[i]-param->sn->egy_loss_egy[i-1];
  return (t-param->sn->egy_loss_t[i-1]) * dy/dx + param->sn->egy_loss_egy[i-1];
}

REAL get_egy(struct RUNPARAMS *param, struct PART *curp, REAL aexp, REAL dt,int level){
  REAL cur_t = param->cosmo->tphy;
  REAL age =  cur_t - curp->age;
  REAL dt_in_yr = dt * aexp*aexp * param->unit.unit_t /(365*24*3600);

  REAL cur_egy  = interp_egy(param,age);
  REAL next_egy = interp_egy(param,age+dt_in_yr);

  REAL de = next_egy-cur_egy; // j/kg
  de*= curp->mass/POW(0.5,3*level) * POW(aexp,2.)/POW(param->unit.unit_v,2.); // j/m3 in code unit

  return de<0?0:de;
}

REAL interp_mass(struct RUNPARAMS *param,REAL t){
  int i;
  for(i=0;i<10000;i++){
    if (param->sn->mass_loss_t[i]>t){
      break;
    }
  }

  REAL dx = param->sn->mass_loss_t[i]    -param->sn->mass_loss_t[i-1];
  REAL dy = param->sn->mass_loss_mass[i] -param->sn->mass_loss_mass[i-1];
  return (t-param->sn->mass_loss_t[i-1]) * dy/dx + param->sn->mass_loss_mass[i-1];
}

REAL get_mass(struct RUNPARAMS *param, struct PART *curp, REAL aexp, REAL dt){
  REAL age =  param->cosmo->tphy - curp->age;
  REAL dt_in_yr = dt * aexp*aexp * param->unit.unit_t /(365*24*3600);

  REAL cur_mass  = interp_mass(param,age);
  REAL next_mass = interp_mass(param,age+dt_in_yr);

  REAL dm = curp->mass*(next_mass-cur_mass);

  return dm<0?0:dm;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedbackCell(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, int level, REAL E){
// ----------------------------------------------------------
/// Inject an energy "E" in the cell "cell" on thermal form.
//----------------------------------------------------------
#ifdef SNTEST
#endif // SNTEST

  printf("injecting Energy in thermal form within a cell\n");

  //REAL mtot_feedback = curp->mass* param->sn->ejecta_proportion;
  //REAL dv = POW(0.5,3.*level);
  //REAL d_e = mtot_feedback/dv ;

  //curp->mass -= mtot_feedback;
  //cell->field.d += d_e;


  cell->field.E += E;
  cell->field.p += E*(GAMMA-1.);

  //Energy conservation
#ifdef DUAL_E
//    struct Utype U; // conservative field structure
//    W2U(&cell->field, &U); // primitive to conservative
//    U.eint*=1.+d_e/cell->field.d; // compute new internal energy
//    U2W(&U, &cell->field); // back to primitive
#else
//    cell->field.p*=1.+d_e/cell->field.d; // compute new internal energy
#endif

    getE(&cell->field); //compute new total energy
    cell->field.p=FMAX(cell->field.p,PMIN);
    cell->field.a=SQRT(GAMMA*cell->field.p/cell->field.d); // compute new sound speed
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalFeedbackOct(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, int level, REAL E){
// ----------------------------------------------------------
/// Inject an energy "E" in all cells of the oct contening
/// the cell "cell" uniformly on thermal form.
// ----------------------------------------------------------
#ifdef SNTEST
  printf("injecting Energy in thermal form within a oct\n");
#endif // SNTEST

  //REAL mtot_feedback = curp->mass* param->sn->ejecta_proportion;
  //REAL dv = POW(0.5,3.*level);
  //REAL d_e = mtot_feedback/dv ;

  //curp->mass -= mtot_feedback;

  struct OCT* oct = cell2oct(cell);
  int i;
  for(i=0;i<8;i++){
      struct CELL* curcell = &oct->cell[i];

      //cell->field.d += d_e/8.;

      REAL e = E/8.;
      curcell->field.E += e;
      curcell->field.p += e*(GAMMA-1.);

      //Energy conservation
#ifdef DUAL_E
//    struct Utype U; // conservative field structure
//    W2U(&cell->field, &U); // primitive to conservative
//    U.eint*=1.+d_e/cell->field.d; // compute new internal energy
//    U2W(&U, &cell->field); // back to primitive
#else
//    cell->field.p*=1.+d_e/cell->field.d; // compute new internal energy
#endif

    getE(&cell->field); //compute new total energy
    cell->field.p=FMAX(cell->field.p,PMIN);
    cell->field.a=SQRT(GAMMA*cell->field.p/cell->field.d); // compute new sound speed

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
/* #ifdef DUAL_E */
/*     struct Utype U; // conservative field structure */
/*     W2U(&curcell->field, &U); // primitive to conservative */
/*     U.eint*=1.+rho_e/curcell->field.d; // compute new internal energy */
/*     U2W(&U, &curcell->field); // back to primitive */
/* #else */
/*     curcell->field.p*=1.+rho_e/curcell->field.d; // compute new internal energy */
/* #endif */

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

void kineticFeedback_simple(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, REAL aexp, int level, REAL E){

  REAL mtot_feedback = curp->mass* param->sn->ejecta_proportion;
  curp->mass -= mtot_feedback;

  struct OCT* oct = cell2oct(cell);
  int i;
  for(i=0;i<8;i++){
    struct CELL* curcell = &oct->cell[i];
    REAL dv = POW(0.5,3*level);// curent volume

    REAL dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};// diagonal projection
    REAL dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
    REAL dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

    REAL E_e = E/8.;
    REAL m_e = mtot_feedback/8.;
    REAL d_e = m_e/dv; // ejecta density

    REAL v_e = SQRT(2.*E_e/curcell->field.d);

    curcell->field.d += m_e/dv;
    curcell->field.u += v_e*dir_x[i]/SQRT(3.);
    curcell->field.v += v_e*dir_y[i]/SQRT(3.);
    curcell->field.w += v_e*dir_z[i]/SQRT(3.);

    //Energy conservation
/* #ifdef DUAL_E */
/*     struct Utype U; // conservative field structure */
/*     W2U(&cell->field, &U); // primitive to conservative */
/*     U.eint*=1.+d_e/cell->field.d; // compute new internal energy */
/*     U2W(&U, &cell->field); // back to primitive */
/* #else */
/*     cell->field.p*=1.+d_e/cell->field.d; // compute new internal energy */
/* #endif */

    getE(&curcell->field); //compute new total energy
    curcell->field.p=FMAX(curcell->field.p,PMIN);
    curcell->field.a=SQRT(GAMMA*curcell->field.p/curcell->field.d); // compute new sound speed
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kineticFeedback_impulsion(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, REAL aexp, int level, REAL E){
// ----------------------------------------------------------//
/// Inject an energy "E" in all cells of the oct contening
/// the cell "cell" on kinetic form, radially to the center
/// of the oct and uniformly in all cells
// ----------------------------------------------------------//

  REAL mtot_feedback = curp->mass* param->sn->ejecta_proportion;
  curp->mass -= mtot_feedback;
  struct OCT* oct = cell2oct(cell);

  int i;
//#define LOAD_FACTOR
#ifdef LOAD_FACTOR

  REAL load_factor=10;

  REAL local_rho_min=FLT_MAX;

  for(i=0;i<8;i++){
    struct CELL* curcell = &oct->cell[i];
    local_rho_min = FMIN(curcell->field.d,local_rho_min);
  }

  REAL lf = FMIN(load_factor*curp->mass,local_rho_min);

  mtot_feedback += lf;
  for(i=0;i<8;i++){
    struct CELL* curcell = &oct->cell[i];
    curcell->field.d -= lf/8.;
  }
#endif // LOAD_FACTOR

  for(i=0;i<8;i++){
    struct CELL* curcell = &oct->cell[i];
    REAL dv = POW(0.5,3*level);// curent volume

    REAL dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};// diagonal projection
    REAL dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
    REAL dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

    REAL E_e = E/8.; // uniform distribution
    REAL m_e = mtot_feedback/8.;

    REAL d_e = m_e/dv; // ejecta density
    REAL v_e = SQRT(2.*E_e/d_e);// ejecta velocity / particle
    REAL vxe = curp->vx + v_e*dir_x[i]/SQRT(3.); // ejecta velocity /grid
    REAL vye = curp->vy + v_e*dir_y[i]/SQRT(3.);
    REAL vze = curp->vz + v_e*dir_z[i]/SQRT(3.);

    REAL d_i = curcell->field.d; // initial density
    REAL vxi = curcell->field.u; // initial velocity
    REAL vyi = curcell->field.v;
    REAL vzi = curcell->field.w;


#ifdef WRADHYD
    REAL xion=curcell->field.dX/d_i;
#endif

    curcell->field.d += d_e; //new density
    curcell->field.u = (vxi*d_i + vxe*d_e)/(d_i+d_e); //new velocity
    curcell->field.v = (vyi*d_i + vye*d_e)/(d_i+d_e);
    curcell->field.w = (vzi*d_i + vze*d_e)/(d_i+d_e);

#ifdef WRADHYD
     curcell->field.dX=curcell->field.d*xion;
#endif

    //Energy conservation
/* #ifdef DUAL_E */
/*     struct Utype U; // conservative field structure */
/*     W2U(&curcell->field, &U); // primitive to conservative */
/*     U.eint*=1.+d_e/curcell->field.d; // compute new internal energy */
/*     U2W(&U, &curcell->field); // back to primitive */
/* #else */
/*     curcell->field.p*=1.+d_e/curcell->field.d; // compute new internal energy */
/* #endif */

    getE(&curcell->field); //compute new total energy
    curcell->field.p=FMAX(curcell->field.p,PMIN);
    curcell->field.a=SQRT(GAMMA*curcell->field.p/curcell->field.d); // compute new sound speed
  }
}

int kineticFeedback_mixt(struct RUNPARAMS *param, struct CELL *cell,struct PART *curp, REAL aexp, int level, REAL E, REAL dt){
// ----------------------------------------------------------//
/// Inject an energy "E" in all cells of the oct contening
/// the cell "cell" on kinetic form, radially to the center
/// of the oct and uniformly in all cells
// ----------------------------------------------------------//

  //REAL mtot_feedback = curp->mass* param->sn->ejecta_proportion;
  //printf("injecting_old m=%e, e=%e\t",curp->mass* param->sn->ejecta_proportion,E);

  REAL mtot_feedback = get_mass(param,curp,aexp,dt);

  if (mtot_feedback==0){
    printf("WARNING null mass in kinetic feedback!\n");
    return 1;
  }

  printf("injecting_new m=%e, e=%e\n",mtot_feedback ,E);

  int i;
{
//#define LOAD_FACTOR
#ifdef LOAD_FACTOR

  REAL load_factor=10;
  REAL local_rho_min=FLT_MAX;

  for(i=0;i<8;i++){
    struct OCT* oct = cell2oct(cell);
    struct CELL* curcell = &oct->cell[i];
    local_rho_min = FMIN(curcell->field.d,local_rho_min);
  }

  REAL lf = FMIN(load_factor*curp->mass,local_rho_min);
  mtot_feedback += lf;

  for(i=0;i<8;i++){
    struct OCT* oct = cell2oct(cell);
    struct CELL* curcell = &oct->cell[i];
    curcell->field.d -= lf/8.;
  }
#endif // LOAD_FACTOR
}

  for(i=0;i<8;i++){
    struct OCT* oct = cell2oct(cell);
    struct CELL* curcell = &oct->cell[i];
    REAL dv = POW(0.5,3*level);// curent volume

    REAL dir_x[]={-1., 1.,-1., 1.,-1., 1.,-1., 1.};// diagonal projection
    REAL dir_y[]={-1.,-1., 1., 1.,-1.,-1., 1., 1.};
    REAL dir_z[]={-1.,-1.,-1.,-1., 1., 1., 1., 1.};

    REAL fact = 1.0;

    REAL m_e = mtot_feedback/8.;
    REAL d_e = m_e/dv; // ejecta density
//--------------------------------------------------------------------------//

    // kinetic energy injection
    REAL E_e = E/8. *fact;
    REAL v_e = SQRT(2.*E_e/curcell->field.d);

    //printf("field.d=%e\n",curcell->field.d);
    REAL sqrt3 = SQRT(3.);
    curcell->field.u += v_e*dir_x[i]/sqrt3;
    curcell->field.v += v_e*dir_y[i]/sqrt3;
    curcell->field.w += v_e*dir_z[i]/sqrt3;
//--------------------------------------------------------------------------//
/*
    // momentum injection
    E_e = E/8. *(1.-fact) ; // uniform distribution
    v_e = SQRT(2.*E_e/d_e);// ejecta velocity / particle

    REAL vxe = curp->vx + v_e*dir_x[i]/sqrt3; // ejecta velocity /grid
    REAL vye = curp->vy + v_e*dir_y[i]/sqrt3;
    REAL vze = curp->vz + v_e*dir_z[i]/sqrt3;

    REAL d_i = curcell->field.d; // initial density
    REAL vxi = curcell->field.u; // initial velocity
    REAL vyi = curcell->field.v;
    REAL vzi = curcell->field.w;

    curcell->field.u = (vxi*d_i + vxe*d_e)/(d_i+d_e); //new velocity
    curcell->field.v = (vyi*d_i + vye*d_e)/(d_i+d_e);
    curcell->field.w = (vzi*d_i + vze*d_e)/(d_i+d_e);
*/
//--------------------------------------------------------------------------//

    // mass conservation
    curp->mass       -= m_e;
    curcell->field.d += d_e; //new density

    //Energy conservation
/* #ifdef DUAL_E */
/*     struct Utype U; // conservative field structure */
/*     W2U(&curcell->field, &U); // primitive to conservative */
/*     U.eint*=1.+d_e/curcell->field.d; // compute new internal energy */
/*     U2W(&U, &curcell->field); // back to primitive */
/* #else */
/*     curcell->field.p*=1.+d_e/curcell->field.d; // compute new internal energy */
/* #endif */

    getE(&curcell->field); //compute new total energy
    curcell->field.p=FMAX(curcell->field.p,PMIN);
    curcell->field.a=SQRT(GAMMA*curcell->field.p/curcell->field.d); // compute new sound speed
  }
  return 0;
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


  static REAL total_sn_egy=0;
  int Nsn = 0;

  struct PART *nexp=cell->phead;

  if(nexp==NULL) return 0;
  do{
    struct PART *curp=nexp;
    nexp=curp->next;

//#define CONTINUOUS_SN
#ifdef CONTINUOUS_SN
    REAL age = param->cosmo->tphy - curp->age;
    if (curp->isStar && age>0 && age<5e7){ // if curp is in SN state
      REAL E= get_egy(param,curp,aexp,dt,level);
      
#else
      if (curp->isStar==5 || curp->isStar==7){ // if curp is in SN state
	REAL E = computeFeedbackEnergy(param, aexp, level, curp->mass);
	printf("mass =%e %d\n",curp->mass,curp->isStar);
#endif // CONTINUOUS_SN

      total_sn_egy+=E;

      //printf("Energy injected =%e mass=%e aexp=%e cellE=%e on RANK %d\n",E,curp->mass,aexp, cell->field.E,cpu->rank);

      if(param->sn->feedback_frac){
         /* oct feedback
          * if kinetic feedback
          * the energy is injected in a oct
          */

        thermalFeedbackOct(param,cell,curp,level, E*(1.-param->sn->feedback_frac));
        //kineticFeedback_simple(param, cell,curp,aexp,level, E*param->sn->feedback_frac);
        //kineticFeedback_impulsion(param, cell,curp,aexp,level, E*param->sn->feedback_frac);

        if (param->sn->ejecta_proportion){
          kineticFeedback_impulsion(param, cell,curp,aexp,level, E*param->sn->feedback_frac);
          //kineticFeedback_mixt(param, cell,curp,aexp,level, E*param->sn->feedback_frac, dt);
        }else{
          kineticFeedback_simple(param, cell,curp,aexp,level, E*param->sn->feedback_frac);
        }

      }else{
         /* cell feedback
          * if thermal feedback only
          * the energy can be injected in just a cell
          */

        thermalFeedbackCell(param,cell,curp,level, E*(1.-param->sn->feedback_frac));

      }

      Nsn++;

      if(curp->isStar==5) curp->isStar=4;
      if(curp->isStar==7) curp->isStar=8;

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
      thermalFeedbackOct(cell, E *0.3);
      kineticFeedback_test(param, cell,NULL,aexp,level, E *0.7);

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

    //    printf("%d\tActive SN on rank %d\n",Nsn,cpu->rank);

#ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&Nsn,1,MPI_INT,MPI_SUM,cpu->comm);
    MPI_Allreduce(MPI_IN_PLACE,&param->sn->trig_sn,1,MPI_INT,MPI_SUM,cpu->comm);
#endif

    if(cpu->rank==RANK_DISP){
      if (Nsn){
        if(!(param->sn->trig_sn)){
          printf("FIRST_SN at z=%e\n",1./aexp-1.);
          param->sn->trig_sn=1;
        }
        printf("%d\tActive SN\n",Nsn);
      }
    }

  }
}

#endif // SUPERNOVAE
