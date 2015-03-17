/**
  * \file convert.c
  * \brief Contain unity convertion fonctions
  *
  */

#include "prototypes.h"

void setScale(struct RUNPARAMS *param, REAL aexp){
/**
  * Compute all the scaling factors between physical unit and code unit
  * according to Martel & Shapiro 1998
  *
  * Need to be called each time aexp change.
  */

#ifdef TESTCOSMO
  param->scale.l = 1./(aexp * param->unit.unit_l);
  param->scale.d = POW(aexp,3) /(param->unit.unit_d);
  param->scale.v = aexp / param->unit.unit_v;
  param->scale.t = 1./(POW(aexp,2) * param->unit.unit_t);
  param->scale.p = POW(aexp,5) / (param->unit.unit_d * POW(param->unit.unit_v,2));
  param->scale.E = POW(aexp,2) / POW(param->unit.unit_v,2);
  param->scale.mass = 1. /(param->unit.unit_d*POW(param->unit.unit_l,3));
  param->scale.n = 1.;
  param->scale.N = 1.;
#else
  // TODO consider the non cosmological case.
#endif
}

// ----------------------------------------------------------//

REAL l2code(struct RUNPARAMS *param, REAL l){
  return l * param->scale.l;
}
REAL code2l(struct RUNPARAMS *param, REAL l){
  return l / param->scale.l;
}

// ----------------------------------------------------------//

REAL d2code(struct RUNPARAMS *param, REAL d){
  return d * param->scale.d;
}
REAL code2d(struct RUNPARAMS *param, REAL d){
  return d / param->scale.d;
}

// ----------------------------------------------------------//

REAL v2code(struct RUNPARAMS *param, REAL v){
  return v * param->scale.v;
}
REAL code2v(struct RUNPARAMS *param, REAL v){
  return v / param->scale.v;
}

// ----------------------------------------------------------//

REAL t2code(struct RUNPARAMS *param, REAL t){
  return t * param->scale.t;
}
REAL code2t(struct RUNPARAMS *param, REAL t){
  return t / param->scale.t;
}

// ----------------------------------------------------------//

REAL p2code(struct RUNPARAMS *param, REAL p){
  return p * param->scale.p;
}
REAL code2p(struct RUNPARAMS *param, REAL p){
  return p / param->scale.p;
}

// ----------------------------------------------------------//

REAL E2code(struct RUNPARAMS *param, REAL E){
  return E * param->scale.E;
}
REAL code2E(struct RUNPARAMS *param, REAL E){
  return E / param->scale.E;
}

// ----------------------------------------------------------//

REAL mass2code(struct RUNPARAMS *param, REAL mass){
  return mass * param->scale.mass;
}
REAL code2mass(struct RUNPARAMS *param, REAL mass){
  return mass / param->scale.mass;
}

