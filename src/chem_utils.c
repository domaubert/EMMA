/**
  * \file chem_utils.c
  */

#ifdef WRAD

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>

#include "prototypes.h"
#include "oct.h"

#define idloc 0 // KEPT FROM CUDATON FOR SIMPLICITY
#define FSCHAYE 1.0


//================================================================================
void E2T(struct Rtype *R, REAL aexp,struct RUNPARAMS *param){

  REAL tloc;
  REAL eint=R->eint;
  REAL nH=R->nh;
  REAL x=R->nhplus/R->nh;
  REAL pstar=param->unit.unit_n*param->unit.unit_d*POW(param->unit.unit_v,2);
  nH=nH/POW(aexp,3)*param->unit.unit_N;
  eint=eint/POW(aexp,5)*pstar;
  tloc=eint/(1.5*nH*KBOLTZ*(1.+x)*(1.+yHE));
  R->temp=tloc;
}

#ifdef WCHEM

// ============================================================================
REAL cucompute_alpha_b(REAL temp, REAL unit_number, REAL aexp)
{
  // CASE B recombination rate m**3 s*-1
  // temperature should be given in Kelvin

  REAL alpha_b,lambda;
  lambda=2e0*157807e0/temp;
  alpha_b=2.753e-14*POW(lambda,1.5)/POW(1e0+POW(lambda/2.740,0.407),2.242); //cm3/s
#ifdef TESTCOSMO
  alpha_b=alpha_b*1e-6*unit_number;//(aexp*aexp*aexp); //m3/s
#else
  alpha_b=alpha_b*1e-6*unit_number; //m3/s
#endif
  return alpha_b;
}

//=========================================================
//=========================================================

REAL cucompute_alpha_a(REAL temp, REAL unit_number, REAL aexp)
{
  // CASE A recombination rate m**3 s*-1
  // temperature should be given in Kelvin

  REAL alpha_a,lambda;
  lambda=2e0*157807e0/temp;
  alpha_a=1.269e-13*POW(lambda,1.503)/POW(1e0+POW(lambda/0.522,0.470),1.923); //cm3/s
#ifdef TESTCOSMO
  alpha_a=alpha_a*1e-6*unit_number;//(aexp*aexp*aexp); //m3/s
#else
  alpha_a=alpha_a*1e-6*unit_number; //m3/s
#endif
  return alpha_a;
}

//=========================================================
//=========================================================

REAL cucompute_beta(REAL temp, REAL unit_number, REAL aexp)
{
  // Collizional ionization rate m**3 s*-1
  // temperature in Kelvin
  REAL beta,T5;
  T5=temp/1e5;
  beta=1.17e-10*SQRT(temp)/(1+SQRT(T5))*exp(-(157809e0/temp)); //cm3/s
#ifdef TESTCOSMO
  beta=beta*1e-6*unit_number;//(aexp*aexp*aexp); // !m3/s
#else
  beta=beta*1e-6*unit_number; // !m3/s
#endif
  return beta;
}

//**********************************************************************************
//**********************************************************************************
void cuCompCooling(REAL temp, REAL x, REAL nH, REAL *lambda, REAL *tcool, REAL aexp,REAL CLUMPF,int src)
{

  REAL c1,c2,c3,c4,c5,c6;
  REAL unsurtc;
  REAL nh2;


  nh2=nH*1e-6;// ! m-3 ==> cm-3


  // Collisional Ionization Cooling

  c1=EXP(-157809.1e0/temp)*2.54e-21*SQRT(temp)/(1.+SQRT(temp/1e5))*x*(1.-x)*nh2*nh2*CLUMPF*(1+yHE);

  // Case A Recombination Cooling

  c2=0.;
#ifndef OTSA
  c2=1.778e-29*temp*POW(2e0*157807e0/temp,1.965e0)/POW(1.+POW(2e0*157807e0/temp/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2*CLUMPF*(1+yHE);
#ifdef NORECSOURCE
  if(src) c2*=NORECSOURCE_FACT;
#endif
#endif

  // Case B Recombination Cooling
  c6=0.;
#ifdef OTSA
  c6=3.435e-30*temp*POW(2e0*157807e0/temp,1.970e0)/POW(1.+(POW(2e0*157807e0/temp/2.250e0,0.376e0)),3.720e0)*x*x*nh2*nh2*CLUMPF;
#ifdef NORECSOURCE
  if(src) c6*=NORECSOURCE_FACT;
#endif
#endif



  // Collisional excitation cooling

  c3=EXP(-118348e0/temp)*7.5e-19/(1+SQRT(temp/1e5))*x*(1.-x)*nh2*nh2*CLUMPF*(1+yHE);


  // Bremmsstrahlung

  c4=1.42e-27*1.5e0*SQRT(temp)*x*(1.+yHE)*nh2*x*nh2*(1.+yHE)*CLUMPF;

  // Compton Cooling

  c5=0;
#ifndef WRADTEST
  c5=5.406e-24*(temp-2.727/aexp)/POW(aexp/0.001,4)*x*nh2;
  REAL Ta=2.727/aexp; c5=5.406e-36*(temp-Ta)/(aexp*aexp*aexp*aexp)*x*nh2*(1.+yHE);
#endif
  // Overall Cooling

  *lambda=c1+c2+c3+c4+c5+c6;// ! erg*cm-3*s-1


#ifdef HESIMPLE
  REAL c7,c8,c9,c10;
  //Ionisation Cooling
  c7=1.88e-21*SQRT(temp)*EXP(-285335.4/temp)/(1.+SQRT(temp/1e5))*nh2*x*(1.+yHE)*nh2*(1.-x)*yHE*CLUMPF;

  //Recombination Cooling
  c8=1.55e-26*POW(temp,0.3647)*nh2*x*(1.+yHE)*nh2*x*yHE*CLUMPF;

  // Dielectric Recombination Cooling
  c9=1.24e-13/POW(temp,1.5)*EXP(-470000./temp)*(1+0.3*EXP(-94000./temp))*nh2*x*(1.+yHE)*nh2*x*yHE*CLUMPF;

  // Coolisional Excitation cooling
  c10=5.54e-17/POW(temp,0.397)*EXP(-473638./temp)/(1.+SQRT(temp/1e5))*nh2*x*(1.+yHE)*nh2*x*yHE*CLUMPF;

  *lambda+=c7+c8+c9+c10;
#endif

  // Unit Conversion

  *lambda=(*lambda)*1e-7*1e6;// ! J*m-3*s-1

  // cooling times

  unsurtc=FMAX(c1,c2);
  unsurtc=FMAX(unsurtc,c3);
  unsurtc=FMAX(unsurtc,c4);
  unsurtc=FMAX(unsurtc,FABS(c5));

#ifdef HESIMPLE
  unsurtc=FMAX(unsurtc,c7);
  unsurtc=FMAX(unsurtc,c8);
  unsurtc=FMAX(unsurtc,c9);
  unsurtc=FMAX(unsurtc,c10);
#endif

  unsurtc=FMAX(unsurtc,c6)*1e-7;// ==> J/cm3/s

  *tcool=1.5e0*nh2*(1.+x)*(1+yHE)*KBOLTZ*temp/unsurtc; //Myr
}


//REAL _exp(REAL EARG){
//
//  assert(EARG!=0);
//
//  REAL EE;
//  if(EARG<1e-3){
//    REAL DL=-EARG+0.5*EARG*EARG-EARG*EARG*EARG/6.;
//    EE=1.+DL;
//  }else{
//    EE=EXP(-EARG);
//  }
//}

void chemrad(struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew, struct RUNPARAMS *param, REAL aexporg, int chemonly)
{
  int i,icell,igrp;
  //int idloc=0;
  int nitcool=0;

  REAL hnu0=13.6*1.6022e-19,
    Cool,
    tcool,
    dtcool,
    tcool1,
    currentcool_t=0.,
    beta,
    tloc,
    xt,
    eintt,
    ai_tmp1=0.,
    et[NGRP],
    p[NGRP];

    REAL alpha,alphab;

    REAL aexp;
    REAL ebkg[NGRP];
    REAL z=1./aexporg-1.;


    REAL hnu[NGRP];
    REAL alphae[NGRP];
    REAL alphai[NGRP];
    REAL factgrp[NGRP];

    for(igrp=0;igrp<NGRP;igrp++) {
      hnu[igrp]=param->atomic.hnu[igrp];
      alphae[igrp]=param->atomic.alphae[igrp];
      alphai[igrp]=param->atomic.alphai[igrp];
      factgrp[igrp]=param->atomic.factgrp[igrp];
    }

#define BLOCKCOOL 1 // KEPT FROM CUDATON FOR SIMPLICITY
#define idloc3 0 // KEPT FROM CUDATON FOR SIMPLICITY

#define FUDGEFACT 10

  REAL egyloc[BLOCKCOOL*NGRP];
  REAL  floc[3*BLOCKCOOL*NGRP];
  REAL  srcloc[BLOCKCOOL*NGRP];
  REAL x0[BLOCKCOOL];
  REAL nH[BLOCKCOOL];
  REAL eint[BLOCKCOOL];

  REAL fxt[NGRP],fyt[NGRP],fzt[NGRP];

  REAL dt=dtnew*param->unit.unit_t*POW(aexporg,2);

  REAL emin;
  struct Rtype R;
  REAL fudgecool=param->fudgecool;
  int ncvgcool=param->ncvgcool;
  REAL E0;
  REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS)*param->unit.unit_d*(1.-YHE)*(1.+yHE);
  REAL xorg;
  int compcool; // do we need to compute the cooling ?

  for(i=0;i<nread;i++){  // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells

      if(stencil[i].oct[6].cell[icell].split) continue; // we dont treat split cells

      memcpy(&R,&stencil[i].New.cell[icell].rfieldnew,sizeof(struct Rtype));// We get the local physical quantities after transport update
      for(igrp=0;igrp<NGRP;igrp++) ebkg[igrp]=0.;

      for (igrp=0;igrp<NGRP;igrp++)
      {
        egyloc[idloc+igrp*BLOCKCOOL]   =R.e[igrp]/(aexporg*aexporg*aexporg)*param->unit.unit_N;
        floc[0+idloc3+igrp*BLOCKCOOL*3]=R.fx[igrp]/POW(aexporg,4)*param->unit.unit_l/param->unit.unit_t*param->unit.unit_N;
        floc[1+idloc3+igrp*BLOCKCOOL*3]=R.fy[igrp]/POW(aexporg,4)*param->unit.unit_l/param->unit.unit_t*param->unit.unit_N;
        floc[2+idloc3+igrp*BLOCKCOOL*3]=R.fz[igrp]/POW(aexporg,4)*param->unit.unit_l/param->unit.unit_t*param->unit.unit_N;
      }

      x0[idloc]=R.nhplus/R.nh;
      xorg= x0[idloc];

      nH[idloc]=R.nh/(aexporg*aexporg*aexporg)*param->unit.unit_N;

      eint[idloc]=R.eint/POW(aexporg,5)*param->unit.unit_n*param->unit.unit_d*POW(param->unit.unit_v,2);
      emin=PMIN/(GAMMA-1.)/POW(aexporg,5)*param->unit.unit_n*param->unit.unit_d*POW(param->unit.unit_v,2); // physical minimal pressure

      for (igrp=0;igrp<NGRP;igrp++){
        srcloc[idloc+igrp*BLOCKCOOL]=R.src[igrp]*param->unit.unit_N/param->unit.unit_t/POW(aexporg,5); //phot/s/dv (physique)
        egyloc[idloc+igrp*BLOCKCOOL]+=srcloc[idloc+igrp*BLOCKCOOL]*factgrp[igrp]*dt*(!chemonly);
      }

      REAL eorg=eint[idloc];
      REAL etorg=egyloc[idloc+0*BLOCKCOOL];
      REAL Nfree=(1.+x0[idloc])*nH[idloc];
      REAL torg=eint[idloc]/(1.5*Nfree*KBOLTZ);

      if(etorg<EMIN) {
       	printf("ERROR : %e %e %e\n",R.e[1],aexporg,param->unit.unit_N);
       	abort();
      }

      /// local cooling loop -------------------------------
      aexp=aexporg;
      fudgecool=param->fudgecool;
      currentcool_t=0.;
      nitcool=0.;
      REAL dtcool_tmp=0;

      while(currentcool_t<dt){

        Nfree=(1.+x0[idloc])*nH[idloc];
        tloc=eint[idloc]/(1.5*Nfree*KBOLTZ);
        REAL tmin=emin/(1.5*Nfree*KBOLTZ);

        cuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,1,(srcloc[idloc]>0.));

        if(fudgecool<1e-20) abort();

        ai_tmp1=0.;
        for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += (alphae[igrp]*hnu[igrp]-alphai[igrp]*hnu0)*egyloc[idloc+igrp*BLOCKCOOL]*nH[idloc]*(1.0-x0[idloc]);
        tcool=FABS(eint[idloc]/(ai_tmp1*(!chemonly)-Cool));
        dtcool=FMIN(fudgecool*tcool,dt-currentcool_t);

        alpha=cucompute_alpha_a(tloc,1.,1.);
        alphab=cucompute_alpha_b(tloc,1.,1.);
        beta=cucompute_beta(tloc,1.,1.);

        // ABSORPTION
        int test = 0;
        for(igrp=0;igrp<NGRP;igrp++){

          if(chemonly){
            et[igrp]=egyloc[idloc+igrp*BLOCKCOOL];
          }else{

            REAL EE=EXP(-dtcool*alphai[igrp]*(1.-x0[idloc])*nH[idloc]);

            et[igrp]=egyloc[idloc+igrp*BLOCKCOOL]*EE;

            //ai_tmp1 = alphai[igrp]*(1.-x0[idloc]);
            //et[igrp]+=(srcloc[idloc+igrp*BLOCKCOOL]*factgrp[igrp])/(ai_tmp1*nH[idloc]+(ai_tmp1==0.))*(1.-EE);

            fxt[igrp]=floc[0+idloc3+igrp*BLOCKCOOL*3]*EE;
            fyt[igrp]=floc[1+idloc3+igrp*BLOCKCOOL*3]*EE;
            fzt[igrp]=floc[2+idloc3+igrp*BLOCKCOOL*3]*EE;
          }
          if((et[igrp]<0)||(isnan(et[igrp]))){
            test=1;
            printf("eint=%e nH=%e x0=%e T=%e N=%e %e %e %e (%e)\n",eint[idloc],nH[idloc],x0[idloc],tloc,et[0],et[1],et[2],etorg,egyloc[idloc+0*BLOCKCOOL]);
          }
        }

        if(test){
          fudgecool=fudgecool/FUDGEFACT;
          continue;
        }

        // IONISATION
        ai_tmp1=0.;
        for(igrp=0;igrp<NGRP;igrp++) ai_tmp1 += alphai[igrp]*et[igrp];
        REAL EE = EXP(-dtcool*ai_tmp1*(!chemonly));
        xt=1.-(1.0 -x0[idloc])*EE;

        REAL deltaX=FABS(xt/x0[idloc]-1.);
        if( (xt>1.) || (xt<0.) || isnan(xt) || (deltaX>FRAC_VAR) ){
          fudgecool/=FUDGEFACT;
          continue;
        }else{
          fudgecool=FMIN(fudgecool*(FUDGEFACT-1.)/2.,param->fudgecool);
        }


        cuCompCooling(tloc,xt,nH[idloc],&Cool,&tcool1,aexp,1,(srcloc[idloc]>0.));

        ai_tmp1=0.;
        for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += (alphae[igrp]*hnu[igrp]-alphai[igrp]*hnu0)*et[igrp]*nH[idloc]*(1.0-xt);
        eintt=(eint[idloc]+ dtcool*((ai_tmp1*(!chemonly))-Cool));

        if(eintt<0.){
          fudgecool=fudgecool/FUDGEFACT;
          continue;
        }

        REAL deltaE=FABS(eintt/eint[idloc]-1.);
        if(deltaE>FRAC_VAR){
          fudgecool=fudgecool/FUDGEFACT;
          continue;
        }else{
          fudgecool=FMIN(fudgecool*(FUDGEFACT-1.)/2.,param->fudgecool);
        }

        /// Cosmological Adiabatic expansion effects ==============
        REAL aold=aexp;
        REAL hubblet=param->cosmo->H0*SQRT(param->cosmo->om/aexp+param->cosmo->ov*(aexp*aexp))/aexp*(1e3/(1e6*PARSEC)); // s-1 // SOMETHING TO CHECK HERE
        REAL da=hubblet*dtcool*aexp;
        aexp+=da;

        for(igrp =0;igrp<NGRP;igrp++){
          egyloc[idloc+igrp*BLOCKCOOL]=et[igrp]*POW(aold/aexp,3); // pourquoi egy tout le temps et les flux seulement si pas chemonly
          if(!chemonly){
            floc[0+idloc3+igrp*BLOCKCOOL*3]=fxt[igrp]*POW(aold/aexp,4);
            floc[1+idloc3+igrp*BLOCKCOOL*3]=fyt[igrp]*POW(aold/aexp,4);
            floc[2+idloc3+igrp*BLOCKCOOL*3]=fzt[igrp]*POW(aold/aexp,4);
          }
        }

        x0[idloc]=xt;

        eintt=FMAX(emin,eintt);
        eint[idloc]=eintt*POW(aold/aexp,5);

        currentcool_t+=dtcool;
        fudgecool=param->fudgecool;
        nitcool++;
        if((nitcool==ncvgcool)&&(ncvgcool!=0)) break;

      }/// ====================== End of the cooling loop


      if(!chemonly){
        for(igrp=0;igrp<NGRP;igrp++){
          R.e[igrp]=FMAX(egyloc[idloc+igrp*BLOCKCOOL]*POW(aexp,3)/param->unit.unit_N,EMIN*factgrp[igrp]/param->unit.unit_N);
          R.fx[igrp]=floc[0+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
          R.fy[igrp]=floc[1+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
          R.fz[igrp]=floc[2+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
        }
      }

      R.nhplus=x0[idloc]*R.nh;
      R.eint=eint[idloc]*POW(aexp,5)/param->unit.unit_n/param->unit.unit_d/POW(param->unit.unit_v,2);
      E2T(&R,aexp,param);
      memcpy(&stencil[i].New.cell[icell].rfieldnew,&R,sizeof(struct Rtype));

    }
  }
}
#endif // WCHEM
#endif // WRAD

