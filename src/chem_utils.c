/**
  * \file chem_utils.c
  */

#ifdef WRAD

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

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

# if 1
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



#ifdef HESIMPLE
  REAL hnu0HE=24.6*1.6022e-19;
  REAL hnu0HE2=54.4*1.6022e-19;
#endif

  REAL aexp;
  REAL ebkg[NGRP];
  REAL z=1./aexporg-1.;

  REAL c=param->clightorg*LIGHT_SPEED_IN_M_PER_S; 			// switch back to physical velocity m/s

  REAL hnu[NGRP];
  REAL alphae[NGRP];
  REAL alphai[NGRP];
  REAL factgrp[NGRP];

#ifdef HESIMPLE
  REAL alphaeHE[NGRP];
  REAL alphaiHE[NGRP];
  REAL alphaeHE2[NGRP];
  REAL alphaiHE2[NGRP];
#endif

  for(igrp=0;igrp<NGRP;igrp++) {
    hnu[igrp]=param->atomic.hnu[igrp];
    alphae[igrp]=param->atomic.alphae[igrp];
    alphai[igrp]=param->atomic.alphai[igrp];
    factgrp[igrp]=param->atomic.factgrp[igrp];
#ifdef HESIMPLE
    alphaeHE[igrp]=param->atomic.alphaeHE[igrp];
    alphaiHE[igrp]=param->atomic.alphaiHE[igrp];

    alphaeHE2[igrp]=param->atomic.alphaeHE2[igrp];
    alphaiHE2[igrp]=param->atomic.alphaiHE2[igrp];
#endif
  }

#ifdef S_X
  REAL E0overI[NGRP];
  REAL N2[NGRP];
  REAL F2[NGRP];
#endif

#define BLOCKCOOL 1 // KEPT FROM CUDATON FOR SIMPLICITY
#define idloc3 0 // KEPT FROM CUDATON FOR SIMPLICITY

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
  REAL deltaE;
  int compcool; // do we need to compute the cooling ?

  for(i=0;i<nread;i++){  // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells

      if(stencil[i].oct[6].cell[icell].split) continue; // we dont treat split cells

      memcpy(&R,&stencil[i].New.cell[icell].rfieldnew,sizeof(struct Rtype));// We get the local physical quantities after transport update

#ifdef HOMOSOURCE
      // we override the value with the homogeneous source density
      R.src=param->bkg;
#endif


	  //if(eint[idloc]!=E0) printf("1!\n");
	  /// ==================== UV Background
#ifdef UVBKG
	  if(NGRP>1) printf("WARNING BAD BEHAVIOR FOR BKG with NGRP>1 !\n");
	  //for(igrp=0;igrp<NGRP;igrp++) ebkg[igrp]=3.6*(z<3?1.:4./(1+z))  ;  // Katz simple model

	  // Poor FIT to Haardt & MAdau 2012
  /*
	  for(igrp=0;igrp<NGRP;igrp++){
	    REAL amp=1.2e-16,sig=1.,zavg=2,mz=1e-18,pz=1.2e-17;
	    ebkg[igrp]=amp/(sig*SQRT(2*M_PI))*exp(-POW((z-zavg),2)/(2.*POW(sig,2)))+mz*z+pz; // comoving photons/s/m3
	  }
  */

#else
	  for(igrp=0;igrp<NGRP;igrp++) ebkg[igrp]=0.;
#endif

      // switch to physical units, chemistry remains unchanged with and without cosmo
      for (igrp=0;igrp<NGRP;igrp++)
	{
	  egyloc[idloc+igrp*BLOCKCOOL]   =R.e[igrp]/(aexporg*aexporg*aexporg)*param->unit.unit_N;//+ebkg[igrp];
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
	srcloc[idloc+igrp*BLOCKCOOL]=(R.src[igrp]*param->unit.unit_N/param->unit.unit_t/(aexporg*aexporg))/POW(aexporg,3); //phot/s/dv (physique)
      }

      // R.src phot/unit_t/unit_dv (comobile)
      REAL eorg=eint[idloc];
      REAL etorg=egyloc[idloc+0*BLOCKCOOL];



      // ========================== POLYTROP ================================

      compcool=1;
      REAL eintschaye;
      int onpoly=0;

#ifdef POLYTROP
      if((nH[idloc]>1e6)&&(R.nh>0)){
	eintschaye=(1.08e9*KBOLTZ)*POW(nH[idloc]/1e5,4./3.)/(GAMMA-1)/FSCHAYE; // polytropic EOS
	if(eint[idloc]<eintschaye){
	  eint[idloc]=eintschaye;
	  onpoly=1;
	}
      }
      else if(srcloc[idloc]>0.){
	eintschaye=(1.08e9*KBOLTZ)*POW(nH[idloc]/1e5,4./3.)/(GAMMA-1)/FSCHAYE; // polytropic EOS
	if(eint[idloc]<eintschaye){
	  eint[idloc]=eintschaye;
	  onpoly=1;
	}

      }
#endif

      // ========================== END POLYTROP ================================




      //REAL torg=eint[idloc]/(1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc]));

      REAL Nfree;
#ifdef HESIMPLE
      Nfree=(1.+x0[idloc])*nH[idloc]*(1.+yHE);
#else
      Nfree=(1.+x0[idloc])*nH[idloc];
#endif
      REAL torg=eint[idloc]/(1.5*Nfree*KBOLTZ);

       if(etorg<EMIN) {
       	printf("ERROR : %e %e %e\n",R.e[1],aexporg,param->unit.unit_N);
       	abort();
       }



      /*  if(xorg>0.01){ */
      /* 	printf("SN HEAT xorg=%e nhp=%e nh=%e t=%e src=%e\n",xorg,R.nhplus,R.nh,torg,srcloc[idloc]); */
      /* } */

      //if(srcloc[0]>0) 	printf("nh=%e %e %e %e\n",R.nh,R.e[0],eint[idloc],3[idloc]);

      // at this stage we are ready to do the calculations

      // DEALING WITH CLUMPING ----------------------
#ifdef WCLUMP
      REAL CLUMPF2=FMIN(FMAX(POW(nH[idloc]/6.,0.7),1.),40.);
      REAL CLUMPI=1.;
#else
      REAL CLUMPF2=1.;
      REAL CLUMPI=1.;
#endif


      for(igrp=0;igrp<NGRP;igrp++)
	{
	  alphai[igrp] *= CLUMPI;
	  alphae[igrp] *= CLUMPI;
#ifdef HESIMPLE
	  alphaiHE[igrp] *= CLUMPI;
	  alphaeHE[igrp] *= CLUMPI;
	  alphaiHE2[igrp] *= CLUMPI;
	  alphaeHE2[igrp] *= CLUMPI;
#endif
	}



      // -------------------------------------------------

      /// local cooling loop -------------------------------
      aexp=aexporg;
      fudgecool=param->fudgecool;
      currentcool_t=0.;
      nitcool=0.;
      REAL da;
      //printf("cpu=%d fudge=%e ncv=%d currentcool_t=%e dt=%e\n",cpu->rank,param->fudgecool,ncvgcool,currentcool_t,dt);

      // local cooling loop -------------------------------
      while(currentcool_t<dt)
	{


	  /// Cosmological Adiabatic expansion effects ==============
#ifdef TESTCOSMO
	  REAL hubblet=param->cosmo->H0*SQRT(param->cosmo->om/aexp+param->cosmo->ov*(aexp*aexp))/aexp*(1e3/(1e6*PARSEC)); // s-1 // SOMETHING TO CHECK HERE
#else
	  REAL hubblet=0.;
#endif


	  //tloc=eint[idloc]/(1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc]));
#ifdef HESIMPLE
	  Nfree=(1.+x0[idloc])*nH[idloc]*(1.+yHE);
#else
	  Nfree=(1.+x0[idloc])*nH[idloc];
#endif

	  tloc=eint[idloc]/(1.5*Nfree*KBOLTZ);
	  REAL tmin=emin/(1.5*Nfree*KBOLTZ);
	  //== Getting a timestep
	  cuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2,(srcloc[idloc]>0.));
	  /* if((srcloc[idloc+igrp*BLOCKCOOL]>0.)){ */
	  /*   Cool*=1e-5; */
	  /* } */

	  if(fudgecool<1e-20){
	    printf("ERROR : eint=%e(%e<%e) eint_temp=%e (delta=%e)  nH=%e x0=%e(%e) T=%e(%e<%e) N=%e %e %e (%e)\n",eint[idloc],eorg,emin,eintt,deltaE,nH[idloc],x0[idloc],xorg,tloc,torg,tmin,et[0],et[1],et[2],etorg);
	    if(fudgecool<1e-20) abort();
	  }

	  ai_tmp1=0.;
#ifdef HESIMPLE
	  for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += ((((alphae[igrp])*hnu[igrp]-(alphai[igrp])*hnu0)
						     +yHE*((alphaeHE[igrp])*hnu[igrp]-(alphaiHE[igrp])*hnu0HE))*(1.0-x0[idloc])
						    +yHE*x0[idloc]*(alphaeHE2[igrp]*hnu[igrp]-alphaiHE2[igrp]*hnu0HE2))*egyloc[idloc+igrp*BLOCKCOOL]*nH[idloc];
#else
	  for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += (alphae[igrp]*hnu[igrp]-alphai[igrp]*hnu0)*egyloc[idloc+igrp*BLOCKCOOL]*nH[idloc]*(1.0-x0[idloc]);
#endif

	  tcool=FABS(eint[idloc]/(ai_tmp1*(!chemonly)-Cool));
	  //tcool=FABS(eint[idloc]/(Cool));
	  ai_tmp1=0.;
	  dtcool=FMIN(fudgecool*tcool,dt-currentcool_t);

	  alpha=cucompute_alpha_a(tloc,1.,1.)*CLUMPF2;
	  alphab=cucompute_alpha_b(tloc,1.,1.)*CLUMPF2;
	  beta=cucompute_beta(tloc,1.,1.)*CLUMPF2;

	  //== Update

	  // ABSORPTION
	  int test = 0;
	  REAL factotsa[NGRP];
	  for(igrp=0;igrp<NGRP;igrp++){

#ifdef OTSA
	  factotsa[igrp]=0;
	  alpha=alphab; // recombination is limited to non ground state levels
#else
	  factotsa[igrp]=(igrp==0);
#endif

#ifdef HESIMPLE
	  ai_tmp1 = (alphai[igrp]+yHE*alphaiHE[igrp])*(1.-x0[idloc])+(alphaiHE2[igrp])*(x0[idloc])*yHE;
#else
	  ai_tmp1 = alphai[igrp]*(1.-x0[idloc]);

#endif

	  if(chemonly){
	  et[igrp]=egyloc[idloc+igrp*BLOCKCOOL];
	}
	  else{

	    REAL EARG=dtcool*ai_tmp1*nH[idloc];
	    REAL EE;
	    REAL ONEMINUSEE;
	    if(EARG<1e-3){
	      REAL DL=-EARG+0.5*EARG*EARG-EARG*EARG*EARG/6.;
	      EE=1.+DL;
	      ONEMINUSEE=-DL;
	    }
	    else{
	      EE=EXP(-EARG);
	      ONEMINUSEE=1.-EE;
	    }

#ifdef HESIMPLE
	    //et[igrp]=(egyloc[idloc+igrp*BLOCKCOOL]+srcloc[idloc+igrp*BLOCKCOOL]*dtcool*factgrp[igrp])/(1.+dtcool*(ai_tmp1*nH[idloc]));
	    et[igrp]=egyloc[idloc+igrp*BLOCKCOOL]*+(srcloc[idloc+igrp*BLOCKCOOL]*factgrp[igrp])/(ai_tmp1*nH[idloc]+(ai_tmp1==0.))*(ONEMINUSEE);
#else

#ifndef OTSA
	    et[igrp]=egyloc[idloc+igrp*BLOCKCOOL]*EE+(srcloc[idloc+igrp*BLOCKCOOL]*factgrp[igrp]+(alpha-alphab)*x0[idloc]*x0[idloc]*nH[idloc]*nH[idloc]*factotsa[igrp])/(ai_tmp1*nH[idloc]+(ai_tmp1==0.))*(ONEMINUSEE);
#else
	    et[igrp]=egyloc[idloc+igrp*BLOCKCOOL]*EE+(srcloc[idloc+igrp*BLOCKCOOL]*factgrp[igrp])/(ai_tmp1*nH[idloc]+(ai_tmp1==0.))*(ONEMINUSEE);
#endif

	    //et[igrp]=((alpha-alphab)*x0[idloc]*x0[idloc]*nH[idloc]*nH[idloc]*dtcool*factotsa[igrp]+egyloc[idloc+igrp*BLOCKCOOL]+srcloc[idloc+igrp*BLOCKCOOL]*dtcool*factgrp[igrp])/(1.+dtcool*(ai_tmp1*nH[idloc]));
#endif

	    fxt[igrp]=floc[0+idloc3+igrp*BLOCKCOOL*3]*EE;
	    fyt[igrp]=floc[1+idloc3+igrp*BLOCKCOOL*3]*EE;
	    fzt[igrp]=floc[2+idloc3+igrp*BLOCKCOOL*3]*EE;

	  }

	  if((et[igrp]<0)||(isnan(et[igrp]))){
	    test=1;
	    printf("eint=%e nH=%e x0=%e T=%e N=%e %e %e %e (%e)\n",eint[idloc],nH[idloc],x0[idloc],tloc,et[0],et[1],et[2],etorg,egyloc[idloc+0*BLOCKCOOL]);
	  }

	  //p[igrp]=(1.+(ai_tmp1*nH[idloc])*dtcool);

	  }

	  if(test)
	    {
	      fudgecool=fudgecool/10.;
	      continue;
	    }

	  // IONISATION

	  ai_tmp1=0.;
#ifndef S_X
#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*et[igrp]*(!chemonly);}
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*egyloc[idloc+igrp*BLOCKCOOL]*(!chemonly);}
#endif
#else
	  N2[0]=1.0;
	  REAL pp=(1.-POW(x0[idloc],0.4092));
	  if(pp<0.) pp=0.;

	  for(igrp=1;igrp<NGRP;igrp++){
	    N2[igrp]=1.0+0.3908*POW(pp,1.7592)*E0overI[igrp];
	    if(N2[igrp]<1.0) N2[igrp]=1.0;
	  }
#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*et[igrp]*N2[igrp]*(!chemonly);}
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*egyloc[idloc+igrp*BLOCKCOOL]*N2[igrp]*(!chemonly);}
#endif
#endif

	  REAL RECT=alpha*x0[idloc]*x0[idloc]*nH[idloc]*dtcool;
#ifdef NORECSOURCE
	  if(srcloc[idloc]>0.) RECT*=NORECSOURCE_FACT;
#endif
	  xt=1.-(RECT+(1.0 -x0[idloc]))/(1.+dtcool*(beta*x0[idloc]*nH[idloc]+ai_tmp1));

	  /* if(srcloc[idloc]>0.){ */
	  /*   printf("xt=%e x0=%e xorg=%e RECT=%e\n",xt,x0[idloc],xorg,RECT); */
	  /*   abort(); */
	  /* } */
	  ai_tmp1=0.;


	  if(((xt>1.)||(xt<0.))||(isnan(xt)))
 	    {
	      //printf("XION ERR eintt=%e xt=%e et=%e %e %e\n",eintt,xt,et[0],et[1],et[2]);
	      fudgecool/=10.;
	      continue;
	    }

#ifdef SEMI_IMPLICIT
	  cuCompCooling(tloc,xt,nH[idloc],&Cool,&tcool1,aexp,CLUMPF2,(srcloc[idloc]>0.));
#else
	  cuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2,(srcloc[idloc]>0.));
#endif

	  /* if((srcloc[idloc+igrp*BLOCKCOOL]>0.)){ */
	  /*   Cool*=1e-5; */
	  /* } */


#ifdef COOLING
	  // HEATING + COOLING


#ifdef SCHAYE
	  if((nH[idloc]>1e6)&&(R.nh>(param->stars->overdensity_cond*navg))){
	    REAL tlocs;
	    REAL eintschaye;
	    tlocs=eintt/(1.5*nH[idloc]*KBOLTZ*(1.+xt)*(1.+yHE));
	    eintschaye=(1.08e9*KBOLTZ)*POW(nH[idloc]/1e5,4./3.)/(GAMMA-1)/FSCHAYE; // polytropic EOS
	    if(eintt<eintschaye){
	      eintt=eintschaye;
	      compcool=0.; // cancel cooling calculation
	      fudgecool=FMIN(fudgecool*1.5,param->fudgecool);
	    }
	  }
#endif // SCHAYE





	  if(compcool){

#ifndef S_X
#ifdef SEMI_IMPLICIT


	    ai_tmp1=0.;
#ifdef HESIMPLE
	    for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += (((alphae[igrp]*hnu[igrp]-alphai[igrp]*hnu0)+yHE*(alphaeHE[igrp]*hnu[igrp]-alphaiHE[igrp]*hnu0HE))*(1.0-xt)+yHE*(alphaeHE2[igrp]*hnu[igrp]-alphaiHE2[igrp]*hnu0HE2)*xt)*et[igrp]*nH[idloc];
	  //for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += ((alphae[igrp])*hnu[igrp]-(alphai[igrp])*hnu0)*et*nH[idloc]*(1.0-xt);
#else
	  for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += (alphae[igrp]*hnu[igrp]-alphai[igrp]*hnu0)*et[igrp]*nH[idloc]*(1.0-xt);
#endif
	  //for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*(!chemonly);}
	  //eintt=(eint[idloc]+ dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool+SN));

	  /* double HminusC; */
	  /* HminusC=(double)(ai_tmp1*(!chemonly))-(double)(Cool); */


	  eintt=(eint[idloc]+ dtcool*((ai_tmp1*(!chemonly))-Cool));
	  //eintt=(eint[idloc]+ dtcool*(HminusC+SN));
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-x0[idloc])*(ai_tmp1)-Cool));
#endif //SEMI


#else
	  //===================================== X RAYS ==============================
	  REAL pp2;
	  F2[0]=1.0;

	  //if(eint[idloc]!=E0) printf("7!\n");

#ifdef SEMI_IMPLICIT
	  pp2=1.0-POW(xt,0.2663);
#else
	  pp2=1.0-POW(x0[idloc],0.2663);
#endif
	  if(pp2<0.) pp2=0.;
	  for(igrp=1;igrp<NGRP;igrp++){
	    F2[igrp]=1.0;
	    F2[igrp]=0.9971*(1.0-POW(pp2,1.3163));

	    if(F2[igrp]>1.0) F2[igrp]=1.0;
	    if(F2[igrp]<0.0) F2[igrp]=0.0;
	  }

	  ai_tmp1=0.;
#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp]*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool+SN));
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp]*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-x0[idloc])*(ai_tmp1)-Cool+SN));
#endif
	  //================================================================================
#endif //S_X

	  if(onpoly){
	    if(eintt<eintschaye) eintt=eintschaye;
	  }

	  if(eintt<0.)
 	    {
	    //printf("E NEG eintt=%e xt=%e et=%e eorg=%e torg=%e\n",eintt,xt,et[0],eorg,torg);
	      fudgecool=fudgecool/10.;
	      continue;
	    }

	  deltaE=FABS(eintt/eint[idloc]-1.);
	  if(deltaE>FRAC_VAR)
	    {
	    //printf("DELTA E eintt=%e xt=%e eint0=%e eorg=%e torg=%e delta=%e dtcool=%e ft=%e\n",eintt,xt,eint[idloc],eorg,torg,deltaE,dtcool,fudgecool*tcool);
	    fudgecool=fudgecool/10.;
	    continue;
	    //}
	  }
  	  else{
 	    fudgecool=FMIN(fudgecool*1.5,param->fudgecool);
	  }

	  ai_tmp1=0;

	  eintt=FMAX(emin,eintt);

 	  }

#else
	  eintt=eint[idloc];
#endif

	  // inner update
	  REAL aold=aexp;
#ifdef TESTCOSMO
	  da=hubblet*dtcool*aexp;
	  aexp+=da;
#endif

	  for(igrp =0;igrp<NGRP;igrp++)
	    {
	      egyloc[idloc+igrp*BLOCKCOOL]=et[igrp]*POW(aold/aexp,3);
	      if(!chemonly){
		/* floc[0+idloc3+igrp*BLOCKCOOL*3]=floc[0+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4); */
		/* floc[1+idloc3+igrp*BLOCKCOOL*3]=floc[1+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4); */
		/* floc[2+idloc3+igrp*BLOCKCOOL*3]=floc[2+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4); */

		floc[0+idloc3+igrp*BLOCKCOOL*3]=fxt[igrp]*POW(aold/aexp,4);
		floc[1+idloc3+igrp*BLOCKCOOL*3]=fyt[igrp]*POW(aold/aexp,4);
		floc[2+idloc3+igrp*BLOCKCOOL*3]=fzt[igrp]*POW(aold/aexp,4);
	      }
	    }

	  x0[idloc]=xt;
	  //printf("xt=%e\n",xt);
#ifdef COOLING
	  eint[idloc]=eintt*POW(aold/aexp,5);
#endif

	  currentcool_t+=dtcool;
	  fudgecool=param->fudgecool;
	  nitcool++;
	  if((nitcool==ncvgcool)&&(ncvgcool!=0)) break;
	}

      /// ====================== End of the cooling loop

      //aexp=aexporg;
      // FIlling the rad structure to send it back

      if(!chemonly){
	for(igrp=0;igrp<NGRP;igrp++)
	  {
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
#endif // 1
#endif // WCHEM
#endif // WRAD

