
#ifdef WRAD
#ifdef WCHEM

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include <mpi.h>
#include "atomic_data/Atomic.h"


#define FRAC_VAR (0.1)

//================================================================================
void E2T(struct Rtype *R, REAL aexp,struct RUNPARAMS *param){

  REAL tloc;
  REAL eint=R->eint;
  REAL nH=R->nh;
  REAL x=R->xion;
  REAL pstar=param->unit.unit_n*param->unit.unit_mass/pow(param->unit.unit_l,3)*pow(param->unit.unit_v,2);

  nH=nH/pow(aexp,3)/pow(param->unit.unit_l,3)*param->unit.unit_n;
  eint=eint/pow(aexp,5)*pstar;
  tloc=eint/(1.5*nH*KBOLTZ*(1.+x));
  R->temp=tloc;
}


// ============================================================================
REAL cucompute_alpha_b(REAL temp, REAL unit_number, REAL aexp)
{
  // CASE B recombination rate m**3 s*-1
  // temperature should be given in Kelvin
  
  REAL alpha_b,lambda;
  lambda=2e0*157807e0/temp;
  alpha_b=2.753e-14*pow(lambda,1.5)/pow(1e0+pow(lambda/2.740,0.407),2.242); //cm3/s
#ifdef TESTCOSMO
  alpha_b=alpha_b*1e-6*unit_number;///(aexp*aexp*aexp); //m3/s
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
  alpha_a=1.269e-13*pow(lambda,1.503)/pow(1e0+pow(lambda/0.522,0.470),1.923); //cm3/s
#ifdef TESTCOSMO
  alpha_a=alpha_a*1e-6*unit_number;///(aexp*aexp*aexp); //m3/s
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
  beta=5.85e-11*sqrt(temp)/(1+sqrt(T5))*expf(-(157809e0/temp)); //cm3/s
#ifdef TESTCOSMO
  beta=beta*1e-6*unit_number;///(aexp*aexp*aexp); // !m3/s
#else
  beta=beta*1e-6*unit_number; // !m3/s
#endif
  return beta;
}

//**********************************************************************************
//**********************************************************************************
void cuCompCooling(REAL temp, REAL x, REAL nH, REAL *lambda, REAL *tcool, REAL aexp,REAL CLUMPF)
{

  REAL c1,c2,c3,c4,c5,c6;
  REAL unsurtc;
  REAL nh2;

  nh2=nH*1e-6;// ! m-3 ==> cm-3
  

  // Collisional Ionization Cooling

  c1=expf(-157809.1e0/temp)*1.27e-21*sqrt(temp)/(1.f+sqrt(temp/1e5))*x*(1.f-x)*nh2*nh2*CLUMPF;
  

  // Case A Recombination Cooling

  c2=1.778e-29*temp*pow(2e0*157807e0/temp,1.965e0)/pow(1.f+pow(2e0*157807e0/temp/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2*CLUMPF;
  
  
  // Case B Recombination Cooling

  c6=3.435e-30*temp*pow(2e0*157807e0/temp,1.970e0)/pow(1.f+(pow(2e0*157807e0/temp/2.250e0,0.376e0)),3.720e0)*x*x*nh2*nh2*CLUMPF;
  c6=0.;

  // Collisional excitation cooling

  c3=expf(-118348e0/temp)*7.5e-19/(1+sqrt(temp/1e5))*x*(1.f-x)*nh2*nh2*CLUMPF;
  
  
  // Bremmsstrahlung

  c4=1.42e-27*1.5e0*sqrt(temp)*x*x*nh2*nh2*CLUMPF;
  
  // Compton Cooling
  
  c5=1.017e-37*pow(2.727/aexp,4)*(temp-2.727/aexp)*nh2*x;
  
  // Overall Cooling
  
  *lambda=c1+c2+c3+c4+c5+c6;// ! erg*cm-3*s-1
  

  // Unit Conversion

  *lambda=(*lambda)*1e-7*1e6;// ! J*m-3*s-1

  // cooling times

  unsurtc=fmax(c1,c2);
  unsurtc=fmax(unsurtc,c3);
  unsurtc=fmax(unsurtc,c4);
  unsurtc=fmax(unsurtc,fabs(c5));
  unsurtc=fmax(unsurtc,c6)*1e-7;// ==> J/cm3/s

  *tcool=1.5e0*nh2*(1.+x)*KBOLTZ*temp/unsurtc; //Myr
}

// ===========================================================================================================================

void chemrad(struct OCT *octstart, struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew, struct RUNPARAMS *param, REAL aexp)
{
  int i,icell,igrp;
  int idloc;
  int nitcool=0;

  REAL hnu0=13.6*1.6022e-19,
    Cool,
    tcool,
    dtcool,
    tcool1,
    currentcool_t=0.,
    alpha,
    alphab,
    beta,
    tloc,
    xt,
    eintt,
    ai_tmp1=0.,
    hnu[NGRP],		// ! Average Photon Energy (J)
    factgrp[NGRP],		
    alphae[NGRP],
    alphai[NGRP],		
    et[NGRP],
    p[NGRP];

#ifdef TESTCOSMO
  REAL hubblet=0.;//;param->cosmo->H0*sqrtf(param->cosmo->om/aexp+param->cosmo->ov*(aexp*aexp))/aexp*(1e3/(1e6*PARSEC)); // s-1

#else
  REAL hubblet=0.;
#endif

  REAL ebkg[NGRP];
  REAL z=1./aexp-1.;
#ifdef UVBKG
  for(igrp=0;igrp<NGRP;igrp++) ebkg[igrp]=3.6*(z<3?1.:4./(1+z))  ;  // Katz simple model
#else
  for(igrp=0;igrp<NGRP;igrp++) ebkg[igrp]=0.;
#endif

  REAL c=param->clight*LIGHT_SPEED_IN_M_PER_S; 			// switch back to physical velocity m/s
  SECTION_EFFICACE; // defined in Atomic.h
  FACTGRP; //defined in Atomic.h

#define BLOCKCOOL 1 // KEPT FROM CUDATON FOR SIMPLICITY
#define idloc 0 // KEPT FROM CUDATON FOR SIMPLICITY
#define idloc3 0 // KEPT FROM CUDATON FOR SIMPLICITY

  REAL
    egyloc[BLOCKCOOL*NGRP],
    floc[3*BLOCKCOOL*NGRP],
    x0[BLOCKCOOL],
    nH[BLOCKCOOL],
    eint[BLOCKCOOL],
    srcloc[BLOCKCOOL];
  
  
  REAL dt=dtnew*param->unit.unit_t*pow(aexp,2);
  REAL emin;
  struct Rtype R;
  REAL fudgecool=param->fudgecool;
  int ncvgcool=param->ncvgcool;

  for(i=0;i<nread;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells

      fudgecool=param->fudgecool;
      if(stencil[i].oct[6].cell[icell].split) continue; // we dont treat split cells
      memcpy(&R,&stencil[i].New.cell[icell].rfieldnew,sizeof(struct Rtype));// We get the local physical quantities after transport update
      
      // switch to physical units, chemistry remains unchanged with and without cosmo
      for (igrp=0;igrp<NGRP;igrp++)
	{			
	  egyloc[idloc+igrp*BLOCKCOOL]   =R.e[igrp]/(aexp*aexp*aexp)/pow(param->unit.unit_l,3)*param->unit.unit_n+ebkg[igrp]; 
	  floc[0+idloc3+igrp*BLOCKCOOL*3]=R.fx[igrp]/pow(aexp,4)/pow(param->unit.unit_l,2)/param->unit.unit_t*param->unit.unit_n;
	  floc[1+idloc3+igrp*BLOCKCOOL*3]=R.fy[igrp]/pow(aexp,4)/pow(param->unit.unit_l,2)/param->unit.unit_t*param->unit.unit_n;
	  floc[2+idloc3+igrp*BLOCKCOOL*3]=R.fz[igrp]/pow(aexp,4)/pow(param->unit.unit_l,2)/param->unit.unit_t*param->unit.unit_n;
	}


      x0[idloc]=R.xion;
      nH[idloc]=R.nh/(aexp*aexp*aexp)/pow(param->unit.unit_l,3)*param->unit.unit_n;
      eint[idloc]=R.eint/pow(aexp,5)/pow(param->unit.unit_l,3)*param->unit.unit_n*param->unit.unit_mass*pow(param->unit.unit_v,2);
      emin=PMIN/(GAMMA-1.)/pow(aexp,5)/pow(param->unit.unit_l,3)*param->unit.unit_n*param->unit.unit_mass*pow(param->unit.unit_v,2); // physical minimal pressure
      srcloc[idloc]=R.src/pow(param->unit.unit_l,3)*param->unit.unit_n/param->unit.unit_t/(aexp*aexp)/pow(aexp,3); 
      
      /* if(((isnan(eint[idloc]))||(isnan(x0[idloc])))||(eint[idloc]==0.)){ */
      /* 	printf("start with nans or ZErO egy %e\n",eint[idloc]); */
      /* 	abort(); */
      /* } */

      // at this stage we are ready to do the calculations

      // DEALING WITH CLUMPING ----------------------
#ifdef WCLUMP
      REAL CLUMPF2=fmin(fmax(pow(nH[idloc]/6.,0.7),1.),40.);
      REAL CLUMPI=1.;
#else
      REAL CLUMPF2=1.;
      REAL CLUMPI=1.;
#endif


      for(igrp=0;igrp<NGRP;igrp++)
	{
	  alphai[igrp] *= CLUMPI;
	  alphae[igrp] *= CLUMPI;
	}

      // -------------------------------------------------


      // local cooling loop -------------------------------

      currentcool_t=0.;
      nitcool=0.;
      
      while(currentcool_t<dt)
	{
	  nitcool++;

	  tloc=eint[idloc]/(1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc]));
	  
	  //== Getting a timestep
	  cuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
	  
	  ai_tmp1=0.;
	  for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += ((alphae[igrp])*hnu[igrp]-(alphai[igrp])*hnu0)*egyloc[idloc+igrp*BLOCKCOOL];
	  
	  tcool=fabs(eint[idloc]/(nH[idloc]*(1.0-x0[idloc])*ai_tmp1-Cool));
	  ai_tmp1=0.;
	  dtcool=fmin(fudgecool*tcool,dt-currentcool_t);
	  
	  alpha=cucompute_alpha_a(tloc,1.,1.)*CLUMPF2;
	  alphab=cucompute_alpha_b(tloc,1.,1.)*CLUMPF2;
	  beta=cucompute_beta(tloc,1.,1.)*CLUMPF2;
      
	  //== Update
	  
	  // ABSORPTION
	  int test = 0;
	  REAL factotsa[NGRP];
	  for(igrp=0;igrp<NGRP;igrp++)
	    {
#ifdef OTSA
	      factotsa[igrp]=0;
#else
	      factotsa[igrp]=(igrp==0);
#endif

	      ai_tmp1 = alphai[igrp];
	      et[igrp]=((alpha-alphab)*x0[idloc]*x0[idloc]*nH[idloc]*nH[idloc]*dtcool*factotsa[igrp]+egyloc[idloc+igrp*BLOCKCOOL]+srcloc[idloc]*dtcool*factgrp[igrp])/(1.+dtcool*(ai_tmp1*(1.-x0[idloc])*nH[idloc]+3*hubblet));
	      
	      if(et[igrp]<0){
		test=1;
	      }
	      p[igrp]=(1.+(alphai[igrp]*nH[idloc]*(1-x0[idloc])+2*hubblet)*dtcool);
	    }
	  ai_tmp1=0.;

	  
	  if(test) 
	    {
	      fudgecool=fudgecool/10.; 
	      continue;	
	    } 

	  /* if((et[0]<0.)||(isnan(et[0]))){ */
	  /*   printf("bnd err ? x0=%e eint=%e nH=%e src=%e xt=%e et=%e eintt=%e %e %d\n",x0[idloc],eint[idloc],nH[idloc],srcloc[idloc],xt,et[0],eintt,fudgecool,nitcool); */
	  /* } */
	  
	  // IONISATION
#ifndef S_X
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*et[igrp];}
	  //for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*egyloc[idloc+igrp*BLOCKCOOL];}
#else
	  N2[0]=1.0;
	  REAL pp=(1.-pow(x0[idloc],0.4092)); 
	  if(pp<0.) pp=0.; 
    
	  N2[1]=1.0+0.3908*pow(pp,1.7592)*E0overI[1]; 
	  if(N2[1]<1.0) N2[1]=1.0; 
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*et[igrp]*N2[igrp];}
#endif
	  
	  xt=1.-(alpha*x0[idloc]*x0[idloc]*nH[idloc]*dtcool+(1. -x0[idloc]))/(1.+dtcool*(beta*x0[idloc]*nH[idloc]+ai_tmp1));
	  ai_tmp1=0.;

	  if(((xt>1.)||(xt<0.))||(isnan(xt))) 
 	    {
	      fudgecool/=10.; 
 	      /* if(nitcool>10000){ */
	      /* 	printf("loopin xion g x0=%e eint=%e nH=%e src=%e xt=%e et=%e eintt=%e %e\n",x0[idloc],eint[idloc],nH[idloc],srcloc[idloc],xt,et[0],eintt,fudgecool); */
	      /* } */
	      continue;	
	    } 

	  cuCompCooling(tloc,xt,nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
	  //cuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);

#ifdef COOLING
	  // HEATING
#ifndef S_X
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0));}
	  //for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0));}
#else
	  REAL pp2;
	  F2[0]=1.0;
	  F2[1]=1.0;
	  pp2=1.0-pow(xt,0.2663); 
	  if(pp2<0.) pp2=0.; 
	  F2[1]=0.9971*(1.0-pow(pp2,1.3163)); 
	  
	  if(F2[1]>1.0) F2[1]=1.0; 
	  if(F2[1]<0.0) F2[1]=0.0; 
	  
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp];}
#endif
	  
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool))/(1.+3*hubblet*dtcool);
	  ai_tmp1=0;


	  if(eintt<0.)
 	    {
	      fudgecool=fudgecool/10.;
	      /* if(nitcool>10000){ */
	      /* 	printf("loopin g e<0  x0=%e eint=%e nH=%e src=%e xt=%e et=%e eintt=%e %e\n",x0[idloc],eint[idloc],nH[idloc],srcloc[idloc],xt,et[0],eintt,fudgecool); */
	      /* } */

	      continue;
	    }

	  if(fabs(eintt-eint[idloc])>FRAC_VAR*eint[idloc])
	    {
	      if(srcloc[idloc]==0.){
		fudgecool=fudgecool/10.;
	      /* if(nitcool>10000){ */
	      /* 	printf("loopin fracvar g x0=%e eint=%e nH=%e src=%e xt=%e et=%e eintt=%e f=%e\n",x0[idloc],eint[idloc],nH[idloc],srcloc[idloc],xt,et[0],eintt,fudgecool); */
	      /* } */

		continue;
	      }
	    }
  	  else{
 	    fudgecool=fmin(fudgecool*1.5,param->fudgecool);
	      /* if(nitcool>10000){ */
	      /* 	printf(" going up loopin g x0=%e eint=%e nH=%e src=%e xt=%e et=%e eintt=%e %e\n",x0[idloc],eint[idloc],nH[idloc],srcloc[idloc],xt,et[0],eintt,fudgecool); */
	      /* } */
	  }

	  eintt=fmax(emin,eintt);

	  
#else
	  eintt=eint[idloc];
#endif
	  
	  for(igrp =0;igrp<NGRP;igrp++)
	    {
	      egyloc[idloc+igrp*BLOCKCOOL]=et[igrp];
	      floc[0+idloc3+igrp*BLOCKCOOL*3]=floc[0+idloc3+igrp*BLOCKCOOL*3]/p[igrp];
	      floc[1+idloc3+igrp*BLOCKCOOL*3]=floc[1+idloc3+igrp*BLOCKCOOL*3]/p[igrp];
	      floc[2+idloc3+igrp*BLOCKCOOL*3]=floc[2+idloc3+igrp*BLOCKCOOL*3]/p[igrp];	
	    }
	  
	  x0[idloc]=xt;

#ifdef COOLING
	  eint[idloc]=eintt;
#endif
	  currentcool_t+=dtcool;
	  fudgecool=param->fudgecool;
	  if((nitcool==ncvgcool)&&(ncvgcool!=0)) break;
	}

      // ====================== End of the cooling loop

      // FIlling the rad structure to send it back
      for(igrp=0;igrp<NGRP;igrp++)
	{
	  R.e[igrp]=fmax(egyloc[idloc+igrp*BLOCKCOOL]*aexp*aexp*aexp,EMIN*factgrp[igrp])*pow(param->unit.unit_l,3)/param->unit.unit_n;
	  R.fx[igrp]=floc[0+idloc3+igrp*BLOCKCOOL*3]*pow(aexp,4)*pow(param->unit.unit_l,2)*param->unit.unit_t/param->unit.unit_n;
	  R.fy[igrp]=floc[1+idloc3+igrp*BLOCKCOOL*3]*pow(aexp,4)*pow(param->unit.unit_l,2)*param->unit.unit_t/param->unit.unit_n;
	  R.fz[igrp]=floc[2+idloc3+igrp*BLOCKCOOL*3]*pow(aexp,4)*pow(param->unit.unit_l,2)*param->unit.unit_t/param->unit.unit_n;
	}
       
      /* if(R.e[0]==0.){ */
      /* 	printf("nil energy x0=%e eint=%e nH=%e src=%e xt=%e et=%e eintt=%e f=%e\n",R.xion,R.eint,nH[idloc],srcloc[idloc],xt,et[0],eintt,fudgecool); */
      /* 	printf("nil energy, abort()\n"); */
      /* 	abort(); */
      /* } */
      //R.temp=tloc;
      R.xion=x0[idloc];

      R.eint=eint[idloc]*pow(aexp,5)*pow(param->unit.unit_l,3)/param->unit.unit_n/param->unit.unit_mass/pow(param->unit.unit_v,2);
      
      memcpy(&stencil[i].New.cell[icell].rfieldnew,&R,sizeof(struct Rtype));

    }
  }

}

#endif
#endif
