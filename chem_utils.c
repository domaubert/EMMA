
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
  REAL x=R->nhplus/R->nh;
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
  
  //c5=1.017e-37*pow(2.727/aexp,4)*(temp-2.727/aexp)*nh2*x;
  c5=5.406e-36*(temp-2.727/aexp)/pow(aexp,4)*x/(1.+x);
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

void chemrad(struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew, struct RUNPARAMS *param, REAL aexporg)
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

  REAL aexp;
  REAL ebkg[NGRP];
  REAL z=1./aexporg-1.;

  REAL c=param->clight*LIGHT_SPEED_IN_M_PER_S; 			// switch back to physical velocity m/s
  

#ifdef S_X
  REAL E0overI[NGRP];
  REAL N2[NGRP];
  REAL F2[NGRP];
#endif

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
  
  
  REAL dt=dtnew*param->unit.unit_t*pow(aexporg,2);

  REAL emin;
  struct Rtype R;
  REAL fudgecool=param->fudgecool;
  int ncvgcool=param->ncvgcool;
  REAL E0;

  for(i=0;i<nread;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells

     
      if(stencil[i].oct[6].cell[icell].split) continue; // we dont treat split cells
      memcpy(&R,&stencil[i].New.cell[icell].rfieldnew,sizeof(struct Rtype));// We get the local physical quantities after transport update
      
      // switch to physical units, chemistry remains unchanged with and without cosmo
      for (igrp=0;igrp<NGRP;igrp++)
	{			
	  egyloc[idloc+igrp*BLOCKCOOL]   =R.e[igrp]/(aexporg*aexporg*aexporg)/pow(param->unit.unit_l,3)*param->unit.unit_n;//+ebkg[igrp]; 
	  floc[0+idloc3+igrp*BLOCKCOOL*3]=R.fx[igrp]/pow(aexporg,4)/pow(param->unit.unit_l,2)/param->unit.unit_t*param->unit.unit_n;
	  floc[1+idloc3+igrp*BLOCKCOOL*3]=R.fy[igrp]/pow(aexporg,4)/pow(param->unit.unit_l,2)/param->unit.unit_t*param->unit.unit_n;
	  floc[2+idloc3+igrp*BLOCKCOOL*3]=R.fz[igrp]/pow(aexporg,4)/pow(param->unit.unit_l,2)/param->unit.unit_t*param->unit.unit_n;
	}

      x0[idloc]=R.nhplus/R.nh;
      nH[idloc]=R.nh/(aexporg*aexporg*aexporg)/pow(param->unit.unit_l,3)*param->unit.unit_n;


	/* if (nH[idloc] <= 0) { */
	/* 	printf("densité negative");		 */
	/* 	abort();	 */
	/* } */

      eint[idloc]=R.eint/pow(aexporg,5)/pow(param->unit.unit_l,3)*param->unit.unit_n*param->unit.unit_mass*pow(param->unit.unit_v,2);
      E0=eint[idloc];
      emin=PMIN/(GAMMA-1.)/pow(aexporg,5)/pow(param->unit.unit_l,3)*param->unit.unit_n*param->unit.unit_mass*pow(param->unit.unit_v,2); // physical minimal pressure
      srcloc[idloc]=(R.src/pow(param->unit.unit_l,3)*param->unit.unit_n/param->unit.unit_t/(aexporg*aexporg)+ebkg[0])/pow(aexporg,3); 


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


      aexp=aexporg;
      fudgecool=param->fudgecool;
      currentcool_t=0.;
      nitcool=0.;
      REAL da;
      //printf("fudge=%e ncv=%d currentcool_t=%e dt=%e\n",param->fudgecool,ncvgcool,currentcool_t,dt);

      // local cooling loop -------------------------------
      while(currentcool_t<dt)
	{

	  z=1./aexp-1.;
	  //if(eint[idloc]!=E0) printf("1!\n");
	  // ==================== UV Background
#ifdef UVBKG
	  if(NGRP>1) printf("WARNING BAD BEHAVIOR FOR BKG with NGRP>1 !\n");
	  //for(igrp=0;igrp<NGRP;igrp++) ebkg[igrp]=3.6*(z<3?1.:4./(1+z))  ;  // Katz simple model
	  
	  // Poor FIT to Haardt & MAdau 2012
	  for(igrp=0;igrp<NGRP;igrp++){
	    REAL amp=1.2e-16,sig=1.,zavg=2,mz=1e-18,pz=1.2e-17;
	    ebkg[igrp]=amp/(sig*sqrt(2*M_PI))*exp(-pow((z-zavg),2)/(2.*pow(sig,2)))+mz*z+pz; // comoving photons/s/m3
	  }
#else
	  for(igrp=0;igrp<NGRP;igrp++) ebkg[igrp]=0.;
#endif
	  // Cosmological Adiabatic expansion effects ==============
#ifdef TESTCOSMO
	  REAL hubblet=param->cosmo->H0*sqrtf(param->cosmo->om/aexp+param->cosmo->ov*(aexp*aexp))/aexp*(1e3/(1e6*PARSEC))*0.; // s-1 // SOMETHING TO CHECK HERE
#else
	  REAL hubblet=0.;
#endif
	  
	  //if(eint[idloc]!=E0) printf("2!\n");
	  tloc=eint[idloc]/(1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc]));
	  
	  //== Getting a timestep
	  cuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
	  
	  ai_tmp1=0.;

	  //if(eint[idloc]!=E0) printf("3!\n");


	  for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += ((alphae[igrp])*hnu[igrp]-(alphai[igrp])*hnu0)*egyloc[idloc+igrp*BLOCKCOOL];
	  
	  tcool=fabs(eint[idloc]/(nH[idloc]*(1.0-x0[idloc])*ai_tmp1-Cool));
	  ai_tmp1=0.;
	  if(fudgecool<1e-4){
	    /* printf("eint[idloc]=%e (emin=%e E0=%e) nH[idloc]=%e X0[idloc=%e ai=%e Cool=%e tcool=%e, t2=%e fudge=%e\n",eint[idloc],emin,E0,nH[idloc],x0[idloc],ai_tmp1,Cool,fudgecool*tcool,dt-currentcool_t,fudgecool); */
	    if(fudgecool<1e-25){
	      printf("lack of convergence abort() fudgecool=%e\n",fudgecool);
	      abort();
	    }
	  }
	  //if(eint[idloc]!=E0) printf("4!\n");

	  dtcool=fmin(fudgecool*tcool,dt-currentcool_t);
	  
	  alpha=cucompute_alpha_a(tloc,1.,1.)*CLUMPF2;
	  alphab=cucompute_alpha_a(tloc,1.,1.)*CLUMPF2;// cancelling recomb rad
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
	      
	      if((et[igrp]<0)||(isnan(et[igrp]))){
		test=1;
	      }
	      p[igrp]=(1.+(alphai[igrp]*nH[idloc]*(1-x0[idloc])+4.*hubblet)*dtcool);
	    }
	  
	  ai_tmp1=0.;

	  
	  if(test) 
	    {
	      fudgecool=fudgecool/10.; 
	      //printf("loop 1 %e %e %e %e %e\n",et[0],egyloc[idloc],nH[idloc],srcloc[0],dtcool);
	      continue;	
	    } 

	  //if(eint[idloc]!=E0) printf("5!\n");

	  // IONISATION
#ifndef S_X
#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*et[igrp];}
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*egyloc[idloc+igrp*BLOCKCOOL];}
#endif
#else
	  N2[0]=1.0;
	  REAL pp=(1.-pow(x0[idloc],0.4092)); 
	  if(pp<0.) pp=0.; 

	  for(igrp=1;igrp<NGRP;igrp++){
	    N2[igrp]=1.0+0.3908*pow(pp,1.7592)*E0overI[igrp]; 
	    if(N2[igrp]<1.0) N2[igrp]=1.0; 
	  }
#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*et[igrp]*N2[igrp];}
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += alphai[igrp]*egyloc[idloc+igrp*BLOCKCOOL]*N2[igrp];}
#endif
#endif
	  
	  xt=1.-(alpha*x0[idloc]*x0[idloc]*nH[idloc]*dtcool+(1. -x0[idloc]))/(1.+dtcool*(beta*x0[idloc]*nH[idloc]+ai_tmp1+3.*hubblet));
	  ai_tmp1=0.;

	  /* printf("xt=%e alphat=%e betat=%e betat2=%e\n",xt,alpha*x0[idloc]*x0[idloc]*nH[idloc]*dtcool,dtcool*(beta*x0[idloc]*nH[idloc]+ai_tmp1+3.*hubblet),dtcool*(beta*x0[idloc]*nH[idloc]+ai_tmp1)); */
	  /* abort(); */
	  //if(eint[idloc]!=E0) printf("!6\n");

	  if(((xt>1.)||(xt<0.))||(isnan(xt))) 
 	    {
	      fudgecool/=10.; 
	      /* printf("loop 2 x0=%e xt=%e e=%e\n",x0[idloc],xt,et[0]);   */
	      /* printf("loop 3 %e %e %e %e %e\n",et[0],egyloc[idloc],nH[idloc],srcloc[0],dtcool);  */

	      continue;	
	    } 

#ifdef SEMI_IMPLICIT
	  cuCompCooling(tloc,xt,nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
#else
	  cuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
#endif
	  Cool=0.;

#ifdef COOLING
	  // HEATING


#ifdef STARS
		REAL SN 	 = R.snfb;
	  	if (R.snfb) Cool = 0;
		//if (R.snfb) printf("dE\t%e\tE0\t%e\n",R.snfb*dtcool,eintt);
#else
		REAL SN = 0;
#endif


#ifndef S_X
	#ifdef SEMI_IMPLICIT
		for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0));}
		eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool+SN))/(1.+5.*hubblet*dtcool);
	#else
		for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0));}
		eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-x0[idloc])*(ai_tmp1)-Cool+SN))/(1.+5*hubblet*dtcool);
	#endif
#else
	REAL pp2;
	F2[0]=1.0;

	//if(eint[idloc]!=E0) printf("7!\n");

	#ifdef SEMI_IMPLICIT
		  pp2=1.0-pow(xt,0.2663); 
	#else
		  pp2=1.0-pow(x0[idloc],0.2663); 
	#endif
	  if(pp2<0.) pp2=0.; 
	  for(igrp=1;igrp<NGRP;igrp++){
	    F2[igrp]=1.0;
	    F2[igrp]=0.9971*(1.0-pow(pp2,1.3163)); 
	  
	    if(F2[igrp]>1.0) F2[igrp]=1.0; 
	    if(F2[igrp]<0.0) F2[igrp]=0.0; 
	  }

	#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp];}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool+SN))/(1.+5.*hubblet*dtcool);
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp];}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-x0[idloc])*(ai_tmp1)-Cool+SN))/(1.+5*hubblet*dtcool);
	#endif
#endif


	  //if(eint[idloc]!=E0) printf("8!\n");

	  if(eintt<0.)
 	    {
	      fudgecool=fudgecool/10.;
	      //printf("loop 3.5 neg e eorg=%e enew=%e looping\n",eint[idloc],eintt);
	      continue;
	    }

	  if(fabs(eintt-eint[idloc])>FRAC_VAR*eint[idloc])
	    {
	      fudgecool=fudgecool/10.;
	      //printf("iread=%d icell=%d\n",i,icell);
	      /* printf("loop 2 x0=%e xt=%e e=%e\n",x0[idloc],xt,et[0]); */
	      /* printf("loop 3 %e %e %e %e %e\n",et[0],egyloc[idloc],nH[idloc],srcloc[0],dtcool); */
	      //printf("loop 4 %e eint[idloc]=%e tloc=%e dt=%e E0=%e fcool=%e a1=%e\n",eintt,eint[idloc],tloc,dt,E0,fudgecool,ai_tmp1);
	      continue;
	    }
  	  else{

 	    fudgecool=fmin(fudgecool*1.5,param->fudgecool);
	  }

	  ai_tmp1=0;

	  


	  //if(eint[idloc]!=E0) printf("9!\n");

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
	  //printf("xt=%e\n",xt);
#ifdef COOLING
	  eint[idloc]=eintt;
#endif
	  
#ifdef TESTCOSMO
	  REAL da=hubblet*dtcool*aexp;
	  aexp+=da;
#endif
	  currentcool_t+=dtcool;
	  fudgecool=param->fudgecool;
	  nitcool++;

	  //printf("nitcool=%d %e %e %e\n",nitcool,aexp,xt,eintt);
	  if((nitcool==ncvgcool)&&(ncvgcool!=0)) break;
	}
      // ====================== End of the cooling loop

      aexp=aexporg;
      // FIlling the rad structure to send it back
      for(igrp=0;igrp<NGRP;igrp++)
	{
	  R.e[igrp]=fmax(egyloc[idloc+igrp*BLOCKCOOL]*aexp*aexp*aexp,EMIN*factgrp[igrp])*pow(param->unit.unit_l,3)/param->unit.unit_n;
	  R.fx[igrp]=floc[0+idloc3+igrp*BLOCKCOOL*3]*pow(aexp,4)*pow(param->unit.unit_l,2)*param->unit.unit_t/param->unit.unit_n;
	  R.fy[igrp]=floc[1+idloc3+igrp*BLOCKCOOL*3]*pow(aexp,4)*pow(param->unit.unit_l,2)*param->unit.unit_t/param->unit.unit_n;
	  R.fz[igrp]=floc[2+idloc3+igrp*BLOCKCOOL*3]*pow(aexp,4)*pow(param->unit.unit_l,2)*param->unit.unit_t/param->unit.unit_n;
	}
       
      R.nhplus=x0[idloc]*R.nh;

      R.eint=eint[idloc]*pow(aexp,5)*pow(param->unit.unit_l,3)/param->unit.unit_n/param->unit.unit_mass/pow(param->unit.unit_v,2);
      
      memcpy(&stencil[i].New.cell[icell].rfieldnew,&R,sizeof(struct Rtype));

    }
  }

}

#endif
#endif
