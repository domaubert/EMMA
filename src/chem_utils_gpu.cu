
#ifdef WRAD
#ifdef WCHEM

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include <mpi.h>
//#include "atomic_data/Atomic.h"
#include "gpu_type.h"


#define FRAC_VAR (0.1)
#define idloc 0

//================================================================================
__device__ void dE2T(struct Rtype *R, REAL aexp,struct RUNPARAMS *param){

  REAL tloc;
  REAL eint=R->eint;
  REAL nH=R->nh;
  REAL x=R->nhplus/R->nh;
  REAL pstar=param->unit.unit_n*param->unit.unit_d*POW(param->unit.unit_v,2);
  nH=nH/POW(aexp,3)*param->unit.unit_N;
  eint=eint/POW(aexp,5)*pstar;
  tloc=eint/(1.5*nH*KBOLTZ*(1.+x));
  R->temp=tloc;
}


// ============================================================================
__device__ REAL dcucompute_alpha_b(REAL temp, REAL unit_number, REAL aexp)
{
  // CASE B recombination rate m**3 s*-1
  // temperature should be given in Kelvin
  
  REAL alpha_b,lambda;
  lambda=2e0*157807e0/temp;
  alpha_b=2.753e-14*POW(lambda,1.5)/POW(1e0+POW(lambda/2.740,0.407),2.242); //cm3/s
#ifdef TESTCOSMO
  alpha_b=alpha_b*1e-6*unit_number;///(aexp*aexp*aexp); //m3/s
#else
  alpha_b=alpha_b*1e-6*unit_number; //m3/s
#endif
  return alpha_b;
}

//=========================================================
//=========================================================

__device__ REAL dcucompute_alpha_a(REAL temp, REAL unit_number, REAL aexp)
{
  // CASE A recombination rate m**3 s*-1
  // temperature should be given in Kelvin
  
  REAL alpha_a,lambda;
  lambda=2e0*157807e0/temp;
  alpha_a=1.269e-13*POW(lambda,1.503)/POW(1e0+POW(lambda/0.522,0.470),1.923); //cm3/s
#ifdef TESTCOSMO
  alpha_a=alpha_a*1e-6*unit_number;///(aexp*aexp*aexp); //m3/s
#else
  alpha_a=alpha_a*1e-6*unit_number; //m3/s
#endif
  return alpha_a;
}

//=========================================================
//=========================================================

__device__ REAL dcucompute_beta(REAL temp, REAL unit_number, REAL aexp)
{
  // Collizional ionization rate m**3 s*-1
  // temperature in Kelvin
  REAL beta,T5;
  T5=temp/1e5;
  beta=5.85e-11*SQRT(temp)/(1+SQRT(T5))*EXP(-(157809e0/temp)); //cm3/s
#ifdef TESTCOSMO
  beta=beta*1e-6*unit_number;///(aexp*aexp*aexp); // !m3/s
#else
  beta=beta*1e-6*unit_number; // !m3/s
#endif
  return beta;
}

//**********************************************************************************
//**********************************************************************************
__device__ void dcuCompCooling(REAL temp, REAL x, REAL nH, REAL *lambda, REAL *tcool, REAL aexp,REAL CLUMPF)
{

  REAL c1,c2,c3,c4,c5,c6;
  REAL unsurtc;
  REAL nh2;

  nh2=nH*1e-6;// ! m-3 ==> cm-3
  

  // Collisional Ionization Cooling

  c1=EXP(-157809.1e0/temp)*1.27e-21*SQRT(temp)/(1.+SQRT(temp/1e5))*x*(1.-x)*nh2*nh2*CLUMPF;
  

  // Case A Recombination Cooling

  c2=1.778e-29*temp*POW(2e0*157807e0/temp,1.965e0)/POW(1.+POW(2e0*157807e0/temp/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2*CLUMPF;
  
  
  // Case B Recombination Cooling

  c6=3.435e-30*temp*POW(2e0*157807e0/temp,1.970e0)/POW(1.+(POW(2e0*157807e0/temp/2.250e0,0.376e0)),3.720e0)*x*x*nh2*nh2*CLUMPF;
  c6=0.;

  // Collisional excitation cooling

  c3=EXP(-118348e0/temp)*7.5e-19/(1+SQRT(temp/1e5))*x*(1.-x)*nh2*nh2*CLUMPF;
  
  
  // Bremmsstrahlung

  c4=1.42e-27*1.5e0*SQRT(temp)*x*x*nh2*nh2*CLUMPF;
  
  // Compton Cooling
  
  /* c5=1.017e-37*POW(2.727/aexp,4)*(temp-2.727/aexp)*nh2*x; */
   c5=0.;
#ifndef WRADTEST
   c5=5.406e-24*(temp-2.727/aexp)/POW(aexp/0.001,4)*x*nh2;
  REAL Ta=2.727/aexp; c5=5.406e-36*(temp-Ta)/(aexp*aexp*aexp*aexp)*x*nh2;
#endif
  // Overall Cooling
  
  *lambda=c1+c2+c3+c4+c5+c6;// ! erg*cm-3*s-1
  

  // Unit Conversion

  *lambda=(*lambda)*1e-7*1e6;// ! J*m-3*s-1

  // cooling times

  unsurtc=FMAX(c1,c2);
  unsurtc=FMAX(unsurtc,c3);
  unsurtc=FMAX(unsurtc,c4);
  unsurtc=FMAX(unsurtc,FABS(c5));
  unsurtc=FMAX(unsurtc,c6)*1e-7;// ==> J/cm3/s

  *tcool=1.5e0*nh2*(1.+x)*KBOLTZ*temp/unsurtc; //Myr
}


// ===========================================================================================================================

__global__ void dchemrad(struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew, struct RUNPARAMS *param, REAL aexporg, int chemonly)
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
    alpha,
    alphab,
    beta,
    tloc,
    xt,
    eintt,
    ai_tmp1=0.,
    et[NGRP],
    p[NGRP];

  REAL aexp;
  REAL ebkg[NGRP];
  REAL z=1./aexporg-1.;

  REAL c=param->clightorg*LIGHT_SPEED_IN_M_PER_S; 			// switch back to physical velocity m/s

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


#ifdef S_X
  REAL E0overI[NGRP];
  REAL N2[NGRP];
  REAL F2[NGRP];
#endif

#define BLOCKCOOL 1 // KEPT FROM CUDATON FOR SIMPLICITY
#define idloc3 0 // KEPT FROM CUDATON FOR SIMPLICITY

  REAL
    egyloc[BLOCKCOOL*NGRP],
    floc[3*BLOCKCOOL*NGRP],
    srcloc[BLOCKCOOL*NGRP],
    x0[BLOCKCOOL],
    nH[BLOCKCOOL],
    eint[BLOCKCOOL];


  REAL dt=dtnew*param->unit.unit_t*POW(aexporg,2);

  REAL emin;
  struct Rtype R;
  REAL fudgecool=param->fudgecool;
  int ncvgcool=param->ncvgcool;
  REAL E0;
#ifdef SCHAYE
  REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS*MOLECULAR_MU)*param->unit.unit_d;
#endif
  REAL xorg;

  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  i=bx*blockDim.x+tx;


  if(i<nread){
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
      REAL etorg=egyloc[idloc];
      REAL torg=eint[idloc]/(1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc]));

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


	  //if(eint[idloc]!=E0) printf("2!\n");
	  tloc=eint[idloc]/(1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc]));

	  //== Getting a timestep
	  dcuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
	  ai_tmp1=0.;

	  //if(eint[idloc]!=E0) printf("3!\n");

	  if(fudgecool<1e-20){
	    printf("eint=%e(%e<%e) nH=%e x0=%e(%e) T=%e(%e) N=%e(%e)\n",eint[idloc],eorg,emin,nH[idloc],x0[idloc],xorg,tloc,torg,et[0],etorg);
	    //if(fudgecool<1e-20) abort();
	  }

	  for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += ((alphae[igrp])*hnu[igrp]-(alphai[igrp])*hnu0)*egyloc[idloc+igrp*BLOCKCOOL];

	  tcool=FABS(eint[idloc]/(nH[idloc]*(1.0-x0[idloc])*ai_tmp1*(!chemonly)-Cool));
	  ai_tmp1=0.;
	  dtcool=FMIN(fudgecool*tcool,dt-currentcool_t);

	  alpha=dcucompute_alpha_a(tloc,1.,1.)*CLUMPF2;
	  alphab=dcucompute_alpha_b(tloc,1.,1.)*CLUMPF2;
	  beta=dcucompute_beta(tloc,1.,1.)*CLUMPF2;

	  //== Update

	  // ABSORPTION
	  int test = 0;
	  REAL factotsa[NGRP];
	  for(igrp=0;igrp<NGRP;igrp++)
	      {
#ifdef OTSA
		factotsa[igrp]=0;
		alpha=alphab; // recombination is limited to non ground state levels
#else
		factotsa[igrp]=(igrp==0);
#endif

		ai_tmp1 = alphai[igrp];
		if(chemonly){
		  et[igrp]=egyloc[idloc+igrp*BLOCKCOOL];
		}
		else{
		  et[igrp]=((alpha-alphab)*x0[idloc]*x0[idloc]*nH[idloc]*nH[idloc]*dtcool*factotsa[igrp]+egyloc[idloc+igrp*BLOCKCOOL]+srcloc[idloc+igrp*BLOCKCOOL]*dtcool*factgrp[igrp])/(1.+dtcool*(ai_tmp1*(1.-x0[idloc])*nH[idloc]));
		}

		if((et[igrp]<0)||(isnan(et[igrp]))){
		  test=1;
		  //printf("eint=%e nH=%e x0=%e T=%e N=%e\n",eint[idloc],nH[idloc],x0[idloc],tloc,et[0]);
		}
		p[igrp]=(1.+(alphai[igrp]*nH[idloc]*(1-x0[idloc]))*dtcool);
	      }

	  ai_tmp1=0.;


	  if(test)
	    {
	      fudgecool=fudgecool/10.;
	      continue;
	    }

	  // IONISATION
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

	  xt=1.-(alpha*x0[idloc]*x0[idloc]*nH[idloc]*dtcool+(1. -x0[idloc]))/(1.+dtcool*(beta*x0[idloc]*nH[idloc]+ai_tmp1));
	  ai_tmp1=0.;


	  if(((xt>1.)||(xt<0.))||(isnan(xt)))
 	    {
	      //printf("XION ERR eintt=%e xt=%e et=%e\n",eintt,xt,et[0]);
	      fudgecool/=10.;
	      continue;
	    }

#ifdef SEMI_IMPLICIT
	  dcuCompCooling(tloc,xt,nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
#else
	  dcuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
#endif

#ifdef COOLING
	  // HEATING + COOLING

	  int compcool=1; // do we need to compute the cooling ?

#ifdef SCHAYE
	  if((nH[idloc]>1e6)&&(R.nh>(param->stars->overdensity_cond*navg))){
	    REAL tlocs;
	    tlocs=eintt/(1.5*nH[idloc]*KBOLTZ*(1.+xt));
	    if(tlocs<1e5){
	      eintt=(1.08e9*KBOLTZ)*POW(nH[idloc]/1e5,4./3.)/(GAMMA-1)/FSCHAYE; // polytropic EOS
	      compcool=0.; // cancel cooling calculation
	      fudgecool=FMIN(fudgecool*1.5,param->fudgecool);
	    }
	  }
#endif // SCHAYE

	  if(compcool){
	    REAL SN = 0;
#ifdef SUPERNOVAE
	    SN 	 = R.snfb;
	    if (R.snfb) Cool = 0; // Stop the cooling if supernovae
	    if (R.snfb) printf("dE\t%e\tE0\t%e\tdtcool\t%e\t",R.snfb*dtcool,eintt, dtcool);
#endif

#ifndef S_X
#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*(!chemonly);}
	  eintt=(eint[idloc]+ dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool+SN));
// 	  if (R.snfb) printf("E0\t%e\n",eintt);
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-x0[idloc])*(ai_tmp1)-Cool+SN));
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

#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp]*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool+SN));
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp]*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-x0[idloc])*(ai_tmp1)-Cool+SN));
#endif
	  //================================================================================
#endif //S_X

	  if(eintt<0.)
 	    {
	      //printf("E NEG eintt=%e xt=%e et=%e\n",eintt,xt,et[0]);
	      fudgecool=fudgecool/10.;
	      continue;
	    }

	  if(FABS(eintt-eint[idloc])>FRAC_VAR*eint[idloc])
	    {
	      //	      if(srcloc[idloc]==0.){
	      //printf("DELTA E eintt=%e xt=%e et=%e\n",eintt,xt,et[0]);
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
		floc[0+idloc3+igrp*BLOCKCOOL*3]=floc[0+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4);
		floc[1+idloc3+igrp*BLOCKCOOL*3]=floc[1+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4);
		floc[2+idloc3+igrp*BLOCKCOOL*3]=floc[2+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4);
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
	    R.e[igrp]=FMAX(egyloc[idloc+igrp*BLOCKCOOL]*POW(aexp,3),EMIN*factgrp[igrp])/param->unit.unit_N;
	    R.fx[igrp]=floc[0+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
	    R.fy[igrp]=floc[1+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
	    R.fz[igrp]=floc[2+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
	  }
      }

      R.nhplus=x0[idloc]*R.nh;
      R.eint=eint[idloc]*POW(aexp,5)/param->unit.unit_n/param->unit.unit_d/POW(param->unit.unit_v,2);
      dE2T(&R,aexp,param);
      memcpy(&stencil[i].New.cell[icell].rfieldnew,&R,sizeof(struct Rtype));

    }
  }
}


#if 0
__global__ void dchemrad(struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew, struct RUNPARAMS *param, REAL aexporg, int chemonly)
{
  int i,icell,igrp;
  //int idloc;
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
    et[NGRP],
    p[NGRP];

  REAL aexp;
  REAL ebkg[NGRP];
  REAL z=1./aexporg-1.;

  REAL c=param->clightorg*LIGHT_SPEED_IN_M_PER_S; 			// switch back to physical velocity m/s

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


#ifdef S_X
  REAL E0overI[NGRP];
  REAL N2[NGRP];
  REAL F2[NGRP];
#endif


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
  
  
  REAL dt=dtnew*param->unit.unit_t*POW(aexporg,2);

  REAL emin;
  struct Rtype R;
  REAL fudgecool=param->fudgecool;
  int ncvgcool=param->ncvgcool;
  REAL E0;
#ifdef SCHAYE
  REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS*MOLECULAR_MU)*param->unit.unit_d;
#endif
  REAL xorg;

  unsigned int bx=blockIdx.x;
  unsigned int tx=threadIdx.x;
  i=bx*blockDim.x+tx;
  
  if(i<nread){
    for(icell=0;icell<8;icell++){ // we scan the cells

     
      if(stencil[i].oct[6].cell[icell].split) continue; // we dont treat split cells
      memcpy(&R,&stencil[i].New.cell[icell].rfieldnew,sizeof(struct Rtype));// We get the local physical quantities after transport update
      
#ifdef HOMOSOURCE
      // we override the value with the homogeneous source density
      R.src=param->bkg;
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
      //srcloc[idloc]=(R.src*param->unit.unit_N/param->unit.unit_t/(aexporg*aexporg)+ebkg[0])/POW(aexporg,3); 
      for (igrp=0;igrp<NGRP;igrp++){
      srcloc[idloc+igrp*BLOCKCOOL]=(R.src[igrp]*param->unit.unit_N/param->unit.unit_t/(aexporg*aexporg))/POW(aexporg,3); //phot/s/dv (physique)
      }

      /// BELOW THE FULL EXPRESSION OF E in natural units
      //emin=PMIN/(GAMMA-1.)/POW(aexporg,5)/POW(param->unit.unit_l,3)*param->unit.unit_n*param->unit.unit_d*POW(param->unit.unit_v,2); // physical minimal pressure


      /* if(((isnan(eint[idloc]))||(isnan(x0[idloc])))||(eint[idloc]==0.)){ */
      /* 	printf("start with nans or ZErO egy %e\n",eint[idloc]); */
      /* 	abort(); */
      /* } */

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
	}

      // -------------------------------------------------


      // local cooling loop -------------------------------
      aexp=aexporg;
      fudgecool=param->fudgecool;
      currentcool_t=0.;
      nitcool=0.;
      REAL da;
      //printf("fudge=%e ncv=%d currentcool_t=%e dt=%e\n",fudgecool,ncvgcool,currentcool_t,dt);
#if 1
      while(currentcool_t<dt)
	{

	  // Cosmological Adiabatic expansion effects ==============
#ifdef TESTCOSMO
	  REAL hubblet=param->cosmo->H0*SQRT(param->cosmo->om/aexp+param->cosmo->ov*(aexp*aexp))/aexp*(1e3/(1e6*PARSEC)); // s-1 // SOMETHING TO CHECK HERE
#else
	  REAL hubblet=0.;
#endif
	  tloc=eint[idloc]/(1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc]));

	  //== Getting a timestep
	  dcuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
	  ai_tmp1=0.;

	  for (igrp=0;igrp<NGRP;igrp++) ai_tmp1 += ((alphae[igrp])*hnu[igrp]-(alphai[igrp])*hnu0)*egyloc[idloc+igrp*BLOCKCOOL];
	  
	  tcool=FABS(eint[idloc]/(nH[idloc]*(1.0-x0[idloc])*ai_tmp1-Cool));
	  ai_tmp1=0.;
	  dtcool=FMIN(fudgecool*tcool,dt-currentcool_t);
	  
	  alpha=dcucompute_alpha_a(tloc,1.,1.)*CLUMPF2;
	  alphab=dcucompute_alpha_b(tloc,1.,1.)*CLUMPF2;
	  beta=dcucompute_beta(tloc,1.,1.)*CLUMPF2;
      
	  //== Update
	  
	  // ABSORPTION
	  int test = 0;
	  REAL factotsa[NGRP];
	  for(igrp=0;igrp<NGRP;igrp++)
	    {
#ifdef OTSA
	      factotsa[igrp]=0;
	      alpha=alphab; // recombination is limited to non ground state levels
#else
	      factotsa[igrp]=(igrp==0);
#endif

	      ai_tmp1 = alphai[igrp];
	      if(chemonly){
		et[igrp]=egyloc[idloc+igrp*BLOCKCOOL];
		}
	      else{
		et[igrp]=((alpha-alphab)*x0[idloc]*x0[idloc]*nH[idloc]*nH[idloc]*dtcool*factotsa[igrp]+egyloc[idloc+igrp*BLOCKCOOL]+srcloc[idloc+igrp*BLOCKCOOL]*dtcool*factgrp[igrp])/(1.+dtcool*(ai_tmp1*(1.-x0[idloc])*nH[idloc]));
	      }

	      if((et[igrp]<0)||(isnan(et[igrp]))){
		test=1;
	      }
	      p[igrp]=(1.+(alphai[igrp]*nH[idloc]*(1-x0[idloc]))*dtcool);
	    }
	  
	  ai_tmp1=0.;

	  
	  if(test) 
	    {
	      fudgecool=fudgecool/10.; 
	      continue;	
	    } 

	  
	  // IONISATION
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
	  
	  xt=1.-(alpha*x0[idloc]*x0[idloc]*nH[idloc]*dtcool+(1. -x0[idloc]))/(1.+dtcool*(beta*x0[idloc]*nH[idloc]+ai_tmp1));
	  ai_tmp1=0.;


	  if(((xt>1.)||(xt<0.))||(isnan(xt))) 
 	    {
	      fudgecool/=10.; 
	      continue;	
	    } 

#ifdef SEMI_IMPLICIT
	  dcuCompCooling(tloc,xt,nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
#else
	  dcuCompCooling(tloc,x0[idloc],nH[idloc],&Cool,&tcool1,aexp,CLUMPF2);
#endif

#ifdef COOLING
	  // HEATING
	  
	  int compcool=1; // do we need to compute the cooling ?

#ifdef SCHAYE
	  if((nH[idloc]>1e5)&&(R.nh>(57.7*navg))){
	    REAL tlocs;
	    tlocs=eintt/(1.5*nH[idloc]*KBOLTZ*(1.+xt));
	    if(tlocs<1e5){
	      eintt=(1.08e9*KBOLTZ)*POW(nH[idloc]/1e5,4./3.)/(GAMMA-1); // polytropic EOS
	      compcool=0.; // cancel cooling calculation
	      fudgecool=FMIN(fudgecool*1.5,param->fudgecool);
	    }
	  }
#endif

	  if(compcool){ 
	  REAL SN=0.;

#ifndef S_X
#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool));
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-x0[idloc])*(ai_tmp1)-Cool));
#endif
#else
	  //===================================== X RAYS ==============================
	  REAL pp2;
	  F2[0]=1.0;
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
#ifdef SEMI_IMPLICIT
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += et[igrp]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp]*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-xt)*(ai_tmp1)-Cool));
#else
	  for(igrp=0;igrp<NGRP;igrp++) {ai_tmp1 += egyloc[idloc+igrp*BLOCKCOOL]*(alphae[igrp]*hnu[igrp]-(alphai[igrp]*hnu0))*F2[igrp]*(!chemonly);}
	  eintt=(eint[idloc]+dtcool*(nH[idloc]*(1.-x0[idloc])*(ai_tmp1)-Cool));
#endif
#endif

	  if(eintt<0.)
 	    {
	      fudgecool=fudgecool/10.;

	      continue;
	    }

	  if(FABS(eintt-eint[idloc])>FRAC_VAR*eint[idloc])
	    {
	      fudgecool=fudgecool/10.;
	      continue;
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
	  REAL da=hubblet*dtcool*aexp;
	  aexp+=da;
#endif

	  for(igrp =0;igrp<NGRP;igrp++)
	    {
	      egyloc[idloc+igrp*BLOCKCOOL]=et[igrp]*POW(aold/aexp,3);
	      if(!chemonly){
		floc[0+idloc3+igrp*BLOCKCOOL*3]=floc[0+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4);
		floc[1+idloc3+igrp*BLOCKCOOL*3]=floc[1+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4);
		floc[2+idloc3+igrp*BLOCKCOOL*3]=floc[2+idloc3+igrp*BLOCKCOOL*3]/p[igrp]*POW(aold/aexp,4);
	      }	
	    }
	  
	  x0[idloc]=xt;
#ifdef COOLING
	  eint[idloc]=eintt*POW(aold/aexp,5);
#endif
	  currentcool_t+=dtcool;
	  fudgecool=param->fudgecool;
	  nitcool++;
	  if((nitcool==ncvgcool)&&(ncvgcool!=0)) break;
	}
#endif
      // ====================== End of the cooling loop

      // FIlling the rad structure to send it back
       if(!chemonly){
	 for(igrp=0;igrp<NGRP;igrp++)
	   {
	     R.e[igrp]=FMAX(egyloc[idloc+igrp*BLOCKCOOL]*aexp*aexp*aexp,EMIN*factgrp[igrp])/param->unit.unit_N;
	     R.fx[igrp]=floc[0+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
	     R.fy[igrp]=floc[1+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
	     R.fz[igrp]=floc[2+idloc3+igrp*BLOCKCOOL*3]*POW(aexp,4)/param->unit.unit_l*param->unit.unit_t/param->unit.unit_N;
	   }
       }
      R.nhplus=x0[idloc]*R.nh;
      R.eint=eint[idloc]*POW(aexp,5)/param->unit.unit_n/param->unit.unit_d/POW(param->unit.unit_v,2);
      dE2T(&R,aexp,param);
      memcpy(&stencil[i].New.cell[icell].rfieldnew,&R,sizeof(struct Rtype));

    }
  }

}
#endif
#endif
#endif
