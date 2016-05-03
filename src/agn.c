
#ifdef AGN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h> 

#ifdef WMPI
#include <mpi.h>
#endif

#include "prototypes.h"
#include "oct.h" //cell2oct
#include "hydro_utils.h" // W2U and U2W
#include "convert.h"
#include "stars.h"
#include "tools.h"
#include "particle.h" //findlastpart



// Added by Bacher MOKHTARI
// ==================================
//BH Properties
struct bhprop{
  REAL mass; 	//BH mass
  REAL lumi; 	//Number of photons emited per second by the BH
  REAL aexp; 	//Scale factor
  REAL MBHdot; 	//Accretion ratio
  REAL gasdens; //Density of gas
  REAL gastemp; //Temperature of gas
  REAL Gaindens;//Dens gained at each time step
  REAL Xion;  	//Ionisation Ratio
  REAL MEddot; 	//Eddington accretion ratio
  int level;	//Level of the cell where the BH is
  REAL x;	//x position of the BH
  REAL y;	//y position of the BH
  REAL z;	//z position of the BH
};
// ==================================


//Added by Bacher MOKHTARI
/////////////////////////////////////////////////////////////////////////////////////////
////////  __  __  ___ ___ ////////  __   __   __   //////// 			 		 ////////
//////// |_  |__|  |   |  //////// |__| |__| |__/  ////////   Bacher MOKHTARI	 ////////
//////// |   |  | _|_  |  //////// |    |  | |  \  ////////			 			 ////////
/////////////////////////////////////////////////////////////////////////////////////////


void conservation(struct Wtype *field, struct RUNPARAMS *param, struct PART *BH, int level, REAL aexp ,REAL drho){
// ----------------------------------------------------------//
// compute the conservation equations after the BH is created
// ----------------------------------------------------------//

	struct Utype U;
	struct Wtype W;

	memcpy(&W,field,sizeof(struct Wtype));

	W2U(&W, &U); //Transform into conservative quantities

	//Conservation of the Mass Density
	U.d    -= drho;

#ifdef WRADHYD
	REAL xion=W.dX/W.d;
	U.dX 	= U.d*xion;
#endif

	//Conservation of Momentum Density
	U.du -= BH->vx * drho;
	U.dv -= BH->vy * drho;
	U.dw -= BH->vz * drho;

	//Conservation of Internal Energy
#ifdef DUAL_E
	U.eint=U.eint*(1.-drho/W.d); // assuming T and x remain constant
#endif

	U2W(&U, &W);

	//Conservation of Total Energy
	getE(&(W));
	W.a=SQRT(GAMMA*W.p/W.d);
	W.p=FMAX(W.p,PMIN);
	memcpy(field,&W,sizeof(struct Wtype));

	if(isnan(U.du)){
	  printf("Density suppressed= %e ,BH x velocity= %e\n",drho,BH->vx);
	}
}

//=================================================================================================
//=================================================================================================

REAL computeAGNfeedback(struct CELL *cell,struct PART *curp,struct RUNPARAMS *param, REAL aexp, int level, REAL dt, struct bhprop *BH){

  REAL epsilonr=0.1; //According to Dubois 2011
  REAL epsilonf=0.15; //According to Dubois
  REAL sigmaT=6.6524E-29; //Thomson cross-section
  REAL alpha;
  REAL threshold=(0.1*PROTON_MASS/1E-6);  //Threshold density in Dubois & al. 2011 in SI (kg/m³)
  
  //Print quantities depending of the density
  BH->gasdens=cell->field.d*param->unit.unit_d/pow(aexp,3);//The ratio is the gas density absorbed / the gas density in the cell BEFORE absorption
  //Computing the Boost factor						
  if((cell->field.d)*param->unit.unit_d/pow(aexp,3)>threshold){
    alpha=pow(cell->field.d*param->unit.unit_d/pow(aexp,3)/threshold,2);
  }
  else{
    alpha=1;
  }
  printf("Boost factor= %f\n",alpha);
  //Computing the Speed of sound and the speed of gas relatively to BH
  REAL VRelgaz2=(pow(cell->field.u-curp->vx,2)+pow(cell->field.v-curp->vy,2)+pow(cell->field.w-curp->vz,2))
    *pow(param->unit.unit_v/aexp,2); 		//In (m/s)² physical
  REAL Cs2=pow(cell->field.a,2)*pow(param->unit.unit_v/aexp,2);	//In (m/s)² physical
  //Dubois' criterion : the speed of gas must'nt be over the average gas velocity dispersion in the ISM (10km/s)
  REAL Vlim=10; //In km/s physical
  if(VRelgaz2>pow(Vlim*1000,2)){
    VRelgaz2=pow(Vlim*1000,2);
  }
  //Computing the Accretion ratio
  REAL MBHdot= (1-epsilonr*epsilonf)
    *(4*M_PI*alpha*pow(NEWTON_G,2)*pow(curp->mass*param->unit.unit_mass,2)*(cell->field.d*param->unit.unit_d/pow(aexp,3)))
    /pow(VRelgaz2+Cs2,1.5); //in kg/s
  //Computing Eddington ratio
  BH->MEddot=(4*M_PI*NEWTON_G*curp->mass*param->unit.unit_mass*PROTON_MASS)/(epsilonr*sigmaT*LIGHT_SPEED_IN_M_PER_S); //In kg/s
  //The accretion ratio must'nt be over the Eddington limit			
  if(MBHdot>BH->MEddot){
    MBHdot=BH->MEddot;
  }
  printf("Accretion ratio in solar mass per yr : %0.9E\n",MBHdot/SOLAR_MASS*3.154E7);
  //The gas accreted mustnt be over 90% of the gas density in the cell			
  if(MBHdot*(pow(aexp,2)*param->unit.unit_t/param->unit.unit_mass)*dt/pow(0.5,3*level)>=0.9*cell->field.d){
    BH->Gaindens=0.9*cell->field.d*param->unit.unit_d/pow(aexp,3); //At time t
    conservation(&(cell->field),param,curp,level,aexp,0.9*cell->field.d);
    curp->mass+=0.9*cell->field.d*pow(0.5,3*level);
  }
  else{
    BH->Gaindens=MBHdot*(pow(aexp,2)*param->unit.unit_t/param->unit.unit_mass)*dt/pow(0.5,3*level)*param->unit.unit_d/pow(aexp,3);
    conservation(&(cell->field),param,curp,level,aexp,MBHdot*(pow(aexp,2)*param->unit.unit_t/param->unit.unit_mass)*dt/pow(0.5,3*level));
    curp->mass+=MBHdot*(pow(aexp,2)*param->unit.unit_t/param->unit.unit_mass)*dt;
  }
  //Computing the Emissivity of the BH in photons/s
  REAL E_FEEDBACK=epsilonr*epsilonf/(1-epsilonf*epsilonr)*MBHdot*pow(LIGHT_SPEED_IN_M_PER_S,2); // J/sec
  REAL Ngammadot=epsilonr*epsilonf/(1-epsilonf*epsilonr)*MBHdot*pow(LIGHT_SPEED_IN_M_PER_S,2)/(curp->mass*param->unit.unit_mass);
  Ngammadot/=(20.26*ELECTRONVOLT); //In photons/s/kg
  printf("Number of photons emited per second per kg :%0.6E Efeedback=%e\n",Ngammadot,E_FEEDBACK);
	    //Print BH Mass in a file
  curp->ekin=E_FEEDBACK;

  
  BH->mass=curp->mass*param->unit.unit_mass/SOLAR_MASS; //After modification at t+1
  BH->lumi=Ngammadot;
  BH->aexp=aexp;
  BH->MBHdot=MBHdot/SOLAR_MASS*3.154E7; //In solar mass/yr
  BH->gastemp=cell->rfield.temp; //At time t+1
  BH->Xion=cell->field.dX/cell->field.d;//The ionisation ratio is the ration of the quantity of the gas ionised / density of the cell at time t+1
  BH->level=level;
  BH->x=curp->x;
  BH->y=curp->y;
  BH->z=curp->z;

  return E_FEEDBACK;

}


// ============================================================================

int agn_grow(struct CELL *cell, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL aexp, int level, REAL dt, struct bhprop *BH){
// ----------------------------------------------------------//
/// Scan all particles of cell "cell" and look for AGN
/// If SN found, compute energy and inject it
// ----------------------------------------------------------//


  static REAL total_sn_egy=0;
  int Nsn = 0;

  struct PART *nexp=cell->phead;

  if(nexp==NULL) return 0;
  do{
    struct PART *curp=nexp;
    nexp=curp->next;

    if (curp->isStar==100){ // if curp is an AGN
      
      REAL E= computeAGNfeedback(cell,curp,param,aexp,level,dt, BH);
      
      
#ifdef AGN_MECH_FEEDBACK
      // TO be implemented later, by inspiration from supernovae
      if(param->agn->feedback_frac){
         /* oct feedback
          * if kinetic feedback
          * the energy is injected in a oct
          */
        thermalFeedbackOct(param,cell,curp,level, E*(1.-param->agn->feedback_frac));

        if (param->agn->ejecta_proportion){
          kineticFeedback_impulsion(param, cell,curp,aexp,level, E*param->sn->feedback_frac);
        }else{
          kineticFeedback_simple(param, cell,curp,aexp,level, E*param->sn->feedback_frac);
        }
      }
#endif
      
      Nsn++;
    }
  }while(nexp!=NULL);
}


// ============================================================================

void initagn(struct Wtype *field, struct PART *agn, struct RUNPARAMS *param, int level, REAL xc, REAL yc, REAL zc,int idx, REAL aexp, int is, REAL dt,REAL dx, REAL magn) {
// ----------------------------------------------------------//
/// compute the initial state of agn particle
// ----------------------------------------------------------//

  // some parameters
	agn->next = NULL;
	agn->idx = -1000;
	agn->level = level;
	agn->is = is;
	agn->isStar = 100; // 100 is for AGNs
	agn->rhocell = field->d;

  // set agn position to cell center
	agn->x = xc ;
	agn->y = yc ;
	agn->z = zc ;

  // set agn velocity to fluid velocity
 	agn->vx = field->u;
	agn->vy = field->v;
	agn->vz = field->w;

#define RDM_STARS
#ifdef RDM_STARS
  // random position
	agn->x += rdm(-0.5,0.5)*dx;
	agn->y += rdm(-0.5,0.5)*dx;
	agn->z += rdm(-0.5,0.5)*dx;

  // compute random velocity component
	REAL r = rdm(0,1) * field->a;
	REAL theta  = acos(rdm(-1,1));
	REAL phi = rdm(0,2*M_PI);

  // add random velocity component
	agn->vx += r * sin(theta) * cos(phi);
	agn->vy += r * sin(theta) * sin(phi);
	agn->vz += r * cos(theta) ;
#endif // RDM_AGNS

  //mass
	agn->mass = magn;

  //energy
	agn->epot = 0.0;
	agn->ekin = 0.5 * agn->mass * (POW(agn->vx,2) + POW(agn->vy,2) + POW(agn->vz,2));

  // age
#ifdef TESTCOSMO
	agn->age = param->cosmo->tphy;
#else
/// TODO fix age of star ifndef TESTCOSMO
#endif
}

// ============================================================================

void addagn(struct CELL *cell, struct Wtype *init_field, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dttilde, struct RUNPARAMS *param, REAL aexp, int is,  int nagn, REAL magn){
// ----------------------------------------------------------//
/// Add a stellar particle in the double linked list\n
/// Call the initialisation and the conservation functions
// ----------------------------------------------------------//

#ifdef PIC
  struct PART *agn = cpu->freepart;
  
  if (agn==NULL){
    printf("\n");
    printf("----------------------------\n");
    printf("No more memory for particles\n");
    printf("----------------------------\n");
    //	exit(0);
  }else{
    
    cpu->freepart = cpu->freepart->next;
    cpu->freepart->prev = NULL;
    
    if (cell->phead==NULL){
      cell->phead = agn;
      agn->prev = NULL;
      agn->next =NULL;
    }else{
      struct PART *lasp = findlastpart(cell->phead);
      lasp->next = agn;
      agn->prev = lasp;
      agn->next = NULL;
    }

    cpu->npart[level-1]++;
    cpu->nstar[level-1]++; // agn counted as star

    REAL dx = POW(2.0,-level);
    
    initagn(init_field, agn, param, level, xc, yc, zc, -1, aexp, is, dttilde, dx,magn);

    // for conserveField we use the same function as stars.c
    conserveField(&cell->field   , param, agn,  dx, aexp,magn);
    conserveField(&cell->fieldnew, param, agn,  dx, aexp,magn);


    printf("AGN mass=%e %d \n",agn->mass,agn->isStar);
    //printPart(star);
  }
#endif
}


 void agn(struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt, REAL aexp, int level,int is){
  // ------------------------------------------------------- //
// Call the agn function for all cells of the grid  //
// ------------------------------------------------------- //

	struct bhprop BH;
	BH.mass=0.;
	BH.lumi=0.;
	BH.aexp=0.;
	BH.MBHdot=0.;
	BH.gasdens=0.;
	BH.gastemp=0.;
	BH.Gaindens=0.;
	BH.Xion=0.;
	BH.MEddot=0.;
	BH.level=0;
	BH.x=0.;
	BH.y=0.;
	BH.z=0.;


	int n_part_agn = 0;
	int n_unit_agn = 0;
	REAL mass_agn=1E5*SOLAR_MASS/param->unit.unit_mass;

	int iOct;
	for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
	  struct OCT *curoct=cpu->octList[level-1][iOct];
	  int icell;
	  for(icell=0;icell<8;icell++) {
	    struct CELL *curcell = &curoct->cell[icell];
	    // if conditions are satisfied we create an agn
	    // we use the same fonction as the one in stars.c

	    if((testCond(curcell, param, aexp, level))&&(param->stars->nagn==0)) {

	      REAL dx = POW(2.0,-level);
	      REAL xc=curoct->x+( icell    & 1)*dx+dx*0.5;
	      REAL yc=curoct->y+((icell>>1)& 1)*dx+dx*0.5;
	      REAL zc=curoct->z+( icell>>2    )*dx+dx*0.5;

	      addagn(curcell, &curcell->field, level, xc, yc, zc, cpu, dt, param, aexp, is,n_part_agn++, mass_agn);
	      n_unit_agn++;
	    }

	    // otherwise we check if an agn exists and make it grow if necessary
	    agn_grow(curcell, param, cpu, aexp, level, dt, &BH);
	  }
	}



#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&n_unit_agn,1,MPI_INT,MPI_SUM,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&n_part_agn,1,MPI_INT,MPI_SUM,cpu->comm);
#endif

	param->stars->n += n_part_agn; // NOTE AGNs are counted as stars
	param->stars->nagn += n_part_agn; // NOTE AGNs are counted as stars

	if(cpu->rank==RANK_DISP) {
	  if (n_unit_agn){
	    printf("%d AGN added on level %d\n",n_unit_agn,level);
	  }
	}




	double dummy;
	dummy=BH.mass;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.mass=dummy;

	dummy=BH.lumi;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.lumi=dummy;

	dummy=BH.aexp;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.aexp=dummy;

	dummy=BH.MBHdot;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.MBHdot=dummy;

	dummy=BH.gasdens;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.gasdens=dummy;

	dummy=BH.gastemp;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.gastemp=dummy;

	dummy=BH.Gaindens;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.Gaindens=dummy;

	dummy=BH.Xion;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.Xion=dummy;

	dummy=BH.level;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.level=dummy;
	
	dummy=BH.x;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.x=dummy;

	dummy=BH.y;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.y=dummy;

	dummy=BH.z;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.z=dummy;

	dummy=BH.MEddot;
	MPI_Allreduce(MPI_IN_PLACE,&dummy,1,MPI_DOUBLE,MPI_MAX,cpu->comm);
	BH.MEddot=dummy;

	if(cpu->rank==0){
		if(BH.aexp>0){
			FILE *fp;
			fp=fopen("BH_mass_1E5_Vlim_10_NotEdd.txt","a+");
			//printf("%0.10E %0.6E %0.6E\n",1./BH.aexp-1,BH.mass,BH.lumi); 
			fprintf(fp,
				"%0.10E %0.6E %0.6E %0.6E %0.6E %0.6E %0.6E %0.6E %0.6E %i %0.6E %0.6E %0.6E\n",
				1./BH.aexp-1,
				BH.mass,
				BH.lumi,
				BH.MBHdot,
				BH.MEddot,
				BH.gasdens,
				BH.gastemp,
				BH.Gaindens/BH.gasdens,
				BH.Xion,
				BH.level,
				BH.x,
				BH.y,
				BH.z
			);
			fclose(fp);
		}
	}
}

	/////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// .__ .  . .__  /////////////////////////////////////
	///////////////////////////////////// |__ |\ | |  \ /////////////////////////////////////
	///////////////////////////////////// |__ | \| |__/ ///////////////////////////////////// 
	/////////////////////////////////////////////////////////////////////////////////////////

#endif // AGN
