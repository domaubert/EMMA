#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "friedmann.h"
#include "segment.h"
#include <string.h>
#include "stars.h"


void cell2lcell(struct CELL *cell, struct LCELL *lcell){

  lcell->marked=cell->marked;
  lcell->child=(cell->child!=NULL);
#ifdef PIC
  lcell->density=cell->density;
#endif

#ifdef WGRAV
  lcell->den=cell->gdata.d;
  lcell->pot=cell->gdata.p;
  lcell->res=cell->res;
 
  lcell->f[0]=cell->f[0];
  lcell->f[1]=cell->f[1];
  lcell->f[2]=cell->f[2];
#endif

#ifdef WHYDRO2
  lcell->d=cell->field.d;
  lcell->u=cell->field.u;
  lcell->v=cell->field.v;
  lcell->w=cell->field.w;
  lcell->p=cell->field.p;
#ifdef WRADHYD
  lcell->dX=cell->field.dX;
#endif
#endif

#ifdef WRAD
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    lcell->e[igrp]=cell->rfield.e[igrp];
    lcell->fx[igrp]=cell->rfield.fx[igrp];
    lcell->fy[igrp]=cell->rfield.fy[igrp];
    lcell->fz[igrp]=cell->rfield.fz[igrp];
  }
  lcell->src=cell->rfield.src;
#ifdef STARS
  lcell->snfb=cell->rfield.snfb;
#endif
  lcell->xion=cell->rfield.nhplus/cell->rfield.nh;
  lcell->temp=cell->rfield.temp;
#endif

}

void oct2loct(struct OCT *oct, struct LOCT *loct){
  int icell;

  for(icell=0;icell<8;icell++){ 
    cell2lcell(&oct->cell[icell],&loct->cell[icell]);
    //    memcpy(&loct->cell[icell],&oct->cell[icell],sizeof(struct CELL)); 
  } 

  loct->x=oct->x; 
  loct->y=oct->y; 
  loct->z=oct->z; 

  loct->cpu=oct->cpu;
  loct->level=oct->level;
}


void dumpHeaderOnScreen(struct RUNPARAMS *param, struct CPUINFO *cpu){


printf( "SINGLEPRECISION\t");
#ifdef SINGLEPRECISION
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "PIC\t");
#ifdef PIC
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WHYDRO2\t");
#ifdef WHYDRO2
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WGRAV\t");
#ifdef WGRAV 
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WRAD\t");
#ifdef WRAD
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WRADHYD\t");
#ifdef WRADHYD
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TESTCOSMO\t");
#ifdef TESTCOSMO
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WDBG\t");
#ifdef WDBG
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "STARS\t");
#ifdef STARS
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WMPI\t");
#ifdef WMPI
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "FLOORDT\t");
#ifdef FLOORDT
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WCUDA_ERR\t");
#ifdef WCUDA_ERR
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "NOCOMP\t");
#ifdef NOCOMP
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "GRAFIC\t");
#ifdef GRAFIC
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "ZELDOVICH\t");
#ifdef ZELDOVICH
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "EVRARD\t");
#ifdef EVRARD
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "EDBERT\t");
#ifdef EDBERT
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TUBE\t");
#ifdef TUBE
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "PARTN\t");
#ifdef PARTN
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "PART2\t");
#ifdef PART2
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WRADTEST\t");
#ifdef WRADTEST  
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TESTCLUMP\t");
#ifdef TESTCLUMP 
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "PART_EGY\t");
#ifdef PART_EGY
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "PERFECT\t");
#ifdef PERFECT
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "FASTGRAV\t");
#ifdef FASTGRAV 
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "ONFLYRED\t");
#ifdef ONFLYRED
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "RIEMANN_HLLC\t");
#ifdef RIEMANN_HLLC
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "RIEMANN_EXACT\t");
#ifdef RIEMANN_EXACT
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "PRIMITIVE\t");
#ifdef PRIMITIVE
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "DUAL_E\t");
#ifdef DUAL_E
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "WCHEM\t");
#ifdef WCHEM 
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "S_100000\t");
#ifdef S_100000
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "COOLING\t");
#ifdef COOLING
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "ACCEL_RAD_STAR\t");
#ifdef ACCEL_RAD_STAR
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "SCHAYE\t");
#ifdef SCHAYE
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "OTSA\t");
#ifdef OTSA
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "COARSERAD\t");
#ifdef COARSERAD
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "UVBKG\t");
#ifdef UVBKG
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TRANSZM\t");
#ifdef TRANSZM
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TRANSZP\t");
#ifdef TRANSZP
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TRANSYM\t");
#ifdef TRANSYM
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TRANSYP\t");
#ifdef TRANSYP
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TRANSXM\t");
#ifdef TRANSXM
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "TRANSXP\t");
#ifdef TRANSXP
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "REFXM\t");
#ifdef REFXM 
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "REFYM\t");
#ifdef REFYM 
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif

printf( "REFZM\t");
#ifdef REFZM
  printf( "%d\n", 1);
#else
  printf( "%d\n", 0);
#endif


  printf( "npartmax\t%d\n",(param->npartmax)); 			// the max particles number (per process)
  printf( "ngridmax\t%d\n",(param->ngridmax) ); 		// the max oct numbers (per process)
  printf( "nbuff\t%d\n",(param->nbuff) ); 			// the mpi buffer size
  printf( "ndumps\t%d\n",(param->ndumps) ); 			// the frequency of outputs
  printf( "nsteps\t%d\n",(param->nsteps) ); 			// the maximal number of timesteps

  printf( "lcoarse\t%d\n",(param->lcoarse) ); 			// the coarse level
  printf( "lmax\t%d\n",(param->lmax) ); 			// the max level of refinement

  printf( "niter\t%d\n",(param->niter) ); 			// the maximal number of iterations for the Poisson solver
  		
  printf( "gstride\t%d\n",(param->gstride) ); 			// the size of the stencil for vector based computations (gravity)
  printf( "hstride\t%d\n",(param->hstride) ); 			// the size of the stencil for vector based computations (hydro)

  printf( "dt\t%e\n",(param->dt) ); 				// the timsestep
  printf( "tmax\t%e\n",(param->tmax) ); 			// the simulation stops at tmax : corresponds to amax in cosmo
  printf( "time_max\t%e\n",(param->time_max) ); 		// for cosmo only : contains the time equivalent to amax (contained in tmax, yeah its obfuscated)

  printf( "maxhash\t%d\n",(param->maxhash) ); 			// the hash table size between hilbert keys and oct adress (should be typically = to (2^levelmax-1)^3
  
  printf( "amrthresh\t%e\n",(param->amrthresh) ); 		// the refinement criterion (refine if mcell>amrthresh)
  printf( "nsmooth\t%d\n",(param->nsmooth) ); 			// the number of neighbour refinement steps

  printf( "poissonacc\t%e\n",(param->poissonacc) ); 		// relaxation accuracy for Poisson equation
  printf( "mgridlmin\t%d\n",(param->mgridlmin) ); 		// coarsest level for multigrid relaxation
  printf( "nvcycles\t%d\n",(param->nvcycles) ); 		// number of vcycles for multigrid relaxation
  printf( "nrelax\t%d\n",(param->nrelax) ); 			// number of smoothing cycles
	
  printf( "nrestart\t%d\n",(param->nrestart) ); 		// the restart snapshot
  printf( "nsubcycles\t%d\n",(param->nsubcycles) ); 		// number of subcyles in AMR advance procedure

  printf( "nthread\t%d\n",(param->nthread) );			// number of thread
  printf( "nstream\t%d\n",(param->nstream) );			// number of stream

  printf( "egy_rhs\t%e\n",(param->egy_rhs) ); 			// the right hand side of the energy conservation equation (0 in non cosmological case);
  printf( "egy_0\t%e\n",(param->egy_0) ); 			// the initial energy
  printf( "egy_last\t%e\n",(param->egy_last) ); 		// the last integrand for the energy equation (used for trapezoidal rule)
  printf( "egy_timelast\t%e\n",(param->egy_timelast) ); 	// the last time for the integrand (used for trapezoidal rule)
  printf( "egy_totlast\t%e\n",(param->egy_totlast) ); 

#ifdef WRAD
  printf( "unit_l\t%e\n",(param->unit.unit_l) );		// comoving length size of the box [meters]
  printf( "unit_v\t%e\n",(param->unit.unit_v) );		// unit velocity
  printf( "unit_t\t%e\n",(param->unit.unit_t) );		// unit time [seconds]
  printf( "unit_\t%e\n",(param->unit.unit_n) );			// unit number [moles typically]
  printf( "unit_mass\t%e\n",(param->unit.unit_mass) );		// unit mass [in kg, total mass is equal to one in unit codes]

  printf( "clight\t%e\n",(param->clightorg) ); 			// speed of light in units of the real one
  printf( "fudgecool\t%e\n",(param->fudgecool) ); 		// cooling fraction
  printf( "ncvgcool\t%d\n",(param->ncvgcool) ); 		// cooling max iterations
  
  printf( "denthresh\t%e\n",(param->denthresh) );		// density threshold to turn the sources on
  printf( "tmpthresh\t%e\n",(param->tmpthresh) );		// temperature threshold to turn the sources on
  printf( "srcint\t%e\n",(param->srcint) );			// intensity of the sources
#endif

#ifdef TESTCOSMO
  printf( "om\t%e\n",(param->cosmo->om) );			// Omega matter
  printf( "ov\t%e\n",(param->cosmo->ov) );			// Omega vacuum
  printf( "ob\t%e\n",(param->cosmo->ob) );			// Omega baryon
  printf( "H0\t%e\n",(param->cosmo->H0) );			// Hubble constant
#endif

#ifdef STARS 
  printf( "overdensity_cond\t%e\n",(param->stars->overdensity_cond) );	// need overdensity_cond times the mean density to begin star formation
  printf( "density_cond\t%e\n",(param->stars->density_cond) );		// minimum Hydrogen density [cm-3]
  printf( "tcar\t%e\n",(param->stars->tcar) );				// carateristic time [yr]
  printf( "tlife\t%e\n",(param->stars->tlife) );			// radiative life time of a stellar particle [yr]
  printf( "feedback_eff\t%e\n",(param->stars->feedback_eff) );		// SN feedback efficiency
  printf( "feedback_frac\t%e\n",(param->stars->feedback_frac) );	// fraction of thermal feedback (the other part goes to kinetic feedback) 

#endif

}

void dumpHeader(struct RUNPARAMS *param, struct CPUINFO *cpu){
  printf("Dumping parameters file\n");

  FILE *fp; 
  fp=fopen("data/param.00000.p00000","w");

  int t = 1;
  int f = 0;


fprintf(fp, "PIC\t");
#ifdef PIC
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WHYDRO2\t");
#ifdef WHYDRO2
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WGRAV\t");
#ifdef WGRAV 
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WRAD\t");
#ifdef WRAD
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WRADHYD\t");
#ifdef WRADHYD
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TESTCOSMO\t");
#ifdef TESTCOSMO
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WDBG\t");
#ifdef WDBG
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "STARS\t");
#ifdef STARS
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WMPI\t");
#ifdef WMPI
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "FLOORDT\t");
#ifdef FLOORDT
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WCUDA_ERR\t");
#ifdef WCUDA_ERR
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "NOCOMP\t");
#ifdef NOCOMP
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "GRAFIC\t");
#ifdef GRAFIC
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "ZELDOVICH\t");
#ifdef ZELDOVICH
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "EVRARD\t");
#ifdef EVRARD
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "EDBERT\t");
#ifdef EDBERT
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TUBE\t");
#ifdef TUBE
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "PARTN\t");
#ifdef PARTN
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "PART2\t");
#ifdef PART2
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WRADTEST\t");
#ifdef WRADTEST  
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TESTCLUMP\t");
#ifdef TESTCLUMP 
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "PART_EGY\t");
#ifdef PART_EGY
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "PERFECT\t");
#ifdef PERFECT
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "FASTGRAV\t");
#ifdef FASTGRAV 
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "ONFLYRED\t");
#ifdef ONFLYRED
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "RIEMANN_HLLC\t");
#ifdef RIEMANN_HLLC
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "RIEMANN_EXACT\t");
#ifdef RIEMANN_EXACT
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "PRIMITIVE\t");
#ifdef PRIMITIVE
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "DUAL_E\t");
#ifdef DUAL_E
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "WCHEM\t");
#ifdef WCHEM 
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "S_100000\t");
#ifdef S_100000
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "COOLING\t");
#ifdef COOLING
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "UVBKG\t");
#ifdef UVBKG
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TRANSZM\t");
#ifdef TRANSZM
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TRANSZP\t");
#ifdef TRANSZP
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TRANSYM\t");
#ifdef TRANSYM
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TRANSYP\t");
#ifdef TRANSYP
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TRANSXM\t");
#ifdef TRANSXM
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "TRANSXP\t");
#ifdef TRANSXP
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "REFXM\t");
#ifdef REFXM 
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "REFYM\t");
#ifdef REFYM 
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif

fprintf(fp, "REFZM\t");
#ifdef REFZM
  fprintf(fp, "%d\n", 1);
#else
  fprintf(fp, "%d\n", 0);
#endif


  fprintf(fp, "npartmax\t%d\n",(param->npartmax)); 		// the max particles number (per process)
  fprintf(fp, "ngridmax\t%d\n",(param->ngridmax) ); 		// the max oct numbers (per process)
  fprintf(fp, "nbuff\t%d\n",(param->nbuff) ); 			// the mpi buffer size
  fprintf(fp, "ndumps\t%d\n",(param->ndumps) ); 		// the frequency of outputs
  fprintf(fp, "nsteps\t%d\n",(param->nsteps) ); 		// the maximal number of timesteps

  fprintf(fp, "lcoarse\t%d\n",(param->lcoarse) ); 		// the coarse level
  fprintf(fp, "lmax\t%d\n",(param->lmax) ); 			// the max level of refinement

  fprintf(fp, "niter\t%d\n",(param->niter) ); 			// the maximal number of iterations for the Poisson solver
  
  fprintf(fp, "gstride\t%d\n",(param->gstride) ); 		// the size of the stencil for vector based computations (gravity)
  fprintf(fp, "hstride\t%d\n",(param->hstride) ); 		// the size of the stencil for vector based computations (hydro)

  fprintf(fp, "dt\t%e\n",(param->dt) ); 			// the timsestep
  fprintf(fp, "tmax\t%e\n",(param->tmax) ); 			// the simulation stops at tmax : corresponds to amax in cosmo
  fprintf(fp, "time_max\t%e\n",(param->time_max) ); 		// for cosmo only : contains the time equivalent to amax (contained in tmax, yeah its obfuscated)

  fprintf(fp, "maxhash\t%d\n",(param->maxhash) ); 		// the hash table size between hilbert keys and oct adress (should be typically = to (2^levelmax-1)^3
  
  fprintf(fp, "amrthresh\t%e\n",(param->amrthresh) ); 		// the refinement criterion (refine if mcell>amrthresh)
  fprintf(fp, "nsmooth\t%d\n",(param->nsmooth) ); 		// the number of neighbour refinement steps

  fprintf(fp, "poissonacc\t%e\n",(param->poissonacc) ); 	// relaxation accuracy for Poisson equation
  fprintf(fp, "mgridlmin\t%d\n",(param->mgridlmin) ); 		// coarsest level for multigrid relaxation
  fprintf(fp, "nvcycles\t%d\n",(param->nvcycles) ); 		// number of vcycles for multigrid relaxation
  fprintf(fp, "nrelax\t%d\n",(param->nrelax) ); 		// number of smoothing cycles

  fprintf(fp, "nrestart\t%d\n",(param->nrestart) ); 		// the restart snapshot
  fprintf(fp, "nsubcycles\t%d\n",(param->nsubcycles) ); 	// number of subcyles in AMR advance procedure

  fprintf(fp, "nthread\t%d\n",(param->nthread) );		// number of thread
  fprintf(fp, "nstream\t%d\n",(param->nstream) );		// number of stream

  fprintf(fp, "egy_rhs\t%e\n",(param->egy_rhs) ); 		// the right hand side of the energy conservation equation (0 in non cosmological case);
  fprintf(fp, "egy_0\t%e\n",(param->egy_0) ); 			// the initial energy
  fprintf(fp, "egy_last\t%e\n",(param->egy_last) ); 		// the last integrand for the energy equation (used for trapezoidal rule)
  fprintf(fp, "egy_timelast\t%e\n",(param->egy_timelast) ); 	// the last time for the integrand (used for trapezoidal rule)
  fprintf(fp, "egy_totlast\t%e\n",(param->egy_totlast) ); 

//  printf("%lf\n",param->egy_0);

#ifdef WRAD
  fprintf(fp, "unit_l\t%e\n",(param->unit.unit_l) );		// comoving length size of the box [meters]
  fprintf(fp, "unit_v\t%e\n",(param->unit.unit_v) );		// unit velocity
  fprintf(fp, "unit_t\t%e\n",(param->unit.unit_t) );		// unit time [seconds]
  fprintf(fp, "unit_\t%e\n",(param->unit.unit_n) );		// unit number [moles typically]
  fprintf(fp, "unit_mass\t%e\n",(param->unit.unit_mass) );	// unit mass [in kg, total mass is equal to one in unit codes]

  fprintf(fp, "clight\t%e\n",(param->clightorg) ); 		// speed of light in units of the real one
  fprintf(fp, "fudgecool\t%e\n",(param->fudgecool) ); 		// cooling fraction
  fprintf(fp, "ncvgcool\t%d\n",(param->ncvgcool) ); 		// cooling max iterations
  
  fprintf(fp, "denthresh\t%e\n",(param->denthresh) );		// density threshold to turn the sources on
  fprintf(fp, "tmpthresh\t%e\n",(param->tmpthresh) );		// temperature threshold to turn the sources on
  fprintf(fp, "srcint\t%e\n",(param->srcint) );			// intensity of the sources
#endif

#ifdef TESTCOSMO
  fprintf(fp, "om\t%e\n",(param->cosmo->om) );			// Omega matter
  fprintf(fp, "ov\t%e\n",(param->cosmo->ov) );			// Omega vacuum
  fprintf(fp, "ob\t%e\n",(param->cosmo->ob) );			// Omega baryon
  fprintf(fp, "H0\t%e\n",(param->cosmo->H0) );			// Hubble constant
#endif

#ifdef STARS 
  fprintf(fp, "overdensity_cond\t%e\n",(param->stars->overdensity_cond) );	// need overdensity_cond times the mean density to begin star formation
  fprintf(fp, "density_cond\t%e\n",(param->stars->density_cond) );		// minimum Hydrogen density [cm-3]
  fprintf(fp, "tcar\t%e\n",(param->stars->tcar) );				// carateristic time [yr]
  fprintf(fp, "tlife\t%e\n",(param->stars->tlife) );				// radiative life time of a stellar particle [yr]
  fprintf(fp, "feedback_eff\t%e\n",(param->stars->feedback_eff) );		// SN feedback efficiency
  fprintf(fp, "feedback_frac\t%e\n",(param->stars->feedback_frac) );		// fraction of thermal feedback (the other part goes to kinetic feedback) 

#endif

  fclose(fp);


  dumpHeaderOnScreen(param,cpu);
  //abort();
}

//====================================================================================================
void save_amr(char filename[], struct OCT **firstoct,REAL tsim, REAL tinit,int nsteps, int ndumps, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt){
  
  FILE *fp;
  int level;
  struct OCT * nextoct;
  struct OCT * root;
  struct CELL * rootcell;
  struct PART * rootpart;
  struct OCT *oct;
  int ic=0;


  root=firstoct[0];
  rootcell=&(root->cell[0]);
  rootpart=proot;

  fp=fopen(filename,"wb");
  
  fwrite(&tsim,sizeof(REAL),1,fp); 
  fwrite(&tinit,sizeof(REAL),1,fp); 
  fwrite(&nsteps,sizeof(int),1,fp); 
  fwrite(&ndumps,sizeof(int),1,fp); 
  fwrite(&(param->lcoarse),sizeof(int),1,fp);
  fwrite(&(param->lmax),sizeof(int),1,fp);

  // dumping timesteps
  for(level=1;level<=param->lmax;level++){
    fwrite(adt+level-1,sizeof(REAL),1,fp);
  }
  
#ifdef WRAD
  fwrite(&(param->unit.unit_l),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_v),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_t),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_n),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_mass),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_d),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_N),sizeof(REAL),1,fp);
#endif


#ifdef TESTCOSMO
  fwrite(&(param->cosmo->om),sizeof(REAL),1,fp);
  fwrite(&(param->cosmo->ov),sizeof(REAL),1,fp);
  fwrite(&(param->cosmo->ob),sizeof(REAL),1,fp);
  fwrite(&(param->cosmo->H0),sizeof(REAL),1,fp);
#endif



  // writing pointer informations
  fwrite(&root,sizeof(struct OCT*),1,fp);
  fwrite(&rootcell,sizeof(struct CELL*),1,fp);
  fwrite(&rootpart,sizeof(struct PART*),1,fp);

  
  for(level=1;level<=param->lmax;level++) // looping over octs
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];

      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=nextoct;
	  nextoct=oct->next;
	  
	  fwrite(&oct,sizeof(struct OCT*),1,fp);
	  fwrite(oct,sizeof(struct OCT),1,fp);
	  ic++;
	}while(nextoct!=NULL);
    }

  
  //printf("%d octs dumped by proc %d\n",ic,cpu->rank);


  fclose(fp);

}



//====================================================================================================
struct OCT * restore_amr(char filename[], struct OCT **firstoct,struct OCT **lastoct, REAL *tsim, REAL *tinit, int *nsteps, int *ndumps,struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt, struct CELL *root){
  
  FILE *fp;
  int level,lcoarse,lmax;
  struct OCT * curoct;
  struct OCT * nextoct;

  struct OCT * root_sna;
  struct CELL * rootcell_sna;
  struct PART * rootpart_sna;

  struct OCT * root_mem;
  struct CELL * rootcell_mem;
  struct PART * rootpart_mem;


  struct OCT oct;
  struct OCT * oct_ad;
  int ic,ioct=0;

  struct OCT *freeoct=NULL;
  


  root_mem=firstoct[0]; // the root oct of the grid in memory
  rootcell_mem=&(root_mem->cell[0]); // the root cell of the grid in memory
  rootpart_mem=proot; // the root cell of the grid in memory

  // reset global pointers

  for(level=1;level<=param->lmax;level++){
    firstoct[level-1]=NULL;
    lastoct[level-1]=NULL;
  }


  // opening the file
  fp=fopen(filename,"rb");
  
  size_t outf;

  // reading snapshot time
  outf=fread(tsim,sizeof(REAL),1,fp); 
  outf=fread(tinit,sizeof(REAL),1,fp); 
  outf=fread(nsteps,sizeof(int),1,fp); 
  outf=fread(ndumps,sizeof(int),1,fp); 
  outf=fread(&lcoarse,sizeof(int),1,fp); 
  outf=fread(&lmax,sizeof(int),1,fp); 

  for(level=1;level<=lmax;level++){
    outf=fread(adt+level-1,sizeof(REAL),1,fp); 
  }
  
#ifdef WRAD
  // reading units
  outf=fread(&(param->unit.unit_l),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_v),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_t),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_n),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_mass),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_d),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_N),sizeof(REAL),1,fp);
#endif
  //printf("UNIT L=%e\n",param->unit.unit_l);

#ifdef TESTCOSMO
  outf=fread(&(param->cosmo->om),sizeof(REAL),1,fp);
  outf=fread(&(param->cosmo->ov),sizeof(REAL),1,fp);
  outf=fread(&(param->cosmo->ob),sizeof(REAL),1,fp);
  outf=fread(&(param->cosmo->H0),sizeof(REAL),1,fp);
#endif

  // reading pointer informations in the snapshot
  outf=fread(&root_sna,sizeof(struct OCT*),1,fp);
  outf=fread(&rootcell_sna,sizeof(struct CELL*),1,fp);
  outf=fread(&rootpart_sna,sizeof(struct PART*),1,fp);
  
  
  if(cpu->rank==RANK_DISP) printf(" STARTING OCT READ%p %p %p\n",root_sna,rootcell_sna,rootpart_sna);
  // reading the octs sequence
   
  outf=fread(&oct_ad,sizeof(struct OCT *),1,fp);
  outf=fread(&oct,sizeof(struct OCT),1,fp);

  while(!feof(fp)){

    // do stuff
    ioct++;
    // 1 copy the content of the oct at the right location
    curoct=root_mem+(oct_ad-root_sna);
    //printf("cpu=%d ioct=%d curoct=%p  oct-ad=%p root_sna=%p dif=%ld lev=%d chi=%p\n",cpu->rank,ioct,curoct,oct_ad,root_sna,(unsigned long int) (oct_ad-root_sna),oct.level,oct.cell[0].child);
    //if(curoct-root_mem>param->ngridmax) printf("ERROR BKP\n");
    memcpy(curoct,&oct,sizeof(struct OCT));

    //if(cpu->rank==RANK_DISP) printf("le=%d ic=%d\n",curoct->level,ic);
    // 2.a modify the oct pointers within curoct
 
    curoct->next=(curoct->next==NULL?NULL:(curoct->next-root_sna)+root_mem);
    curoct->prev=(curoct->prev==NULL?NULL:(curoct->prev-root_sna)+root_mem);
    curoct->nexthash=(curoct->nexthash==NULL?NULL:(curoct->nexthash-root_sna)+root_mem);
    
    // 2.c modifity the oct pointers within cells
    for(ic=0;ic<8;ic++){
      curoct->cell[ic].child=(curoct->cell[ic].child==NULL?NULL:(curoct->cell[ic].child-root_sna)+root_mem);
    }
 
    
    // 2.b modify the cell pointers within curoct
    
    for(ic=0;ic<6;ic++){

      if(curoct->nei[ic]!=NULL){

	curoct->nei[ic]=(struct CELL *)(((unsigned long long int)(curoct->nei[ic])-(unsigned long long int)rootcell_sna)+(unsigned long long int)rootcell_mem);
      }
      

    }

    if(curoct->parent!=NULL){

	struct CELL *co;
	if(curoct->level>1){
	  co=(struct CELL *)((unsigned long long int)curoct->parent-(unsigned long long int) rootcell_sna+(unsigned long long int) rootcell_mem);
	  curoct->parent=co;
	}
	else{
	  curoct->parent=root;
	}

    }
    

#ifdef PIC
    // 2.c modifity the particle pointers within cells
    
    for(ic=0;ic<8;ic++){
      curoct->cell[ic].phead=(curoct->cell[ic].phead==NULL?NULL:(curoct->cell[ic].phead-rootpart_sna)+rootpart_mem);
    }
#endif

    // Overall pointer management
    if(firstoct[curoct->level-1]==NULL){
      if(curoct->prev==NULL){
	firstoct[curoct->level-1]=curoct;
      }
    }

    if(lastoct[curoct->level-1]==NULL){
      if(curoct->next==NULL){
	lastoct[curoct->level-1]=curoct;
      }
    }

    // read next oct

    outf=fread(&oct_ad,sizeof(struct OCT *),1,fp);
    outf=fread(&oct,sizeof(struct OCT),1,fp);

  }

  //printf("%d octs recovered by proc %d with root=%p firstoct=%p\n",ioct,cpu->rank,root_mem,firstoct[0]);

  // done
  fclose(fp);


  // Building back the freeoct chained list

  struct OCT * loct;
  for(curoct=root_mem;curoct<root_mem+param->ngridmax;curoct++){
    if(curoct->level==0){
      if(freeoct==NULL){
	freeoct=curoct;
	curoct->prev=NULL;
	curoct->next=NULL;
	loct=curoct;
      }
      else{
	curoct->next=NULL;
	curoct->prev=loct;
	loct->next=curoct;
	loct=curoct;
      }
    }
  }

  //abort();

  /* // check and count the number of freeocts */
  
  /* nextoct=freeoct; */
  /* int ifree=0; */
  /* while(nextoct!=NULL){ */
  /*   curoct=nextoct; */
  /*   ifree++; */
  /*   nextoct=curoct->next; */
  /* } */

  /* printf("ifree=%d\n",ifree); */


  
  return freeoct;

}


//====================================================================================================

void dumpgrid(int levelmax,struct OCT **firstoct, char filename[],REAL tsim, struct RUNPARAMS *param)
{

  int icur,ii,jj,kk;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  struct LOCT loct;
  FILE *fp;
  int noct=0;

  fp=fopen(filename,"wb");

  //printf("tsim=%f\n",tsim);
  fwrite(&tsim,sizeof(REAL),1,fp); 

#ifdef WRAD
  fwrite(&(param->unit.unit_l),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_v),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_t),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_n),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_d),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_N),sizeof(REAL),1,fp);
#endif

  fwrite(&(firstoct[0]),sizeof(struct OCT*),1,fp);


  for(level=param->lcoarse;level<=levelmax;level++) // looping over octs
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];

      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  nextoct=oct.next;

	  oct2loct(&oct,&loct);
	  fwrite(&loct,sizeof(struct LOCT),1,fp);
	  noct++;
/* 	  fwrite(&oct,sizeof(struct OCT),1,fp); */
	}while(nextoct!=NULL);
    }

  //printf("noct=%d\n",noct);
  fwrite(&noct,sizeof(int),1,fp); 
  
  fclose(fp);
}


//====================================================================================================
//=================================================================================================

  //------------------------------------------------------------------------

  //------------------------------------------------------------------------



#ifdef PIC
void dumppart(struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu){

  FILE *fp;
  float val;
  int vali;
  int ipart=0;
  int level;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  struct PART *nexp;
  struct PART *curp;
  int icell;
  float tsimf=tsim;

  int npart=0; 

#ifdef STARS
  int nstar=0; 

  char filenamestar[128];							char filenamepart[128];	
  FILE *fstar;									FILE *fpart;
  sprintf(filenamestar,"data/star.%05d.p%05d",*(cpu->ndumps),cpu->rank);	sprintf(filenamepart,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
  fstar=fopen(filenamestar,"wb");						fpart=fopen(filenamepart,"wb");
  fwrite(&nstar,1,sizeof(int)  ,fstar);						fwrite(&npart,1,sizeof(int)  ,fpart);
  fwrite(&tsimf,1,sizeof(float),fstar);						fwrite(&tsimf,1,sizeof(float),fpart);



#else
  fp=fopen(filename,"wb");
  fwrite(&npart,1,sizeof(int)  ,fp);
  fwrite(&tsimf,1,sizeof(float),fp);
#endif

  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      //printf("level=%d\n",level);
      // setting the first oct
      
      nextoct=firstoct[level-1];
      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
//	  dxcur=1./POW(2,oct.level);
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
		do{
		  curp=nexp;
		  nexp=curp->next;

#ifdef STARS
		  if(curp->isStar) 	{	fp=fstar;	nstar++;	}
		  else 			{	fp=fpart;	npart++;	}
#else
		  npart++;
#endif	
		  val=curp->x;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->y;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->z;			fwrite(&val,1,sizeof(float),fp);
#ifndef PARTN
#ifdef PART_EGY
		  val=curp->ekin+curp->epot;	fwrite(&val,1,sizeof(float),fp);
		  val=curp->fx;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;			fwrite(&val,1,sizeof(float),fp);
#else
		  val=curp->vx;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->vy;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->vz;			fwrite(&val,1,sizeof(float),fp);
#endif
#else
		  val=curp->fx;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->fz;			fwrite(&val,1,sizeof(float),fp);
#endif
		  val=(float)(curp->idx);	fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->mass);	fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->epot);	fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->ekin);	fwrite(&val,1,sizeof(float),fp);
#ifdef STARS
		  if(curp->isStar) {
		    val = curp->age;		fwrite(&val,1,sizeof(float),fp);
		  }
#endif
		  ipart++;

		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }

#ifdef STARS
  rewind(fpart);	fwrite(&npart,1,sizeof(int)  ,fpart);	fclose(fpart);
  rewind(fstar);	fwrite(&nstar,1,sizeof(int)  ,fstar);	fclose(fstar);
#else
  rewind(fp);		fwrite(&npart,1,sizeof(int)  ,fp);	fclose(fp);
#endif

  //printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}

//=======================================================================================================

void save_part(char filename[],struct OCT **firstoct, int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct PART* proot){

  FILE *fp;
  int vali;
  int ipart=0;
  int level;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  struct PART *nexp;
  struct PART *curp;
  int icell;


  int npart=0; 
  for(level=levelcoarse;level<=levelmax;level++) npart+=cpu->npart[level-1];

  fp=fopen(filename,"wb");
  fwrite(&npart,1,sizeof(int),fp);		
  fwrite(&tsim,1,sizeof(REAL),fp);		
  fwrite(&proot,1,sizeof(struct PART *),fp);	

//	printf("cpu %d \t%d\t%e\t%p \n", cpu->rank, npart, tsim,proot );

  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      //printf("level=%d\n",level);
      // setting the first oct
      
      nextoct=firstoct[level-1];
      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./POW(2,oct.level);
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
		do{ 
		  curp=nexp; 
		  nexp=curp->next; 

		  fwrite(&curp,1,sizeof(struct PART *),fp);
		  fwrite(curp,1,sizeof(struct PART),fp);

		  ipart++;
		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }
  	
  fclose(fp);
  //printf("wrote %d particles (%d expected) in %s on proc %d\n",ipart,npart,filename,cpu->rank);

}




// ===================================================================================================
// ===================================================================================================
// ===================================================================================================

struct PART * restore_part(char filename[], struct OCT **firstoct, REAL *tsim, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot){
  
  FILE *fp;
  int level;

  struct PART * rootpart_sna;
  struct PART * rootpart_mem;

  rootpart_mem=proot;

  int ic,ipart=0;
  int npart;
  int sp;

  struct PART part;
  struct PART *part_ad;
  struct PART *curp;

#ifdef STARS
  int nstar=0;
#endif

  rootpart_mem=proot; // the root cell of the grid in memory

  // opening the file
  fp=fopen(filename,"rb");

  // reading snapshot time
  size_t outf;
  outf=fread(&npart,1,sizeof(int),fp);
  outf=fread(tsim,1,sizeof(REAL),fp);
  outf=fread(&rootpart_sna,1,sizeof(struct PART *),fp);

//	printf("cpu %d \t%d\t%e\t%p \n", cpu->rank, npart, *tsim,rootpart_sna);

  // reading the particle sequence
   
  outf=fread(&part_ad,sizeof(struct PART*),1,fp);
  outf=fread(&part,sizeof(struct PART),1,fp);


 // printf("debut de lecture %ld\t%p\t%p\t%p\n ", (unsigned long long int)(part_ad-rootpart_sna), part_ad,rootpart_sna,rootpart_mem );

  while(!feof(fp)){

    /* if(cpu->rank==RANK_DISP){ */
    /*   printf("ipart=%d\n",ipart); */
    /* } */
    // do stuff
    ipart++;

    // 1 copy the content of the particle at the right location

    curp=(part_ad-rootpart_sna)+rootpart_mem;

    memcpy(curp,&part,sizeof(struct PART));

#ifdef STARS
    if(curp->isStar)  nstar++;
#endif

 //   printf("memcpy OK \n");

    // 2.a modify the particle pointers 
    curp->next=(curp->next==NULL?NULL:(curp->next-rootpart_sna)+rootpart_mem);
    curp->prev=(curp->prev==NULL?NULL:(curp->prev-rootpart_sna)+rootpart_mem);
    
  //  printf("*part ok \n");

    // read next particle
    outf=fread(&part_ad,sizeof(struct PART*),1,fp);
    outf=fread(&part,sizeof(struct PART),1,fp);

//    printf("%d\t ", ipart);
  }

//  printf("READ OK \n ");
  // Building back the freepart chained list

  struct PART * lpart;
  struct PART *freepart;
  freepart=NULL;

  for(curp=rootpart_mem;curp<rootpart_mem+param->npartmax;curp++){
    if(curp->mass==-1.){ // flag empty particles
      if(freepart==NULL){
	freepart=curp;
	curp->prev=NULL;
	curp->next=NULL;
	lpart=curp;
      }
      else{
	curp->next=NULL;
	curp->prev=lpart;
	lpart->next=curp;
	lpart=curp;
      }
    }
    /* else if(cpu->rank==RANK_DISP){ */
    /*   printf("%e\n",curp->mass); */
    /* } */
  }

  //printf("%d/%d part recovered by proc %d with freepart=%p\n",ipart,param->npartmax,cpu->rank,freepart);

  // done
  fclose(fp);

#ifdef STARS
#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&nstar,1,MPI_INT,   MPI_SUM,cpu->comm);
#endif
	//printf("nstar=%d\n",nstar);
  param->stars->n=nstar;
#endif

  return freepart;
}

#endif

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

//------------------------------------------------------------------------

void GetParameters(char *fparam, struct RUNPARAMS *param)
{
  FILE *buf; 
  char stream[256];
  size_t rstat;
  float dummyf;


  buf=fopen(fparam,"r");
  if(buf==NULL)
    {
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
    }
  else
    {
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ngridmax);
      rstat=fscanf(buf,"%s %d",stream,&param->npartmax);
      rstat=fscanf(buf,"%s %d",stream,&param->nbuff);
      
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ndumps);
      rstat=fscanf(buf,"%s %d",stream,&param->nsteps);
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->dt=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->tmax=(REAL)dummyf;

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->lcoarse);
      rstat=fscanf(buf,"%s %d",stream,&param->lmax);
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->amrthresh=(REAL)dummyf;
      rstat=fscanf(buf,"%s %d",stream,&param->nsmooth);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->niter);
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->poissonacc=(REAL)dummyf;
      rstat=fscanf(buf,"%s %d",stream,&param->mgridlmin);
      rstat=fscanf(buf,"%s %d",stream,&param->nvcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nrelax);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->nrestart);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->gstride);
      rstat=fscanf(buf,"%s %d",stream,&param->hstride);
      rstat=fscanf(buf,"%s %d",stream,&param->nsubcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nthread);
      rstat=fscanf(buf,"%s %d",stream,&param->nstream);
      rstat=fscanf(buf,"%s %d",stream,&param->ompthread);

#ifdef WRAD
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->clight=(REAL)dummyf;param->clightorg=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->denthresh=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->tmpthresh=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->srcint=(REAL)dummyf;
      param->fudgecool=1.0;
      param->ncvgcool=0;
#else
	int i;
				rstat=fscanf(buf,"%s",stream);
	for (i=0; i<4; i++)	rstat=fscanf(buf,"%s %f",stream,&dummyf);
#endif

#ifdef STARS
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->stars->overdensity_cond	=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->stars->density_cond			=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->stars->tcar							=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->stars->tlife							=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->stars->feedback_eff			=(REAL)dummyf;
      rstat=fscanf(buf,"%s %f",stream,&dummyf);param->stars->feedback_frac			=(REAL)dummyf;
#endif

#ifdef MOVIE
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d", stream,&param->movie->lmap);
      rstat=fscanf(buf,"%s %f",stream,&param->movie->xmin);
      rstat=fscanf(buf,"%s %f",stream,&param->movie->xmax);
      rstat=fscanf(buf,"%s %f",stream,&param->movie->ymin);
      rstat=fscanf(buf,"%s %f",stream,&param->movie->ymax);
      rstat=fscanf(buf,"%s %f",stream,&param->movie->zmin);
      rstat=fscanf(buf,"%s %f",stream,&param->movie->zmax);

#endif

      fclose(buf);
    }




  // computing the maxhash
  int val=(POW(2,param->lmax-1)<256?POW(2,param->lmax-1):256); // limit to 2097152 octs in hash table i.e. 16e6 cells
  param->maxhash=POW(val,3);
  //printf("maxhash=%d\n",param->maxhash);

  // ====================== some checks

  // stencil/streams conformity
#ifdef GPUAXL
  if(param->hstride<(param->nthread*param->nstream)){
    printf(" Stream Thread granulosity too high : nt=%d ns=%d stencil=%d\n",param->hstride,param->nthread,param->nstream);
    abort();
  }
#endif

#ifdef STARS
    param->stars->n		= 0;
#endif

}

//==================================================================================
//==================================================================================
#ifdef PIC
#ifdef GRAFIC
struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param,int level)
{
  FILE *fx, *fy, *fz;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
  int dummy;
  struct PART *lastpart;
  int ip;
  char filename[256];

  if(cpu->rank==0){

    sprintf(filename,"./level_%03d/ic_velbx",level); 
    fx=fopen(filename,"rb");
    sprintf(filename,"./level_%03d/ic_velby",level); 
    fy=fopen(filename,"rb");
    sprintf(filename,"./level_%03d/ic_velbz",level); 
    fz=fopen(filename,"rb");
    

    // reading the headers
    size_t outf;

    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&np1,1,4,fx);
    outf=fread(&np2,1,4,fx);
    outf=fread(&np3,1,4,fx);
    outf=fread(&dx,1,4,fx);
    printf("DX=%e\n",dx);
    outf=fread(&x1o,1,4,fx);
    outf=fread(&x2o,1,4,fx);
    outf=fread(&x3o,1,4,fx);
    outf=fread(&astart,1,4,fx);
    outf=fread(&om,1,4,fx);
    outf=fread(&ov,1,4,fx);
    outf=fread(&h0,1,4,fx);
    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&dummy,1,sizeof(dummy),fy);
    outf=fread(&np1,1,4,fy);
    outf=fread(&np2,1,4,fy);
    outf=fread(&np3,1,4,fy);
    outf=fread(&dx,1,4,fy);
    outf=fread(&x1o,1,4,fy);
    outf=fread(&x2o,1,4,fy);
    outf=fread(&x3o,1,4,fy);
    outf=fread(&astart,1,4,fy);
    outf=fread(&om,1,4,fy);
    outf=fread(&ov,1,4,fy);
    outf=fread(&h0,1,4,fy);
    outf=fread(&dummy,1,sizeof(dummy),fy);

    outf=fread(&dummy,1,sizeof(dummy),fz);
    outf=fread(&np1,1,4,fz);
    outf=fread(&np2,1,4,fz);
    outf=fread(&np3,1,4,fz);
    outf=fread(&dx,1,4,fz);
    outf=fread(&x1o,1,4,fz);
    outf=fread(&x2o,1,4,fz);
    outf=fread(&x3o,1,4,fz);
    outf=fread(&astart,1,4,fz);
    outf=fread(&om,1,4,fz);
    outf=fread(&ov,1,4,fz);
    outf=fread(&h0,1,4,fz);
    outf=fread(&dummy,1,sizeof(dummy),fz);
  }

#ifdef WMPI
  // Massive broadcast of grafic header
  MPI_Bcast(&np1,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np2,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np3,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&dx,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x1o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x2o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x3o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&astart,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&om,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&ov,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&h0,1,MPI_FLOAT,0,cpu->comm);
  MPI_Barrier(cpu->comm);
#endif

  if(cpu->rank==0){
    printf("============================================\n");
    printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
    printf("om=%f ov=%f h0=%f\n",om,ov,h0);
    printf("dx=%f np1*dx=%f\n",dx,np1*dx);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }

  if(level==param->lcoarse){
  if((np1*np2*np3)/(cpu->nproc)*1.>param->npartmax){
    printf("Error : Number of particles greater than npartmax Np(grafic, est. )=%d Npartmax=%d!\n",(np1*np2*np3)/(cpu->nproc),param->npartmax);
    abort();
  }
  }
  //setting omegab

  ob=OMEGAB;


  // computing Zeldovich displacement quantities
  
  double vfact;
  vfact=fomega(astart,om,ov)*h0*dladt(astart,om,ov)/astart;
  if(cpu->rank==0) printf("vfact=%f\n",vfact);
  // reading the grafic planes

  float *velx;
  float *vely;
  float *velz;

  double x;
  double y;
  double z;

  double vx;
  double vy;
  double vz;

  double x0,y0,z0;

  int i1,i2,i3;
  int offidx=0;
  int keep;
  double mass;

#ifdef WHYDRO2
  mass=(1.-ob/om)/(np1*np2*np3);
#else
  mass=1./(np1*np2*np3);
#endif

#ifdef ZOOM
  if(level>param->lcoarse){
    offidx=POW(2,3*(level-param->lcoarse)); // for ids of particles
  }
#endif

  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);

  //REAL rmin=2.;

  ip=0;
  size_t outf;
  for(i3=1;i3<=np3;i3++){

    if(cpu->rank==0){
      printf("\r %f percent done",(i3*1.0)/np3*100.);
      outf=fread(&dummy,1,sizeof(dummy),fx);
      outf=fread(velx,np1*np2,sizeof(float),fx);
      outf=fread(&dummy,1,sizeof(dummy),fx);
      
      outf=fread(&dummy,1,sizeof(dummy),fy);
      outf=fread(vely,np1*np2,sizeof(float),fy);
      outf=fread(&dummy,1,sizeof(dummy),fy);
      
      outf=fread(&dummy,1,sizeof(dummy),fz);
      outf=fread(velz,np1*np2,sizeof(float),fz);
      outf=fread(&dummy,1,sizeof(dummy),fz);
    }


#ifdef WMPI
    MPI_Barrier(cpu->comm);
    MPI_Bcast(velx,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(vely,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velz,np1*np2,MPI_FLOAT,0,cpu->comm);
#endif
    

    z0=(i3-0.5)*dx;
    for(i2=1;i2<=np2;i2++){
      y0=(i2-0.5)*dx;
      for(i1=1;i1<=np1;i1++){
	x0=(i1-0.5)*dx;
	// computing the displacements
	x=(x0+velx[(i1-1)+(i2-1)*np1]/vfact)/(np1*dx);
	y=(y0+vely[(i1-1)+(i2-1)*np1]/vfact)/(np2*dx);
	z=(z0+velz[(i1-1)+(i2-1)*np1]/vfact)/(np3*dx);

	// periodic boundary conditions

	x+=(x<=0.)*((int)(-x)+1.)-(x>1.)*((int)x); 
	y+=(y<=0.)*((int)(-y)+1.)-(y>1.)*((int)y); 
	z+=(z<=0.)*((int)(-z)+1.)-(z>1.)*((int)z); 

	// computing the velocities
	vx=velx[(i1-1)+(i2-1)*np1]*astart/(np1*dx*h0)/(sqrt(om)*0.5);
	vy=vely[(i1-1)+(i2-1)*np1]*astart/(np2*dx*h0)/(sqrt(om)*0.5);
	vz=velz[(i1-1)+(i2-1)*np1]*astart/(np3*dx*h0)/(sqrt(om)*0.5);

	// if it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,cpu,cpu->levelcoarse)){

	  keep=1;
#ifdef ZOOM
	  // is the current particle at the correct level?
	  int lzoom;
	  lzoom=pos2levelzoom(x,y,z,param);
	  REAL rloc;
	  rloc=sqrt(((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)));
	  
	  if(lzoom!=level){
	    keep=0;
	  }
 #endif

	  if(keep) {
 	    part[ip].x=x;
	    part[ip].y=y;
	    part[ip].z=z;
	    
 	    //rmin=(rloc<rmin?rloc:rmin);

	    part[ip].vx=vx;
	    part[ip].vy=vy;
	    part[ip].vz=vz;
	    part[ip].level=level;
	    
	    part[ip].mass=mass;
	    part[ip].idx=(i1-1)+(i2-1)*np1+(i3-1)*np1*np2+offidx;
	    lastpart=part+ip;
	    ip++; 
	  }
	}
      }
    }
  }

  free(velx);
  free(vely);
  free(velz);

  if(cpu->rank==0){
    fclose(fx);
    fclose(fy);
    fclose(fz);
  }

  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=ip;
#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif

  if(cpu->rank==RANK_DISP){
    printf("Grafic Particle Read ok\n");
  }
  return lastpart;
}
#endif

// ========================== ZELDOVICH
// ====================================
// ====================================

#ifdef ZELDOVICH
struct PART * read_zeldovich_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param, struct OCT **firstoct)
{
  FILE *fd;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
  int dummy;
  struct PART *lastpart;
  int ip;
  
  int nploc;
  float munit_z;
  float lbox;
  float ainit_z;
  REAL mass;
  size_t outf;

  fd=fopen("utils/grafic_src/ZEL.PM.0","rb");

  outf=fread(&dummy,sizeof(dummy),1,fd);
  outf=fread(&nploc,sizeof(int),1,fd);	 
  outf=fread(&munit_z,sizeof(float),1,fd); 
  outf=fread(&ainit_z,sizeof(float),1,fd); 
  outf=fread(&lbox,sizeof(float),1,fd);	 
  outf=fread(&om,sizeof(float),1,fd);
  outf=fread(&ov,sizeof(float),1,fd);
  outf=fread(&h0,sizeof(float),1,fd);
  outf=fread(&dummy,sizeof(dummy),1,fd);  

  astart=ainit_z;


  if(cpu->rank==0){
    printf("============================================\n");
    printf("ntot%d\n",nploc);
    printf("om=%f ov=%f h0=%f\n",om,ov,h0);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }

  if((nploc)/(cpu->nproc)*1.>param->npartmax){
    printf("Error : Number of particles greater than npartmax Np(grafic, est. )=%d Npartmax=%d!\n",(nploc)/(cpu->nproc),param->npartmax);
    abort();
  }

  //setting omegab

  ob=OMEGAB;

  
#ifdef WHYDRO2
  mass=(1.-ob/om)/(nploc);
#else
  mass=1./(nploc);
#endif


  // reading the grafic planes
  float *pos;
  float *vel;
  int nread=nploc; // quick fixes
  int npatch=1.; // quick fixes
  int ipatch;
  int i;
  pos=(float *)malloc(sizeof(REAL)*3*nread);
  vel=(float *)malloc(sizeof(REAL)*3*nread);

  double x;
  double y;
  double z;

  double vx;
  double vy;
  double vz;
  size_t rstat;

  int pstart=ftell(fd);

  ip=0.;
  for(ipatch=0;ipatch<npatch;ipatch++) {
    //    rstat=outf=fread(&dummy,sizeof(dummy),1,fd); 
    //    fseek(fd,pstart,SEEK_SET);
    fseek(fd,pstart+(0*nploc+ipatch*nread)*sizeof(float)+1*sizeof(dummy),SEEK_SET);
    outf=fread(pos,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(1*nploc+ipatch*nread)*sizeof(float)+3*sizeof(dummy),SEEK_SET);
    outf=fread(pos+nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(2*nploc+ipatch*nread)*sizeof(float)+5*sizeof(dummy),SEEK_SET);
    outf=fread(pos+2*nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);
  
    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(3*nploc+ipatch*nread)*sizeof(float)+7*sizeof(dummy),SEEK_SET);
    outf=fread(vel,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(4*nploc+ipatch*nread)*sizeof(float)+9*sizeof(dummy),SEEK_SET);
    outf=fread(vel+nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(5*nploc+ipatch*nread)*sizeof(float)+11*sizeof(dummy),SEEK_SET);
    outf=fread(vel+2*nread,sizeof(float),nread,fd);
    //outf=fread(&dummy,sizeof(dummy),1,fd);
    

    for(i=0;i<nread;i++)
      {
	x=pos[i];
	y=pos[i+nread];
	z=pos[i+2*nread];
	vx=vel[i];
	vy=vel[i+nread];
	vz=vel[i+2*nread];
	// periodic boundary conditions
	
	x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
	y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
	z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 
	// it it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,cpu,cpu->levelcoarse)){

	  part[ip].x=x;
	  part[ip].y=y;
	  part[ip].z=z;
	  
	  part[ip].vx=vx;
	  part[ip].vy=vy;
	  part[ip].vz=vz;
	  
	  part[ip].mass=mass;
	  part[ip].idx=i;
	  lastpart=part+ip;
	  ip++;
	}
      }
  }


  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=ip;

  if(cpu->rank==RANK_DISP){
    printf("Zeldovich Particle Read ok\n");
  }

#ifdef WHYDRO2
  
  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;
  REAL ZC=10; // hard coded collapse of redshift
  REAL ZI=1./astart-1.;

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
 		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		W.d=(1.+(1.+ZC)/(1.+ZI)*cos(2.*M_PI*(xc-0.5)))*ob/om;
		W.p=PMIN;
		W.u=-(1.+ZC)/POW(1.+ZI,1.5)*sin(2.*M_PI*(xc-0.5))/(M_PI); // for omegam=1. only
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);
		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

#endif

  return lastpart;
}
#endif



// =====================================================================================
// =====================================================================================

#ifdef EDBERT
struct PART * read_edbert_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param, struct OCT **firstoct)
{
  float astart,om,ov,h0,ob;
  int dummy;

  struct PART *lastpart;
  int ip,iploc;
  int i,j,k;

  float lbox;
  float ainit_z;
  REAL mass;
  REAL x,y,z,r;
  REAL delta;
  REAL lsphere;
  printf("Start EDBERT\n");

  lsphere=0.05;
  om=0.99999;
  ov=0.00001;
  ob=OMEGAB;
  delta=0.2;
  h0=70.;
  lbox=1.;
  astart=1e-3;

  int n1d=POW(2,param->lcoarse);
  REAL dx=1./n1d;
  iploc=0;
  int nin=0;
  int nout=0;
  REAL m;

#ifdef PIC

  REAL mout,mint;
  REAL vsphere=0.;
  for(k=0;k<n1d;k++)
    {
      for(j=0;j<n1d;j++)
	{
	  for(i=0;i<n1d;i++)
	    {
	      x=(i+0.5)*dx;
	      y=(j+0.5)*dx;
	      z=(k+0.5)*dx;
	      
	      r=sqrt(POW(x-0.5,2)+POW(y-0.5,2)+POW(z-0.5,2));
	      if(r<lsphere){
		nin++;
		mout=-1;
		vsphere+=POW(dx,3);
	      }
	      else{
		nout++;
		mout=1;
	      }
	      
	      if(segment_part(x,y,z,cpu,cpu->levelcoarse)){
		part[iploc].x=x;
		part[iploc].y=y;
		part[iploc].z=z;
		
		part[iploc].vx=0;
		part[iploc].vy=0;
		part[iploc].vz=0;
	
		part[iploc].mass=mout;
		part[iploc].idx=-nin;
		lastpart=part+iploc;
		iploc++;
	      }
	    }
	}
    }

  mint=(om-ob)*(1.+delta)*vsphere/nin;
  mout=((om-ob)-mint*nin)/nout;
  printf("mint=%e mout=%e\n",mint,mout);
  
  for(i=0;i<iploc;i++){
    part[i].mass=(part[i].mass<0?mint:mout);
  }

#endif
  
  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=iploc;

  if(cpu->rank==RANK_DISP){
    printf("Edbert Particle Read ok\n");
  }

#ifdef WHYDRO2

  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		rad=sqrt((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5));
		if(rad<lsphere){
		  W.d=(1.+delta)*ob;
		}
		else{
		  W.d=ob*(1.-(1.+delta)*vsphere)/(1.-vsphere);
		}
		W.p=PMIN;
		W.u=0.; // vstar is expressed in m/s and grafic vel are in km/s
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);

		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

#endif

  return lastpart;
}
#endif


#endif
//==================================================================================
//==================================================================================

#ifdef WHYDRO2

#ifdef TUBE
// =====================================================================================================
// =====================================================================================================

void read_shocktube(struct CPUINFO *cpu, REAL *ainit, struct RUNPARAMS *param, struct OCT **firstoct)
{
  FILE *fd;
  
  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;


  struct Wtype WL, WR;
  if(cpu->rank==RANK_DISP) printf("Init Hydro\n");
  
  /* /\*  /\\* // TEST 1 *\\/ *\/ */
  
  WL.d=1.;
  WL.u=0.;
  WL.v=0.;
  WL.w=0.;
  WL.p=1.0;
  WL.a=sqrt(GAMMA*WL.p/WL.d);
  getE(&WL);
  
  WR.d=0.125;
  WR.u=0.;
  WR.v=0.;
  WR.w=0.;
  WR.p=0.1;
  WR.a=sqrt(GAMMA*WR.p/WR.d);
  getE(&WR);

  REAL X0=0.3125; 

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		if(xc<X0){
		  memcpy(&(curoct->cell[icell].field),&WL,sizeof(struct Wtype)); 
		}
		else{
		  memcpy(&(curoct->cell[icell].field),&WR,sizeof(struct Wtype)); 
		}
	      }
	  }while(nextoct!=NULL);
      }
}
#endif

#ifdef EVRARD
int read_evrard_hydro(struct CPUINFO *cpu,struct OCT **firstoct, struct RUNPARAMS *param){
  
  int level;
  REAL dxcur;
  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  struct Wtype W;
  REAL rad;
  REAL xc,yc,zc;
  
  //==== parameters of evrard sphere
  REAL R=0.35;
  REAL M=1.;
  REAL rhostar=M/(4./3.*M_PI*R*R*R);
  REAL estar=M/R; //assuming G=1
  REAL pstar=rhostar*estar;
  REAL tstar=sqrt(M_PI*M_PI/8.)*POW(R,1.5)/POW(M,0.5);
  if(cpu->rank==RANK_DISP) printf("Generating Evrard Test Case ts=%e, rhostar=%e\n",tstar,rhostar);

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		rad=sqrt((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5))/R;
		if(rad<1.){
		  W.d=rhostar/rad;
		  W.p=pstar/rad*0.05;
		}
		else{
		  W.d=1e-3;
		  W.p=1e-5;
		}

		W.u=0.; // vstar is expressed in m/s and grafic vel are in km/s
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);

		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

  }
#endif


#ifdef TESTCOSMO
#ifdef GRAFIC
 int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, struct RUNPARAMS *param,int level){
  
  FILE *fx;
  FILE *fy;
  FILE *fz;
  FILE *fdx;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,ob,h0;
  int dummy;
  int ip;
  struct Wtype W;
  size_t outf;
  char filename[256];

  // Note only the rank 0 reads the file.
 
  if(cpu->rank==0){

    sprintf(filename,"./level_%03d/ic_deltab",level); 
    fdx=fopen(filename,"rb");
    sprintf(filename,"./level_%03d/ic_velbx",level); 
    fx=fopen(filename,"rb");
    sprintf(filename,"./level_%03d/ic_velby",level); 
    fy=fopen(filename,"rb");
    sprintf(filename,"./level_%03d/ic_velbz",level); 
    fz=fopen(filename,"rb");

    /* fdx=fopen("./ic_deltab","rb"); */
    /* fx=fopen("./ic_velcx","rb"); */
    /* fy=fopen("./ic_velcy","rb"); */
    /* fz=fopen("./ic_velcz","rb"); */
  
    // reading the headers

    outf=fread(&dummy,1,sizeof(dummy),fdx);
    outf=fread(&np1,1,4,fdx);
    outf=fread(&np2,1,4,fdx);
    outf=fread(&np3,1,4,fdx);
    outf=fread(&dx,1,4,fdx);
    outf=fread(&x1o,1,4,fdx);
    outf=fread(&x2o,1,4,fdx);
    outf=fread(&x3o,1,4,fdx);
    outf=fread(&astart,1,4,fdx);
    outf=fread(&om,1,4,fdx);
    outf=fread(&ov,1,4,fdx);
    outf=fread(&h0,1,4,fdx);
    outf=fread(&dummy,1,sizeof(dummy),fdx);

    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&np1,1,4,fx);
    outf=fread(&np2,1,4,fx);
    outf=fread(&np3,1,4,fx);
    outf=fread(&dx,1,4,fx);
    outf=fread(&x1o,1,4,fx);
    outf=fread(&x2o,1,4,fx);
    outf=fread(&x3o,1,4,fx);
    outf=fread(&astart,1,4,fx);
    outf=fread(&om,1,4,fx);
    outf=fread(&ov,1,4,fx);
    outf=fread(&h0,1,4,fx);
    outf=fread(&dummy,1,sizeof(dummy),fx);

    outf=fread(&dummy,1,sizeof(dummy),fy);
    outf=fread(&np1,1,4,fy);
    outf=fread(&np2,1,4,fy);
    outf=fread(&np3,1,4,fy);
    outf=fread(&dx,1,4,fy);
    outf=fread(&x1o,1,4,fy);
    outf=fread(&x2o,1,4,fy);
    outf=fread(&x3o,1,4,fy);
    outf=fread(&astart,1,4,fy);
    outf=fread(&om,1,4,fy);
    outf=fread(&ov,1,4,fy);
    outf=fread(&h0,1,4,fy);
    outf=fread(&dummy,1,sizeof(dummy),fy);

    outf=fread(&dummy,1,sizeof(dummy),fz);
    outf=fread(&np1,1,4,fz);
    outf=fread(&np2,1,4,fz);
    outf=fread(&np3,1,4,fz);
    outf=fread(&dx,1,4,fz);
    outf=fread(&x1o,1,4,fz);
    outf=fread(&x2o,1,4,fz);
    outf=fread(&x3o,1,4,fz);
    outf=fread(&astart,1,4,fz);
    outf=fread(&om,1,4,fz);
    outf=fread(&ov,1,4,fz);
    outf=fread(&h0,1,4,fz);
    outf=fread(&dummy,1,sizeof(dummy),fz);
  }

  // setting baryon density parameter
  ob=OMEGAB;

#ifdef WMPI
  // Massive broadcast of grafic header
  MPI_Bcast(&np1,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np2,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np3,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&dx,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x1o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x2o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x3o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&astart,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&om,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&ov,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&h0,1,MPI_FLOAT,0,cpu->comm);
  MPI_Barrier(cpu->comm);

#endif
  
  if(cpu->rank==0){
    printf("============================================\n");
    printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
    printf("om=%f ov=%f ob=%f h0=%f\n",om,ov,ob,h0);
    printf("dx=%f np1*dx=%f\n",dx,np1*dx);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }



  if(np1!=(int)POW(2,level)){
    printf("ERROR !ABORT! Grafic hydro  file not compliant with parameter file : ngrafic=%d nquartz=%d\n",np1,(int)POW(2,level));
    abort();
  }


  // reading the grafic planes

  float *deltab;
  float *velz;
  float *vely;
  float *velx;

  int i1,i2,i3;
  int icx,icy,icz,icell;
  unsigned long long key;
  struct OCT *curoct;
  struct OCT *nextoct;
  unsigned long hidx;
  int found;
  float z0,y0,x0;
  int ifound=0;
  

  deltab=(float*)malloc(sizeof(float)*np1*np2);
  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);

  double rhob,pressure;
  double H0=h0*1e3/3.08568025e22; // Hubble constant (s-1)
  double rhoc=3.*H0*H0/(8.*M_PI*NEWTON_G); // comoving critical density (kg/m3)
  double zstart=1./astart-1.;

  // ---------- Setting the initial temperature ---- //
  double temp;
  //  double temp=550.*((1.0+zstart)*(1.0+zstart)); // baryon temperature (to check) in K
  //double temp=2.7*(1+zstart);
  //double temp=1e4;
  //double temp=170.*(1.+zstart)*(1.+zstart)/10000.;

  //double temp=0.0874545+0.0302621*zstart+0.00675076*zstart*zstart; // recfast ob fit
#ifdef COOLING
  if(om==1.) {
    temp=33.64/POW(41.,2)*POW(1.+zstart,2);
    if(cpu->rank==RANK_DISP) printf("WARNING: YOU ARE USING SCDM COSMOLOGY\n");
  }
  else{
    if(cpu->rank==RANK_DISP) printf("No temperature law for cosmologies other than SCDM -> F** it\n");
    temp=33.64/POW(41.,2)*POW(1.+zstart,2);
    //    abort();
  }
#else
  temp=1e4;
#endif

  // supercomoving unit values
  double rhostar;
  double rstar;
  double vstar;
  double tstar;
  double tstar2;
  double pstar;
  double mpc=3.08568025e22; // Mpc in m

  rstar= np1*dx*mpc; // box size in m
  rhostar=rhoc*om;
  tstar=2./H0/sqrt(om); // sec
  tstar2=2./h0/sqrt(om); // Mpc sec / km
  vstar=rstar/tstar; //m/s
  pstar=rhostar*vstar*vstar;
  
  if(cpu->rank==RANK_DISP) printf("rhoc=%e temperature=%lf rstar=%e(%e) pstar=%e tstar=%e vstar=%e rhostar=%e\n",rhoc,temp,rstar,np1*dx,pstar,tstar,vstar,rhostar);


  for(i3=0;i3<np3;i3++){

    if(cpu->rank==0){
      printf("\r %f percent done",(i3*1.0)/np3*100.);

      outf=fread(&dummy,1,sizeof(dummy),fdx);
      outf=fread(deltab,np1*np2,sizeof(float),fdx);
      outf=fread(&dummy,1,sizeof(dummy),fdx);

      outf=fread(&dummy,1,sizeof(dummy),fx);
      outf=fread(velx,np1*np2,sizeof(float),fx);
      outf=fread(&dummy,1,sizeof(dummy),fx);

      outf=fread(&dummy,1,sizeof(dummy),fy);
      outf=fread(vely,np1*np2,sizeof(float),fy);
      outf=fread(&dummy,1,sizeof(dummy),fy);

      outf=fread(&dummy,1,sizeof(dummy),fz);
      outf=fread(velz,np1*np2,sizeof(float),fz);
      outf=fread(&dummy,1,sizeof(dummy),fz);
    }

#ifdef WMPI
    MPI_Barrier(cpu->comm);
    MPI_Bcast(deltab,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velx,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(vely,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velz,np1*np2,MPI_FLOAT,0,cpu->comm);
#endif

    z0=(i3*1.0)/(np3);
    for(i2=0;i2<np2;i2++){
      y0=(i2*1.0)/(np2);
      for(i1=0;i1<np1;i1++){
	x0=(i1*1.0)/(np1);
	
	key=pos2key(x0,y0,z0,level);

	// first we compute the adress from the hashfunction
	hidx=hfun(key,cpu->maxhash);
	nextoct=cpu->htable[hidx];

	// looking for the oct
	found=0;
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,level)==key)&&(curoct->level==level));
	  }while((nextoct!=NULL)&&(!found));
	}

	// filling the cell

	if(found){
	  icx=i1%2;
	  icy=i2%2;
	  icz=i3%2;
	  
	  icell=icx+icy*2+icz*4;
	
	  rhob=(deltab[i1+i2*np1]+1.0)*ob*rhoc/POW(astart,3); // physical baryon density in kg/m3
	  pressure=(GAMMA-1.0)*1.5*(rhob/(PROTON_MASS*MOLECULAR_MU))*KBOLTZ*temp; // physical pressure
	  
	  //printf("pres=%e\n",pressure);
	  // filling the cells using supercomoving values
	  
	  //abort();

	  W.d=(deltab[i1+i2*np1]+1.0)*ob/om;
	  W.u=(velx[i1+i2*np1]*1e3)*astart/vstar; // vstar is expressed in m/s and grafic vel are in km/s
	  W.v=(vely[i1+i2*np1]*1e3)*astart/vstar;
	  W.w=(velz[i1+i2*np1]*1e3)*astart/vstar;
	  W.p=pressure/pstar*POW(astart,5);
	  W.a=sqrt(GAMMA*W.p/W.d);
	  getE(&W);

#ifdef WRADHYD
	  // Testing ADVECTION
	  //W.X=(i1/6)%2+((i2+1)/6)%2;
	  W.dX=0.2e-3*W.d;
#endif
	  memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	  ifound++;
	}
	else{

	  // this branch corresponds to cell out of the domain
	  //	  printf("euh pas trouve! hidx=%d %p",hidx,cpu->htable[hidx]);
	  //	  abort();
	}
      }
    }
  }

  if(cpu->rank==0){
    fclose(fdx);
    fclose(fx);
    fclose(fy);
    fclose(fz);
  }


  free(deltab);
  free(velx);
  free(vely);
  free(velz);

  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  param->cosmo->H0=h0;
  param->cosmo->unit_l=rstar;

#ifdef WRAD
  param->unit.unit_l=rstar;
  param->unit.unit_v=vstar;
  param->unit.unit_t=param->unit.unit_l/param->unit.unit_v;
  param->unit.unit_n=1.;//(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS)*rhostar; // 
  param->unit.unit_mass=rhostar*POW(param->unit.unit_l,3);
  param->unit.unit_d=rhostar; // kg/m3
  param->unit.unit_N=rhostar/PROTON_MASS; // atom/m3
#endif

//  REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS*MOLECULAR_MU/param->unit.unit_mass)/POW(param->unit.unit_l,3);
 /*   REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS*MOLECULAR_MU)*rhostar; */

/* if(cpu->rank==RANK_DISP) printf("navg=%e \n",navg); */
  
#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif
  if(cpu->rank==RANK_DISP) printf("Grafic hydro read ok\n");
  return ifound;
}
#endif
#endif
#endif

// ======================================================================
// ======================================================================

void dumpIO(REAL tsim, struct RUNPARAMS *param,struct CPUINFO *cpu, struct OCT **firstoct, REAL *adt, int pdump){

  REAL tdump,adump;
  char filename[128]; 
  int idir=cpu->rank%8;

#ifndef TESTCOSMO
#ifdef WRAD
	tdump=(tsim)*param->unit.unit_t/MYR;
#else
	tdump=(tsim);
#endif
	adump=tdump;
#else
	tdump=interp_aexp(tsim,(double*)param->cosmo->tab_aexp,(double*)param->cosmo->tab_ttilde);
	adump=tdump;
#endif

	if(pdump){
	  // === particle dump
#ifdef PIC
	  sprintf(filename,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
	  if(cpu->rank==RANK_DISP){
	    printf("Dumping .......");
	    printf("%s %p\n",filename,cpu->part);
	  }
	  dumppart(firstoct,filename,param->lcoarse,param->lmax,tdump,cpu);

#endif
	}
	else{
	  // === Hydro dump
    
	  sprintf(filename,"data/grid.%05d.p%05d",*(cpu->ndumps),cpu->rank); 
	  if(cpu->rank==RANK_DISP){
	    printf("Dumping .......");
	    printf("%s\n",filename);
	  }
	  dumpgrid(param->lmax,firstoct,filename,adump,param); 

	  // backups for restart

	  if(*(cpu->ndumps)%FBKP==0){

	    if(cpu->rank==RANK_DISP){
	      printf("BACKUP .......#%d\n",*cpu->ndumps%2);
	    }

	    
	    sprintf(filename,"bkp/grid.%05d.p%05d",*(cpu->ndumps)%2+1,cpu->rank); 
	    save_amr(filename,firstoct,tdump,cpu->tinit,cpu->nsteps,*(cpu->ndumps),param,cpu,cpu->firstpart,adt);
	    
#ifdef PIC
	    // backups for restart
	    sprintf(filename,"bkp/part.%05d.p%05d",*(cpu->ndumps)%2+1,cpu->rank); 
	    save_part(filename,firstoct,param->lcoarse,param->lmax,tdump,cpu,cpu->firstpart);
#endif
	  }
	}	
}
