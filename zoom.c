#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"
#include "amr.h"
#include "hydro_utils.h"
#include "tools.h"
#ifdef WRAD
#include "rad_utils.h"
#include "src_utils.h"
#endif
#include "friedmann.h"
#include "cic.h"
#include "particle.h"


#include "communication.h"

#ifdef STARS
#include "stars.h"
#endif

#ifdef ZOOM

int queryzoom(struct OCT *curoct, int icell, REAL dxcur, REAL Rin) {
  
  REAL xc,yc,zc;
  int res=0.;

  // we assume that the zoom region is at the center
  xc=curoct->x+( icell&1)*dxcur+dxcur*0.5    -0.5;
  yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5-0.5;
  zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5  -0.5;

  if((xc*xc+yc*yc+zc*zc)<(Rin*Rin)){
    res=1;
  }
  return res;  
}

// ===========================================================
// ===========================================================

int pos2levelzoom(REAL xc, REAL yc, REAL zc, struct RUNPARAMS *param){

  int res=param->lmaxzoom;

  // we assume that the zoom region is at the center
  xc-=0.5;
  yc-=0.5;
  zc-=0.5;

  REAL rc=sqrt(xc*xc+yc*yc+zc*zc);
  REAL rcur=param->rzoom;
  
  while((rc>rcur)&&(res>param->lcoarse)){
    //      printf("res=%d lcoarse=%d\n",res,param->lcoarse);
      rcur*=param->fzoom;
      res=res-1;
      //   abort();
  }
  
  return res;

}




void zoom_level(int level, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct){
 
#ifdef TESTCOSMO
  struct COSMOPARAM *cosmo;
  cosmo=param->cosmo;
#endif
  struct OCT *curoct;
  REAL dtnew;
  REAL dt=0.;
  REAL dtold;
  REAL dtvel;
  REAL dtfine;
  int npart=0;
  int mtot;
  int is;
  REAL aexp;
  int nsource;
  int hstride=param->hstride;
  int gstride=param->gstride;
  int nsub=param->nsubcycles;
  int ip;

  int ptot[2];
  int deltan[2];
  ptot[0]=0;
  ptot[1]=0;


#ifdef TESTCOSMO
  aexp=cosmo->aexp;
#else
  aexp=1.0;
#endif

  if(cpu->rank==RANK_DISP){
    printf("\n === zooming level =%d\n",level);
  }

  // ==================================== Check the number of particles and octs
  ptot[0]=0; for(ip=1;ip<=param->lmax;ip++){
    ptot[0]+=cpu->npart[ip-1]; // total of local particles
  }

#ifdef STARS
  ptot[1]=0; for(ip=1;ip<=param->lmax;ip++) ptot[1]+=cpu->nstar[ip-1];
#endif


  mtot=multicheck(firstoct,ptot,param->lcoarse,param->lmax,cpu->rank,cpu,param,0);

  if(level==param->lcoarse){
    // =========== Grid Census ========================================
    grid_census(param,cpu);
  }
  
  is=0;

#ifdef PIC
  //reset substep index of the particles
  L_reset_is_part(level,firstoct);
#endif
 
  
  // ================= I we refine the current level
  if((param->lmax!=param->lcoarse)&&(level<param->lmaxzoom)){

    // enforcing the 2 levels rule
      L_check_rule(level,param,firstoct,cpu);
      
#ifdef WMPI
      mpi_cic_correct(cpu, cpu->sendbuffer, cpu->recvbuffer, 3);
      mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,3,1,level); // propagate the rule check
#ifdef WHYDRO2
      mpi_exchange_hydro(cpu,cpu->hsendbuffer,cpu->hrecvbuffer,0); // propagate hydro for refinement
#endif
#ifdef WRAD
      mpi_exchange_rad_level(cpu,cpu->Rsendbuffer,cpu->Rrecvbuffer,0,level); // propagate rad for refinement
#endif
#endif
      // refining (and destroying) octs
      curoct=L_refine_cells(level,param,firstoct,lastoct,cpu->freeoct,cpu,firstoct[0]+param->ngridmax,0.);
      cpu->freeoct=curoct;
    
      // ==================================== Check the number of particles and octs
      ptot[0]=0; for(ip=1;ip<=param->lmax;ip++){
	ptot[0]+=cpu->npart[ip-1]; // total of local particles
      }

#ifdef STARS
      ptot[1]=0; for(ip=1;ip<=param->lmax;ip++) ptot[1]+=cpu->nstar[ip-1];
#endif


      mtot=multicheck(firstoct,ptot,param->lcoarse,param->lmax,cpu->rank,cpu,param,1);

      ptot[0]=0; 
      for(ip=1;ip<=param->lmax;ip++){
	ptot[0]+=cpu->npart[ip-1]; // total of local particles
      }

#ifdef WMPI
      //reset the setup in case of refinement
      setup_mpi(cpu,firstoct,param->lmax,param->lcoarse,param->ngridmax,1); // out of WMPI to compute the hash table
      MPI_Barrier(cpu->comm);
#endif
    }
    
    // =============================== cleaning 
#ifdef WHYDRO2
    clean_new_hydro(level,param,firstoct,cpu);
#endif


#ifdef PIC
    L_clean_dens(level,param,firstoct,cpu);
#endif

#ifdef WRAD
    clean_new_rad(level,param,firstoct,cpu,0);
#endif

    
   // ================= III Recursive call to finer level

  double tt2,tt1;
    MPI_Barrier(cpu->comm);
    tt1=MPI_Wtime();

    int nlevel=cpu->noct[level]; // number at next level
#ifdef WMPI
    MPI_Allreduce(MPI_IN_PLACE,&nlevel,1,MPI_INT,MPI_SUM,cpu->comm);
#endif


    if(level<param->lmaxzoom){
      if(nlevel>0){
	zoom_level(level+1,cpu,param,firstoct,lastoct);
      }
    }

    MPI_Barrier(cpu->comm);
    tt2=MPI_Wtime();
    
    // ==================================== Check the number of particles and octs
    ptot[0]=0; for(ip=1;ip<=param->lmax;ip++) ptot[0]+=cpu->npart[ip-1]; // total of local particles
#ifdef STARS
    ptot[1]=0; for(ip=1;ip<=param->lmax;ip++) ptot[1]+=cpu->nstar[ip-1];
#endif
    mtot=multicheck(firstoct,ptot,param->lcoarse,param->lmax,cpu->rank,cpu,param,2);

    // ================= V Computing the new refinement map
    

    if((param->lmax!=param->lcoarse)&&(level<param->lmax)){
      // cleaning the marks
      L_clean_marks(level,firstoct);
      // marking the cells of the current level
      L_mark_cells(level,param,firstoct,1,param->amrthresh,cpu,cpu->sendbuffer,cpu->recvbuffer);
    }
  
  
  if(cpu->rank==RANK_DISP){
    printf("--\n");
    printf("exiting level =%d\n",level);
  }

}

#endif
