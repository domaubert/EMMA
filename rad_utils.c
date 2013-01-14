
#ifdef WRAD
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include "vector.h"

#include <mpi.h>

REAL Eddington(REAL fx, REAL fy, REAL fz, REAL ee, REAL c,int i,int j)
{
  REAL c2e=ee*c*c; // 2 flop
  REAL ff=0.;
  REAL arg,chi,res=0.;
  REAL n[3];
  n[0]=0.;n[1]=0.;n[2]=0.;

  if(ee>0)
    {
      ff=sqrtf(fx*fx+fy*fy+fz*fz); // 6 flop
      if(ff>0)
	{
	  n[0]=fx/ff; 
	  n[1]=fy/ff;
	  n[2]=fz/ff; 
	}
      ff=ff/(c*ee); // 2flop
    }
  
  arg=fmaxf(4.-3.*ff*ff,0.); // 4 flop
  chi=(3.+4.*ff*ff)/(5.+2.*sqrtf(arg)); // 7 flops

  if(i==j) res=(1.-chi)/2.*c2e; // 1 flops on average
  arg=(3.*chi-1.)/2.*c2e;
  res+=arg*n[i]*n[j];

  return res;
}



// =============================================================================================================

int radM_sweepX(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];
  struct Rtype Rold;

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fact;
  REAL fp,fm;
  REAL up,um;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
      
    for(i=0;i<nread;i++){ // we scan the octs
      
      memset(FL,0,sizeof(REAL)*NVAR_R);
      memset(FR,0,sizeof(REAL)*NVAR_R);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);
      
      /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);
      for(iface=0;iface<2;iface++){
	memcpy(RC+iface,curcell,sizeof(struct Rtype));
      }

      // Neighbor "MUSCL-LIKE" reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface;
	memcpy(RN+iface,stencil[i].oct[ioct[vnei[inei]]].cell[vell[icell]].rfield,sizeof(struct Rtype));

	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }




      // X DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems LEFT


      
      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RC[0].e[igrp];
	um=RN[0].e[igrp];

	fp=RC[0].fx[igrp];
	fm=RN[0].fx[igrp];

	FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RC[0].fx[igrp];
	um=RN[0].fx[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,0);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,0);

	FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RC[0].fy[igrp];
	um=RN[0].fy[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,0);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,0);

	FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RC[0].fz[igrp];
	um=RN[0].fz[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,0);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,0);

	FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      

	// ===========================================


      

      // --------- solving the Riemann Problems RIGHT

      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RN[1].e[igrp];
	um=RC[1].e[igrp];

	fp=RN[1].fx[igrp];
	fm=RC[1].fx[igrp];

	FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RN[1].fx[igrp];
	um=RC[1].fx[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,0);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,0);

	FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RN[1].fy[igrp];
	um=RC[1].fy[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,0);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,0);

	FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RN[1].fz[igrp];
	um=RC[1].fz[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,0);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,0);

	FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      
      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours


      for(igrp=0;igrp<NGRP;igrp++){
	for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0]; 
	for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1]; 
      }

      memcpy(stencil[i].New.cell[icell].rflux+0*NVAR_R*NGRP,FL,sizeof(REAL)*NVAR_R*NGRP);
      memcpy(stencil[i].New.cell[icell].rflux+1*NVAR_R*NGRP,FR,sizeof(REAL)*NVAR_R*NGRP);


      // ready for the next cell
    }
    //ready for the next oct
  }

  return 0;
}



// =============================================================================================================

int radM_sweepY(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];
  struct Rtype Rold;

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fact;
  REAL fp,fm;
  REAL up,um;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
      
    for(i=0;i<nread;i++){ // we scan the octs
      
      memset(FL,0,sizeof(REAL)*NVAR_R);
      memset(FR,0,sizeof(REAL)*NVAR_R);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);
      
      /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);
      for(iface=0;iface<2;iface++){
	memcpy(RC+iface,curcell,sizeof(struct Rtype));
      }

      // Neighbor "MUSCL-LIKE" reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+2;
	memcpy(RN+iface,stencil[i].oct[ioct[vnei[inei]]].cell[vell[icell]].rfield,sizeof(struct Rtype));

	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }




      // Y DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems FRONT


      
      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RC[0].e[igrp];
	um=RN[0].e[igrp];

	fp=RC[0].fx[igrp];
	fm=RN[0].fx[igrp];

	FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RC[0].fx[igrp];
	um=RN[0].fx[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,1);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,1);

	FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RC[0].fy[igrp];
	um=RN[0].fy[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,1);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,1);

	FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RC[0].fz[igrp];
	um=RN[0].fz[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,1);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,1);

	FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      

	// ===========================================


      

      // --------- solving the Riemann Problems BACK

      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RN[1].e[igrp];
	um=RC[1].e[igrp];

	fp=RN[1].fx[igrp];
	fm=RC[1].fx[igrp];

	FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RN[1].fx[igrp];
	um=RC[1].fx[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,1);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,1);

	FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RN[1].fy[igrp];
	um=RC[1].fy[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,1);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,1);

	FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RN[1].fz[igrp];
	um=RC[1].fz[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,1);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,1);

	FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      
      
      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours

      for(igrp=0;igrp<NGRP;igrp++){
	for(iface=0;iface<NVAR_R;iface++) FL[iface+igrp*NVAR_R]*=ffact[0]; 
	for(iface=0;iface<NVAR_R;iface++) FR[iface+igrp*NVAR_R]*=ffact[1]; 
      }

      memcpy(stencil[i].New.cell[icell].rflux+2*NVAR_R*NGRP,FL,sizeof(REAL)*NVAR_R*NGRP);
      memcpy(stencil[i].New.cell[icell].rflux+3*NVAR_R*NGRP,FR,sizeof(REAL)*NVAR_R*NGRP);


      // ready for the next cell
    }
    //ready for the next oct
  }

  return 0;
}


// ===================================================================================================

int radM_sweepZ(struct RGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt, REAL c){

  int inei,icell,iface;
  int i;
  int vnei[6],vcell[6];

  REAL FL[NVAR_R*NGRP],FR[NVAR_R*NGRP];
  struct Rtype Rold;

  struct Rtype RC[2];
  struct Rtype RN[2];

  int ioct[7]={0,1,2,3,4,5,6};

  struct Rtype *curcell;

  int ffact[2]={0,0};
  REAL fact;
  REAL fp,fm;
  REAL up,um;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
      
    for(i=0;i<nread;i++){ // we scan the octs
      
      memset(FL,0,sizeof(REAL)*NVAR_R);
      memset(FR,0,sizeof(REAL)*NVAR_R);

      // Getting the original state ===========================
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].rfield);
      
      /* // "MUSCL-LIKE" STATE RECONSTRUCTION */
      memset(ffact,0,sizeof(int)*2);
      for(iface=0;iface<2;iface++){
	memcpy(RC+iface,curcell,sizeof(struct Rtype));
      }

      // Neighbor "MUSCL-LIKE" reconstruction
      for(iface=0;iface<2;iface++){
	inei=iface+4;
	memcpy(RN+iface,stencil[i].oct[ioct[vnei[inei]]].cell[vell[icell]].rfield,sizeof(struct Rtype));

	if(!stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].split){
	  ffact[iface]=1; // we cancel the contriubtion of split neighbors
	}

      }




      // Y DIRECTION =========================================================================
      
      // --------- solving the Riemann Problems BOTTOM


      
      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RC[0].e[igrp];
	um=RN[0].e[igrp];

	fp=RC[0].fx[igrp];
	fm=RN[0].fx[igrp];

	FL[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RC[0].fx[igrp];
	um=RN[0].fx[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,0,2);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,0,2);

	FL[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RC[0].fy[igrp];
	um=RN[0].fy[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,1,2);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,1,2);

	FL[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RC[0].fz[igrp];
	um=RN[0].fz[igrp];

	fp=Eddington(RC[0].fx[igrp],RC[0].fy[igrp],RC[0].fz[igrp],RC[0].e[igrp],c,2,2);
	fm=Eddington(RN[0].fx[igrp],RN[0].fy[igrp],RN[0].fz[igrp],RN[0].e[igrp],c,2,2);

	FL[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      

	// ===========================================


      

      // --------- solving the Riemann Problems TOP

      for(igrp=0;igrp<NGRP;igrp++){
	
	// E
	up=RN[1].e[igrp];
	um=RC[1].e[igrp];

	fp=RN[1].fx[igrp];
	fm=RC[1].fx[igrp];

	FR[0+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);
	
	//FX

	up=RN[1].fx[igrp];
	um=RC[1].fx[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,0,2);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,0,2);

	FR[1+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FY

	up=RN[1].fy[igrp];
	um=RC[1].fy[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,1,2);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,1,2);

	FR[2+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

	//FX

	up=RN[1].fz[igrp];
	um=RC[1].fz[igrp];

	fp=Eddington(RN[1].fx[igrp],RN[1].fy[igrp],RN[1].fz[igrp],RN[1].e[igrp],c,2,2);
	fm=Eddington(RC[1].fx[igrp],RC[1].fy[igrp],RC[1].fz[igrp],RC[1].e[igrp],c,2,2);

	FR[3+igrp*NVAR_R]=0.5*(fp+fm)+0.5*c*(um-up);

      }
      
      
      //========================= copy the fluxes
      // Cancelling the fluxes from splitted neighbours

      for(iface=0;iface<NVAR;iface++) FL[iface]*=ffact[0]; 
      for(iface=0;iface<NVAR;iface++) FR[iface]*=ffact[1]; 

      memcpy(stencil[i].New.cell[icell].rflux+4*NVAR_R,FL,sizeof(REAL)*NVAR_R);
      memcpy(stencil[i].New.cell[icell].rflux+5*NVAR_R,FR,sizeof(REAL)*NVAR_R);

      // ready for the next cell
    }
    //ready for the next oct
  }

  return 0;
}



//==================================================================================
//==================================================================================

void recursive_neighbor_gather_oct(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct RGRID *stencil,char *visit){

  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;
  static int face[8]={0,1,2,3,4,5,6,7};
  REAL dxcur;

  struct Rtype Wi[8];
  char child[8];

  struct OCT *oct;
  struct OCT *neioct;
  struct CELL *neicell;

  
  if(cell->child!=NULL){
    // the oct at the right level exists
    neicell=cell->child->nei[ineiloc];
  }
  else{
    getcellnei(cell->idx, vnei, vcell); // we get the neighbors
    oct=cell2oct(cell);
    if(vnei[ineiloc]==6){
      neicell=&(oct->cell[vcell[ineiloc]]);
    }
    else{
      if(oct->nei[ineiloc]->child!=NULL){
	neicell=&(oct->nei[ineiloc]->child->cell[vcell[ineiloc]]);
      }
      else{
	printf("big problem\n");
	abort();
      }
    }
  }


  if(neicell->child!=NULL){
    // optimal case
    for(icell=0;icell<8;icell++){
      memcpy(Ri+icell,&(neicell->child->cell[icell].field),sizeof(struct Rtype));
      child[icell]=(neicell->child->cell[icell].child!=NULL);
    }

  }
  else{
    coarse2fine_rad2(neicell,Ri);
    for(icell=0;icell<8;icell++){
      child[icell]=0;
    }
  }




  for(icell=0;icell<8;icell++){
    memcpy(&(stencil->oct[ioct].cell[icell].rfield),Ri+face[icell],sizeof(struct Wtype)); //
    
    stencil->oct[ioct].cell[icell].split=child[icell];

  }
  
}

// ===================================================================================================
// ===================================================================================================

struct OCT *gatherstencilrad(struct OCT *octstart, struct RGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  char visit[27]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int ioct;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;

      //memset(visit,0,27*sizeof(char));
      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(&(stencil[iread].oct[6].cell[icell].rfield),&(curoct->cell[icell].rfield),sizeof(struct Rtype)); // for calculations
	stencil[iread].oct[6].cell[icell].split=(curoct->cell[icell].child!=NULL);
      }
      visit[13]=1;

      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  visit[ioct]=1;
	  recursive_neighbor_gather_oct_rad(ioct, inei, -1, -1, 1, cell, stencil+iread,visit);
	}
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}


// ===================================================================================================
// ===================================================================================================

void updatefieldrad(struct OCT *octstart, struct RGRID *stencil, int nread, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew)
{
  int i,icell,igrp;
  struct Rtype R;
  REAL one;
  int flx;
  REAL dtsurdx=dtnew/dxcur;
  REAL F[NFLUX_R];

  for(i=0;i<nread;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells
      
      if(stencil[i].oct[6].cell[icell].split) continue;
      memcpy(F,stencil[i].New.cell[icell].rflux,sizeof(REAL)*NFLUX_R);// New fluxes from the stencil
      
      
      // ==== updating

      // actually we compute and store the delta U only
      one=1.;
      memset(&R,0,sizeof(struct Rtype)); // setting delta U
      for(flx=0;flx<6;flx++){
	for(igrp=0;igrp<NGRP;igrp++){
	  R.e[igrp]  +=F[0+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fx[igrp] +=F[1+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fy[igrp] +=F[2+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;
	  R.fz[igrp] +=F[3+igrp*NVAR_R+flx*NVAR_R*NGRP]*dtsurdx*one;

	}
	one*=-1.;
      }
      // scatter back the delta Uwithin the stencil

      memcpy(&(stencil[i].New.cell[icell].deltaR),&R,sizeof(struct Rtype));
    }
  }

}

// ===========================================================================================================

struct OCT *scatterstencilrad(struct OCT *octstart, struct RGRID *stencil, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int ioct[7]={0,1,2,3,4,5,6};
  int vnei[6],vcell[6],inei;
  
  nextoct=octstart;
  iread=0;

  struct Rtype R,deltaR;
  REAL dtsurdx=dtnew/dxcur;


  int one;
  REAL F[NVAR_R*NGRP];

  //printf("let's scatter\n");
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	//we scatter the values in the central cell
	
	memcpy(&deltaR,&(stencil[iread].New.cell[icell].deltaR),sizeof(struct Rtype)); // getting the delta U back


	for(igrp=0;igrp<NGRP;igrp++){
	  R.e[igrp]   +=deltaR.e[igrp];
	  R.fx[igrp]  +=deltaR.fx[igrp];
	  R.fy[igrp]  +=deltaR.fy[igrp];
	  R.fz[igrp]  +=deltaR.fz[igrp];
	}
	
	memcpy(&(curoct->cell[icell].rfieldnew),&R,sizeof(struc Rtype)); // at this stage the central cell has been updated
	
	// let us now deal with coarser neighbors
	getcellnei(icell, vnei, vcell); // we get the neighbors
	
	for(inei=0;inei<6;inei++){


	  if(vnei[inei]!=6){
	    if(curoct->nei[vnei[inei]]->child==NULL){
	      // the neighbor cell is unsplit we update its value with fluxes
	   
	      // initial data from the new value
	      memcpy(&R,&(curoct->nei[vnei[inei]]->rfieldnew),sizeof(struct Rtype));

	      // getting the flux

	      /// ==== > A FIXER A PARTIR D'ICI AVEC IGRP SUR LES GROUPES  (NVAR ? NFLUX_R ?)

	      memcpy(F,stencil[iread].New.cell[icell].rflux+inei*NVAR,sizeof(REAL)*NVAR_R*NGRP);
	      
	      for(igrp=0;igrp<NGRP;igrp++){
		// update
		one=pow(-1.,inei+1);
		R.e[igrp] += F[0+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fx[igrp]+= F[1+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fy[igrp]+= F[2+igrp*NVAR_R]*dtsurdx*one*0.125;
		R.fz[igrp]+= F[3+igrp*NVAR_R]*dtsurdx*one*0.125;
	      }
 	      
	      memcpy(&(curoct->nei[vnei[inei]]->rfieldnew),&R,sizeof(struct Rtype));
	      
	    }
	  }
	}

      }

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}

// ====================================================================================================================

int advancerad(struct OCT **firstoct, int level, struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dxcur, REAL dtnew){

  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot,nread;
  double t[10];
  double tg=0.,th=0.,tu=0.,ts=0.,tfu=0.,ttot=0.;
  
  
  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  nreadtot=0;
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next; 

      t[0]=MPI_Wtime();
  
      // ------------ gathering the stencil value values
      nextoct= gatherstencilrad(curoct,stencil,stride,cpu, &nread);
      
      t[2]=MPI_Wtime();
      // ------------ solving the hydro
	    
      rad_sweepX(stencil,level,cpu->rank,nread,stride,dxcur,dtnew);   
      rad_sweepY(stencil,level,cpu->rank,nread,stride,dxcur,dtnew); 
      rad_sweepZ(stencil,level,cpu->rank,nread,stride,dxcur,dtnew); 
      
      // ------------ updating values within the stencil

      t[4]=MPI_Wtime();

      updatefieldrad(curoct,stencil,nread,stride,cpu,dxcur,dtnew);
      
      // ------------ scatter back the FLUXES
      
      t[6]=MPI_Wtime();
   
      nextoct=scatterstencilrad(curoct,stencil, nread, cpu,dxcur,dtnew);


      t[8]=MPI_Wtime();

      nreadtot+=nread;


      ts+=(t[8]-t[6]);
      tu+=(t[6]-t[4]);
      th+=(t[4]-t[2]);
      tg+=(t[2]-t[0]);

    }while(nextoct!=NULL);
  }
  
  printf("CPU | tgat=%e tcal=%e tup=%e tscat=%e\n",tg,th,tu,ts);

  return nreadtot;
}

// =================================================================================================
// =================================================================================================


void RadSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dtnew){
  
  int nread,nreadtot;;
  struct OCT *curoct;
  struct OCT *nextoct;
  
  REAL dxcur=pow(0.5,level);
  REAL one;
  struct Rtype R;
  int icell;
  int nocthydro=cpu->noct[level-1];
  double t[10];

  t[0]=MPI_Wtime();
  if(cpu->rank==0) printf("Start Radiation on %d octs with dt=%e on level %d with stride=%d\n",nocthydro,dtnew,level,stride);

  // ===== COMPUTING THE FLUXES
  
#ifndef GPUAXL
  nreadtot=advancerad(firstoct,level,cpu,stencil,stride,dxcur,dtnew);
#else
  nreadtot=advancerad(firstoct,level,cpu,stencil,stride,dxcur,dtnew);
  //  nreadtot=advanceradGPU(firstoct,level,cpu,stencil,stride,dxcur,dtnew);
#endif

  // FINAL UPDATE OF THE VALUES
  if(nreadtot>0){
    nextoct=firstoct[level-1];
    do {
      curoct=nextoct;
      nextoct=curoct->next; 
      for(icell=0;icell<8;icell++){
	if(curoct->cell[icell].child==NULL){
	  // unsplit case
	  memcpy(&(curoct->cell[icell].rfield),&(curoct->cell[icell].rfieldnew),sizeof(struct Rtype));
	}
	else{
	  // split case : -> quantities are averaged
	  struct OCT *child;
	  int i;
	  child=curoct->cell[icell].child;
	  memset(&R,0,sizeof(struct Rtype));
	  for(i=0;i<8;i++){
	    for(igrp=0;igrp<NGRP;igrp++){
	      R.e[igrp] +=child->cell[i].rfield.e[igrp]*0.125;
	      R.fx[igrp]+=child->cell[i].rfield.fx[igrp]*0.125;
	      R.fy[igrp]+=child->cell[i].rfield.fy[igrp]*0.125;
	      R.fz[igrp]+=child->cell[i].rfield.fz[igrp]*0.125;
	    }
	  }
	  memcpy(&curoct->cell[icell].rfield,&R,sizeof(struct Rtype));
	}
      }
    }while(nextoct!=NULL);
  }

  t[9]=MPI_Wtime();
  
#ifndef GPUAXL
  printf("==== CPU RAD TOTAL TIME =%e\n",t[9]-t[0]);
#else
  //printf(" === GPU RAD TOTAL TIME =%e\n",t[9]-t[0]);
  printf(" === GPU (DISABLED) RAD TOTAL TIME =%e\n",t[9]-t[0]);
#endif
  
}


// ==============================================================
void clean_new_rad(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu){

  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next; 

      //for(icell=0;icell<8;icell++) memset(&(curoct->cell[icell].fieldnew),0,sizeof(struct Wtype));
      for(icell=0;icell<8;icell++) memcpy(&(curoct->cell[icell].rfieldnew),&(curoct->cell[icell].rfield),sizeof(struct Rtype));
      
    }while(nextoct!=NULL);
  }
}




#endif
