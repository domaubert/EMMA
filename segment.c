
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include "hilbert.h"
#include "oct.h"
#include "cic.h"

#ifdef WMPI
#include <mpi.h>
#endif

#include "communication.h"

int segment_cell(struct OCT *curoct, int icell, struct CPUINFO *cpu, int levelcoarse)
{

  int res=-1;
  float xc,yc,zc;
  float xc0,yc0,zc0;
  int ix,iy,iz;
  int ii,jj,kk;
  int i,j,k;
  float dxcur;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  char first=1;
  int max,min;
  unsigned keyloc;

  // First we compute the current cell position (lower left corner)
  // Note : this is equivalenet to the position of the next level oct
  dxcur=pow(0.5,curoct->level); 
  xc0=curoct->x+( icell   %2)*dxcur; 
  yc0=curoct->y+((icell/2)%2)*dxcur;
  zc0=curoct->z+( icell/4   )*dxcur; 

  //if(curoct->level==1) printf("--\n");

  // we test the 27 cells centered around the current position
  for(ii=-1;ii<=1;ii++){
    for(jj=-1;jj<=1;jj++){
      for(kk=-1;kk<=1;kk++){
	// offset the current cell (for neighbor search)
	
	xc=xc0+ii*dxcur;
	yc=yc0+jj*dxcur;
	zc=zc0+kk*dxcur;

	// Periodic boundary conditions
	
	if(xc<0.){
	  xc+=1.;
	}
	else if(xc>=1.){
	  xc-=1.;
	}

	if(yc<0.){
	  yc+=1.;
	}
	else if(yc>=1.){
	  yc-=1.;
	}

	if(zc<0.){
	  zc+=1.;
	}
	else if(zc>=1.){
	  zc-=1.;
	}
	
	
	// we convert it in dxcoarse unit
	
	ix=(int)(xc/pow(0.5,levelcoarse-1));
	iy=(int)(yc/pow(0.5,levelcoarse-1));
	iz=(int)(zc/pow(0.5,levelcoarse-1));
	
	//if((curoct->level==1)) printf("--------------- %d %d %d\n",ix,iy,iz);
	// we compute the keys of the 8 corners and extract the min/max key
	first=1;
	for(i=0;i<2;i++){
	  for(j=0;j<2;j++){
	    for(k=0;k<2;k++){
	      
	      c[0]=ix+(int)(i*pow(0.5,curoct->level-levelcoarse+1))-(i==1);
	      c[1]=iy+(int)(j*pow(0.5,curoct->level-levelcoarse+1))-(j==1);
	      c[2]=iz+(int)(k*pow(0.5,curoct->level-levelcoarse+1))-(k==1);
	      keyloc=(unsigned)(hilbert_c2i(3,levelcoarse,c));
	      
	      //if((curoct->level==1)) printf("i=%d j=%d k=%d|| %d %d %d || %d || %d %d\n",i,j,k,(int)(c[0]),(int)(c[1]),(int)(c[2]),keyloc,cpu->kmin,cpu->kmax);
	      if(first){
		min=keyloc;
		max=keyloc;
		first=0;
	      }
	      else{
		min=(min>keyloc?keyloc:min);
		max=(max<keyloc?keyloc:max);
	      }
	    }
	  }
	}
	
	// we check if these keys intersect the cpu domain
	if((max<cpu->kmin)||(min>cpu->kmax)){
	  res=0;
	}
	else
	  {
	      res=1;
	  }
	// No need to go furhter if the intersection exists
	if(res==1) break;
      }
      if(res==1) break;
    }
    if(res==1) break;
  }

  return res;

}


//------------------------------------------------------------------------

unsigned oct2key(struct OCT *curoct,int level){
  
  float xc,yc,zc;
  int ix,iy,iz;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  unsigned keyloc;

  xc=curoct->x;
  yc=curoct->y;
  zc=curoct->z;
	
  // we convert it in dxcoarse octs unit
  
  ix=(int)(xc/pow(0.5,level-1));
  iy=(int)(yc/pow(0.5,level-1));
  iz=(int)(zc/pow(0.5,level-1));
	
  // we compute the keys of the lowest corners 
  c[0]=ix;
  c[1]=iy;
  c[2]=iz;
  keyloc=(unsigned)(hilbert_c2i(3,level,c));

  return keyloc;

}


//------------------------------------------------------------------------
void assigncpu2coarseoct(struct OCT *curoct, struct CPUINFO *cpu, int levelcoarse)
{

  float xc,yc,zc;
  int ix,iy,iz;
  float dxcur;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  unsigned keyloc;
  int cpuloc;

  xc=curoct->x;
  yc=curoct->y;
  zc=curoct->z;
	
  // we convert it in dxcoarse octs unit
  
  ix=(int)(xc/pow(0.5,levelcoarse-1));
  iy=(int)(yc/pow(0.5,levelcoarse-1));
  iz=(int)(zc/pow(0.5,levelcoarse-1));
	
  // we compute the keys of the lowest corners 
  c[0]=ix;
  c[1]=iy;
  c[2]=iz;
  keyloc=(unsigned)(hilbert_c2i(3,levelcoarse,c));
  cpuloc=keyloc/cpu->nkeys;
  curoct->cpu=(cpuloc>(cpu->nproc-1)?cpu->nproc-1:cpuloc);
}

//------------------------------------------------------------------------

int segment_part(float xc,float yc,float zc, struct CPUINFO *cpu, int levelcoarse)
{

  int ix,iy,iz;
  float dxcur;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  unsigned keyloc;
  int cpuloc;
	
  // we convert it in dxcoarse octs unit
  
  ix=(int)(xc/pow(0.5,levelcoarse-1));
  iy=(int)(yc/pow(0.5,levelcoarse-1));
  iz=(int)(zc/pow(0.5,levelcoarse-1));
	
  // we compute the keys of the lowest corners 
  c[0]=ix;
  c[1]=iy;
  c[2]=iz;
  keyloc=(unsigned)(hilbert_c2i(3,levelcoarse,c));

  if((keyloc>=cpu->kmin)&&(keyloc<=cpu->kmax)){
    return 1;
  }
  else{
    return 0;
  }

}

 //------------------------------------------------------------------------

 //------------------------------------------------------------------------


 //------------------------------------------------------------------------

 // the hash function
int hfun(unsigned key){
  return key>>6;
}


 //------------------------------------------------------------------------
void load_balance(int levelcoarse,struct CPUINFO *cpu){

  int keymax=pow(2,3*(levelcoarse-1))-1; // the maximal key along the Hilbert curve
  
  cpu->kmin=((keymax+1)/cpu->nproc)*cpu->rank; // the min key of the current cpu
  cpu->nkeys=((keymax+1)/cpu->nproc); // the number of keys per cpu

  if(cpu->rank!=(cpu->nproc-1)){
    cpu->kmax=((keymax+1)/cpu->nproc)*(cpu->rank+1)-1; // the max key of the current cpu
  }
  else{
    cpu->kmax=keymax; // the last proc should go until the end of the chain
  }
    
  printf("proc %d cpu min=%d cpu max=%d delta=%d\n",cpu->rank,cpu->kmin,cpu->kmax,(keymax+1)/cpu->nproc);
}
