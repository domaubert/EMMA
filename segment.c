
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include "hilbert.h"
#include "oct.h"
#include "cic.h"
#include "particle.h"

#ifdef WMPI
#include <mpi.h>
#endif

#include "communication.h"
#define LCO 0


int segment_cell(struct OCT *curoct, int icell, struct CPUINFO *cpu, int levelcoarse)
{

  int res=-1;
  REAL xc,yc,zc;
  REAL xc0,yc0,zc0;
  int ix,iy,iz;
  int ii,jj,kk;
  int i,j,k;
  REAL dxcur;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  char first=1;
  int max,min;
  unsigned long long keyloc;


  // First we compute the current cell position (lower left corner)
  // Note : this is equivalenet to the position of the next level oct
  dxcur=POW(0.5,curoct->level); 
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
	
	ix=(int)(xc/POW(0.5,levelcoarse-1));
	iy=(int)(yc/POW(0.5,levelcoarse-1));
	iz=(int)(zc/POW(0.5,levelcoarse-1));
	
	//if((curoct->level==1)) printf("--------------- %d %d %d\n",ix,iy,iz);
	// we compute the keys of the 8 corners and extract the min/max key
	first=1;
	for(i=0;i<2;i++){
	  for(j=0;j<2;j++){
	    for(k=0;k<2;k++){
	      
	      c[0]=ix+(int)(i*POW(0.5,curoct->level-levelcoarse+1))-(i==1);
	      c[1]=iy+(int)(j*POW(0.5,curoct->level-levelcoarse+1))-(j==1);
	      c[2]=iz+(int)(k*POW(0.5,curoct->level-levelcoarse+1))-(k==1);
#ifndef TUBE
	      keyloc=(unsigned long long)(hilbert_c2i(3,levelcoarse,c));
#else
	      keyloc=(unsigned long long)(c[0]+c[1]*POW(2,levelcoarse-1)+c[2]*POW(2,levelcoarse-1)*POW(2,levelcoarse-1));
#endif	      
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

unsigned long long oct2key(struct OCT *curoct,int level){
  
  REAL xc,yc,zc;
  int ix,iy,iz;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  unsigned long long keyloc;

  xc=curoct->x;
  yc=curoct->y;
  zc=curoct->z;
	
  // we convert it in dxcoarse octs unit
  
  ix=(int)(xc/POW(0.5,level-1));
  iy=(int)(yc/POW(0.5,level-1));
  iz=(int)(zc/POW(0.5,level-1));
	
  // we compute the keys of the lowest corners 
  c[0]=ix;
  c[1]=iy;
  c[2]=iz;
#ifndef TUBE
  keyloc=(unsigned long long)(hilbert_c2i(3,level,c));
#else
  keyloc=(unsigned long long)(ix+iy*POW(2,level-1)+iz*POW(2,level-1)*POW(2,level-1));
#endif


  return keyloc;

}

//------------------------------------------------------------------------
unsigned long long pos2key(REAL xc, REAL yc, REAL zc,int level){
  
  int ix,iy,iz;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  unsigned long long keyloc;

	
  // we convert it in dxcoarse octs unit
  
  ix=(int)(xc/POW(0.5,level-1));
  iy=(int)(yc/POW(0.5,level-1));
  iz=(int)(zc/POW(0.5,level-1));
	
  // we compute the keys of the lowest corners 
  c[0]=ix;
  c[1]=iy;
  c[2]=iz;
#ifndef TUBE
  keyloc=(unsigned long long)(hilbert_c2i(3,level,c));
#else
  keyloc=(unsigned long long)(ix+iy*POW(2,level-1)+iz*POW(2,level-1)*POW(2,level-1));
#endif

  return keyloc;

}


//------------------------------------------------------------------------
void assigncpu2coarseoct(struct OCT *curoct, struct CPUINFO *cpu, int levelcoarse)
{

  REAL xc,yc,zc;
  int ix,iy,iz;
  REAL dxcur;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  unsigned long long keyloc;
  int cpuloc;

  if(curoct->level<=LCO){
    curoct->cpu=0;
  }
  else{

    xc=curoct->x;
    yc=curoct->y;
    zc=curoct->z;
	
    // we convert it in dxcoarse octs unit
  
    ix=(int)(xc/POW(0.5,levelcoarse-1));
    iy=(int)(yc/POW(0.5,levelcoarse-1));
    iz=(int)(zc/POW(0.5,levelcoarse-1));
	
    // we compute the keys of the lowest corners 
    c[0]=ix;
    c[1]=iy;
    c[2]=iz;
#ifndef TUBE
    keyloc=(unsigned long long)(hilbert_c2i(3,levelcoarse,c));
#else
    keyloc=(unsigned long long)(ix+iy*POW(2,levelcoarse-1)+iz*POW(2,levelcoarse-1)*POW(2,levelcoarse-1));
#endif
    cpuloc=keyloc/cpu->nkeys;
    curoct->cpu=(cpuloc>(cpu->nproc-1)?cpu->nproc-1:cpuloc);

    /* if(xc==0.53125) */
    /* if(yc==0.53125) */
    /*   if(zc==0.53125){ */
    /* 	printf("TESTO rank=%d\n",curoct->cpu); */
    /* 	abort(); */
    /*   } */
  }
}

//------------------------------------------------------------------------

int segment_part(REAL xc,REAL yc,REAL zc, struct CPUINFO *cpu, int levelcoarse)
{

  int ix,iy,iz;
  REAL dxcur;
  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
  unsigned long long keyloc;
  int cpuloc;
	
  // we convert it in dxcoarse octs unit
  
  ix=(int)(xc/POW(0.5,levelcoarse-1));
  iy=(int)(yc/POW(0.5,levelcoarse-1));
  iz=(int)(zc/POW(0.5,levelcoarse-1));
	
  // we compute the keys of the lowest corners 
  c[0]=ix;
  c[1]=iy;
  c[2]=iz;
#ifndef TUBE
  keyloc=(unsigned long long)(hilbert_c2i(3,levelcoarse,c));
#else
  keyloc=(unsigned long long)(ix+iy*POW(2,levelcoarse-1)+iz*POW(2,levelcoarse-1)*POW(2,levelcoarse-1));
#endif
  if((keyloc>=cpu->kmin)&&(keyloc<=cpu->kmax)){
    return 1;
  }
  else{
    return 0;
  }

}

 //------------------------------------------------------------------------


 //------------------------------------------------------------------------

 // the hash function
unsigned long long hfun(unsigned long long key, unsigned long long maxval){
  // warnign maxval must be a power of two
  //return key>>6;
  return key&(maxval-1);
}

 //------------------------------------------------------------------------

#ifdef PIC
void part2grid(struct PART *part, struct CPUINFO *cpu,int npart){

  int i;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct PART *curp;

  for(i=0;i<npart;i++){
    unsigned long long key;
    unsigned long long hidx;
    
    int found=0;
    key=pos2key(part[i].x,part[i].y,part[i].z,part[i].level);
    hidx=hfun(key,cpu->maxhash);
    nextoct=cpu->htable[hidx];
    if(nextoct!=NULL){
      do{ // resolving collisions
	curoct=nextoct;
	nextoct=curoct->nexthash;
	found=((oct2key(curoct,curoct->level)==key)&&(part[i].level==curoct->level));
      }while((nextoct!=NULL)&&(!found));
    }
    
    if(found){ // the reception oct has been found
      
      REAL dxcur=POW(0.5,part[i].level);
      
      int xp=(int)((part[i].x-curoct->x)/dxcur);//xp=(xp>1?1:xp);xp=(xp<0?0:xp);
      int yp=(int)((part[i].y-curoct->y)/dxcur);//yp=(yp>1?1:yp);yp=(yp<0?0:yp);
      int zp=(int)((part[i].z-curoct->z)/dxcur);//zp=(zp>1?1:zp);zp=(zp<0?0:zp);
      int ip=xp+yp*2+zp*4;
      
      
      // the reception cell 
      struct CELL *newcell=&(curoct->cell[ip]);
      struct OCT *newoct;
      REAL dxcur2;

      // if refined we assign the particle to a refined cell
      if(newcell->child!=NULL){
	newoct=newcell->child;
	dxcur2=1./POW(2.,newoct->level);

	// new cell coordinates
	xp=(int)(DFACT*(part[i].x-newoct->x)/dxcur2);

	/* if(xp>1){ */
	/*   printf("WEIRDA !!\n"); */
	/*   printf("px=%e ox=%e dxcur2=%e xp=%d\n",part[i].x,newoct->x,dxcur2,xp); */
	/*   abort(); */
	/* } */
	xp=(xp>1?1:xp);
	xp=(xp<0?0:xp);
	yp=(int)(DFACT*(part[i].y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
	zp=(int)(DFACT*(part[i].z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
	ip=xp+yp*2+zp*4;
	
	newcell=&(newoct->cell[ip]);

      }
      // ready to assign part to cell

      curp=findlastpart(newcell->phead);
      if(curp!=NULL){
	curp->next=part+i;
	part[i].next=NULL;
	part[i].prev=curp;
      }
      else{
	newcell->phead=part+i;
	part[i].next=NULL;
	part[i].prev=NULL;
      }
    }
  }
}
#endif

 //------------------------------------------------------------------------
void load_balance(int levelcoarse,struct CPUINFO *cpu){

  unsigned long long keymax=POW(2,3*(levelcoarse-1))-1; // the maximal key along the Hilbert curve
  
  cpu->kmin=((keymax+1)/cpu->nproc)*cpu->rank; // the min key of the current cpu
  cpu->nkeys=((keymax+1)/cpu->nproc); // the number of keys per cpu

  if(cpu->rank!=(cpu->nproc-1)){
    cpu->kmax=((keymax+1)/cpu->nproc)*(cpu->rank+1)-1; // the max key of the current cpu
  }
  else{
    cpu->kmax=keymax; // the last proc should go until the end of the chain
  }
    
  //printf("proc %d cpu min=%d cpu max=%d delta=%d\n",cpu->rank,cpu->kmin,cpu->kmax,(keymax+1)/cpu->nproc);
}
