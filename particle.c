#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include "oct.h"
#include "cic.h"
#include "vector.h"

#ifdef WMPI
#include <mpi.h>
#include "communication.h"
#endif


#ifdef PIC
//------------------------------------------------------------------------

struct PART* findlastpart(struct PART* phead)
{
  struct PART* curp;
  struct PART* nexp;

  curp=NULL;
  nexp=phead; //sweeping the particles of the current cell */
  if(nexp!=NULL){ 
    do{  
      curp=nexp; 
      nexp=curp->next; 
    }while(nexp!=NULL); 
  }
  return curp;
}


  
//------------------------------------------------------------------------

int countpart(struct PART* phead)
{
  struct PART* curp;
  struct PART* nexp;
  int npart=0;

  curp=NULL;
  nexp=phead; //sweeping the particles of the current cell */
  if(nexp!=NULL){ 
    do{  
      curp=nexp; 
      nexp=curp->next; 
      npart++;
    }while(nexp!=NULL); 
  }

  return npart;
}
  // ==========================================

int countpartDM(struct CELL* cell, int *npart)
{
  struct PART* curp;
  struct PART* nexp;
  
  if(cell->child==NULL){ // tree leaf -> explore particle
    curp=NULL;
    nexp=cell->phead; //sweeping the particles of the current cell */
    if(nexp!=NULL){ 
      do{  
	curp=nexp; 
	nexp=curp->next; 
#ifdef STARS
	(*npart)+=(curp->isStar!=1);
#else
	(*npart)++;
#endif
      }while(nexp!=NULL); 
    }
  }
  else{
    int icell;
    for(icell=0;icell<8;icell++){
      countpartDM(&cell->child->cell[icell],npart);
    }
  }
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

struct PART* modifpospart(struct PART* phead, REAL len, int dir)
{
  struct PART* curp;
  struct PART* nexp;

  curp=NULL;
  nexp=phead; //sweeping the particles of the current cell */
  if(nexp!=NULL){ 
    do{  
      curp=nexp; 
      nexp=curp->next; 
      switch(dir){
      case 0:curp->x+=len;break;
      case 1:curp->y+=len;break;
      case 2:curp->z+=len;break;
      }
    }while(nexp!=NULL); 
  }
  return curp;
}


//=====================================================================
//=====================================================================

REAL comptstep(int levelcoarse,int levelmax,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax){
  
  int level;
  REAL vmax,lmdisp;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  int icell;
  struct PART *nexp;
  struct PART *curp;
  REAL va;
  REAL aa;
  REAL dtlev,dtnew=1e9;
  REAL dt;

  // Computing new timestep
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      dtlev=0.;
      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      if(nextoct==NULL) continue; // in case the level is empty
      do // sweeping through the octs of level
	{
	  oct=(*nextoct);
	  dxcur=1./POW(2,oct.level+1); // +1 to protect level change
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {

#ifdef PIC
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){

		aa=SQRT(oct.cell[icell].f[0]*oct.cell[icell].f[0]+oct.cell[icell].f[1]*oct.cell[icell].f[1]+oct.cell[icell].f[2]*oct.cell[icell].f[2])*fa2;
		
		do{ 
		  curp=nexp; 
		  nexp=curp->next; 

 /* #ifdef PART2 */
/* 		  if(curp->idx==0) continue; */
/* #endif */
		  // particle velocit
		  va=SQRT(curp->vx*curp->vx+curp->vy*curp->vy+curp->vz*curp->vz)*fa;

		  // loc dt
		  if((va>0.)&&(aa>0.)){
		    REAL EPS=2.0*(FRACDX*dxcur)*aa/(va*va);
		    if(EPS<1e-1){
		      dtlev=(FRACDX*dxcur)/va;
		    }
		    else if(EPS>10.){
		      dtlev=SQRT(2.0*(FRACDX*dxcur)/aa);
		    }
		    else{
		      dtlev=va/aa*(SQRT(1.+2.*aa*(FRACDX*dxcur)/(va*va))-1.);
		    }
		  }		 
		  else if((aa==0.)||(va==0.)){
		    dtlev=1e9;
		  } 
#ifdef PART2
		  printf("mass=%e aa=%e va=%e dt=%e tmax=%e\n",curp->mass,aa,va,dtlev,tmax);
#endif
		  if(dtnew>dtlev) dtnew=dtlev;
		}while(nexp!=NULL);
	      }
#endif
	    }
	}while(nextoct!=NULL);

      /* if(vmax>0.){ */
      /* 	dtlev=FRACDX*dxcur/vmax;///fa; */
      /* 	dtnew=(dtlev<dtnew?dtlev:dtnew); */
      /* } */
      /* else{ */
      /* 	dtlev=1e9; */
      /* 	dtnew=(dtlev<dtnew?dtlev:dtnew); */
      /* } */
      //printf("level=%d vmax=%e dtlev=%e\n",level,vmax,dtlev);

    }

  // new tstep
  //printf("original dt=%f chosen dt=%f\n",dt,dtnew);
  dt=dtnew;


#ifdef WMPI
  // reducing by taking the smallest time step
  MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_REEL,MPI_MIN,cpu->comm);
#endif  

  dt=(dt>tmax?tmax:dt);
  return dt;
}


// =================================================================================================================
// =================================================================================================================

REAL L_comptstep(int level,struct RUNPARAMS *param,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax){
  
  REAL vmax,lmdisp;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  int icell;
  struct PART *nexp;
  struct PART *curp;
  REAL va;
  REAL aa;
  REAL dtlev,dtnew=1e9;
  REAL dt;

  // Computing new timestep
  dtlev=0.;
  // setting the first oct
  
  nextoct=firstoct[level-1];
  
  if(nextoct!=NULL){
    do // sweeping through the octs of level
      {
	oct=(*nextoct);
	dxcur=1./POW(2,oct.level+1); // +1 to protect level change
	nextoct=oct.next;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
#ifdef PIC
	    nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	    if(nexp!=NULL){
	    
	      aa=SQRT(oct.cell[icell].f[0]*oct.cell[icell].f[0]+oct.cell[icell].f[1]*oct.cell[icell].f[1]+oct.cell[icell].f[2]*oct.cell[icell].f[2])*fa2;
	      do{ 
		curp=nexp; 
		nexp=curp->next; 
	      
		// particle velocity
		va=SQRT(curp->vx*curp->vx+curp->vy*curp->vy+curp->vz*curp->vz)*fa;
		if(va!=0){
		  dtlev=(FRACDX*dxcur)/va;
		}		 
		else{
		  dtlev=1e9;
		} 
#ifdef PART2
		printf("mass=%e aa=%e va=%e dt=%e tmax=%e\n",curp->mass,aa,va,dtlev,tmax);
#endif
		if(dtnew>dtlev) dtnew=dtlev;
	      }while(nexp!=NULL);
	    }
#endif
	  }
      }while(nextoct!=NULL);
  }  
  // new tstep
  
  dt=dtnew;


/* #ifdef WMPI */
/*   // reducing by taking the smallest time step */
/*   MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_REEL,MPI_MIN,cpu->comm); */
/* #endif   */

  dt=(dt>tmax?tmax:dt);
  return dt;
}


// ================================================================================================

REAL L_movepart(int level,struct OCT** firstoct, REAL*adt, int is, struct CPUINFO* cpu){
  
  REAL mdisp,lmdisp;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  int icell;
  struct PART *nexp;
  struct PART *curp;
  REAL disp;
  REAL dt;
  // === Moving particles

  // setting the first oct
  mdisp=0.;
  nextoct=firstoct[level-1];
  
  dxcur=1./POW(2,level);

  if(nextoct!=NULL){
  do // sweeping through the octs of level
    {
      oct=(*nextoct);
      nextoct=oct.next;

      //if(oct.cpu!=cpu->rank) continue;
      for(icell=0;icell<8;icell++) // looping over cells in oct
	{
	  nexp=oct.cell[icell].phead; //sweeping the particles of the current cell

	  if(nexp!=NULL){
	    do{ 
	      curp=nexp; 
	      nexp=curp->next; 
	      
		// particle displacement
	      if((curp->is==0)||(curp->level==level)){
		dt=adt[curp->level-1];
		disp=SQRT(curp->vx*dt*curp->vx*dt+curp->vy*dt*curp->vy*dt+curp->vz*dt*curp->vz*dt);
		curp->x+=curp->vx*dt;
		curp->y+=curp->vy*dt;
		curp->z+=curp->vz*dt;
		if(disp>mdisp){
		  mdisp=disp;
		  //printf("%e %e %e %e\n",curp->vx,curp->vy,curp->vz,disp);
		}
	      }
	    }while(nexp!=NULL);
	  }
	}
    }while(nextoct!=NULL);
  }

#ifdef WMPI
  REAL mmdisp;
  MPI_Allreduce(&mdisp,&mmdisp,1,MPI_REEL,MPI_MAX,cpu->comm);
  //  mdisp=mmdisp;
#endif

  if(cpu->rank==RANK_DISP) printf("level=%d maxdisp=%e or %e dx\n",level,mmdisp,mmdisp/dxcur);

  return dt;
}

//===============================================
//===============================================

void L_levpart(int level,struct OCT** firstoct, int is){
  
  struct OCT *nextoct;
  struct OCT oct;
  int icell;
  struct PART *nexp;
  struct PART *curp;

  
  // setting the first oct
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
  do // sweeping through the octs of level
    {
      oct=(*nextoct);
      nextoct=oct.next;

      for(icell=0;icell<8;icell++) // looping over cells in oct
	{
	  nexp=oct.cell[icell].phead; //sweeping the particles of the current cell

	  if(nexp!=NULL){
	    do{ 
	      curp=nexp; 
	      nexp=curp->next; 
	      curp->level=level;
	      curp->is=is;
	    }while(nexp!=NULL);
	  }
	}
    }while(nextoct!=NULL);
  }
}





void L_reset_is_part(int level,struct OCT** firstoct){
  
  struct OCT *nextoct;
  struct OCT oct;
  int icell;
  struct PART *nexp;
  struct PART *curp;

  
  // setting the first oct
  nextoct=firstoct[level-1];
  
  if(nextoct!=NULL){
  do // sweeping through the octs of level
    {
      oct=(*nextoct);
      nextoct=oct.next;

      for(icell=0;icell<8;icell++) // looping over cells in oct
	{
	  nexp=oct.cell[icell].phead; //sweeping the particles of the current cell

	  if(nexp!=NULL){
	    do{ 
	      curp=nexp; 
	      nexp=curp->next; 
	      
	      curp->is=0;
	    }while(nexp!=NULL);
	  }
	}
    }while(nextoct!=NULL);
  }
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void partcellreorg(int levelcoarse,int levelmax,struct OCT **firstoct){

  int dir;
  char out;
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *newoct;
  REAL dxcur,dxcur2;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  struct PART *part;
  REAL xc,yc,zc;
  int xp,yp,zp;
  int vnei[6],vcell[6];
  struct CELL *newcell;
  int ip;

  for(dir=0;dir<3;dir++) 
    { 
      for(level=levelcoarse;level<=levelmax;level++) // looping over levels
	{
	  //printf("dir=%d level=%d\n",dir,level);
	  // setting the first oct
	
	  nextoct=firstoct[level-1];
	  if(nextoct==NULL) continue; // in case the level is empty

	  do{ // sweeping through the octs of level
	    curoct=nextoct;
	    dxcur=1./POW(2,curoct->level);
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
		{
		  nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell
		  newoct=NULL;
		  newcell=NULL;
		  // A few informations about the current cell
		  getcellnei(icell, vnei, vcell);
		
		  // sweeping the particles
		  if(nexp!=NULL){
		    do{ 
		      curp=nexp; 
		      nexp=curp->next; 
		      out=0;

		      switch(dir){

		      case 0: // ======================================   x displacement============

			xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
			out=((curp->x-xc)>dxcur*0.5)-((curp->x-xc)<-dxcur*0.5);
			if(out==1){
			  if(vnei[1]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[1]]);
			    
			    // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
				
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 
				
				newcell=&(newoct->cell[ip]);

			      }

			  }
			  else{ // the particle drifts to a new oct
			    // periodic boundaries
			    curp->x+=(curp->x<0.)-(curp->x>1.);

			    if(curoct->nei[vnei[1]]->child==NULL){ // the particle will go to level-1
			      newcell=curoct->nei[vnei[1]];
			    }
			    else{ // the particle will remain at the same level or more
			      newcell=&(curoct->nei[vnei[1]]->child->cell[vcell[1]]);
			    
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
				
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);

			      }
			    }
			  }
			}
			else if(out==-1){
			  if(vnei[0]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[0]]);

			    // if refined we assign the particle to a refined cell
			    if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
			      

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);
			      }
			  }
			  else{
			    // periodic boundaries
			    curp->x+=(curp->x<0.)-(curp->x>1.);
			    
			    if(curoct->nei[vnei[0]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[0]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[0]]->child->cell[vcell[0]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);
				
				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
			      

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);
			      }

			    }
			  }
			}
			break;

		    case 1: // ======================================   y displacement============
		      yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
		
		      out=((curp->y-yc)>dxcur*0.5)-((curp->y-yc)<-dxcur*0.5);
		      // getting the new cell // TO FIX PERIODICITY IS ILL DEFINED HERE
		      if(out==1){
			  if(vnei[3]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[3]]);
			    // if refined we assign the particle to a refined cell at level+1
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				//if((ip>7)||(ip<0)) ip=0; 
				newcell=&(newoct->cell[ip]);
			      }

			  }
			  else{ // the particle moves to another oct
			    // periodic boundaries
			    curp->y+=(curp->y<0.)-(curp->y>1.);
			    
			    if(curoct->nei[vnei[3]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[3]]);
			      //if(curp->idx==222044) printf("coucou hello");
			    }
			    else{ // the particle will remain at the same level or level+1
			      newcell=&(curoct->nei[vnei[3]]->child->cell[vcell[3]]);
			      
			     
			      // if refined we assign the particle to a refined cell at level+1
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				//if((ip>7)||(ip<0)) ip=0; 
				newcell=&(newoct->cell[ip]);
			      }
			    }
			  }
			}
			else if(out==-1){
			  if(vnei[2]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[2]]);
			    // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);
				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 
				newcell=&(newoct->cell[ip]);
			      }

			  }
			  else{
			    // periodic boundaries
			    curp->y+=(curp->y<0.)-(curp->y>1.);

			    if(curoct->nei[vnei[2]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[2]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[2]]->child->cell[vcell[2]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);
				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 
				newcell=&(newoct->cell[ip]);
			      }

			    }
			  }
			}

			break;

		      case 2: // ======================================   z displacement============
			zc=curoct->z+(icell/4)*dxcur+dxcur/2;
			
			out=((curp->z-zc)>dxcur*0.5)-((curp->z-zc)<-dxcur*0.5);
			// getting the new cell
			if(out==1){
			  if(vnei[5]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[5]]);
			      // if refined we assign the particle to a refined cell
			    if(newcell->child!=NULL){
			      //if(curp->idx==82279) printf("to level +1");
			      
			      newoct=newcell->child;
			      dxcur2=1./POW(2.,newoct->level);
			      
			      xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);
			      //printf("xp=%d\n",xp);
			      xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			      yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			      zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			      ip=xp+yp*2+zp*4;


			      
			      //ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
			      // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			      //if((ip>7)||(ip<0)) ip=0; 
			      
			      newcell=&(newoct->cell[ip]);
			    }
			  }
			  else{
			    // periodic boundaries
			    curp->z+=(curp->z<0.)-(curp->z>1.);

			    if(curoct->nei[vnei[5]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[5]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[5]]->child->cell[vcell[5]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				//if(curp->idx==82279) printf("to level +1");

				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);
				xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);
				yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);
				zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;

				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 

				

				newcell=&(newoct->cell[ip]);
			      }

			    }
			  }
			}
			else if(out==-1){
			  if(vnei[4]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[4]]);
			    // if refined we assign the particle to a refined cell
			    if(newcell->child!=NULL){
			      newoct=newcell->child;
			      dxcur2=1./POW(2.,newoct->level);
			      
			      xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			      yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			      zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			      ip=xp+yp*2+zp*4;
			      
			      //				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
			      // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			      //if((ip>7)||(ip<0)) ip=0; 
			      
			      newcell=&(newoct->cell[ip]);
			    }

			  }
			  else{
			    // periodic boundaries
			    curp->z+=(curp->z<0.)-(curp->z>1.);

			    if(curoct->nei[vnei[4]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[4]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[4]]->child->cell[vcell[4]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);
				
				xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;

				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//				if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);
			      }
			    }
			  }
			}
			break;
		      }
		      
		      if(out!=0){
			// ========= assigning the particle to the newcell + pointer management
			
			

			// removing the cell from its old list
			if(curp->prev !=NULL){
			  curp->prev->next=curp->next;
			}
			else{
			  curoct->cell[icell].phead=curp->next; 
			}
		      
			if(curp->next!=NULL){
			  curp->next->prev=curp->prev;
			}
		      
			// update the linking of the current particle
		      
			// adding the particle to the new cell
			if(newcell->child!=NULL){
			  printf("erro in part displacement\n");
			  abort();
			}


			if(newcell->phead!=NULL){
			  part=findlastpart(newcell->phead);
			  part->next=curp;
			  curp->prev=part;
			}
			else{
			  newcell->phead=curp;
			  curp->prev=NULL;
			}
		      
			curp->next=NULL;
		      }

		    }while(nexp!=NULL);
		  }
		}
	  }while(nextoct!=NULL);
	}	
      
    }

}    

// =======================================================
// =======================================================

void L_partcellreorg(int level,struct OCT **firstoct){

  int dir;
  char out;
 
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *newoct;
  REAL dxcur,dxcur2;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  struct PART *part;
  REAL xc,yc,zc;
  int xp,yp,zp;
  int vnei[6],vcell[6];
  struct CELL *newcell;
  int ip;

  for(dir=0;dir<3;dir++) 
    { 
	  //printf("dir=%d level=%d\n",dir,level);
	  // setting the first oct
	
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue; // in case the level is empty

      do{ // sweeping through the octs of level
	curoct=nextoct;
	dxcur=1./POW(2,curoct->level);
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell
	    newoct=NULL;
	    newcell=NULL;
	    // A few informations about the current cell
	    getcellnei(icell, vnei, vcell);
		
	    // sweeping the particles
	    if(nexp!=NULL){
	      do{ 
		curp=nexp; 
		nexp=curp->next; 
		out=0;
		switch(dir){

		case 0: // ======================================   x displacement============

		  xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
		  out=((curp->x-xc)>dxcur*0.5)-((curp->x-xc)<-dxcur*0.5);
		  if(out==1){
		    if(vnei[1]==6){ // the particle will remain in the same oct
		      newcell=&(curoct->cell[vcell[1]]);
			    
		      // if refined we assign the particle to a refined cell
		      if(newcell->child!=NULL){
			newoct=newcell->child;
			dxcur2=1./POW(2.,newoct->level);

			// new cell coordinates
			//ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			ip=xp+yp*2+zp*4;
				
			// some particles will experience more than one displacement (along diagonals) we store them in cell 0
			//if((ip>7)||(ip<0)) ip=0; 
				
			newcell=&(newoct->cell[ip]);

		      }

		    }
		    else{ // the particle drifts to a new oct
		      // periodic boundaries
		      curp->x+=(curp->x<0.)-(curp->x>1.);

		      if(curoct->nei[vnei[1]]->child==NULL){ // the particle will go to level-1
			newcell=curoct->nei[vnei[1]];
		      }
		      else{ // the particle will remain at the same level or more
			newcell=&(curoct->nei[vnei[1]]->child->cell[vcell[1]]);
			    
			// if refined we assign the particle to a refined cell
			if(newcell->child!=NULL){
			  newoct=newcell->child;
			  dxcur2=1./POW(2.,newoct->level);

			  // new cell coordinates
			  //ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			  xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			  yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			  zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			  ip=xp+yp*2+zp*4;
				
			  // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			  //if((ip>7)||(ip<0)) ip=0; 

			  newcell=&(newoct->cell[ip]);

			}
		      }
		    }
		  }
		  else if(out==-1){
		    if(vnei[0]==6){ // the particle will remain in the same oct
		      newcell=&(curoct->cell[vcell[0]]);

		      // if refined we assign the particle to a refined cell
		      if(newcell->child!=NULL){
			newoct=newcell->child;
			dxcur2=1./POW(2.,newoct->level);

			xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			ip=xp+yp*2+zp*4;
			      

			//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
			// some particles will experience more than one displacement (along diagonals) we store them in cell 0
			//if((ip>7)||(ip<0)) ip=0; 

			newcell=&(newoct->cell[ip]);
		      }
		    }
		    else{
		      // periodic boundaries
		      curp->x+=(curp->x<0.)-(curp->x>1.);
			    
		      if(curoct->nei[vnei[0]]->child==NULL){ // the particle will go to level-1
			newcell=(curoct->nei[vnei[0]]);
		      }
		      else{ // the particle will remain at the same level
			newcell=&(curoct->nei[vnei[0]]->child->cell[vcell[0]]);
			// if refined we assign the particle to a refined cell
			if(newcell->child!=NULL){
			  newoct=newcell->child;
			  dxcur2=1./POW(2.,newoct->level);
				
			  xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			  yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			  zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			  ip=xp+yp*2+zp*4;
			      

			  //				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
			  // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			  //if((ip>7)||(ip<0)) ip=0; 

			  newcell=&(newoct->cell[ip]);
			}

		      }
		    }
		  }
		  break;

		case 1: // ======================================   y displacement============
		  yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
		
		  out=((curp->y-yc)>dxcur*0.5)-((curp->y-yc)<-dxcur*0.5);
		  // getting the new cell // TO FIX PERIODICITY IS ILL DEFINED HERE
		  if(out==1){
		    if(vnei[3]==6){ // the particle will remain in the same oct
		      newcell=&(curoct->cell[vcell[3]]);
		      // if refined we assign the particle to a refined cell at level+1
		      if(newcell->child!=NULL){
			newoct=newcell->child;
			dxcur2=1./POW(2.,newoct->level);

			xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			ip=xp+yp*2+zp*4;

			//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			//if((ip>7)||(ip<0)) ip=0; 
			newcell=&(newoct->cell[ip]);
		      }

		    }
		    else{ // the particle moves to another oct
		      // periodic boundaries
		      curp->y+=(curp->y<0.)-(curp->y>1.);
			    
		      if(curoct->nei[vnei[3]]->child==NULL){ // the particle will go to level-1
			newcell=(curoct->nei[vnei[3]]);
			//if(curp->idx==222044) printf("coucou hello");
		      }
		      else{ // the particle will remain at the same level or level+1
			newcell=&(curoct->nei[vnei[3]]->child->cell[vcell[3]]);
			      
			     
			// if refined we assign the particle to a refined cell at level+1
			if(newcell->child!=NULL){
			  newoct=newcell->child;
			  dxcur2=1./POW(2.,newoct->level);

			  xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			  yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			  zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			  ip=xp+yp*2+zp*4;

			  //				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			  //if((ip>7)||(ip<0)) ip=0; 
			  newcell=&(newoct->cell[ip]);
			}
		      }
		    }
		  }
		  else if(out==-1){
		    if(vnei[2]==6){ // the particle will remain in the same oct
		      newcell=&(curoct->cell[vcell[2]]);
		      // if refined we assign the particle to a refined cell
		      if(newcell->child!=NULL){
			newoct=newcell->child;
			dxcur2=1./POW(2.,newoct->level);
			xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			ip=xp+yp*2+zp*4;

			//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			// some particles will experience more than one displacement (along diagonals) we store them in cell 0
			//if((ip>7)||(ip<0)) ip=0; 
			newcell=&(newoct->cell[ip]);
		      }

		    }
		    else{
		      // periodic boundaries
		      curp->y+=(curp->y<0.)-(curp->y>1.);

		      if(curoct->nei[vnei[2]]->child==NULL){ // the particle will go to level-1
			newcell=(curoct->nei[vnei[2]]);
		      }
		      else{ // the particle will remain at the same level
			newcell=&(curoct->nei[vnei[2]]->child->cell[vcell[2]]);
			// if refined we assign the particle to a refined cell
			if(newcell->child!=NULL){
			  newoct=newcell->child;
			  dxcur2=1./POW(2.,newoct->level);
			  xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			  yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			  zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			  ip=xp+yp*2+zp*4;

			  //				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			  // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			  //if((ip>7)||(ip<0)) ip=0; 
			  newcell=&(newoct->cell[ip]);
			}

		      }
		    }
		  }

		  break;

		case 2: // ======================================   z displacement============
		  zc=curoct->z+(icell/4)*dxcur+dxcur/2;
			
		  out=((curp->z-zc)>dxcur*0.5)-((curp->z-zc)<-dxcur*0.5);
		  // getting the new cell
		  if(out==1){
		    if(vnei[5]==6){ // the particle will remain in the same oct
		      newcell=&(curoct->cell[vcell[5]]);
		      // if refined we assign the particle to a refined cell
		      if(newcell->child!=NULL){
			//if(curp->idx==82279) printf("to level +1");
			      
			newoct=newcell->child;
			dxcur2=1./POW(2.,newoct->level);
			      
			xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			ip=xp+yp*2+zp*4;
			      
			//ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
			// some particles will experience more than one displacement (along diagonals) we store them in cell 0
			//if((ip>7)||(ip<0)) ip=0; 
			      
			newcell=&(newoct->cell[ip]);
		      }
		    }
		    else{
		      // periodic boundaries
		      curp->z+=(curp->z<0.)-(curp->z>1.);

		      if(curoct->nei[vnei[5]]->child==NULL){ // the particle will go to level-1
			newcell=(curoct->nei[vnei[5]]);
		      }
		      else{ // the particle will remain at the same level
			newcell=&(curoct->nei[vnei[5]]->child->cell[vcell[5]]);
			// if refined we assign the particle to a refined cell
			if(newcell->child!=NULL){
			  //if(curp->idx==82279) printf("to level +1");

			  newoct=newcell->child;
			  dxcur2=1./POW(2.,newoct->level);

			  xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			  yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			  zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			  ip=xp+yp*2+zp*4;

			  //ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;

			  // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			  //if((ip>7)||(ip<0)) ip=0; 

			  newcell=&(newoct->cell[ip]);
			}

		      }
		    }
		  }
		  else if(out==-1){
		    if(vnei[4]==6){ // the particle will remain in the same oct
		      newcell=&(curoct->cell[vcell[4]]);
		      // if refined we assign the particle to a refined cell
		      if(newcell->child!=NULL){
			newoct=newcell->child;
			dxcur2=1./POW(2.,newoct->level);
			      
			xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			ip=xp+yp*2+zp*4;
			      
			//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
			// some particles will experience more than one displacement (along diagonals) we store them in cell 0
			//if((ip>7)||(ip<0)) ip=0; 
			      
			newcell=&(newoct->cell[ip]);
		      }

		    }
		    else{
		      // periodic boundaries
		      curp->z+=(curp->z<0.)-(curp->z>1.);

		      if(curoct->nei[vnei[4]]->child==NULL){ // the particle will go to level-1
			newcell=(curoct->nei[vnei[4]]);
		      }
		      else{ // the particle will remain at the same level
			newcell=&(curoct->nei[vnei[4]]->child->cell[vcell[4]]);
			// if refined we assign the particle to a refined cell
			if(newcell->child!=NULL){
			  newoct=newcell->child;
			  dxcur2=1./POW(2.,newoct->level);
				
			  xp=(int)(DFACT*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			  yp=(int)(DFACT*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			  zp=(int)(DFACT*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			  ip=xp+yp*2+zp*4;

			  //				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;

			  // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			  //				if((ip>7)||(ip<0)) ip=0; 

			  newcell=&(newoct->cell[ip]);
			}
		      }
		    }
		  }
		  break;
		}
		      
		if(out!=0){
		  // ========= assigning the particle to the newcell + pointer management

		  /* if(curp->idx==128352){ */
		  /*   printf("COUCOUC idx cell=%d pos=%e\n",newcell->idx,curp->x); */
		  /* } */
		  

		  struct OCT *ooct;
		  ooct=cell2oct(newcell);
		  
		  // we change the level of the particle
		  curp->level=ooct->level;

		  // removing the cell from its old list
		  if(curp->prev !=NULL){
		    curp->prev->next=curp->next;
		  }
		  else{
		    curoct->cell[icell].phead=curp->next; 
		  }
		      
		  if(curp->next!=NULL){
		    curp->next->prev=curp->prev;
		  }
		      
		  // update the linking of the current particle
		      
		  // adding the particle to the new cell
		  if(newcell->child!=NULL){
		    printf("erro in part displacement\n");
		    abort();
		  }

		  
		  if(newcell->phead!=NULL){
		    part=findlastpart(newcell->phead);
		    part->next=curp;
		    curp->prev=part;
		  }
		  else{
		    newcell->phead=curp;
		    curp->prev=NULL;
		  }
		      
		  curp->next=NULL;
		  
		  
		}

	      }while(nexp!=NULL);
	    }
	  }
      }while(nextoct!=NULL);
    }

}    

void  partcellreorg_GPU(int levelcoarse,int levelmax,struct OCT **firstoct){

  int dir;
  char out;
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *newoct;
  REAL dxcur,dxcur2;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  struct PART *part;
  REAL xc,yc,zc;
  int xp,yp,zp;
  int vnei[6],vcell[6];
  struct CELL *newcell;
  int ip;

  //if(cpu->rank==RANK_DISP) printf("particles exchange\n");
  for(dir=0;dir<3;dir++) 
    { 
      for(level=levelcoarse;level<=levelmax;level++) // looping over levels
	{
	  //printf("dir=%d level=%d\n",dir,level);
	  // setting the first oct
	
	  nextoct=firstoct[level-1];
	  if(nextoct==NULL) continue; // in case the level is empty

	  do{ // sweeping through the octs of level
	    curoct=nextoct;
	    dxcur=1./POW(2,curoct->level);
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
		{
		  nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell
		  newoct=NULL;
		  newcell=NULL;
		  // A few informations about the current cell
		  getcellnei(icell, vnei, vcell);
		
		  // sweeping the particles
		  if(nexp!=NULL){
		    do{ 
		      curp=nexp; 
		      nexp=curp->next; 
		      out=0;
		      switch(dir){
		      case 0: // ======================================   x displacement============
			xc=curoct->x+( icell   %2)*dxcur+0.5*dxcur; // coordinates of the cell 
			out=((curp->x-xc)>dxcur)-((curp->x-xc)<0.);
			if(out==1){
			  if(vnei[1]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[1]]);
			    
			    // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
				
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				newcell=&(newoct->cell[ip]);

			      }

			  }
			  else{ // the particle drifts to a new oct
			    // periodic boundaries
			    curp->x+=(curp->x<0.)-(curp->x>1.);

			    if(curoct->nei[vnei[1]]->child==NULL){ // the particle will go to level-1
			      newcell=curoct->nei[vnei[1]];
			    }
			    else{ // the particle will remain at the same level or more
			      newcell=&(curoct->nei[vnei[1]]->child->cell[vcell[1]]);
			    
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
				
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);

			      }
			    }
			  }
			}
			else if(out==-1){
			  if(vnei[0]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[0]]);

			    // if refined we assign the particle to a refined cell
			    if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
			      

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);
			      }
			  }
			  else{
			    // periodic boundaries
			    curp->x+=(curp->x<0.)-(curp->x>1.);
			    
			    if(curoct->nei[vnei[0]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[0]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[0]]->child->cell[vcell[0]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);
				
				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
			      

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
			      
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);
			      }

			    }
			  }
			}
			break;
		    case 1: // ======================================   y displacement============
			yc=curoct->y+((icell/2)%2)*dxcur+0.5*dxcur;
			out=((curp->y-yc)>dxcur)-(curp->y<yc);
			// getting the new cell 
			if(out==1){
			  if(vnei[3]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[3]]);
			    // if refined we assign the particle to a refined cell at level+1
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				//if((ip>7)||(ip<0)) ip=0; 
				newcell=&(newoct->cell[ip]);
			      }

			  }
			  else{ // the particle moves to another oct
			    // periodic boundaries
			    curp->y+=(curp->y<0.)-(curp->y>1.);
			    
			    if(curoct->nei[vnei[3]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[3]]);
			      //if(curp->idx==222044) printf("coucou hello");
			    }
			    else{ // the particle will remain at the same level or level+1
			      newcell=&(curoct->nei[vnei[3]]->child->cell[vcell[3]]);
			      
			     
			      // if refined we assign the particle to a refined cell at level+1
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(DFACT*(curp->x-newoct->x)/dxcur)+(int)(DFACT*(curp->y-newoct->y)/dxcur)*2+(int)(DFACT*(curp->z-newoct->z)/dxcur)*4;
				//if((ip>7)||(ip<0)) ip=0; 
				newcell=&(newoct->cell[ip]);
			      }
			    }
			  }
			}
			else if(out==-1){
			  if(vnei[2]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[2]]);
			    // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);
				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 
				newcell=&(newoct->cell[ip]);
			      }

			  }
			  else{
			    // periodic boundaries
			    curp->y+=(curp->y<0.)-(curp->y>1.);

			    if(curoct->nei[vnei[2]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[2]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[2]]->child->cell[vcell[2]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);
				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 
				newcell=&(newoct->cell[ip]);
			      }

			    }
			  }
			}

			break;
		      case 2: // ======================================   z displacement============
			zc=curoct->z+(icell/4)*dxcur+0.5*dxcur;
			out=((curp->z-zc)>dxcur)-((curp->z-zc)<0.);
			// getting the new cell
			if(out==1){
			  if(vnei[5]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[5]]);
			      // if refined we assign the particle to a refined cell
			    if(newcell->child!=NULL){
			      //if(curp->idx==82279) printf("to level +1");
			      
			      newoct=newcell->child;
			      dxcur2=1./POW(2.,newoct->level);
			      
			      xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			      yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			      zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			      ip=xp+yp*2+zp*4;
			      
			      //ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
			      
			      // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			      //if((ip>7)||(ip<0)) ip=0; 
			      
			      newcell=&(newoct->cell[ip]);
			    }
			  }
			  else{
			    // periodic boundaries
			    curp->z+=(curp->z<0.)-(curp->z>1.);

			    if(curoct->nei[vnei[5]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[5]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[5]]->child->cell[vcell[5]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				//if(curp->idx==82279) printf("to level +1");

				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);

				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;

				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);
			      }

			    }
			  }
			}
			else if(out==-1){
			  if(vnei[4]==6){ // the particle will remain in the same oct
			    newcell=&(curoct->cell[vcell[4]]);
			    // if refined we assign the particle to a refined cell
			    if(newcell->child!=NULL){
			      newoct=newcell->child;
			      dxcur2=1./POW(2.,newoct->level);
			      
			      xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			      yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			      zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
			      ip=xp+yp*2+zp*4;
			      
			      //				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
			      
			      // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			      //if((ip>7)||(ip<0)) ip=0; 
			      
			      newcell=&(newoct->cell[ip]);
			    }

			  }
			  else{
			    // periodic boundaries
			    curp->z+=(curp->z<0.)-(curp->z>1.);

			    if(curoct->nei[vnei[4]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[4]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[4]]->child->cell[vcell[4]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./POW(2.,newoct->level);
				
				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;

				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				//				if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);
			      }
			    }
			  }
			}
			break;
		      }
		      
		      if(out!=0){
			// ========= assigning the particle to the newcell + pointer management

			// removing the cell from its old list
			if(curp->prev !=NULL){
			  curp->prev->next=curp->next;
			}
			else{
			  curoct->cell[icell].phead=curp->next; 
			}
		      
			if(curp->next!=NULL){
			  curp->next->prev=curp->prev;
			}
		      
			// update the linking of the current particle
		      
			// adding the particle to the new cell
			if(newcell->child!=NULL){
			  printf("erro in part displacement\n");
			  abort();
			}


			if(newcell->phead!=NULL){
			  part=findlastpart(newcell->phead);
			  part->next=curp;
			  curp->prev=part;
			}
			else{
			  newcell->phead=curp;
			  curp->prev=NULL;
			}
		      
			curp->next=NULL;
		      }
		      
		      
		    }while(nexp!=NULL);
		  }
		}
	  }while(nextoct!=NULL);
	}	
      
    }

}    


//------------------------------------------------------------------------
#ifdef PIC

void L_accelpart(int level,struct OCT **firstoct, REAL *adt, int is, struct CPUINFO *cpu){

  struct OCT *nextoct;
  struct OCT *curoct;
  int icomp,icell;
  struct PART *curp;
  struct PART *nexp;
  
  // ==================================== Computing the Velocities
  // ==================================== performing the INVERSE CIC assignement
  
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){ 
    do // sweeping level
      {
	curoct=nextoct;
	nextoct=curoct->next;
      
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
	    if(nexp!=NULL){ 
	      do{  
		curp=nexp; 
		nexp=curp->next;
		if(((curp->level==level)||(curp->is==0))||(is==-1)) cell2part_cic(curp, curoct, icell,adt[curp->level-1]*0.5); // here 0.5 because of half timestep for midpoint rule
	      }while(nexp!=NULL); 
	    }
	  }
      }while(nextoct!=NULL);
  }
  
}
#endif



//------------------------------------------------------------------------

void egypart(struct CPUINFO *cpu, REAL *ekintot, REAL *epottot,struct RUNPARAMS *param, REAL tsim)
{
  struct PART *curp;
  struct PART *lastp;
  
  REAL ekinloc=0.;
  REAL epotloc=0.;
  lastp=cpu->firstpart+param->npartmax;
  int ipart=0;
  for(curp=cpu->firstpart;curp<lastp;curp++){
    if(curp->mass>0){
      ekinloc+=curp->ekin*curp->mass;
      epotloc+=(0.5*curp->mass*curp->epot*tsim); // here we adjust for the lack of aexp in supercomoving expression we used for Poisson
      ipart++;
    }
  }

  //printf("ipart =%d\n",ipart);
  *ekintot=ekinloc;
  *epottot=epotloc;
}
  

#endif
