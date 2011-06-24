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

#define FRACDX 0.4


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
      npart++;
      curp=nexp; 
      nexp=curp->next; 
    }while(nexp!=NULL); 
  }

  return npart;
}

//------------------------------------------------------------------------

struct PART* modifpospart(struct PART* phead, float len, int dir)
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

float comptstep(int levelcoarse,int levelmax,struct OCT** firstoct, float fa, float fa2, struct CPUINFO* cpu, float tmax){
  
  int level;
  float vmax,lmdisp;
  struct OCT *nextoct;
  struct OCT oct;
  float dxcur;
  int icell;
  struct PART *nexp;
  struct PART *curp;
  float va;
  float aa;
  float dtlev,dtnew=1e9;
  float dt;

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
	  dxcur=1./pow(2,oct.level+1); // +1 to protect level change
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {

#ifdef PIC
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){

#ifdef AXLFORCE
		//cell acceleration
		aa=sqrt(oct.cell[icell].fx*oct.cell[icell].fx+oct.cell[icell].fy*oct.cell[icell].fy+oct.cell[icell].fz*oct.cell[icell].fz)*fa2;
#else
		aa=0.;
#endif
		do{ 
		  curp=nexp; 
		  nexp=curp->next; 

#ifdef PART2
		  if(curp->idx==0) continue;
#endif
		  // particle velocit
		  va=sqrt(curp->vx*curp->vx+curp->vy*curp->vy+curp->vz*curp->vz)*fa;

		  // loc dt
		  if((va>0.)&&(aa>0.)){
		    float EPS=2.0*(FRACDX*dxcur)*aa/(va*va);
		    if(EPS<1e-1){
		      dtlev=(FRACDX*dxcur)/va;
		    }
		    else if(EPS>10.){
		      dtlev=sqrt(2.0*(FRACDX*dxcur)/aa);
		    }
		    else{
		      dtlev=va/aa*(sqrt(1.+2.*aa*(FRACDX*dxcur)/(va*va))-1.);
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
  MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_FLOAT,MPI_MIN,cpu->comm);
#endif  

  dt=(dt>tmax?tmax:dt);
  return dt;
}

//=====================================================================
//=====================================================================

float movepart(int levelcoarse,int levelmax,struct OCT** firstoct, float dt, struct CPUINFO* cpu){
  
  int level;
  float mdisp,lmdisp;
  struct OCT *nextoct;
  struct OCT oct;
  float dxcur;
  int icell;
  struct PART *nexp;
  struct PART *curp;
  float disp;
  float dtlev,dtnew;

  // === Moving particles

  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      // setting the first oct
      mdisp=0.;
      nextoct=firstoct[level-1];
      
      if(nextoct==NULL) continue; // in case the level is empty
      do // sweeping through the octs of level
	{
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
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
		  disp=sqrt(curp->vx*dt*curp->vx*dt+curp->vy*dt*curp->vy*dt+curp->vz*dt*curp->vz*dt);
		  /* if(curp->idx==96571){ */
		  /*   printf("bug disp=%f or %f dx// xinit=%f xinfal=%f\n",disp,disp/dxcur,curp->x,curp->x+curp->vx*dt); */
		  /* } */
		  curp->x+=curp->vx*dt;
		  curp->y+=curp->vy*dt;
		  curp->z+=curp->vz*dt;
		  //if(curp->idx==0) printf("%e %e %e\n",curp->vx,curp->vy,curp->vz);
		  if(disp>mdisp) mdisp=disp;

		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
      if(cpu->rank==0) printf("level=%d maxdisp=%f or %f dx\n",level,mdisp,mdisp/dxcur);
    }

  return dt;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void  partcellreorg(int levelcoarse,int levelmax,struct OCT **firstoct){

  int dir;
  char out;
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *newoct;
  float dxcur,dxcur2;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  struct PART *part;
  float xc,yc,zc;
  int xp,yp,zp;
  int vnei[6],vcell[6];
  struct CELL *newcell;
  int ip;

  //printf("particles exchange\n");
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
	    dxcur=1./pow(2,curoct->level);
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
				dxcur2=1./pow(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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
				dxcur2=1./pow(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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
				dxcur2=1./pow(2.,newoct->level);

				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
			      

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
			      
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
				dxcur2=1./pow(2.,newoct->level);
				
				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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
				dxcur2=1./pow(2.,newoct->level);

				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
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
				dxcur2=1./pow(2.,newoct->level);

				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
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
				dxcur2=1./pow(2.,newoct->level);
				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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
				dxcur2=1./pow(2.,newoct->level);
				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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
			      dxcur2=1./pow(2.,newoct->level);
			      
			      xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			      yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			      zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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
				dxcur2=1./pow(2.,newoct->level);

				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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
			      dxcur2=1./pow(2.,newoct->level);
			      
			      xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			      yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			      zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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
				dxcur2=1./pow(2.,newoct->level);
				
				xp=(int)(2*(curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
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

			/* if(curp->idx==96571){ */
			/*   struct OCT* ooct; */
			/*   if(newcell!=NULL){ */
			/*     ooct=cell2oct(newcell); */
			/*     printf("PPZ %f %f %f %d %d  curoct=%f %f %f newoct=%f %f %f icell=%d dir=%d curoct=%p newoct=%p dxcur=%f final level=%d init level=%d\n",curp->x,curp->y,curp->z,ip,out,curoct->x,curoct->y,curoct->z,ooct->x,ooct->y,ooct->z,icell,dir,curoct,ooct,dxcur,ooct->level,curoct->level); */
			/*   } */
			/* } */


		    }while(nexp!=NULL);
		  }
		}
	  }while(nextoct!=NULL);
	}	
      
    }

}    

void  partcellreorg_GPU(int levelcoarse,int levelmax,struct OCT **firstoct){

  int dir;
  char out;
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *newoct;
  float dxcur,dxcur2;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  struct PART *part;
  float xc,yc,zc;
  int xp,yp,zp;
  int vnei[6],vcell[6];
  struct CELL *newcell;
  int ip;

  printf("particles exchange\n");
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
	    dxcur=1./pow(2,curoct->level);
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
				dxcur2=1./pow(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
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
				dxcur2=1./pow(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
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
				dxcur2=1./pow(2.,newoct->level);

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
			    curp->x+=(curp->x<0.)-(curp->x>1.);
			    
			    if(curoct->nei[vnei[0]]->child==NULL){ // the particle will go to level-1
			      newcell=(curoct->nei[vnei[0]]);
			    }
			    else{ // the particle will remain at the same level
			      newcell=&(curoct->nei[vnei[0]]->child->cell[vcell[0]]);
			      // if refined we assign the particle to a refined cell
			      if(newcell->child!=NULL){
				newoct=newcell->child;
				dxcur2=1./pow(2.,newoct->level);
				
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
				dxcur2=1./pow(2.,newoct->level);

				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
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
				dxcur2=1./pow(2.,newoct->level);

				xp=(int)((curp->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)((curp->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)((curp->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
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
				dxcur2=1./pow(2.,newoct->level);
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
				dxcur2=1./pow(2.,newoct->level);
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
			      dxcur2=1./pow(2.,newoct->level);
			      
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
				dxcur2=1./pow(2.,newoct->level);

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
			      dxcur2=1./pow(2.,newoct->level);
			      
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
				dxcur2=1./pow(2.,newoct->level);
				
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

			/* if(curp->mass>0.5){ */
			/*   struct OCT *ooct; */
			/*   ooct=cell2oct(newcell); */
			/*   switch(dir){ */
			/*   case 0: */
			/*     printf("xp=%e xoorg=%e xo=%e xc=%e icell=%d dir=%d out=%d dif=%e\n",curp->x,curoct->x,ooct->x,xc,icell,dir,out,curp->x-xc); */
			/*     break; */
			/*   case 1: */
			/*     printf("xp=%e xoorg=%e xo=%e xc=%e icell=%d dir=%d out=%d dif=%e\n",curp->y,curoct->y,ooct->y,yc,icell,dir,out,curp->y-yc); */
			/*     break; */
			/*   case 2: */
			/*     printf("xp=%e xoorg=%e xo=%e xc=%e icell=%d dir=%d out=%d dif=%e\n",curp->z,curoct->z,ooct->z,zc,icell,dir,out,curp->z-zc); */
			/*     break; */
			/*   } */
			/* } */

		      }
		      
		      
		    }while(nexp!=NULL);
		  }
		}
	  }while(nextoct!=NULL);
	}	
      
    }

}    

//------------------------------------------------------------------------
//------------------------------------------------------------------------
void forcevel(int levelcoarse,int levelmax,struct OCT **firstoct, float **vcomp,int stride,float dt, struct CPUINFO *cpu, struct PACKET **sendbuffer,struct PACKET **recvbuffer){

  int dir;
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  float dx;
  int icomp,icell;
  struct PART *curp;
  struct PART *nexp;
  int nread;

#ifndef AXLFORCE
  for(dir=0;dir<3;dir++){
#endif
    //printf("Force Start dir=%d\n",dir);
    for(level=levelcoarse;level<=levelmax;level++)
      {
	// COARSE LEVEL TREATMENT ============================================
	//if(level==levelcoarse) fixbound_reg(grid,levelcoarse,NBND);
	nextoct=firstoct[level-1];
	
	dx=pow(0.5,level);
	if(nextoct==NULL){
	  continue;
	}
	else{
	  do{ 
	    curoct=nextoct;
#ifndef AXLFORCE	    
	    // First we gather the potential in all neighbors
	    for(icomp=2*dir;icomp<=2*dir+1;icomp++){
	      memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
	      nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,cpu,&nread);
	    }
	    // Next we perform the finite difference along x
	    grad(vcomp,stride,dx,dir);
	  
	    // Then we scatter back the results in the temp variable
	    nextoct=scattercomp(curoct, vcomp[6], 6, 2, stride,cpu);
#else
	    // First we gather the potential in all neighbors
	    for(icomp=0;icomp<6;icomp++){
	      memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
	      nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,cpu,&nread);
	    }

	    // Next we perform the finite difference along x
	    grad(vcomp,stride,dx,0);
	    // Next we perform the finite difference along y
	    grad(vcomp,stride,dx,1);
	    // Next we perform the finite difference along z
	    grad(vcomp,stride,dx,2);

	    // Then we scatter back the results in the force variable
	    nextoct=scattercomp(curoct, vcomp[6], 6, 3, stride,cpu);
	    nextoct=scattercomp(curoct, vcomp[7], 6, 4, stride,cpu);
	    nextoct=scattercomp(curoct, vcomp[8], 6, 5, stride,cpu);
#endif


	  
	  }while(nextoct!=NULL);
	
	}
      }

#ifdef WMPI
#ifndef AXLFORCE
    mpi_exchange(cpu,sendbuffer,recvbuffer,4,1); // temp field
#else
    mpi_exchange(cpu,sendbuffer,recvbuffer,5,1); // fx field
    mpi_exchange(cpu,sendbuffer,recvbuffer,6,1); // fy field
    mpi_exchange(cpu,sendbuffer,recvbuffer,7,1); // fz field
#endif
#endif
  
#ifdef PIC
    // ==================================== Computing the Velocities
    // ==================================== performing the INVERSE CIC assignement
  
    //printf("start INVERSE CIC\n");
    //start INVERSE CIC
    for(level=levelmax;level>=levelcoarse;level--)
      {
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
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
		    cell2part_cic(curp, curoct, icell,dir,dt); 
		  }while(nexp!=NULL); 
		}
	      }
	  }while(nextoct!=NULL);
      }
#endif
#ifndef AXLFORCE
  }
#endif
}


//------------------------------------------------------------------------

float egypart(int levelcoarse,int levelmax,struct OCT **firstoct, struct CPUINFO *cpu){

  int dir;
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  float dx;
  int icomp,icell;
  struct PART *curp;
  struct PART *nexp;
  int nread;

  float upart, utot=0.; // potential energy
  float tpart, ttot=0.; // kinetic energy
  float etot; // total energy

  float vx,vy,vz;

  // ==================================== performing the INVERSE CIC assignement
  
    //printf("start INVERSE CIC\n");
    //start INVERSE CIC
    for(level=levelmax;level>=levelcoarse;level--)
      {
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	  
	    //==== FIRST WE CONSIDER THE PARTICLES INSIDE THE BOUNDARIES AT LEVEL L
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
		if(nexp!=NULL){ 
		  do{  
		    curp=nexp; 
		    nexp=curp->next; 
		    upart=cell2part_cic_egy(curp, curoct, icell); 
		    utot=utot-upart;
		    
		    vx=curp->vx;
		    vy=curp->vy;
		    vz=curp->vz;

		    tpart=0.5*curp->mass*(vx*vx+vy*vy+vz*vz);
		    //printf("tpart=%e upart=%e\n",tpart,upart);
		    ttot+=tpart;
		  }while(nexp!=NULL); 
		}
	      }
	  }while(nextoct!=NULL);
      }

    //printf("%e %e \n",utot,ttot);
    etot=utot*0.5+ttot;
    return etot;
}
