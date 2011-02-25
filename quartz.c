#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef WMPI
#include <mpi.h>
#endif

#include "hilbert.h"


// TO COMPILE

//mpicc -lm -DWMPI -DNBUFF=4096 -DNEWASSIGN -DTESTPLUM -DNDUMP=1 -DNSTEP=10 -DLCOARSE=6 -DLMAX=6  -DCIC2 -DNITER=4096 -DALLCELL -DSTRIDE=128 -DDT=1e-4 hilbert.c gilgamesh.c

//gcc -g -lm -DNEWASSIGN -DPART2 -DNDUMP=1 -DNSTEP=10 -DLCOARSE=6 -DLMAX=6  -DCIC2 -DNITER=4096 -DALLCELL -DSTRIDE=128 -DDT=1e-4 hilbert.o gilgamesh.c

//mpicc -lm -DWMPI -DNEWASSIGN -DPART2 -DNDUMP=1 -DNSTEP=10 -DLCOARSE=6 -DLMAX=6  -DCIC2 -DNITER=4096 -DALLCELL -DSTRIDE=128 -DDT=1e-4 hilbert.c gilgamesh.c


//=======================================

// this structure exists for MPI communication protocol

struct PACKET{
  float data[8]; // the data to be transfered (8 since we transmit data per octs)
  int key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
};


//=======================================


struct CPUINFO{
  int rank;
  int nproc;

  int kmin;
  int kmax;
  int nkeys;

  float load;

  struct OCT **bndoct; // the list of external boundary octs

  int nebnd; // the number of external boundary octs
  int nnei; // the number of neighbors procs
  
  int *mpinei; // the ranks of the neighbors procs

  int *dict; // a hash table to converts ranks in local neighbour index

  struct OCT **htable; // the hashing table to recover the octs from hilbert keys

  int *allkmin;
  int *allkmax;

  int nbuff; // the max number of buffer cells

#ifdef WMPI
  MPI_Datatype *MPI_PACKET; // the structured type for MPI messages (fields)
  MPI_Datatype *MPI_PART; // the structured type for MPI messages (particles)
  MPI_Comm comm;
#endif
  
};

//=======================================

struct PART
{
  float x;
  float y;
  float z;

  float vx;
  float vy;
  float vz;

  struct PART *next;
  struct PART *prev;

  float mass;

  int idx;

};

struct PART_MPI // For mpi communications
{
  float x;
  float y;
  float z;

  float vx;
  float vy;
  float vz;

  float mass;
  int idx;

  int key; // the destination hilbert key
  int level; // the level of the destination (to remove the key degeneracy)
  int icell; // the cell of destination
};

//=======================================



//-------------------------------------
struct CELL
{
  struct OCT *child;
  float marked; // float for consistency with physical quantities during communications
  char idx; //index of the cell within the oct

  // the head particle
  struct PART * phead;

  // the physical quantities
  float density;
  float pot;
  float temp;

};



//-------------------------------------
//-------------------------------------------------------------------------
struct OCT
{
  // the cell properties
  struct CELL cell[8];

  int level;// level of the cells in the oct
  struct CELL *parent; // parent cell 
  struct CELL *nei[6];// neighbor cells at level - 1
  
  // the next two pointers are required for sweeps through a single level
  struct OCT *next; // next oct on the same level
  struct OCT *prev; // previous oct on the same level

  // the oct position (lowest left corner)
  float x;
  float y;
  float z;

  // parallel data
  int cpu; 
  struct OCT *nexthash;


  
};


void breakmpi()
{
  {
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
      sleep(5);
  }
}

//------------------------------------------------------------------------

// This function returns the parent oct of a given cell

struct OCT* cell2oct(struct CELL* cell)
{
  long int adress;
  struct OCT *poct;
  adress=(long int) cell;
  adress=adress-cell->idx*sizeof(struct CELL);
  poct=(struct OCT*) adress;

  return poct;
}

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


//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

void getcellnei(int cindex, int *neip, int *cell)
{
  switch(cindex){
  case 0:
    neip[0]=0;cell[0]=1;
    neip[1]=6;cell[1]=1;
    neip[2]=2;cell[2]=2;
    neip[3]=6;cell[3]=2;
    neip[4]=4;cell[4]=4;
    neip[5]=6;cell[5]=4;
    break;
  case 1:
    neip[0]=6;cell[0]=0;
    neip[1]=1;cell[1]=0;
    neip[2]=2;cell[2]=3;
    neip[3]=6;cell[3]=3;
    neip[4]=4;cell[4]=5;
    neip[5]=6;cell[5]=5;
    break;
  case 2:
    neip[0]=0;cell[0]=3;
    neip[1]=6;cell[1]=3;
    neip[2]=6;cell[2]=0;
    neip[3]=3;cell[3]=0;
    neip[4]=4;cell[4]=6;
    neip[5]=6;cell[5]=6;
    break;
  case 3:
    neip[0]=6;cell[0]=2;
    neip[1]=1;cell[1]=2;
    neip[2]=6;cell[2]=1;
    neip[3]=3;cell[3]=1;
    neip[4]=4;cell[4]=7;
    neip[5]=6;cell[5]=7;
    break;
  case 4:
    neip[0]=0;cell[0]=5;
    neip[1]=6;cell[1]=5;
    neip[2]=2;cell[2]=6;
    neip[3]=6;cell[3]=6;
    neip[4]=6;cell[4]=0;
    neip[5]=5;cell[5]=0;
    break;
  case 5:
    neip[0]=6;cell[0]=4;
    neip[1]=1;cell[1]=4;
    neip[2]=2;cell[2]=7;
    neip[3]=6;cell[3]=7;
    neip[4]=6;cell[4]=1;
    neip[5]=5;cell[5]=1;
    break;
  case 6:
    neip[0]=0;cell[0]=7;
    neip[1]=6;cell[1]=7;
    neip[2]=6;cell[2]=4;
    neip[3]=3;cell[3]=4;
    neip[4]=6;cell[4]=2;
    neip[5]=5;cell[5]=2;
    break;
  case 7:
    neip[0]=6;cell[0]=6;
    neip[1]=1;cell[1]=6;
    neip[2]=6;cell[2]=5;
    neip[3]=3;cell[3]=5;
    neip[4]=6;cell[4]=3;
    neip[5]=5;cell[5]=3;
    break;
  }

}
//==================================================================
//==================================================================
//==================================================================

void part2cell_cic(struct PART *curp, struct OCT *curoct, int icell, char full)
{
  float xc,yc,zc;
  float dxcur=1./pow(2,curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int neip[3];
  float tx,ty,tz;
  float dx,dy,dz;
  float contrib;
  struct OCT **cicoct;
  char ciccell[8];
  char fullok=0;

  cicoct=(struct OCT **)calloc(8,sizeof(struct OCT *));

  xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
  yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
  zc=curoct->z+( icell/4   )*dxcur+dxcur/2; 

  // getting the neighbors
  getcellnei(icell, vnei, vcell);
		  
  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);
		  
  // here we denote the offset in ZYX
  ciccell[0]=icell;          //cell 000
  ciccell[1]=vcell[neip[0]]; //cell 001
  ciccell[2]=vcell[neip[1]]; //cell 010
  ciccell[4]=vcell[neip[2]]; //cell 100
  cicoct[0]=curoct;

  // the CIC weights
  tx=fabs((curp->x-xc)/dxcur);
  ty=fabs((curp->y-yc)/dxcur);
  tz=fabs((curp->z-zc)/dxcur);
  if(!full){fullok=(tx<1.)*(ty<1.)*(tz<1.);} // necessary for boundary cells
	  
  dx=1.-tx;
  dy=1.-ty;
  dz=1.-tz;

  // contrib to current cell 000 =====
  contrib=dx*dy*dz;
  if((contrib<=1.)&&(contrib>=0.)){
    if(full||fullok){
      curoct->cell[icell].density+=contrib/pow(dxcur,3)*curp->mass;
      //if(!full) printf("%f %f %f contrib=%f\n",tx,ty,tz,contrib);
    }
  }
  

  if(full) // full=0 for level-1 corrections
    {
      // contrib to 100 cell ===========================================================
  
      contrib=tx*dy*dz;
      if((contrib<=1.)&&(contrib>=0)){
	if(vnei[neip[0]]==6){
	  curoct->cell[vcell[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	  cicoct[1]=curoct;
	}
	else if(curoct->nei[vnei[neip[0]]]->child!=NULL){
	  curoct->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	  cicoct[1]=curoct->nei[vnei[neip[0]]]->child;
	}
      }
      // contrib to 010 cell ===========================================================
      contrib=dx*ty*dz;
      if((contrib<=1.)&&(contrib>=0)){
	if(vnei[neip[1]]==6){
	  curoct->cell[vcell[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	  cicoct[2]=curoct;
	}
	else{
	  if(curoct->nei[vnei[neip[1]]]->child!=NULL){
	    curoct->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[2]=curoct->nei[vnei[neip[1]]]->child;
	  }
	}
      }
  
      // contrib to 001 cell ===========================================================
      contrib=dx*dy*tz;
      if((contrib<=1.)&&(contrib>=0)){
	if(vnei[neip[2]]==6){
	  curoct->cell[vcell[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	  cicoct[4]=curoct;
	}
	else{
	  if(curoct->nei[vnei[neip[2]]]->child!=NULL){
	    curoct->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[4]=curoct->nei[vnei[neip[2]]]->child;
	  }
	}
      }
      //contrib to 110 cell 2 paths ======================================================
  
      contrib=tx*ty*dz;
      if((contrib<=1.)&&(contrib>=0)){
	if(cicoct[1]!=NULL){
	  getcellnei(ciccell[1], vnei2, vcell2);
	  ciccell[3]=vcell2[neip[1]];
	  if(vnei2[neip[1]]==6){
	    cicoct[1]->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[3]=cicoct[1];
	  }
	  else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
	    cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
	  }
	}
  
	if(cicoct[3]==NULL){
	  if(cicoct[2]!=NULL){
	    getcellnei(ciccell[2], vnei2, vcell2);
	    ciccell[3]=vcell2[neip[0]];
	    if(vnei2[neip[0]]==6){
	      cicoct[2]->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[3]=cicoct[2];
	    }
	    else if(cicoct[2]->nei[vnei2[neip[0]]]->child!=NULL){
	      cicoct[2]->nei[vnei2[neip[0]]]->child->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[3]=cicoct[2]->nei[vnei2[neip[0]]]->child;
	    }
	  }
	}
      }
  
      //contrib to 101 cell 2 paths ======================================================
      contrib=tx*dy*tz;
      if((contrib<=1.)&&(contrib>=0)){

	if(cicoct[1]!=NULL){
	  getcellnei(ciccell[1], vnei2, vcell2);
	  ciccell[5]=vcell2[neip[2]];
	  if(vnei2[neip[2]]==6){
	    cicoct[1]->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[5]=cicoct[1];
	  }
	  else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
	    cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
	  }
	}
  
	if(cicoct[5]==NULL){
	  if(cicoct[4]!=NULL){
	    getcellnei(ciccell[4], vnei2, vcell2);
	    ciccell[5]=vcell2[neip[0]];
	    if(vnei2[neip[0]]==6){
	      cicoct[4]->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[5]=cicoct[4];
	    }
	    else if(cicoct[4]->nei[vnei2[neip[0]]]->child!=NULL){
	      cicoct[4]->nei[vnei2[neip[0]]]->child->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[5]=cicoct[4]->nei[vnei2[neip[0]]]->child;
	    }
	  }
	}
      }
		  
      //contrib to 011 cell 2 paths ======================================================
      contrib=dx*ty*tz;
      if((contrib<=1.)&&(contrib>=0)){

	if(cicoct[2]!=NULL){
	  getcellnei(ciccell[2], vnei2, vcell2);
	  ciccell[6]=vcell2[neip[2]];
	  if(vnei2[neip[2]]==6){
	    cicoct[2]->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[6]=cicoct[2];
	  }
	  else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
	    cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
	  }
	}
		  
	if(cicoct[6]==NULL){
	  if(cicoct[4]!=NULL){
	    getcellnei(ciccell[4], vnei2, vcell2);
	    ciccell[6]=vcell2[neip[1]];
	    if(vnei2[neip[1]]==6){
	      cicoct[4]->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[6]=cicoct[4];
	    }
	    else if(cicoct[4]->nei[vnei2[neip[1]]]->child!=NULL){
	      cicoct[4]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[6]=cicoct[4]->nei[vnei2[neip[1]]]->child;
	    }
	  }
	}
      }		  
      // contrib to 111 ======================================================
      contrib=tx*ty*tz;
      if((contrib<=1.)&&(contrib>=0)){
		  
	if(cicoct[3]!=NULL){
	  getcellnei(ciccell[3], vnei2, vcell2);
	  ciccell[7]=vcell2[neip[2]];
	  if(vnei2[neip[2]]==6){
	    cicoct[3]->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[7]=cicoct[3];
	  }
	  else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
	    cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
	  }
	}
		  
	if(cicoct[7]==NULL){
	  if(cicoct[5]!=NULL){
	    getcellnei(ciccell[5], vnei2, vcell2);
	    ciccell[7]=vcell2[neip[1]];
	    if(vnei2[neip[1]]==6){
	      cicoct[5]->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[7]=cicoct[5];
	    }
	    else if(cicoct[5]->nei[vnei2[neip[1]]]->child!=NULL){
	      cicoct[5]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[7]=cicoct[5]->nei[vnei2[neip[1]]]->child;
	    }
	  }
	}
		  
	if(cicoct[7]==NULL){
	  if(cicoct[6]!=NULL){
	    getcellnei(ciccell[6], vnei2, vcell2);
	    ciccell[7]=vcell2[neip[0]];
	    if(vnei2[neip[0]]==6){
	      cicoct[6]->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[7]=cicoct[6];
	    }
	    else if(cicoct[6]->nei[vnei2[neip[0]]]->child!=NULL){
	      cicoct[6]->nei[vnei2[neip[0]]]->child->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[7]=cicoct[6]->nei[vnei2[neip[0]]]->child;
	    }
	  }
	}
      }
    }

#if 0
  //check total
  float tot=0.;
  for(ii=0;ii<8;ii++)
    {
      float xcc,ycc,zcc;
      tot+=cicoct[ii]->cell[ciccell[ii]].density;
      xcc=cicoct[ii]->x+( ciccell[ii]   %2)*dxcur+dxcur/2; // coordinates of the cell center 
      ycc=cicoct[ii]->y+((ciccell[ii]/2)%2)*dxcur+dxcur/2;
      zcc=cicoct[ii]->z+( ciccell[ii]/4   )*dxcur+dxcur/2; 
      
      printf("%f %f %f %d %f\n",xcc,ycc,zcc,ciccell[ii],cicoct[ii]->cell[ciccell[ii]].density*pow(dxcur,3));
    }
  printf("tot=%f\n",tot*pow(dxcur,3));
#endif
}



void cell2part_cic(struct PART *curp, struct OCT *curoct, int icell, char dir, float dt)
{
  float xc,yc,zc;
  float dxcur=1./pow(2,curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int neip[3];
  float tx,ty,tz;
  float dx,dy,dz;
  float contrib;
  struct OCT **cicoct;
  struct OCT *curoctlr;
  char ciccell[8];
  char fullok=0;
  float accel=0.;
  char hres=1; // could be switched to hres=0 if particle is not deep enough

  cicoct=(struct OCT **)calloc(8,sizeof(struct OCT *));

  xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
  yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
  zc=curoct->z+( icell/4   )*dxcur+dxcur/2; 

  // getting the neighbors
  getcellnei(icell, vnei, vcell);
		  
  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);
		  
  // here we denote the offset in ZYX
  ciccell[0]=icell;          //cell 000
  ciccell[1]=vcell[neip[0]]; //cell 001
  ciccell[2]=vcell[neip[1]]; //cell 010
  ciccell[4]=vcell[neip[2]]; //cell 100
  cicoct[0]=curoct;

  // the CIC weights
  tx=fabs((curp->x-xc)/dxcur);
  ty=fabs((curp->y-yc)/dxcur);
  tz=fabs((curp->z-zc)/dxcur);
	  
  dx=1.-tx;
  dy=1.-ty;
  dz=1.-tz;

  // contrib from current cell 000 =====
  contrib=dx*dy*dz;
  if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoct->cell[icell].temp*contrib;
      //curoct->cell[icell].density+=contrib/pow(dxcur,3);
  }
  

  // contrib from 100 cell ===========================================================
      
  if(hres==1){
    contrib=tx*dy*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[0]]==6){
	accel+=curoct->cell[vcell[neip[0]]].temp*contrib;
	cicoct[1]=curoct;
      }
      else {
	if(curoct->nei[vnei[neip[0]]]->child!=NULL){
	  accel+=curoct->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].temp*contrib;
	  cicoct[1]=curoct->nei[vnei[neip[0]]]->child;
	}
	else{ // the particle is not deep enough
	  accel=0.;
	  hres=0;
	}
      }
    }
  }
  // contrib from 010 cell ===========================================================
  if(hres==1){
    contrib=dx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[1]]==6){
	accel+=curoct->cell[vcell[neip[1]]].temp*contrib;
	cicoct[2]=curoct;
      }
      else{
	if(curoct->nei[vnei[neip[1]]]->child!=NULL){
	  accel+=curoct->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].temp*contrib;
	  cicoct[2]=curoct->nei[vnei[neip[1]]]->child;
	}
	else{ // the particle is not deep enough
	  accel=0.;
	  hres=0;
	}
      }
    }
  }
  // contrib from 001 cell ===========================================================

  if(hres==1){
    contrib=dx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[2]]==6){
	accel+=curoct->cell[vcell[neip[2]]].temp*contrib;
	cicoct[4]=curoct;
      }
      else{
	if(curoct->nei[vnei[neip[2]]]->child!=NULL){
	  accel+=curoct->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].temp*contrib;
	  cicoct[4]=curoct->nei[vnei[neip[2]]]->child;
	}
	else{ // the particle is not deep enough
	  accel=0.;
	  hres=0;
	}
      }
    }
  }
  //contrib from 110 cell 2 paths ======================================================
  
  contrib=tx*ty*dz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){
      if(cicoct[1]!=NULL){
	getcellnei(ciccell[1], vnei2, vcell2);
	ciccell[3]=vcell2[neip[1]];
	if(vnei2[neip[1]]==6){
	  accel+=cicoct[1]->cell[vcell2[neip[1]]].temp*contrib;
	  cicoct[3]=cicoct[1];
	}
	else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
	  accel+=cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].temp*contrib;
	  cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
	}
	else{ // the particle is not deep enough
	  accel=0.;
	  hres=0;
	}
      }
    }
  }
      

  //contrib from 101 cell 2 paths ======================================================
  contrib=tx*dy*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[1]!=NULL){
	getcellnei(ciccell[1], vnei2, vcell2);
	ciccell[5]=vcell2[neip[2]];
	if(vnei2[neip[2]]==6){
	  accel+=cicoct[1]->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[5]=cicoct[1];
	}
	else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
	  accel+=cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
	}
	else{ // the particle is not deep enough
	  accel=0.;
	  hres=0;
	}
	
      }
    }
  }

  //contrib from 011 cell 2 paths ======================================================
  contrib=dx*ty*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[2]!=NULL){
	getcellnei(ciccell[2], vnei2, vcell2);
	ciccell[6]=vcell2[neip[2]];
	if(vnei2[neip[2]]==6){
	  accel+=cicoct[2]->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[6]=cicoct[2];
	}
	else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
	  accel+=cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
	}
	else{ // the particle is not deep enough
	  accel=0.;
	  hres=0;
	}

      }
    }
  }
      
  // contrib from 111 ======================================================
  contrib=tx*ty*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){
		  
      if(cicoct[3]!=NULL){
	getcellnei(ciccell[3], vnei2, vcell2);
	ciccell[7]=vcell2[neip[2]];
	if(vnei2[neip[2]]==6){
	  accel+=cicoct[3]->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[7]=cicoct[3];
	}
	else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
	  accel+=cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
	}
	else{ // the particle is not deep enough
	  accel=0.;
	  hres=0;
	}
      }
    }
  }
  

  // we must recompute the force from level-1 if the particle is not deep enough LOWRES
  if(hres!=1){
    
    // Getting the new curoct at low resolution
    curoctlr=cell2oct(curoct->parent);
    icell=curoct->parent->idx;
    dxcur=1./pow(2,curoctlr->level);

    // start again
    
    xc=curoctlr->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
    yc=curoctlr->y+((icell/2)%2)*dxcur+dxcur/2;
    zc=curoctlr->z+( icell/4   )*dxcur+dxcur/2; 
    
    // getting the neighbors
    getcellnei(icell, vnei, vcell);
    
    // here we compute the indexes of the direct neighbors which are involved in the cic
    neip[0]=(curp->x<xc?0:1);
    neip[1]=(curp->y<yc?2:3);
    neip[2]=(curp->z<zc?4:5);
    
    // here we denote the offset in ZYX
    ciccell[0]=icell;          //cell 000
    ciccell[1]=vcell[neip[0]]; //cell 001
    ciccell[2]=vcell[neip[1]]; //cell 010
    ciccell[4]=vcell[neip[2]]; //cell 100
    cicoct[0]=curoctlr;

    // the CIC weights
    tx=fabs((curp->x-xc)/dxcur);
    ty=fabs((curp->y-yc)/dxcur);
    tz=fabs((curp->z-zc)/dxcur);
    
    dx=1.-tx;
    dy=1.-ty;
    dz=1.-tz;
    
    // contrib from current cell 000 =====
    contrib=dx*dy*dz;
    if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoctlr->cell[icell].temp*contrib;
      //curoct->cell[icell].density+=contrib/pow(dxcur,3);
    }
  

    // contrib from 100 cell ===========================================================
      
    contrib=tx*dy*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[0]]==6){
	accel+=curoctlr->cell[vcell[neip[0]]].temp*contrib;
	cicoct[1]=curoctlr;
      }
      else {
	if(curoctlr->nei[vnei[neip[0]]]->child!=NULL){
	  accel+=curoctlr->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].temp*contrib;
	  cicoct[1]=curoctlr->nei[vnei[neip[0]]]->child;
	}
      }
    }

    // contrib from 010 cell ===========================================================
    contrib=dx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[1]]==6){
	accel+=curoctlr->cell[vcell[neip[1]]].temp*contrib;
	cicoct[2]=curoctlr;
      }
      else{
	if(curoctlr->nei[vnei[neip[1]]]->child!=NULL){
	  accel+=curoctlr->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].temp*contrib;
	  cicoct[2]=curoctlr->nei[vnei[neip[1]]]->child;
	}
      }
    }
    
    // contrib from 001 cell ===========================================================

    contrib=dx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[2]]==6){
	accel+=curoctlr->cell[vcell[neip[2]]].temp*contrib;
	cicoct[4]=curoctlr;
      }
      else{
	if(curoctlr->nei[vnei[neip[2]]]->child!=NULL){
	  accel+=curoctlr->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].temp*contrib;
	  cicoct[4]=curoctlr->nei[vnei[neip[2]]]->child;
	}
      }
    }

    //contrib from 110 cell 2 paths ======================================================
  
    contrib=tx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(cicoct[1]!=NULL){
	getcellnei(ciccell[1], vnei2, vcell2);
	ciccell[3]=vcell2[neip[1]];
	if(vnei2[neip[1]]==6){
	  accel+=cicoct[1]->cell[vcell2[neip[1]]].temp*contrib;
	  cicoct[3]=cicoct[1];
	}
	else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
	  accel+=cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].temp*contrib;
	  cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
	}
      }
    }
    
    //contrib from 101 cell 2 paths ======================================================
    contrib=tx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[1]!=NULL){
	getcellnei(ciccell[1], vnei2, vcell2);
	ciccell[5]=vcell2[neip[2]];
	if(vnei2[neip[2]]==6){
	  accel+=cicoct[1]->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[5]=cicoct[1];
	}
	else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
	  accel+=cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
	}
      }
    }

    //contrib from 011 cell 2 paths ======================================================
    contrib=dx*ty*tz;
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[2]!=NULL){
	getcellnei(ciccell[2], vnei2, vcell2);
	ciccell[6]=vcell2[neip[2]];
	if(vnei2[neip[2]]==6){
	  accel+=cicoct[2]->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[6]=cicoct[2];
	}
	else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
	  accel+=cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
	}
      }
    }
      
    // contrib from 111 ======================================================
    contrib=tx*ty*tz;
    if((contrib<=1.)&&(contrib>=0)){
		  
      if(cicoct[3]!=NULL){
	getcellnei(ciccell[3], vnei2, vcell2);
	ciccell[7]=vcell2[neip[2]];
	if(vnei2[neip[2]]==6){
	  accel+=cicoct[3]->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[7]=cicoct[3];
	}
	else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
	  accel+=cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
	  cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
	}
      }
    }
  }

  
  // Once we have the acceleration we can compute the velocity
  switch(dir){
  case(0):
    curp->vx+=-accel*dt;
    break;
  case(1):
    curp->vy+=-accel*dt;
    break;
  case(2):
    curp->vz+=-accel*dt;
    break;
  }

}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

struct OCT* gathercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int icell;
  int vnei[6],vcell[6];
  int ioct=0;
  float ominterp=0.2;
  nextoct=octstart;
  int ncell=0;

  (*nread)=0; // nread will return the effective number of CELLS read (important for reductions)

  if(nextoct!=NULL){
    do{// sweeping level
	ioct++;
	curoct=nextoct;
	nextoct=curoct->next;

	if(cpu->rank!=curoct->cpu) continue; // we consider only current cpu octs

	for(icell=0;icell<8;icell++){
	  ncell++;
	  getcellnei(icell, vnei, vcell);
	  switch(var){
	    //=================================
	  case 0: //------>  density
	    if(nei==6){
	      vec[(*nread)]=curoct->cell[icell].density;
	    }
	    else{
	      if(vnei[nei]==6){
		vec[(*nread)]=curoct->cell[vcell[nei]].density;
	      }
	      else{
		if(curoct->nei[vnei[nei]]->child!=NULL){
		  vec[(*nread)]=curoct->nei[vnei[nei]]->child->cell[vcell[nei]].density;
		}
		else{
		  vec[(*nread)]=curoct->nei[vnei[nei]]->density;
		}
	      }
	    }
	    break;
	    //=================================
	  case 1: //------>  pot
	    if(nei==6){
	      vec[(*nread)]=curoct->cell[icell].pot;
	    }
	    else{
	      if(vnei[nei]==6){
		vec[(*nread)]=curoct->cell[vcell[nei]].pot;
	      }
	      else{
		if(curoct->nei[vnei[nei]]->child!=NULL){
		  vec[(*nread)]=curoct->nei[vnei[nei]]->child->cell[vcell[nei]].pot;
		}
		else{
		  vec[(*nread)]=curoct->nei[vnei[nei]]->pot*(1.-ominterp)+curoct->cell[icell].pot*(ominterp);
		}
	      }
	    }
	    break;
	    //=================================
	  case 2: //------>  TEMP
	    if(nei==6){
	      vec[(*nread)]=curoct->cell[icell].temp;
	    }
	    else{
	      if(vnei[nei]==6){
		vec[(*nread)]=curoct->cell[vcell[nei]].temp;
	      }
	      else{
		if(curoct->nei[vnei[nei]]->child!=NULL){
		  vec[(*nread)]=curoct->nei[vnei[nei]]->child->cell[vcell[nei]].temp;
		}
		else{
		  vec[(*nread)]=curoct->nei[vnei[nei]]->temp;
		}
	      }
	    }
	    break;

	  }
	  (*nread)++;
	}
    }while((nextoct!=NULL)&&((*nread)<stride));
  }

  return nextoct;
  
}



struct OCT* scattercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int nread=0;
  int icell;
  int vnei[6],vcell[6];
  int ioct=0;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{// sweeping level
	ioct++;
	curoct=nextoct;
	nextoct=curoct->next;

	if(cpu->rank!=curoct->cpu) continue;

	for(icell=0;icell<8;icell++){

	  getcellnei(icell, vnei, vcell);

	  switch(var){
	    //==============================================
	  case 0: //------>  density
	    if(nei==6){
	      curoct->cell[icell].density=vec[nread];
	    }
	    else{
	      if(vnei[nei]==6){
		curoct->cell[vcell[nei]].density=vec[nread];
	      }
	      else{
		if(curoct->nei[vnei[nei]]->child!=NULL){
		  curoct->nei[vnei[nei]]->child->cell[vcell[nei]].density=vec[nread];
		}
		else{
		  curoct->nei[vnei[nei]]->density=vec[nread];
		}
	      }
	    }
	    break;
	    //==============================================
	  case 1: //------>  potential
	    if(nei==6){
	      curoct->cell[icell].pot=vec[nread];
	    }
	    else{
	      if(vnei[nei]==6){
		curoct->cell[vcell[nei]].pot=vec[nread];
	      }
	      else{
		if(curoct->nei[vnei[nei]]->child!=NULL){
		  curoct->nei[vnei[nei]]->child->cell[vcell[nei]].pot=vec[nread];
		}
		else{
		  curoct->nei[vnei[nei]]->pot=vec[nread];
		}
	      }
	    }
	    break;
	    //==============================================
	  case 2: //------>  temp
	    if(nei==6){
	      curoct->cell[icell].temp=vec[nread];
	    }
	    else{
	      if(vnei[nei]==6){
	  	curoct->cell[vcell[nei]].temp=vec[nread];
	      }
	      else{
	  	if(curoct->nei[vnei[nei]]->child!=NULL){
	  	  curoct->nei[vnei[nei]]->child->cell[vcell[nei]].temp=vec[nread];
	  	}
	  	else{
	  	  curoct->nei[vnei[nei]]->temp=vec[nread];
	  	}
	      }
	    }
	    break;
	  }
	  nread++;
	}
      }while((nextoct!=NULL)&&(nread<stride));
  }
  return nextoct;
}
  

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
 

void laplacian(float **vcomp, int stride, float dx){

  int i;
  float temp;
  float om=1.0;
  for(i=0;i<stride;i++){
    temp=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i])/6.-dx*dx*vcomp[7][i]/6.*4*M_PI;
    vcomp[6][i]=vcomp[6][i]*(1.-om)+om*temp;
  }
}

void grad(float **vcomp, int stride, float dx, int dir){

  int i;
  float temp;
  float om=1.;
  for(i=0;i<stride;i++){
    switch(dir){
    case 0:
      temp=-0.5*(vcomp[0][i]-vcomp[1][i])/dx;
      break;
    case 1:
      temp=-0.5*(vcomp[2][i]-vcomp[3][i])/dx;
      break;
    case 2:
      temp=-0.5*(vcomp[4][i]-vcomp[5][i])/dx;
      break;
    }
    vcomp[6][i]=temp;
  }
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

void remove_avg(float *vcomp, int stride, float avg)
{
  int i;
  for(i=0;i<stride;i++){
    vcomp[i]=vcomp[i]-avg;
  }
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

float square(float *vcomp, int stride)
{
  int i;
  float res=0.;
  for(i=0;i<stride;i++){
    res+=vcomp[i]*vcomp[i];
  }
  return res;
}

float square_res(float **vcomp, int stride, float dx)
{
  int i;
  float res=0.;
  float res2=0.;
  for(i=0;i<stride;i++){
    res2=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i]-6.*vcomp[6][i])/(dx*dx)-vcomp[7][i]*4.*M_PI;
    res+=res2*res2;
  }
  return res;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

void movepart(int levelcoarse,int levelmax,struct OCT** firstoct, float dt){
  
  int level;
  float mdisp;
  struct OCT *nextoct;
  struct OCT oct;
  float dxcur;
  int icell;
  struct PART *nexp;
  struct PART *curp;
  float disp;

  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      printf("level=%d\n",level);
      mdisp=0.;

      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      if(nextoct==NULL) continue; // in case the level is empty
      do // sweeping through the octs of level
	{
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell

	      if(nexp!=NULL){
		do{ 
		  curp=nexp; 
		  nexp=curp->next; 
		 
		  // particle displacement
		  curp->x+=curp->vx*dt;
		  curp->y+=curp->vy*dt;
		  curp->z+=curp->vz*dt;
		  disp=sqrt(curp->vx*dt*curp->vx*dt+curp->vy*dt*curp->vy*dt+curp->vz*dt*curp->vz*dt);
		  if(disp>mdisp) mdisp=disp;

		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
      printf("level=%d maxdisp=%f or %f dx\n",level,mdisp,mdisp/dxcur);
    }
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
  float dxcur;
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
				dxcur=1./pow(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				xp=0;
				yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
				
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);

			      }

			  }
			  else{
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
				dxcur=1./pow(2.,newoct->level);

				// new cell coordinates
				//ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				xp=0;
				yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
				
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				if((ip>7)||(ip<0)) ip=0; 

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
				dxcur=1./pow(2.,newoct->level);

				xp=1;
				yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
			      

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
			      
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				if((ip>7)||(ip<0)) ip=0; 

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
				dxcur=1./pow(2.,newoct->level);

				xp=1;
				yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;
			      

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
			      
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				if((ip>7)||(ip<0)) ip=0; 

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
				dxcur=1./pow(2.,newoct->level);

				xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=0;
				zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				if((ip>7)||(ip<0)) ip=0; 
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
				dxcur=1./pow(2.,newoct->level);

				xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=0;
				zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				if((ip>7)||(ip<0)) ip=0; 
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
				dxcur=1./pow(2.,newoct->level);
				xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=1;
				zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				if((ip>7)||(ip<0)) ip=0; 
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
				dxcur=1./pow(2.,newoct->level);
				xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=1;
				zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				if((ip>7)||(ip<0)) ip=0; 
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
			      dxcur=1./pow(2.,newoct->level);
			      
			      xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			      yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			      zp=0;
			      ip=xp+yp*2+zp*4;
			      
			      //ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
			      
			      // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			      if((ip>7)||(ip<0)) ip=0; 
			      
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
				dxcur=1./pow(2.,newoct->level);

				xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=0;
				ip=xp+yp*2+zp*4;

				//ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;

				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				if((ip>7)||(ip<0)) ip=0; 

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
			      dxcur=1./pow(2.,newoct->level);
			      
			      xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
			      yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
			      zp=1;
			      ip=xp+yp*2+zp*4;
			      
			      //				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
			      
			      // some particles will experience more than one displacement (along diagonals) we store them in cell 0
			      if((ip>7)||(ip<0)) ip=0; 
			      
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
				dxcur=1./pow(2.,newoct->level);
				
				xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
				yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
				zp=1;
				ip=xp+yp*2+zp*4;

				//				ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;

				// some particles will experience more than one displacement (along diagonals) we store them in cell 0
				if((ip>7)||(ip<0)) ip=0; 

				newcell=&(newoct->cell[ip]);
			      }
			    }
			  }
			}
			/* if(curp->idx==82279){ */
			/*   printf("PPZ %f %f %f %d %p %p out=%d\n",curp->x,curp->y,curp->z,ip,curoct,newoct,out); */
			/* } */
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
//------------------------------------------------------------------------


void cic_child(struct OCT* oct,struct OCT* octorg, int icellorg)
{
  struct PART *nexp;
  struct PART *curp;
  int icell;
  
  for(icell=0;icell<8;icell++){
    if(oct->cell[icell].child==NULL){
      nexp=oct->cell[icell].phead; //sweeping the particles of the current cell */
      if(nexp!=NULL){ 
	do{  
	  curp=nexp; 
	  nexp=curp->next; 
	  part2cell_cic(curp, octorg, icellorg,1); 
	}while(nexp!=NULL); 
      }
    }
    else{
      cic_child(oct->cell[icell].child,octorg,icellorg);
    }
  }
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

void fixbound_reg(struct OCT * grid, int levelcoarse, int nbnd){
  
  int level;
  int ci,cj,ck,cur,il;
  int nbnd2;
  int ic;
  nbnd2=2*nbnd;

  for(level=1;level<=levelcoarse;level++){ // sweeping the levels from l=1 to l=levelcoarse
    int nxoct=pow(2,level-1); // number of octs along one direction
    int firstoct_currl=0;
    for(il=1;il<level;il++) firstoct_currl+=pow(pow(2,il-1)+nbnd2,3); // the index of the first oct of current level
    int corg,curorg;
    // z minus
    ck=-1;
    corg=nxoct-1;
    for(cj=0;cj<nxoct;cj++){
      for(ci=0;ci<nxoct;ci++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+corg)*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	memcpy(grid+cur,grid+curorg,sizeof(struct OCT));
	for(ic=0;ic<8;ic++){
	  grid[cur].cell[ic].phead=NULL; // nullifying the particles in order to host properly the new ones
	}
      }
    }

    // z plus
    ck=nxoct;
    corg=0;
    for(cj=0;cj<nxoct;cj++){
      for(ci=0;ci<nxoct;ci++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+corg)*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	memcpy(grid+cur,grid+curorg,sizeof(struct OCT));
	for(ic=0;ic<8;ic++){
	  grid[cur].cell[ic].phead=NULL; // nullifying the particles in order to host properly the new ones
	}
      }
    }

    // y minus
    cj=-1;
    corg=nxoct-1;
    for(ck=0;ck<nxoct;ck++){
      for(ci=0;ci<nxoct;ci++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+ci  )+(nbnd+corg)*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	memcpy(grid+cur,grid+curorg,sizeof(struct OCT));
	for(ic=0;ic<8;ic++){
	  grid[cur].cell[ic].phead=NULL; // nullifying the particles in order to host properly the new ones
	}
      }
    }

    // y plus
    cj=nxoct;
    corg=0;
    for(ck=0;ck<nxoct;ck++){
      for(ci=0;ci<nxoct;ci++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+ci  )+(nbnd+corg)*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	memcpy(grid+cur,grid+curorg,sizeof(struct OCT));
	for(ic=0;ic<8;ic++){
	  grid[cur].cell[ic].phead=NULL; // nullifying the particles in order to host properly the new ones
	}
      }
    }

    // x minus
    ci=-1;
    corg=nxoct-1;
    for(ck=0;ck<nxoct;ck++){
      for(cj=0;cj<nxoct;cj++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+corg)+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	memcpy(grid+cur,grid+curorg,sizeof(struct OCT));
	for(ic=0;ic<8;ic++){
	  grid[cur].cell[ic].phead=NULL; // nullifying the particles in order to host properly the new ones
	}
      }
    }

    // x plus
    ci=nxoct;
    corg=0;
    for(ck=0;ck<nxoct;ck++){
      for(cj=0;cj<nxoct;cj++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+corg)+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	memcpy(grid+cur,grid+curorg,sizeof(struct OCT));
	for(ic=0;ic<8;ic++){
	  grid[cur].cell[ic].phead=NULL; // nullifying the particles in order to host properly the new ones
	}
      }
    }
  }

}

//================================================================================
//================================================================================

void fixbound_cic(struct OCT * grid, int levelcoarse, int nbnd){


  int level;
  int ci,cj,ck,cur,il;
  int nbnd2;
  int ic;
  struct PART *part;
  nbnd2=2*nbnd;

  for(level=1;level<=levelcoarse;level++){ // sweeping the levels from l=1 to l=levelcoarse

    int nxoct=pow(2,level-1); // number of octs along one direction
    
    int firstoct_currl=0;
    for(il=1;il<level;il++) firstoct_currl+=pow(pow(2,il-1)+nbnd2,3); // the index of the first oct of current level
    int corg,curorg;
    // z minus
    ck=-1;
    corg=nxoct-1;
    for(cj=-1;cj<nxoct+1;cj++){
      for(ci=-1;ci<nxoct+1;ci++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+corg)*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	for(ic=0;ic<8;ic++){
	  grid[curorg].cell[ic].density+=grid[cur].cell[ic].density;

	  //dealing with particles (if any)
	  if(grid[cur].cell[ic].phead!=NULL){
	    modifpospart(grid[cur].cell[ic].phead,1.,2); // first we periodize the positions of the particles in boundaries
	    if(grid[curorg].cell[ic].phead!=NULL){
	      part=findlastpart(grid[curorg].cell[ic].phead); // we search for the last part of the cell to be modified
	      part->next=grid[cur].cell[ic].phead;
	      grid[cur].cell[ic].phead->prev=part;
	    }
	    else{
	      grid[curorg].cell[ic].phead=grid[cur].cell[ic].phead; // the new head is the old head (the inner cell gets allt the particles
	    }
	    grid[cur].cell[ic].phead=NULL; // we erase all the particles from the boundaries
	  }
	  // End - particles boundaries
	}
      }
    }

    // z plus
    ck=nxoct;
    corg=0;
    for(cj=-1;cj<nxoct+1;cj++){
      for(ci=-1;ci<nxoct+1;ci++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+corg)*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	for(ic=0;ic<8;ic++){
	  grid[curorg].cell[ic].density+=grid[cur].cell[ic].density;
	 
	  //dealing with particles (if any)
	  if(grid[cur].cell[ic].phead!=NULL){
	    modifpospart(grid[cur].cell[ic].phead,-1.,2); // first we periodize the positions of the particles in boundaries
	    if(grid[curorg].cell[ic].phead!=NULL){
	      part=findlastpart(grid[curorg].cell[ic].phead); // we search for the last part of the cell to be modified
	      part->next=grid[cur].cell[ic].phead;
	      grid[cur].cell[ic].phead->prev=part;
	    }
	    else{
	      grid[curorg].cell[ic].phead=grid[cur].cell[ic].phead; // the new head is the old head (the inner cell gets allt the particles
	    }
	    grid[cur].cell[ic].phead=NULL; // we erase all the particles from the boundaries
	  }
	  // End - particles boundaries
	}
      }
    }

    // y minus
    cj=-1;
    corg=nxoct-1;
    for(ck=-1;ck<nxoct+1;ck++){
      for(ci=-1;ci<nxoct+1;ci++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+ci  )+(nbnd+corg)*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	for(ic=0;ic<8;ic++){
	  grid[curorg].cell[ic].density+=grid[cur].cell[ic].density;

	  //dealing with particles (if any)
	  if(grid[cur].cell[ic].phead!=NULL){
	    modifpospart(grid[cur].cell[ic].phead,1.,1); // first we periodize the positions of the particles in boundaries
	    if(grid[curorg].cell[ic].phead!=NULL){
	      part=findlastpart(grid[curorg].cell[ic].phead); // we search for the last part of the cell to be modified
	      part->next=grid[cur].cell[ic].phead;
	      grid[cur].cell[ic].phead->prev=part;
	    }
	    else{
	      grid[curorg].cell[ic].phead=grid[cur].cell[ic].phead; // the new head is the old head (the inner cell gets allt the particles
	    }
	    grid[cur].cell[ic].phead=NULL; // we erase all the particles from the boundaries
	  }
	  // End - particles boundaries

	}
      }
    }

    // y plus
    cj=nxoct;
    corg=0;
    for(ck=-1;ck<nxoct+1;ck++){
      for(ci=-1;ci<nxoct+1;ci++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+ci  )+(nbnd+corg)*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	for(ic=0;ic<8;ic++){
	  grid[curorg].cell[ic].density+=grid[cur].cell[ic].density;
	  //dealing with particles (if any)
	  if(grid[cur].cell[ic].phead!=NULL){
	    modifpospart(grid[cur].cell[ic].phead,-1.,1); // first we periodize the positions of the particles in boundaries
	    if(grid[curorg].cell[ic].phead!=NULL){
	      part=findlastpart(grid[curorg].cell[ic].phead); // we search for the last part of the cell to be modified
	      part->next=grid[cur].cell[ic].phead;
	      grid[cur].cell[ic].phead->prev=part;
	    }
	    else{
	      grid[curorg].cell[ic].phead=grid[cur].cell[ic].phead; // the new head is the old head (the inner cell gets allt the particles
	    }
	    grid[cur].cell[ic].phead=NULL; // we erase all the particles from the boundaries
	  }
	  // End - particles boundaries

	}
      }
    }

    // x minus
    ci=-1;
    corg=nxoct-1;
    for(ck=-1;ck<nxoct+1;ck++){
      for(cj=-1;cj<nxoct+1;cj++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+corg)+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	for(ic=0;ic<8;ic++){
	  grid[curorg].cell[ic].density+=grid[cur].cell[ic].density;
	  //dealing with particles (if any)
	  if(grid[cur].cell[ic].phead!=NULL){
	    modifpospart(grid[cur].cell[ic].phead,1.,0); // first we periodize the positions of the particles in boundaries
	    if(grid[curorg].cell[ic].phead!=NULL){
	      part=findlastpart(grid[curorg].cell[ic].phead); // we search for the last part of the cell to be modified
	      part->next=grid[cur].cell[ic].phead;
	      grid[cur].cell[ic].phead->prev=part;
	    }
	    else{
	      grid[curorg].cell[ic].phead=grid[cur].cell[ic].phead; // the new head is the old head (the inner cell gets allt the particles
	    }
	    grid[cur].cell[ic].phead=NULL; // we erase all the particles from the boundaries
	  }
	  // End - particles boundaries
	}
      }
    }

    // x plus
    ci=nxoct;
    corg=0;
    for(ck=-1;ck<nxoct+1;ck++){
      for(cj=-1;cj<nxoct+1;cj++){
	cur   =(nbnd+ci  )+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	curorg=(nbnd+corg)+(nbnd+cj  )*(nxoct+nbnd2)+(nbnd+ck  )*(nxoct+nbnd2)*(nxoct+nbnd2)+firstoct_currl;
	for(ic=0;ic<8;ic++){
	  grid[curorg].cell[ic].density+=grid[cur].cell[ic].density;
	  //dealing with particles (if any)
	  if(grid[cur].cell[ic].phead!=NULL){
	    modifpospart(grid[cur].cell[ic].phead,-1.,0); // first we periodize the positions of the particles in boundaries
	    if(grid[curorg].cell[ic].phead!=NULL){
	      part=findlastpart(grid[curorg].cell[ic].phead); // we search for the last part of the cell to be modified
	      part->next=grid[cur].cell[ic].phead;
	      grid[cur].cell[ic].phead->prev=part;
	    }
	    else{
	      grid[curorg].cell[ic].phead=grid[cur].cell[ic].phead; // the new head is the old head (the inner cell gets allt the particles
	    }
	    grid[cur].cell[ic].phead=NULL; // we erase all the particles from the boundaries
	  }
	  // End - particles boundaries
	}
      }
    }
  }

}

//------------------------------------------------------------------------

void dumpmap(int lmap,struct OCT **firstoct,int field,char filename[],float zmin, float zmax)
{
  float *map;
  int imap,jmap;
  int nmap=pow(2,lmap);
  float dxmap=1./nmap,dxcur;
  int icur,ii,jj,ic,icell;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  FILE *fp;
  float xc,yc,zc;

  map=(float *)calloc(nmap*nmap,sizeof(float));

  //printf("==>  start map \n");
  for(level=1;level<=lmap;level++) // looping over octs
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];

      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
	  nextoct=oct.next;
	  //	  printf("%f %f %f\n",oct.x,oct.y,oct.z);
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if((oct.cell[icell].child==NULL)||(oct.level==lmap))
		{
		  xc=oct.x+( icell   %2)*dxcur;//+0.5*dxcur;
		  yc=oct.y+((icell/2)%2)*dxcur;//+0.5*dxcur;
		  zc=oct.z+( icell/4   )*dxcur;//+0.5*dxcur;
		  imap=xc/dxmap;
		  jmap=yc/dxmap;
		  
		  if((zc>zmax)||(zc<zmin)) continue;

		  //if(grid[icur].level==lmap) printf("%d %f %f %d %d %f\n",icell,xc,yc,imap,jmap,grid[icur].dens[icell]);
		  for(jj=0;jj<pow(2,lmap-oct.level);jj++)
		    {
		      for(ii=0;ii<pow(2,lmap-oct.level);ii++)
			{

			  switch(field){
			  case 0:
			    map[(imap+ii)+(jmap+jj)*nmap]=fmax(oct.level,map[(imap+ii)+(jmap+jj)*nmap]);
			    break;
			  case 1:
			    map[(imap+ii)+(jmap+jj)*nmap]+=oct.cell[icell].density*pow(2,lmap-oct.level);
			    break;
			  case 2:
			    map[(imap+ii)+(jmap+jj)*nmap]+=oct.cell[icell].pot*pow(2,lmap-oct.level);
			    break;
			  }
			  //map[(imap+ii)+(jmap+jj)*nmap]=fmin(oct.cell[icell].pot,map[(imap+ii)+(jmap+jj)*nmap]);
			  //map[(imap+ii)+(jmap+jj)*nmap]=fmax(oct.cell[icell].temp*oct.cell[icell].temp,map[(imap+ii)+(jmap+jj)*nmap]);
			  //map[(imap+ii)+(jmap+jj)*nmap]+=oct.cell[icell].pot*pow(2,lmap-oct.level);
			  //map[(imap+ii)+(jmap+jj)*nmap]+=oct.cell[icell].density*pow(2,lmap-oct.level);
			  //map[(imap+ii)+(jmap+jj)*nmap]=fmax(oct.level,map[(imap+ii)+(jmap+jj)*nmap]);
			  //map[(imap+ii)+(jmap+jj)*nmap]=fmax(oct.cell[icell].marked,map[(imap+ii)+(jmap+jj)*nmap]);
			}
		    }
		}
	    }
	}while(nextoct!=NULL);
    }
  
  
  //============= dump

  //printf("dumping %s\n",filename);
  fp=fopen(filename,"wb");
  fwrite(&nmap,1,sizeof(int),fp);
  fwrite(map,nmap*nmap,sizeof(float),fp);
  fclose(fp);
  free(map);

}

void dumpcube(int lmap,struct OCT **firstoct,int field,char filename[])
{
  float *map;
  int imap,jmap,kmap;
  int nmap=pow(2,lmap);
  float dxmap=1./nmap,dxcur;
  int icur,ii,jj,kk,ic,icell;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  FILE *fp;
  float xc,yc,zc;

  map=(float *)calloc(nmap*nmap*nmap,sizeof(float));

  //printf("==> start map \n");
  for(level=1;level<=lmap;level++) // looping over octs
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];

      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
	  nextoct=oct.next;
	  //	  printf("%f %f %f\n",oct.x,oct.y,oct.z);
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if((oct.cell[icell].child==NULL)||(oct.level==lmap))
		{
		  xc=oct.x+( icell   %2)*dxcur;//+0.5*dxcur;
		  yc=oct.y+((icell/2)%2)*dxcur;//+0.5*dxcur;
		  zc=oct.z+( icell/4   )*dxcur;//+0.5*dxcur;
		  imap=xc/dxmap;
		  jmap=yc/dxmap;
		  kmap=zc/dxmap;
		  
		  //if(grid[icur].level==lmap) printf("%d %f %f %d %d %f\n",icell,xc,yc,imap,jmap,grid[icur].dens[icell]);
		  for(kk=0;kk<pow(2,lmap-oct.level);kk++)
		    {
		      for(jj=0;jj<pow(2,lmap-oct.level);jj++)
			{
			  for(ii=0;ii<pow(2,lmap-oct.level);ii++)
			    {
			      
			      switch(field){
			      case 0:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.level;
				break;
			      case 1:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.cell[icell].density;
				break;
			      case 2:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.cell[icell].pot;
				break;
			      case 3:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cpu;
				break;
			      case 4:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].marked;
				break;
			      case 5:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].temp;
				break;
			      }
			      //map[(imap+ii)+(jmap+jj)*nmap]=fmin(oct.cell[icell].pot,map[(imap+ii)+(jmap+jj)*nmap]);
			      //map[(imap+ii)+(jmap+jj)*nmap]=fmax(oct.cell[icell].temp*oct.cell[icell].temp,map[(imap+ii)+(jmap+jj)*nmap]);
			      //map[(imap+ii)+(jmap+jj)*nmap]+=oct.cell[icell].pot*pow(2,lmap-oct.level);
			      //map[(imap+ii)+(jmap+jj)*nmap]+=oct.cell[icell].density*pow(2,lmap-oct.level);
			      //map[(imap+ii)+(jmap+jj)*nmap]=fmax(oct.level,map[(imap+ii)+(jmap+jj)*nmap]);
			      //map[(imap+ii)+(jmap+jj)*nmap]=fmax(oct.cell[icell].marked,map[(imap+ii)+(jmap+jj)*nmap]);
			    }
			}
		    }
		}
	    }
	}while(nextoct!=NULL);
    }
  
  //============= dump
  
  //printf("dumping %s\n",filename);
  fp=fopen(filename,"wb");
  fwrite(&nmap,1,sizeof(int),fp);
  fwrite(map,nmap*nmap*nmap,sizeof(float),fp);
  fclose(fp);
  free(map);

}
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------


void refine_cells(int levelcoarse, int levelmax, struct OCT **firstoct, struct OCT ** lastoct, struct OCT * endoct, struct CPUINFO *cpu)
{
  int nref,ndes;
  struct OCT *newoct;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *desoct;
  int level;
  float dxcur;
  int icell;
  int sump,sum2;
  int ii;
  int vnei[6],vcell[6];
  int ip,xp,yp,zp;
  struct PART *curploc;
  struct PART *nexploc;
  struct PART *nexp;
  struct PART *curp;

  //if(nsteps==1) abort();
  nref=0;
  ndes=0;

  newoct=endoct+1; // the new oct will be the next to the last one

  printf("==> start refining on cpu %d lcoarse=%d lmax=%d\n",cpu->rank,levelcoarse,levelmax);
  //breakmpi();
  //  for(level=levelcoarse;level<levelmax;level++) // looping over levels
  for(level=1;level<levelmax;level++) // looping over levels
    {
      dxcur=1./pow(2,level);
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if(newoct->level!=0) {
		printf("newoct level=%d\n",newoct->level);
		abort();
	      }

	      // destroying octs with level>levelcoarse ==========================
	      if(((curoct->cell[icell].child!=NULL)&&(curoct->cell[icell].marked==0))&&(curoct->level>levelcoarse)){
	      	ndes++;
		
		desoct=curoct->cell[icell].child; // the oct to destroy
		
		if(curoct->cell[icell].phead!=NULL){
	      	  printf("non void split cell !\n");
	      	  abort();
	      	}
		// we remove the child from the refined cell
	      	curoct->cell[icell].child=NULL;
		
	      	// we remove the parent from the oct to be destroyed
	      	desoct->parent=NULL;
		
	      	// we remove the oct from the list
		
	      	if(desoct->prev!=NULL){
	      	  desoct->prev->next=desoct->next;
	      	}
	      	else{
	      	  //the oct is the first one at this level
		  if(firstoct[level]!=desoct){
		    printf("oups 1\n");
		    abort();
		  }
	      	  firstoct[level]=desoct->next;
	      	}

	      	if(desoct->next!=NULL){
	      	  desoct->next->prev=desoct->prev;
	      	}
	      	else{
	      	  //the oct is the last one at this level
		  if(lastoct[level]!=desoct){
		    printf("oups 2\n");
		    abort();
		  }
		  lastoct[level]=desoct->prev;
	      	}

		desoct->level=0;
		
	      	//======== dealing with particles
	      	// we gather the particles of the oct to be destroyed in the parent cell
		
	      	curploc=findlastpart(curoct->cell[icell].phead); // we get the last particle from the parent cell (should be nil in principle)
	      	sump=0;
		sum2=0;
	      	for(ii=0;ii<8;ii++){ // sweeping the cells of desoct
	      	  nexp=desoct->cell[ii].phead;
	      	  sump+=countpart(desoct->cell[ii].phead);
	      	  if(curploc==NULL){
	      	    curoct->cell[icell].phead=nexp;
		    curploc=findlastpart(nexp);
	      	  }
	      	  else if(nexp!=NULL){
	      	    curploc->next=nexp;
	      	    nexp->prev=curploc;
		    curploc=findlastpart(nexp);
	      	  }
	      	}
	      	sum2=countpart(curoct->cell[icell].phead);
	      	if(sum2!=sump){
	      	  printf("sum2=%d sump=%d\n",sum2,sump);
	      	  abort();
	      	}
	      }
	      
	      // creation of a new oct ==================
	      if((curoct->cell[icell].child==NULL)&&(curoct->cell[icell].marked!=0)){
		nref++;
		
		// a new oct is created
		newoct->parent=&(curoct->cell[icell]);
		newoct->level=curoct->level+1;
		newoct->x=curoct->x+( icell   %2)*dxcur;
		newoct->y=curoct->y+((icell/2)%2)*dxcur;
		newoct->z=curoct->z+( icell   /4)*dxcur;
		
		// the new oct is connected to parent
		curoct->cell[icell].child=newoct;
		
		// it shares the same cpu
		newoct->cpu=curoct->cpu;

		//the neighbours
		getcellnei(icell, vnei, vcell);
		for(ii=0;ii<6;ii++){
		  if(vnei[ii]!=6){
		
		    if((curoct->nei[vnei[ii]]->child==NULL)&&(curoct->cpu==cpu->rank)){
		      printf("ouhla rank=%d curoct.cpu=%d\n",cpu->rank,curoct->cpu);
		      abort();
		    }
		    // Note that boundary octs are refined but may have missing neighbors
		    if(curoct->nei[vnei[ii]]->child!=NULL) newoct->nei[ii]=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
		  }else{
		    newoct->nei[ii]=&(curoct->cell[vcell[ii]]);
		  }
		}

		// filling the cells
		for(ii=0;ii<8;ii++){
		  newoct->cell[ii].marked=0;
		  newoct->cell[ii].child=NULL;
		  newoct->cell[ii].density=0.;
		  newoct->cell[ii].idx=ii;
		  newoct->cell[ii].phead=NULL;
		}


		// splitting the particles
#if 1
		nexp=curoct->cell[icell].phead;
		if(nexp!=NULL){
		  do{ //sweeping the particles of the current cell
		    curp=nexp;
		    nexp=curp->next;
		    
		    xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);xp=(xp==2?1:xp);
		    yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);yp=(yp==2?1:yp);
		    zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);zp=(zp==2?1:zp);
		    ip=xp+yp*2+zp*4;

		    //ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;
		    
		    if((ip<0)||(ip>7)){
		      printf("ip=%d idx=%d\n",ip,curp->idx);
		      abort();
		    }
		    
		    if(newoct->cell[ip].phead==NULL)
		      {
			// we create a new particle list in the cell
			newoct->cell[ip].phead=curp;
			//			newoct->cell[ip].density=8./(dxcur*dxcur*dxcur);

			// we remove the particle from the parent cell
			if(curp->prev!=NULL) curp->prev->next=curp->next; 
			if(curp->next!=NULL) curp->next->prev=curp->prev;
			if(curp->prev==NULL) curoct->cell[icell].phead=curp->next; // we removed the head
			curp->prev=NULL;
			curp->next=NULL;
		      }
		    else{
		      //newoct->cell[ip].density+=8./(dxcur*dxcur*dxcur);
		      // we remove the particle from the parent cell
		      if(curp->prev!=NULL) curp->prev->next=curp->next;
		      if(curp->next!=NULL) curp->next->prev=curp->prev;
		      if(curp->prev==NULL) curoct->cell[icell].phead=curp->next; // we removed the head
		      curp->next=NULL;	
		      
		      nexploc=newoct->cell[ip].phead;
		      do{ //sweeping the particles of the new cell
			curploc=nexploc;
			nexploc=curploc->next;
		      }while(nexploc!=NULL);
		      curploc->next=curp;
		      curp->prev=curploc;
		    }
		  }while(nexp!=NULL);
		}

		if(curoct->cell[icell].phead!=NULL){
		  printf("cell not emptied after split !\n");
		  abort();
		};
#endif
		

		// preparing the next creations on level+1
		newoct->next=NULL;

		if(firstoct[level]==NULL){
		  firstoct[level]=newoct;
		  newoct->prev=NULL;
		}
		else{
		  newoct->prev=lastoct[level];
		  lastoct[level]->next=newoct;
		}
		lastoct[level]=newoct;
		newoct++;
	      }

	      
	    }
	  //printf("nextoct=%p endoct=%p\n",nextoct,endoct);
	}while(nextoct!=NULL);
      //printf("level=%d done\n",level);
    }

  printf("octs created by rank %d =%d\n",cpu->rank,nref);
  printf("octs destroyed by rank %d=%d\n",cpu->rank,ndes);

}


//------------------------------------------------------------------------
  //------------------------------------------------------------------------
void dumppart(struct OCT **firstoct,char filename[],int npart, int levelcoarse, int levelmax){

  FILE *fp;
  float val;
  int vali;
  int ipart=0;
  int level;
  struct OCT *nextoct;
  struct OCT oct;
  float dxcur;
  struct PART *nexp;
  struct PART *curp;
  int icell;

  fp=fopen(filename,"wb");
  fwrite(&npart,1,sizeof(int),fp);
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      //printf("level=%d\n",level);
      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
		do{ 
		  curp=nexp; 
		  nexp=curp->next; 
		  val=curp->x;fwrite(&val,1,sizeof(float),fp);
		  val=curp->y;fwrite(&val,1,sizeof(float),fp);
		  val=curp->z;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vy;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vz;fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->idx);fwrite(&val,1,sizeof(float),fp);
		  ipart++;
		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }
  fclose(fp);
  printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

//------------------------------------------------------------------------



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
void  multicheck(struct OCT **firstoct,int npart,int levelmax, int rank){

  int ntot;
  float ntotd;
  float nlevd;
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  float dx;
  int nlev,noct;
  int icell;
  struct PART *nexp;
  struct PART *curp;

  if(rank==0) printf("Check\n");
  ntot=0.;
  ntotd=0.;
  nlevd=0.;

  for(level=1;level<=levelmax;level++)
    {
      nextoct=firstoct[level-1];
      dx=1./pow(2,level);
      nlev=0;
      nlevd=0.;
      noct=0;
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      ntotd+=curoct->cell[icell].density*dx*dx*dx;
	      nlevd+=curoct->cell[icell].density*dx*dx*dx;
	      
	      nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell
	      if((curoct->cell[icell].child!=NULL)&&(curoct->cell[icell].phead!=NULL)){
		printf("check: split cell with particles !\n");
		abort();
	      }
	      if(nexp==NULL)continue;
	      do{ 
		nlev++;
		ntot++;
		curp=nexp;
		nexp=curp->next;
	      }while(nexp!=NULL);
	    }
	  noct++;
	}while(nextoct!=NULL);
      if(rank==0) printf("level=%d npart=%d npartd=%f noct=%d\n",level,nlev,nlevd,noct);
    }
  
  printf("CHECK==> RANK # %d total   npart=%d/%d npartd=%f\n",rank,ntot,npart,ntotd);
  if(ntot!=npart) {
    printf("particles number discrepancy ntot=%d npart=%d\n",ntot,npart);
    abort();
  }
}
 //------------------------------------------------------------------------
 //------------------------------------------------------------------------
#if 1
  // ==================================== performing the CIC assignement
  
void call_cic(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu){

  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  int inei,inei2;
  float dxcur;

  if(cpu->rank==0) printf("==> start CIC\n");

  // First we clean the density
  for(level=levelmax;level>=levelcoarse;level--)
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) curoct->cell[icell].density=0.;
	}while(nextoct!=NULL);
    }

  //second start CIC
  for(level=levelmax;level>=levelcoarse;level--)
    {
      dxcur=1./pow(2,level);
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  
	  // we skip octs which do not belong to the current CPU (they will be considered through mpi)
	  if(curoct->cpu!=cpu->rank) continue;
	  
	  //==== FIRST WE CONSIDER THE PARTICLES INSIDE THE BOUNDARIES AT LEVEL L
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
	      if(nexp!=NULL){ 
		do{  
		  curp=nexp; 
		  nexp=curp->next; 
		  part2cell_cic(curp, curoct, icell,1); 
		}while(nexp!=NULL); 
	      }

	      //==== SECOND WE CONSIDER THE PARTICLES INSIDE THE NEIGHBORS AT LEVEL L-1
	      // first the cartesian neighbors (6)
	      for(inei=0;inei<6;inei++){
		nexp=curoct->nei[inei]->phead; //sweeping the particles of the neighbour cell at level l-1
		if(nexp!=NULL){
		  do{ 
		    curp=nexp;
		    nexp=curp->next;
		    part2cell_cic(curp, curoct, icell,0);
		  }while(nexp!=NULL);
		}
	      }

	      // second the fundamental plane (4)
	      char dir[4]={2,3,1,0};
	      for(inei=0;inei<4;inei++){
		if(curoct->nei[inei]->child!=NULL){
		  nexp=curoct->nei[inei]->child->nei[dir[inei]]->phead; //sweeping the particles of the neighbour cell at level l-1
		  if(nexp!=NULL){
		    do{ 
		      curp=nexp;
		      nexp=curp->next;
		      part2cell_cic(curp, curoct, icell,0);
		    }while(nexp!=NULL);
		  }
		}
	      }

	      // third the top-bottom cross (8)
	      for(inei=5;inei<6;inei++){
		if(curoct->nei[inei]->child!=NULL){
		  for(inei2=0;inei2<4;inei2++){
		    nexp=curoct->nei[inei]->child->nei[inei2]->phead; //sweeping the particles of the neighbour cell at level l-1
		    if(nexp!=NULL){
		      do{ 
			curp=nexp;
			nexp=curp->next;
			part2cell_cic(curp, curoct, icell,0);
		      }while(nexp!=NULL);
		    }
		  }
		}
	      }

	      // fourth the plane of each top/bottom cross returns the corners (8)
	      
	      for(inei=5;inei<6;inei++){
		if(curoct->nei[inei]->child!=NULL){
		  for(inei2=0;inei2<4;inei2++){
		    if(curoct->nei[inei]->child->nei[inei2]->child!=NULL){
		      nexp=curoct->nei[inei]->child->nei[inei2]->child->nei[dir[inei2]]->phead; 
		      if(nexp!=NULL){
			do{ 
			  curp=nexp;
			  nexp=curp->next;
			  part2cell_cic(curp, curoct, icell,0);
			}while(nexp!=NULL);
		      }
		    }
		  }
		}
	      }
	      // THIRD WE LOOK FOR THE PARTICLES IN THE CHILD OCTS
	      if(curoct->cell[icell].child!=NULL){
		cic_child(curoct->cell[icell].child,curoct,icell);
	      }
	    }

	}while(nextoct!=NULL);
    }
  //  printf("great total=%f\n",toto);
#endif
}


 //------------------------------------------------------------------------

void myradixsort(int *a,int n)
{
  int *b;
  int i,m=0,exp=1;
  b=(int*)calloc(n,sizeof(int));
  for(i=0;i<n;i++)
    {
      if(a[i]>m)
	m=a[i];
    }
  
  while(m/exp>0)
    {
      int bucket[10]={0};
      for(i=0;i<n;i++)
	bucket[a[i]/exp%10]++;
      for(i=1;i<10;i++)
	bucket[i]+=bucket[i-1];
      for(i=n-1;i>=0;i--)
	b[--bucket[a[i]/exp%10]]=a[i];
      for(i=0;i<n;i++)
	a[i]=b[i];
      exp*=10;
    }		
  free(b);
}


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

 //------------------------------------------------------------------------
 //------------------------------------------------------------------------
 //------------------------------------------------------------------------

void  setup_mpi(struct CPUINFO *cpu, struct OCT **firstoct, int levelmax, int levelcoarse, int ngridmax){

  int nnei=0;
  int *mpinei; // the COMPLETE list of neighbors cpu (can be redundant)
  int *neicpu; // the reduced list of neighbors CPU (not redundant);
  int hidx;
  int inidx; // counts the number of inner boundaries octs
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *desoct;
  struct OCT *newoct;
  int key;
  int inei;
  int i,j;
  int nbnd;
  int icell;

  mpinei =(int*)calloc(ngridmax,sizeof(int));
  neicpu =(int*)calloc(ngridmax,sizeof(int));

  if(cpu->bndoct!=NULL) free(cpu->bndoct);
  cpu->bndoct=(struct OCT**)calloc(cpu->nbuff,sizeof(struct OCT*));

  // looking for neighbors

  for(level=1;level<=levelmax;level++)
    {
      nextoct=firstoct[level-1];
      if(nextoct!=NULL){
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    if(level>=levelcoarse){

	      if(level==levelcoarse){
		assigncpu2coarseoct(curoct, cpu, levelcoarse);
	      }

	      key=oct2key(curoct,level); // getting the key of the current oct

	      if(curoct->cpu!=cpu->rank){
		// a neighbor has been found
		mpinei[nnei]=curoct->cpu;
		cpu->bndoct[nnei]=curoct; // contains the oct adresses for emission
		nnei++;
	      }

	      //filling the hash table for levelcoarse octs
	      hidx=hfun(key); // getting the hash index from the key
	      
	      if(cpu->htable[hidx]==NULL){ //this bucket is empty
		cpu->htable[hidx]=curoct;
		curoct->nexthash=NULL;
	      }
	      else{ // we solve for collisions by looking for the last oct the chained list
		newoct=cpu->htable[hidx];
		do{
		  desoct=newoct;
		  newoct=desoct->nexthash;
		}while(newoct!=NULL);
		desoct->nexthash=curoct;
		curoct->nexthash=NULL;
	      }
	    }
	    else{
	      curoct->cpu=-1;
	    }
	  }while(nextoct!=NULL);
      }
    }


  // =========================================== SETTING UP THE COMMUNICATIONS BETWEEN NEIGHBORS


  // computing the mpi neighbor list
  myradixsort(mpinei,nnei); // we sort the neighbors list
  neicpu[0]=mpinei[0];
  j=0;
  for(i=1;i<nnei;i++){ // we scan the list
    if(mpinei[i]!=neicpu[j]){
      j++;
      neicpu[j]=mpinei[i];
    }
  }
  nbnd=nnei;
  nnei=j+1;


  free(mpinei);

  cpu->nebnd=nbnd;
  cpu->nnei=nnei;
  cpu->mpinei=(int*)calloc(nnei,sizeof(int)); // we reallocate the array to free some memory
  for(i=0;i<cpu->nnei;i++) cpu->mpinei[i]=neicpu[i];
  free(neicpu);



  // AT THIS STAGE: 
  // nbnd contains the number of boundary octs
  // nnei contains the number of neighbor procs
  // mpinei contains the rank of the neighbor procs and has a size nnei
  // inidx contains the number of innerboundary octs


  // some displays
  printf("Found %d neighbors and %d bnd octs  on rank %d :",cpu->nnei,cpu->nebnd,cpu->rank);
  for(i=0;i<cpu->nnei;i++) printf("%d ",cpu->mpinei[i]);
  printf("\n");
  
  
  // creating a cpu dictionnary to translate from cpu number to inei
  
  cpu->dict=(int*)calloc(cpu->nproc,sizeof(int));
  for(i=0;i<cpu->nproc;i++) cpu->dict[i]=-1;
  for(i=0;i<cpu->nnei;i++){
    cpu->dict[cpu->mpinei[i]]=i;
  }
}

//======================================================================================
void gather_ex(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field){
  
  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
   */
  int *countpacket;
  int icpu;
  struct PACKET *pack;
  int i,ii;

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer
    
    // assigning the values
    pack->level=cpu->bndoct[i]->level;
    pack->key=(int)oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    for(ii=0;ii<8;ii++){
      switch(field){
      case 0:
	pack->data[ii]=cpu->bndoct[i]->cell[ii].density;
	break;
      case 1:
	pack->data[ii]=cpu->bndoct[i]->cell[ii].marked;
	break;
      }
    }
    
    // counting the number of packets for icpu
    countpacket[icpu]++;
  }
  
  //for(i=0;i<cpu->nnei;i++) printf("rank=%d cpu %d nbnd=%d\n",cpu->rank,cpu->mpinei[i],countpacket[i]);
  free(countpacket);

}


//======================================================================================
int gather_ex_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer,struct PART **lastp){
  
  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
   */
  int *countpacket;
  int icpu;
  struct PART_MPI *part;
  struct PART *curp;
  struct PART *nexp;
  int i,ii;
  int nrem=0;

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    for(ii=0;ii<8;ii++){
      nexp=cpu->bndoct[i]->cell[ii].phead;
      // sweeping the particles of the boundary cells
      if(nexp!=NULL){
	do{ 
	  curp=nexp; 
	  nexp=curp->next; 
	  
	  part=psendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer

	  // assigning the values
	  part->level=cpu->bndoct[i]->level;
	  part->key=(int)oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct (is the eventual destination)
	  part->icell=ii;

	  // we assign the data to the communicator
	  part->x=curp->x;
	  part->y=curp->y;
	  part->z=curp->z;

	  part->vx=curp->vx;
	  part->vy=curp->vy;
	  part->vz=curp->vz;

	  part->mass=curp->mass;
	  part->idx=curp->idx;
	  
	  // counting the number of packets for icpu
	  countpacket[icpu]++;

	  // switching the mass to -1 to flag exiting particles
	  curp->mass=-1.;
	  nrem++;

	  // disconnecting the particle from its list

	  curp->prev=NULL;
	  curp->next=NULL;
	  
	  if(nexp==NULL){ // reached the last particle from the cell
	    cpu->bndoct[i]->cell[ii].phead=NULL; // we "clear" the cell
	  }

	  // is it the global last particle ?
	  if(curp==(*lastp)){
	    (*lastp)=curp-1;
	    while((*lastp)->mass<0){ // if true this particle does not exist anymore
	      (*lastp)=(*lastp)-1;
	    }
	  }

	}while(nexp!=NULL);
      }
    }
  }
  free(countpacket);

  return nrem; // we return the number of exiting particles
}



//======================================================================================

void gather_mpi(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field){

  int i,j;
  int found=0;
  int hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  int *countpacket;

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

    for(j=0;j<cpu->nnei;j++){
      for(i=0;i<cpu->nbuff;i++){
	pack=sendbuffer[j]+i; // we assume that the sendbuffer already contains the keys
	if(pack->level!=0){ // we do something
	  countpacket[j]++;

	  // first we compute the adress from the hashfunction
	  hidx=hfun(pack->key);
	  nextoct=cpu->htable[hidx];
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    if(found){ // the reception oct has been found
	      for(icell=0;icell<8;icell++){
		switch(field){
		case 0:
		  pack->data[icell]=curoct->cell[icell].density; // density
		  break;
		case 1:
		  pack->data[icell]=curoct->cell[icell].density; // density again we reproduce the case 1 in order to be consistent with scatter_mpi
		  break;
		case 2:
		  pack->data[icell]=curoct->cell[icell].pot; // potential
		  break;
		case 3:
		  pack->data[icell]=curoct->cell[icell].marked; //refinment mark
		  break;
		case 4:
		  pack->data[icell]=curoct->cell[icell].temp; //temp field for force calculation
		  break;
		}
	      }
	    }
	    else{
	      printf("error no reception oct found !");
	      abort();
	    }
	    
	  }else{
	    printf("error no hash key obtained !!\n");
	    abort();
	  }
	}
      }
    }
    
    free(countpacket);
      
}

 //------------------------------------------------------------------------
void scatter_mpi(struct CPUINFO *cpu, struct PACKET **recvbuffer,  int field){

  int i,j;
  int found=0;
  int hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something
	  

	// first we compute the adress from the hashfunction
	hidx=hfun(pack->key);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      switch(field){
	      case 0:
		curoct->cell[icell].density+=pack->data[icell]; // density += for CIC correction
		break;
	      case 1:
		curoct->cell[icell].density =pack->data[icell]; // density
		break;
	      case 2:
		curoct->cell[icell].pot=pack->data[icell]; // potential
		break;
	      case 3:
		curoct->cell[icell].marked=fmax(pack->data[icell],(float)curoct->cell[icell].marked); // refinement mark
		break;
	      case 4:
		curoct->cell[icell].temp=pack->data[icell]; // temp field for force calculation
		break;
	      }
	    }
	  }
	  else{
	    printf("error no reception oct found !");
	    abort();
	  }
	    
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
      
}


 //------------------------------------------------------------------------
int scatter_mpi_part(struct CPUINFO *cpu, struct PART_MPI **precvbuffer, struct PART **lastp){

  int i,j;
  int found=0;
  int hidx;
  struct PART_MPI *part;
  struct PART *curp;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  int nadd=0;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      part=precvbuffer[j]+i;
      if(part->level!=0){ // we do something
	// first we compute the adress from the hashfunction
	hidx=hfun(part->key);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==part->key)&&(part->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    // we assign the new particle to the global particle list of the client
	    *lastp=*lastp+1; // getting the next slot in the global particle list
	    
	    if((*lastp)->mass>0) {
	      printf("oum\n");
	      abort();
	    }
	    // copying the data
	    (*lastp)->x =part->x;
	    (*lastp)->y =part->y;
	    (*lastp)->z =part->z;

	    (*lastp)->vx=part->vx;
	    (*lastp)->vy=part->vy;
	    (*lastp)->vz=part->vz;
	      
	    (*lastp)->mass=part->mass;
	    (*lastp)->idx=part->idx;
	    nadd++;
  
	    // looking for the tail particle in the destination cell
	    curp=findlastpart(curoct->cell[part->icell].phead);
	    if(curp!=NULL){
	      curp->next=(*lastp);
	      (*lastp)->next=NULL;
	      (*lastp)->prev=curp;
	    }
	    else{
	      curoct->cell[part->icell].phead=(*lastp);
	      (*lastp)->next=NULL;
	      (*lastp)->prev=NULL;
	    }
	  }
	  else{
	    printf("error no reception oct found !");
	    abort();
	  }
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
  return nadd; // returning the new number of particles
}

//------------------------------------------------------------------------
void compute_bndkeys(struct CPUINFO *cpu, struct PACKET **recvbuffer){

  // we create a counter of values for each neighbor
  int *countpacket;
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct PACKET *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    int keyloc;
    int cpuloc;
    int inei;
    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary
    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=keyloc;
    pack->level=cpu->bndoct[i]->level;
    countpacket[inei]++;
  }  

  free(countpacket);
}

//------------------------------------------------------------------------

void  clean_mpibuff(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer){
  int i;
  for(i=0;i<cpu->nnei;i++) {
    memset(sendbuffer[i],0,cpu->nbuff*sizeof(struct PACKET));
    memset(recvbuffer[i],0,cpu->nbuff*sizeof(struct PACKET));
  }
}


void  clean_mpibuff_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer){
  int i;
  for(i=0;i<cpu->nnei;i++) {
    memset(psendbuffer[i],0,cpu->nbuff*sizeof(struct PART_MPI));
    memset(precvbuffer[i],0,cpu->nbuff*sizeof(struct PART_MPI));
  }
}
 //------------------------------------------------------------------------

#ifdef WMPI
void mpi_exchange(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field)
{
  int i;
  MPI_Status *stat;
  MPI_Request *req;
  int mpitag=1;
  

  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff(cpu,sendbuffer,recvbuffer);
  // ----------- I  / we compute the boundary keys and store them in recvbuffer
  compute_bndkeys(cpu,recvbuffer);
  // ----------- II / we send the keys to the server
  /* for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors */
  /*   MPI_Sendrecv(recvbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,sendbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,&stat); */
  /* } */

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors
    mpitag=cpu->rank+cpu->mpinei[i];
    //MPI_Sendrecv(sendbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,recvbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,&stat);
    MPI_Isend(recvbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,req+2*i  );
    MPI_Irecv(sendbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,req+2*i+1);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  // ----------- III/ the server gather the data
  gather_mpi(cpu, sendbuffer, field);
  MPI_Barrier(cpu->comm);

  // ----------- IV / the server send the data back to the client
  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors
    mpitag=cpu->rank+cpu->mpinei[i];
    //MPI_Sendrecv(sendbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,recvbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,&stat);
    MPI_Isend(sendbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,req+2*i  );
    MPI_Irecv(recvbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,req+2*i+1);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);

  // ----------- V  / the client scatter the data back in the oct tree
  scatter_mpi(cpu,recvbuffer,field);
  MPI_Barrier(cpu->comm);

  //
  free(stat);
  free(req);

}

//------------------------------------------------------------------------
void mpi_cic_correct(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field)
{
  int i;
  //MPI_Status stat;
  int mpitag=1;
  MPI_Status *stat;
  MPI_Request *req;


  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff(cpu,sendbuffer,recvbuffer);

  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  gather_ex(cpu,sendbuffer,field);
  MPI_Barrier(cpu->comm);

  // ---------  second we transmit the data through the network
  /* for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors */
  /*   mpitag=cpu->rank+cpu->mpinei[i]; */
  /*   MPI_Sendrecv(sendbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,recvbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,&stat); */
  /* } */

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors
    mpitag=cpu->rank+cpu->mpinei[i];
    //MPI_Sendrecv(sendbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,recvbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,&stat);
    MPI_Isend(sendbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,req+2*i  );
    MPI_Irecv(recvbuffer[i],cpu->nbuff,*cpu->MPI_PACKET,cpu->mpinei[i],mpitag,cpu->comm,req+2*i+1);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);


  // ---------  third we scatter the data back to the INTERNAL boundary octs
  if(field==1) field=3; // fix for an offset of marked scatter
  scatter_mpi(cpu,recvbuffer,field);
  MPI_Barrier(cpu->comm);

  //
  free(stat);
  free(req);


}
#endif


 //------------------------------------------------------------------------
void mark_cells(int levelcoarse,int levelmax,struct OCT **firstoct, int nsmooth, float threshold, struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer){

  int nmark;
  int level;
  int marker;
  int ismooth;
  float dx;
  int pass;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *newoct;
  struct OCT *desoct;
  int icell,ii,il,ip;
  int smark;
  int vnei[6],vcell[6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  struct CELL *newcell;
  struct CELL *newcell2;
  struct CELL *newcell3;
  int ichild;

    printf("==> start marking\n");
    //    for(level=levelmax;level>=levelcoarse;level--) // looping over octs
    for(level=levelmax;level>=1;level--) // looping over octs
      {
	marker=0;
	nmark=0;
	for(ismooth=0;ismooth<nsmooth;ismooth++)
	  {
	    //printf("level=%d ",level);
	    dx=1./pow(2,level);
	    
	    for(pass=0;pass<3;pass++)
	      {
		marker++;
		nextoct=firstoct[level-1];
		if(nextoct==NULL) continue;
		
		do
		  {
		    curoct=nextoct;
		    nextoct=curoct->next;
		    
		    for(icell=0;icell<8;icell++) // looping over cells in oct
		      {
			switch(pass){
			  //=========================================================
			case 0: // marking cell already refined or marked marker=1/4
			  smark=0;
			  if(curoct->cell[icell].child!=NULL){ // CHECK IF ACTING ON ORIGINAL OR COPY
			    newoct=curoct->cell[icell].child;
			    for(ichild=0;ichild<8;ichild++){
			      smark+=((newoct->cell[ichild].marked!=0)||(newoct->cell[ichild].child!=NULL));
			    }
			  }
			  if(smark!=0){
			    curoct->cell[icell].marked=marker;
			    nmark++;
			  }
			  break;
			  //=========================================================
			case 1: //marking neighbors marker=2/5
			  if((curoct->cell[icell].marked<marker)&&(curoct->cell[icell].marked>0)){
			    getcellnei(icell, vnei, vcell);
			    for(ii=0;ii<6;ii++){ // marking the 6 cardinal neighbors
			      if(vnei[ii]==6){
				newcell=&(curoct->cell[vcell[ii]]);
				if(curoct->cell[vcell[ii]].marked==0){
				  curoct->cell[vcell[ii]].marked=marker;
				  nmark++;
				}
			      }
			      else{
				// Note that the neibourgh cell may not exist therefore we have to check
				if(curoct->nei[vnei[ii]]->child!=NULL){
				  newcell=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
				  if(curoct->nei[vnei[ii]]->child->cell[vcell[ii]].marked==0) {
				    curoct->nei[vnei[ii]]->child->cell[vcell[ii]].marked=marker;
				    nmark++;
				  }
				}
				else{
				  newcell=NULL;
				}
			      }

			      // each of the 6 cardinal neighbors will mark 4 side neighbors
			      if(newcell!=NULL){
				newoct=cell2oct(newcell); // we get the parent oct
				getcellnei(newcell->idx, vnei2, vcell2); //we get the neighbors
				for(il=0;il<6;il++){
				  if((il/2)==(ii/2)) continue;
				  if(vnei2[il]==6){
				    newcell2=&(newoct->cell[vcell2[il]]);
				    if(newoct->cell[vcell2[il]].marked==0){
				      newoct->cell[vcell2[il]].marked=marker;
				      nmark++;
				    }
				  }
				  else{
				    if(newoct->nei[vnei2[il]]->child!=NULL){
				      newcell2=&(newoct->nei[vnei2[il]]->child->cell[vcell2[il]]);
				      if(newoct->nei[vnei2[il]]->child->cell[vcell2[il]].marked==0){
					newoct->nei[vnei2[il]]->child->cell[vcell2[il]].marked=marker;
					nmark++;
				      }
				    }
				    else{
				      newcell2=NULL;
				    }
				  }
				  
				  // ecah of the 4 side neighbors will mark 2 corners
				  if(newcell2!=NULL){
				    desoct=cell2oct(newcell2);
				    getcellnei(newcell2->idx, vnei3, vcell3);
				    for(ip=0;ip<6;ip++){
				      if(((ip/2)==(il/2))||((ip/2)==(ii/2))) continue;
				      if(vnei3[ip]==6){
					if(desoct->cell[vcell3[ip]].marked==0){
					  desoct->cell[vcell3[ip]].marked=marker;
					  nmark++;
					}
				      }
				      else{
				  	if(desoct->nei[vnei3[ip]]->child!=NULL){
					  if(desoct->nei[vnei3[ip]]->child->cell[vcell3[ip]].marked==0){
					    desoct->nei[vnei3[ip]]->child->cell[vcell3[ip]].marked=marker;
					    nmark++;
					  }
				  	}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			  break;
			//=========================================================
			case 2: // marking cells satisfying user defined criterion marker=3/6
			  //if(countpart(curoct->cell[icell].phead)>threshold) curoct->cell[icell].marked=marker;
			  //if(curoct->npart>threshold) curoct->cell[icell].marked=marker;
			  if((curoct->cell[icell].density>threshold)&&(curoct->cell[icell].marked==0)) {
			    curoct->cell[icell].marked=marker;
			    nmark++;
			  }
			  break;
			}
		      }
	      
		  }while(nextoct!=NULL);
		//printf("pass=%d nmark=%d\n",pass,nmark);
#ifdef WMPI
		// first we correct from the marker diffusion
		mpi_cic_correct(cpu,sendbuffer,recvbuffer,1);
		
		// second exchange boundaries
		mpi_exchange(cpu,sendbuffer,recvbuffer,3);
#endif
	      }
	    //printf("\n");
	  }
      }

}

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

  for(dir=0;dir<3;dir++){
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
	    
	    // First we gather the potential in all neighbors
	    for(icomp=2*dir;icomp<=2*dir+1;icomp++){
	      memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
	      nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,cpu,&nread);
	    }
	  
	    // Next we perform the finite difference along x
	    grad(vcomp,stride,dx,dir);
	  
	    // Then we scatter back the results in the temp variable
	    nextoct=scattercomp(curoct, vcomp[6], 6, 2, stride,cpu);
	  
	  }while(nextoct!=NULL);
	
	}
      }

#ifdef WMPI
	  mpi_exchange(cpu,sendbuffer,recvbuffer,4);
#endif
  
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
	  
	    //==== FIRST WE CONSIDER THE PARTICLES INSIDE THE BOUNDARIES AT LEVEL L
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
  }
}


//------------------------------------------------------------------------
int mpi_exchange_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer, struct PART **lastpart){

  int nrem,nadd;
  int mpitag=1;
  int i;

  clean_mpibuff_part(cpu,psendbuffer,precvbuffer);
  nrem=gather_ex_part(cpu,psendbuffer,lastpart);
  MPI_Barrier(cpu->comm);

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors
    mpitag=cpu->rank+cpu->mpinei[i];
    MPI_Sendrecv(psendbuffer[i],cpu->nbuff,*cpu->MPI_PART,cpu->mpinei[i],mpitag,precvbuffer[i],cpu->nbuff,*cpu->MPI_PART,cpu->mpinei[i],mpitag,cpu->comm,MPI_STATUS_IGNORE);
  }

  nadd=scatter_mpi_part(cpu,precvbuffer,lastpart);
  MPI_Barrier(cpu->comm);

  // Return delta part

  return nadd-nrem;
}

 //------------------------------------------------------------------------
 // the MAIN CODE
 //------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  struct OCT *grid;
  struct OCT **firstoct;
  struct OCT **lastoct;

  int level,levelcoarse,levelmax;
  int ngridmax,ngrid;
  int npartmax;
  int cur,curnext; // flat indexes with boundaries
  int i,il,ichild,icell,inext,ii,ip,j;
  int xp,yp,zp;
  int NBND=1,NBND2=2*NBND;
  float dx;
  int vnei[6],vcell[6]; // arrays to get neighbors
  int vnei2[6],vcell2[6]; // arrays to get neighbors
  int vnei3[6],vcell3[6]; // arrays to get neighbors
  int neip[7]; // contains the index of the six neighbors of the current parent cell +current parent
  int ci,cj,ck;
  int cinext,cjnext,cknext;
  float threshold;
  struct OCT oct;
  struct OCT* nextoct;
  struct OCT* curoct;
  struct OCT* desoct;
  struct CELL * parcell;
  struct CELL * newcell;
  struct CELL * newcell2;
  int tag;
  float dxcur;
  float *dens;
  int firstoct_currl;
  int nxoct;
  int lmap;
  int npart;
  struct PART *part;
  struct PART *nexploc, *curploc;

  struct OCT *endoct; // the very last oct of all the levels;
  struct OCT *newoct;
  int nref=0,ndes=0;

  float xc,yc,zc;
  int stride;
  float **vcomp;
  int ncomp;
  float acc;
  float dt;
  int ntot=0,nlev,noct;
  float ntotd=0.,nlevd=0.;

  float disp,mdisp;
  
  int dir;

  char filename[128]; 
  FILE *fd;
  struct PART *nexp;
  struct PART *nexp2;
  struct PART *curp;

  struct PART *lastpart;

  int curc;

  int nbnd;

  float x,y,z;
  float vx,vy,vz;
  float mass;
  float idx;
  
  unsigned key;

  struct CPUINFO cpu;

  struct PACKET **sendbuffer; 
  struct PACKET **recvbuffer; 

  struct PART_MPI **psendbuffer; 
  struct PART_MPI **precvbuffer; 

  //=========== some initial calls =============
#ifdef WMPI
  MPI_Status stat;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(cpu.nproc));
  MPI_Comm_rank(MPI_COMM_WORLD,&(cpu.rank));

  //========= creating a PACKET MPI type =======
  MPI_Datatype MPI_PACKET,oldtypes[2]; 
  int          blockcounts[2];
  
  /* MPI_Aint type used to be consistent with syntax of */
  /* MPI_Type_extent routine */
  MPI_Aint    offsets[2], extent;
  
  
  /* Setup description of the 8 MPI_FLOAT fields data */
  offsets[0] = 0;
  oldtypes[0] = MPI_FLOAT;
  blockcounts[0] = 8;
  
  /* Setup description of the 2 MPI_INT fields key, level */
  /* Need to first figure offset by getting size of MPI_FLOAT */
  MPI_Type_extent(MPI_FLOAT, &extent);
  offsets[1] = 8 * extent;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 2;

  /* Now define structured type and commit it */
  MPI_Type_struct(2, blockcounts, offsets, oldtypes, &MPI_PACKET);
  MPI_Type_commit(&MPI_PACKET);

  //========= creating a PART MPI type =======
  MPI_Datatype MPI_PART;

  /* Setup description of the 7 MPI_FLOAT fields x,y,z,vx,vy,vz */
  offsets[0] = 0;
  oldtypes[0] = MPI_FLOAT;
  blockcounts[0] = 7;
  
  /* Setup description of the 4 MPI_INT fields idx key level icell*/
  /* Need to first figure offset by getting size of MPI_FLOAT */
  MPI_Type_extent(MPI_FLOAT, &extent);
  offsets[1] = 7 * extent;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 4;

  /* Now define structured type and commit it */
  MPI_Type_struct(2, blockcounts, offsets, oldtypes, &MPI_PART);
  MPI_Type_commit(&MPI_PART);
  
  //============================================


  cpu.MPI_PACKET=&MPI_PACKET;
  cpu.MPI_PART=&MPI_PART;
  cpu.comm=MPI_COMM_WORLD;
#else
  cpu.rank=0;
  cpu.nproc=1;
#endif

  //=========== assigning values =============
  levelcoarse=LCOARSE;
  levelmax=LMAX;

  ngridmax=1000000;
  npartmax=64*64*64*2;
#ifdef PART2
  npart=2;
#else
  npart=64*64*64;
#endif

  threshold=50;
  lmap=LMAX;
  stride=fmax(8,STRIDE);//pow(2,levelcoarse);
  ncomp=8;
  acc=1e-2;
  dt=DT;

  //breakmpi();
  //========== allocations ===================

  if(cpu.rank==0) printf("Allocating %f GB cell=%f GB part=%f GB",(sizeof(struct OCT)*ngridmax+sizeof(struct PART)*npart)/(1024*1024*1024.),sizeof(struct OCT)*ngridmax/(1024*1024*1024.),sizeof(struct PART)*npart/(1024*1024*1024.));

  grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT));
  firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
  lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
  part=(struct PART*)calloc(npartmax,sizeof(struct PART));
  cpu.htable=(struct OCT**) calloc(pow(2,3*levelmax-3)/64,sizeof(struct OCT *)); //hashtable assunming a hash function >>6


  vcomp=(float **)calloc(ncomp,sizeof(float*));
  for(i=0;i<ncomp;i++)
    {
      vcomp[i]=(float *)calloc(stride,sizeof(float));
    }
    if(cpu.rank==0) printf("Allocations ok\n");


  //========== setting up the parallel topology ===

  // We segment the oct distributions at levelcoarse 
    cpu.bndoct=NULL;
    cpu.nbuff=NBUFF;
    cpu.allkmin=(int*)calloc(cpu.nproc,sizeof(int));
    cpu.allkmax=(int*)calloc(cpu.nproc,sizeof(int));

    load_balance(levelcoarse,&cpu);

#ifdef WMPI
    MPI_Allgather(&cpu.kmin,1,MPI_INT,cpu.allkmin,1,MPI_INT,cpu.comm);
    MPI_Allgather(&cpu.kmax,1,MPI_INT,cpu.allkmax,1,MPI_INT,cpu.comm);
    MPI_Barrier(cpu.comm);
#else
    cpu.allkmin[0]=cpu.kmin;
    cpu.allkmax[0]=cpu.kmax;
#endif    
    
    

    
  //========== building the initial meshes ===

  if(cpu.rank==0) printf("building initial mesh\n");

  //breakmpi();
  // ZERO WE CREATE A ROOT CELL
  
  struct CELL root;
  root.child=grid;
  

  // FIRST WE POPULATE THE ROOT OCT
  grid->x=0.;
  grid->y=0.;
  grid->z=0.;

  grid->parent=NULL;
  grid->level=1;
  for(i=0;i<6;i++) grid->nei[i]=&root;
  grid->prev=NULL;
  grid->next=NULL;

  // setting the densities in the cells and the index
  for(icell=0;icell<8;icell++){ 
    grid->cell[icell].density=0.;
    grid->cell[icell].pot=0.;
    grid->cell[icell].temp=0.;
    grid->cell[icell].idx=icell;
  }

  grid->cpu=-1;

  // start the creation of the initial amr grid from level 1
  firstoct[0]=grid;
  lastoct[0]=grid;
  int noct2;
  int segok;

  newoct=grid+1;
  for(level=1;level<levelcoarse;level++){ // sweeping the levels from l=1 to l=levelcoarse
    dxcur=1./pow(2,level);
    nextoct=firstoct[level-1];
    noct2=0;
    if(nextoct==NULL) continue;
    do // sweeping level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++){ // sweeping the cells

	  segok=segment_cell(curoct,icell,&cpu,levelcoarse);// the current cell will be splitted according to a segmentation condition
	  if(segok==1){ 
	    //if(level==levelcoarse-1) printf(" segok=%d\n",segok);

	    noct2++;
	    // the newoct is connected to its mother cell
	    curoct->cell[icell].child=newoct;
	    
	    // a newoct is created
	    newoct->parent=&(curoct->cell[icell]);
	    newoct->level=curoct->level+1;
	    newoct->x=curoct->x+( icell   %2)*dxcur;
	    newoct->y=curoct->y+((icell/2)%2)*dxcur;
	    newoct->z=curoct->z+( icell   /4)*dxcur;

	    // filling the cells
	    for(ii=0;ii<8;ii++){
	      newoct->cell[ii].marked=0;
	      newoct->cell[ii].child=NULL;
	      newoct->cell[ii].density=0.;
	      newoct->cell[ii].idx=ii;
	      newoct->cell[ii].phead=NULL;
	    }
	    
	    //the neighbours
	    getcellnei(icell, vnei, vcell);
	    for(ii=0;ii<6;ii++){
	      if((vnei[ii]!=6)){ 
		newoct->nei[ii]=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
	      }else{
		newoct->nei[ii]=&(curoct->cell[vcell[ii]]);
	      }
	    }
	    
	    // preparing the next creations on level+1
	    newoct->next=NULL;
	    
	    if(firstoct[level]==NULL){
	      firstoct[level]=newoct;
	      newoct->prev=NULL;
	    }
	    else{
	      newoct->prev=lastoct[level];
	      lastoct[level]->next=newoct;
	    }
	    lastoct[level]=newoct;


	    // next oct ready
	    newoct++; 
	  }
 	}
      }while(nextoct!=NULL);
    if(cpu.rank==0) printf("level=%d noct=%d\n",level,noct2);
  }


 // ==================================== assigning CPU number to levelcoarse OCTS // filling the hash table // Setting up the MPI COMMS

  setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax);

  // allocating the communication buffers

  sendbuffer=(struct PACKET **)(calloc(cpu.nnei,sizeof(struct PACKET*)));
  recvbuffer=(struct PACKET **)(calloc(cpu.nnei,sizeof(struct PACKET*)));
  for(i=0;i<cpu.nnei;i++) {
    sendbuffer[i]=(struct PACKET *) (calloc(cpu.nbuff,sizeof(struct PACKET)));
    recvbuffer[i]=(struct PACKET *) (calloc(cpu.nbuff,sizeof(struct PACKET)));
  }

  psendbuffer=(struct PART_MPI **)(calloc(cpu.nnei,sizeof(struct PART_MPI*)));
  precvbuffer=(struct PART_MPI **)(calloc(cpu.nnei,sizeof(struct PART_MPI*)));
  for(i=0;i<cpu.nnei;i++) {
    psendbuffer[i]=(struct PART_MPI *) (calloc(cpu.nbuff,sizeof(struct PART_MPI)));
    precvbuffer[i]=(struct PART_MPI *) (calloc(cpu.nbuff,sizeof(struct PART_MPI)));
  }


  //===================================================================================================
  
  // ==== some initial dump

  sprintf(filename,"data/levstart.%05d.p%05d",0,cpu.rank);
  dumpcube(lmap,firstoct,0,filename);
  sprintf(filename,"data/cpustart.%05d.p%05d",0,cpu.rank);
  dumpcube(lmap,firstoct,3,filename);

  // =====  computing the memory location of the last oct 

  endoct=lastoct[0];
  for(i=0;i<levelcoarse;i++) {
    if(lastoct[i]>endoct) endoct=lastoct[i];
  }


#if 1  // ==================================== assigning particles to cells

  if(cpu.rank==0) printf("==> starting part\n");
  firstoct_currl=0;
  for(il=1;il<levelcoarse;il++) firstoct_currl+=pow(pow(2,il-1),3); // the index of the first oct of current level
 
  // initialisation of particles
  

#ifdef PART2

  int ir,nr=2;
  ip=0;
  for(ir=0;ir<nr;ir++) {
    // first we read the position etc... (eventually from the file)
    if(ir==0){
      x=0.5;
      y=0.5;
      z=0.5;
      
      vx=0.;
      vy=0.;
      vz=0.;
      
      mass=0.999;
    }
    else if(ir==1){
      x=0.5+0.2;
      y=0.5;
      z=0.5;
      
      vx=0.;
      vy=sqrt(0.999/0.2);
      vz=0.;
      
      mass=0.001;
    }
    
    // periodic boundary conditions
    
    x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
    y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
    z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 
    
    // it it belongs to the current cpu, we proceed and assign the particle to the particle array
    if(segment_part(x,y,z,&cpu,levelcoarse)){
      part[ip].x=x;
      part[ip].y=y;
      part[ip].z=z;
      
      part[ip].vx=vx;
      part[ip].vy=vy;
      part[ip].vz=vz;
      
      part[ip].mass=mass;
      lastpart=part+ip;
      ip++;
    }
  }
  
  npart=ip; // we compute the localnumber of particle

#endif

#ifdef TESTPLUM
  int dummy;
  float dummyf;
  int npartf;

  //breakmpi();
  fd=fopen("utils/data.inp","r");
  if(fd==NULL) {
    printf("Error while reading particle file ABORT\n");
    abort();
  }
  fscanf(fd,"%d",&dummy);
  fscanf(fd,"%d",&npartf);
  fscanf(fd,"%f",&dummyf);

  ip=0.;
  for(i=0;i<npartf;i++)
    {
      //fscanf(fd,"%d %f %f %f %f %f %f %f",&part[i].idx,&part[i].mass,&(part[i].x),&(part[i].y),&(part[i].z),&(part[i].vx),&(part[i].vy),&(part[i].vz));
      fscanf(fd,"%d %f %f %f %f %f %f %f",&dummy,&mass,&x,&y,&z,&vx,&vy,&vz);
      
      x+=0.5;
      y+=0.5;
      z+=0.5;
      // periodic boundary conditions
    
      x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
      y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
      z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 

      // it it belongs to the current cpu, we proceed and assign the particle to the particle array
      if(segment_part(x,y,z,&cpu,levelcoarse)){
	part[ip].x=x;
	part[ip].y=y;
	part[ip].z=z;
	
	part[ip].vx=vx;
	part[ip].vy=vy;
	part[ip].vz=vz;
	
	part[ip].mass=mass;
	lastpart=part+ip;
	ip++;
      }
      
    }
  fclose(fd);
  npart=ip; // we compute the localnumber of particle

#endif  

  // assigning particles to cells in coarse octs (assuming octs are aligned)

  if(cpu.rank==0) printf("start populating coarse grid with particles\n");
  struct PART* lastp[8]; // will contain the last particle of the 8 cells in each oct

  // FIRST WE CONSIDER THE LEVEL 1
  for(ii=0;ii<8;ii++) lastp[ii]=NULL; // we initialise the last part of each sub cell
  dxcur=0.5;
  for(i=0;i<npart;i++)
    {
      curc=(int)((part[i].x-grid[0].x)/dxcur)+(int)((part[i].y-grid[0].y)/dxcur)*2+(int)((part[i].z-grid[0].z)/dxcur)*4;
      
      if(grid[0].cell[curc].phead==NULL){
	grid[0].cell[curc].phead=&part[i];
	lastp[curc]=&part[i];
      }
      else{
	lastp[curc]->next=&part[i];
	part[i].prev=lastp[curc];
	lastp[curc]=&part[i];
      }
    }
  if(cpu.rank==0) printf("Part assigned root level ok\n");
  
  // WE POPULATE THE NEXT LEVELS BY SUBDIVISIONS
  int np=0;
  for(level=1;level<=levelcoarse-1;level++) // we stop at level coarse -1 because it will be assigned from levelcoarse-1
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  dxcur=1./pow(2,level+1); // size of a CELL at level +1
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if(curoct->cell[icell].child!=NULL){ // a child has been detected so we split the particle in 8 cells
		for(ii=0;ii<8;ii++) lastp[ii]=NULL; // we initialise the last part of each sub cell
		newoct=curoct->cell[icell].child; // temp oct for practical reasons
		nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
		if(nexp!=NULL){ 
		  do{  
		    curp=nexp; 
		    nexp=curp->next; 
		    
		    //curc is the index of the new cell at level+1
		    curc=(int)((curp->x-newoct->x)/dxcur)+(int)((curp->y-newoct->y)/dxcur)*2+(int)((curp->z-newoct->z)/dxcur)*4;
		    if(lastp[curc]==NULL){
		      // first particle in the current subcell
		      newoct->cell[curc].phead=curp;
		      curp->prev=NULL;
		      curp->next=NULL;
		      lastp[curc]=curp;
		    }
		    else{
		      // the current particle is linked to the last one in the current part
		      curp->prev=lastp[curc];
		      lastp[curc]->next=curp;
		      curp->next=NULL;
		      lastp[curc]=curp;
		    }
		  }while(nexp!=NULL); 
		  
		  // we empty the mother cell from particles
		  curoct->cell[icell].phead=NULL;
		  
		}
	      }
	    }
	}while(nextoct!=NULL);
    }


#endif
#if 1
  // ==================================== Check the number of particles and octs

  multicheck(firstoct,npart,levelmax,cpu.rank);

  sprintf(filename,"data/parstart.%05d.p%05d",0,cpu.rank);
  dumppart(firstoct,filename,npart,levelcoarse,levelmax);

#endif	

#if 1
  // ==================================== performing the CIC assignement

  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#endif

#ifdef WMPI
    // ==================================== performing the CIC BOUNDARY CORRECTION 

  mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);

  // ======================================= Density boundary mpi update 

  mpi_exchange(&cpu,sendbuffer,recvbuffer,1);

#endif



  sprintf(filename,"data/denstart.%05d.p%05d",0,cpu.rank);
  printf("%s\n",filename);
  dumpcube(lmap,firstoct,1,filename);





  //================================================================================
  //================================================================================
  //================================================================================
  //
  //          AT THIS STAGE THE INITIAL SETUP HAS BEEN COMPLETED
  //
  //================================================================================
  //================================================================================
  //================================================================================


  int nsteps;
  int pass;
  int smark;
  int ismooth,nsmooth=2;
  int marker;

  for(nsteps=0;nsteps<=NSTEP;nsteps++){
    
    if(cpu.rank==0) printf("============== STEP %d ================\n",nsteps);
    //printf("endoct=%p\n",endoct);
#if 1
    // ==================================== marking the cells
    //if(nsteps==1)  breakmpi();
    
    sprintf(filename,"data/dmstart.%05d.p%05d",nsteps+1,cpu.rank);
    printf("%s\n",filename);
    dumpcube(lmap,firstoct,1,filename);

    mark_cells(levelcoarse,levelmax,firstoct,nsmooth,threshold,&cpu,sendbuffer,recvbuffer);

    sprintf(filename,"data/markstart.%05d.p%05d",nsteps+1,cpu.rank);
    printf("%s\n",filename);
    dumpcube(lmap,firstoct,4,filename);
    
    
  
    
    // ==================================== refining (and destroying) the octs

    
    refine_cells(levelcoarse,levelmax,firstoct,lastoct,endoct,&cpu);
    //if(nsteps==1) breakmpi();
    // cleaning the marks
    for(level=1;level<=levelmax;level++) // looping over levels
      {
	float maxd=0.,mind=1e30,avg=0.;
	int ncell=0;
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
/* 	      if(curoct->cell[icell].marked>50){ */
/* 		printf("ouhl %f\n",curoct->cell[icell].marked); */
/* 		//abort(); */
/* 	      } */
	      curoct->cell[icell].marked=0.;
	      ncell++;
	      avg+=curoct->cell[icell].density;
	      if(curoct->cell[icell].density>maxd) maxd=curoct->cell[icell].density;
	      if(curoct->cell[icell].density<mind) mind=curoct->cell[icell].density;
	    }
	}while(nextoct!=NULL);

      printf("level=%d avg=%f mind=%f maxd=%f\n",level,avg/ncell,mind,maxd);
      }
    
  // recomputing the last oct;
  
  printf("==> memory state\n");
  printf("endoct=%p\n",endoct);
  
  endoct=lastoct[0];
  for(i=0;i<levelmax;i++) {
    if(lastoct[i]>endoct) endoct=lastoct[i];
    //printf("i=%d %p %p\n",i+1,firstoct[i],lastoct[i]);
  }
  //printf("endoct=%p\n",endoct);

#endif

  sprintf(filename,"data/levstart.%05d.p%05d",nsteps+1,cpu.rank);
  dumpcube(lmap,firstoct,0,filename);



#ifdef WMPI
  // ==================================== after refinement we should remap the boundary cells
  setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax);
#endif

  sprintf(filename,"data/cpustart.%05d.p%05d",nsteps+1,cpu.rank);
  dumpcube(lmap,firstoct,3,filename);


#if 1
  // ==================================== performing the CIC assignement
  //breakmpi();

  call_cic(levelmax,levelcoarse,firstoct,&cpu);

#ifdef WMPI
    // ==================================== performing the CIC BOUNDARY CORRECTION 

  mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);

  // ======================================= Density boundary mpi update 

  mpi_exchange(&cpu,sendbuffer,recvbuffer,1);
#endif

#endif

  sprintf(filename,"data/denstart.%05d.p%05d",nsteps+1,cpu.rank);
  printf("%s\n",filename);
  dumpcube(lmap,firstoct,1,filename);


#if 1
  // ==================================== Check the number of particles and octs
  multicheck(firstoct,npart,levelmax,cpu.rank);
#endif	 

#if 1
// ==================================== POISSON Testing the jacobi iteration
  printf("==> Poisson Start \n");
  int icomp,iter,niter=NITER;
  float norm_d;
  int nread;
  for(level=levelcoarse;level<=levelmax;level++)
    {
      // COARSE LEVEL TREATMENT ============================================
      if(level==levelcoarse){
	norm_d=0.;
	//fixbound_reg(grid,levelcoarse,NBND);
	for(iter=0;iter<niter;iter++){

	  if((iter%64==0)&&(cpu.rank==0)) printf("iter=%d ",iter);

	  nextoct=firstoct[level-1];
	  if(nextoct!=NULL){
	    dx=pow(0.5,level);
	    do{ 
	      curoct=nextoct;
	    
	      // First we gather the potential in all neighbors
	      for(icomp=0;icomp<=6;icomp++){
		memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
		nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
	      }
	      
	      // Second we gather the local density
	      memset(vcomp[7],0,stride*sizeof(float)); // reset the vcomp;
	      nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);

	      // 2.5 we contrast the density by removing the average density value
	      remove_avg(vcomp[7],stride,1.);

	      // we compute the square of the norm of the density (first iteration only)
	      if(iter==0) norm_d+=square(vcomp[7],nread);
	  
	      // Third we perform the calculation (eventually on GPU)
	      laplacian(vcomp,stride,dx);
	  

	      // Fourth we scatter back the potential estimation to the temp position
	      nextoct=scattercomp(curoct, vcomp[6], 6, 2, stride,&cpu);

	    }while(nextoct!=NULL);
	  }
	  
	  // 4.5 we copy the result in the temp position to the potential
	  nextoct=firstoct[level-1];
	  if(nextoct!=NULL){
	    do{ 
	      curoct=nextoct;
	    
	      memset(vcomp[0],0,stride*sizeof(float)); // reset the vcomp;

	      nextoct=gathercomp(curoct, vcomp[0], 6, 2, stride,&cpu,&nread); // getting the data in the temp field
	      nextoct=scattercomp(curoct, vcomp[0], 6, 1, stride,&cpu);

	    }while(nextoct!=NULL);
	  }

#ifdef WMPI
	  mpi_exchange(&cpu,sendbuffer,recvbuffer,2);
	  if(iter==0) MPI_Allreduce(MPI_IN_PLACE,&norm_d,1,MPI_FLOAT,MPI_SUM,cpu.comm);
#endif

	    // Fifth we compute the residuals

	    nextoct=firstoct[level-1];
	    float res=0.;
	    if(nextoct!=NULL){
	      dx=pow(0.5,level);
	      do{ 
		curoct=nextoct;
	      
		// First we gather the potential in all neighbors + local
		for(icomp=0;icomp<=6;icomp++){
		  memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
		  nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
		}

		// Second we gather the local density
		memset(vcomp[7],0,stride*sizeof(float)); // reset the vcomp;
		nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);

		// 2.5 we contrast the density by removing the average value
		remove_avg(vcomp[7],stride,1.);
	      
		// Third we perform the square of the residual
		res+=square_res(vcomp,nread,dx);

	      }while(nextoct!=NULL);
	    }
	    
#ifdef WMPI
	    // reducing the residuals
	    MPI_Barrier(cpu.comm);
	    float restot;
	    //if(iter%64==0) printf("res = %f on rank=%d\n",res,cpu.rank);
	    MPI_Allreduce(&res,&restot,1,MPI_FLOAT,MPI_SUM,cpu.comm);
	    res=restot;
#endif
	    if((iter%64==0)&&(cpu.rank==0)) printf("dens=%e res=%e relative residual=%e\n ",sqrt(norm_d),sqrt(res),sqrt(res/norm_d));
	    if(sqrt(res/norm_d)<acc) break;
	}
	
      }
      // FINE LEVEL TREATMENT ============================================
      else{
	norm_d=0.;
	
	//initial guess from parent cell
	nextoct=firstoct[level-1];
	//printf("first=%p next=%p\n",nextoct,nextoct->next);
	if(nextoct==NULL){
	  continue; // we skip to next level if the firstoct is empty
	}
	else{
	  do{ 
	    curoct=nextoct;
	    nextoct=curoct->next;
	    if(curoct->cpu!=cpu.rank) continue;
	    for(icell=0;icell<8;icell++){
	      curoct->cell[icell].pot=curoct->parent->pot;
	    }
	  }while(nextoct!=NULL);
	}
#ifdef WMPI
	mpi_exchange(&cpu,sendbuffer,recvbuffer,2);
#endif
	printf("guess ok\n");
#if 1
	// fine level relaxation
	for(iter=0;iter<niter;iter++){
	  nextoct=firstoct[level-1];
	  if((iter%16==0)&&(cpu.rank==0)) printf("level=%d iter=%d ",level,iter);
	  int ncell=0;
	  if(nextoct!=NULL){
	    dx=pow(0.5,level);
	    do{ 
	      ncell++;
	      curoct=nextoct;
	      // First we gather the potential in all neighbors
	      for(icomp=0;icomp<=6;icomp++){
		memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
		nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
	      }

	      // Second we gather the local density
	      memset(vcomp[7],0,stride*sizeof(float)); // reset the vcomp;
	      nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);
	      
	      // 2.5 we contrast the density by removing the average value
	      remove_avg(vcomp[7],stride,1.);
	      
	      // we compute the square of the norm of the density (first iteration only)
	      if(iter==0) norm_d+=square(vcomp[7],nread);
	      
	      // Third we perform the calculation (eventually on GPU)
	      laplacian(vcomp,stride,dx);
	      
	      // Fourth we scatter back the potential estimation
	      nextoct=scattercomp(curoct, vcomp[6], 6, 1, stride,&cpu);

	      
	    }while(nextoct!=NULL);

#ifdef WMPI
	    mpi_exchange(&cpu,sendbuffer,recvbuffer,2);
	    if(iter==0) MPI_Allreduce(MPI_IN_PLACE,&norm_d,1,MPI_FLOAT,MPI_SUM,cpu.comm);
#endif
	    // Fifth we compute the residuals
	    nextoct=firstoct[level-1];
	    float res=0.;
	    if(nextoct!=NULL){
	      dx=pow(0.5,level);
	      do{ 
		curoct=nextoct;
	      
		// First we gather the potential in all neighbors + local
		for(icomp=0;icomp<=6;icomp++){
		  memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
		  nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
		}

		// Second we gather the local density
		memset(vcomp[7],0,stride*sizeof(float)); // reset the vcomp;
		nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);

		// 2.5 we contrast the density by removing the average value
		remove_avg(vcomp[7],nread,1.);
	      
		// Third we perform the square of the residual
		res+=square_res(vcomp,nread,dx);
	      
	      }while(nextoct!=NULL);
	    }

#ifdef WMPI
	    // reducing the residuals
	    MPI_Barrier(cpu.comm);
	    float restot;
	    //if(iter%64==0) printf("res = %f on rank=%d\n",res,cpu.rank);
	    MPI_Allreduce(&res,&restot,1,MPI_FLOAT,MPI_SUM,cpu.comm);
	    res=restot;
#endif
	    if((iter%16==0)&&(cpu.rank==0)) printf("dens=%e res=%e relative residual=%e\n ",sqrt(norm_d),sqrt(res),sqrt(res/norm_d));
	    if(sqrt(res/norm_d)<acc) break;
	  }
	}
#endif
      }
    }
#endif



  sprintf(filename,"data/potstart.%05d.p%05d",nsteps+1,cpu.rank);
  printf("%s\n",filename);
  dumpcube(lmap,firstoct,2,filename);
	    

  // ==================================== Force calculation and velocity update   // corrector step


  if(nsteps!=0){
    forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dt*0.5,&cpu,sendbuffer,recvbuffer);
  }

  // ==================================== DUMP AFTER SYNCHRONIZATION
#if 1
  if(nsteps%NDUMP==0){
    // ===== Casting rays to fill a map

  /*   sprintf(filename,"data/level.%05d",nsteps); */
  /*   dumpmap(lmap,firstoct,0,filename,0.,1.); */
  /*   sprintf(filename,"data/dens.%05d",nsteps); */
  /*   dumpmap(lmap,firstoct,1,filename,0.,1.); */
  /*   sprintf(filename,"data/pot3d.%05d",nsteps); */
  /*   dumpcube(lmap,firstoct,2,filename); */
  /*   sprintf(filename,"data/lev3d.%05d",nsteps); */
  /*   dumpcube(lmap,firstoct,0,filename); */

  
  /* //==== Gathering particles for dump */

  /*   sprintf(filename,"data/part.%05d",nsteps); */
  /*   dumppart(firstoct,filename,npart,levelcoarse,levelmax); */

  }
#endif
  
  // ==================================== Force calculation and velocity update   // predictor step
  
  forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dt*0.5,&cpu,sendbuffer,recvbuffer);
  

  printf("Moving particles\n");



  // ==================================== Moving Particles + Oct management
  

  // Computing displacement (predictor)

  movepart(levelcoarse,levelmax,firstoct,dt);


  // Moving particles through cells (3 passes)

  partcellreorg(levelcoarse,levelmax,firstoct);

#ifdef WMPI

  // Communication of particles
  int deltan;
  deltan=mpi_exchange_part(&cpu,psendbuffer,precvbuffer,&lastpart);

  // Recounting particles
  npart=npart+deltan;
#endif

  //==== Gathering particles for dump

  sprintf(filename,"data/partstart.%05d.p%05d",nsteps+1,cpu.rank);
  dumppart(firstoct,filename,npart,levelcoarse,levelmax);

#if 1
  // ==================================== Check the number of particles and octs
  multicheck(firstoct,npart,levelmax,cpu.rank);
#endif	 

  //==================================== timestep completed, looping

#ifdef WMPI
  if(nsteps==1){
    printf("ABORTING !!\n");
    MPI_Barrier(cpu.comm);
    MPI_Abort(cpu.comm,42);
  }
#else
    abort();
#endif

  }
  return 0;
}
      
