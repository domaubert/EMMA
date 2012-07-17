#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include "oct.h"

//============================================================================
void clean_vec(int levelmax,struct OCT **firstoct)
{
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  for(level=1;level<=levelmax;level++) // looping over levels
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  curoct->vecpos=-1;
	}while(nextoct!=NULL);
      
    }
}


void clean_pot(int levelmax,struct OCT **firstoct)
{
  int level;
  int icell;
  struct OCT *nextoct;
  struct OCT *curoct;
  for(level=1;level<=levelmax;level++) // looping over levels
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) curoct->cell[icell].pot=0.;
	  curoct->vecpos=-1;
	}while(nextoct!=NULL);
      
    }
}


//============================================================================
struct OCT *gathervecnei(struct OCT *octstart, int *vecnei, REAL *vec, int var, int *vecl, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int iread=0;
  int inei;
  int ipos=0;
  int icur;
  int j;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->vecpos<0){
	// the current oct has not been vectorized yet
	curoct->vecpos=ipos;
	ipos++;
      }
      
      // getting the vector element
      icur=curoct->vecpos;

      // we scan the neighbors
      for(inei=0;inei<6;inei++){
	if(curoct->nei[inei]->child!=NULL){ // the neighbor oct exists (always true for levelcoarse)

	  if(curoct->nei[inei]->child->vecpos<0){ // the neighbor does not have a vector position
	    curoct->nei[inei]->child->vecpos=ipos;
	    ipos++;
	  }

	  vecnei[icur+inei*stride]=curoct->nei[inei]->child->vecpos; // we assign a neighbor
	  if(vecnei[icur+inei*stride]>=stride){
	    printf("error vecnei\n");
	    abort();
	  }
	}
	else{ // we need to interpolate from coarser level
	  
	  vecnei[icur+inei*stride]=ipos; // we assign a neighbor
	  vecl[ipos]=curoct->level-1; // we assign a level

	  if(vecnei[icur+inei*stride]>=stride){
	    printf("error vecnei\n");
	    abort();
	  }

	  // the field is directly filled at this stage
	  switch(var){
	  case(0):
	    for(j=0;j<8;j++) vec[ipos+j*stride]=curoct->nei[inei]->density;
	    break;
	  case(1):
	    for(j=0;j<8;j++) vec[ipos+j*stride]=curoct->nei[inei]->pot;
	    break;
	  }
	  ipos++;
	}
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }

  (*nread)=iread;

  return nextoct;
}

//============================================================================
int countvecocts(struct OCT *octstart, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct OCT* coarseoct;
  int iread=0;
  int inei;
  int ipos=0;
  int icur;
  int j;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if(curoct->vecpos<0){
	// the current oct has not been vectorized yet
	curoct->vecpos=ipos;
	ipos++;
      }
      
      // getting the vector element
      icur=curoct->vecpos;

      if(curoct->cpu==cpu->rank){ // we only consider octs which belong to current cpu
	// we scan the neighbors
	for(inei=0;inei<6;inei++){
	  /* if(curoct->nei[inei]==NULL){ */
	  /*   printf("rank=%d cpu=%d level=%d\n",cpu->rank,curoct->cpu,curoct->level); */
	  /*   abort(); */
	  /* } */

	  if(curoct->nei[inei]->child!=NULL){ // the neighbor oct exists (always true for levelcoarse)
	    if(curoct->nei[inei]->child->vecpos<0){ // the neighbor does not have a vector position
	      curoct->nei[inei]->child->vecpos=ipos;
	      ipos++;
	    }
	  }
	  else{ // we need to interpolate from coarser level
	  
	    coarseoct=cell2oct(curoct->nei[inei]);
	    if(coarseoct->vecpos<0){
	      coarseoct->vecpos=ipos;
	      ipos++;
	    }
	  }
	}
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }


  //printf("iread=%d ipos=%d\n",iread,ipos);
  (*nread)=iread; 

  return ipos; // we return the required number of octs
}

//============================================================================
struct OCT *gathervecnei2(struct OCT *octstart, int *vecnei, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct OCT* coarseoct;
  int iread=0;
  int inei;
  int ipos=0;
  int icur;
  int j;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->vecpos<0){
	// the current oct has not been vectorized yet
	curoct->vecpos=ipos;
	ipos++;
      }
      
      if(curoct->cpu==cpu->rank){ // we look for neighbors of current cpu octs only
	// getting the vector element
	icur=curoct->vecpos;

	// we scan the neighbors
	for(inei=0;inei<6;inei++){
	  if(curoct->nei[inei]->child!=NULL){ // the neighbor oct exists (always true for levelcoarse)

	    if(curoct->nei[inei]->child->vecpos<0){ // the neighbor does not have a vector position
	      curoct->nei[inei]->child->vecpos=ipos;
	      ipos++;
	    }

	    vecnei[icur+inei*stride]=curoct->nei[inei]->child->vecpos; // we assign a neighbor

#ifdef TRANSXP
	    if(inei==1){
	      if((curoct->nei[inei]->child->x-curoct->x)<0.){
		// the neighbor is a periodic mirror
		//printf("wouhou\n");
		vecnei[icur+inei*stride]=curoct->vecpos; // the curoct is its own neighbor
	      }
	    }
#endif

#ifdef TRANSXM
	    if(inei==0){
	      if((curoct->nei[inei]->child->x-curoct->x)>0.5){
		// the neighbor is a periodic mirror
		//printf("wouhou\n");
		vecnei[icur+inei*stride]=curoct->vecpos; // the curoct is its own neighbor
	      }
	    }
#endif

#ifdef TRANSYP
	    if(inei==3){
	      if((curoct->nei[inei]->child->y-curoct->y)<0.){
		// the neighbor is a periodic mirror
		//printf("wouhou\n");
		vecnei[icur+inei*stride]=curoct->vecpos; // the curoct is its own neighbor
	      }
	    }
#endif

#ifdef TRANSYM
	    if(inei==2){
	      if((curoct->nei[inei]->child->y-curoct->y)>0.5){
		// the neighbor is a periodic mirror
		//printf("wouhou\n");
		vecnei[icur+inei*stride]=curoct->vecpos; // the curoct is its own neighbor
	      }
	    }
#endif

#ifdef TRANSZP
	    if(inei==5){
	      if((curoct->nei[inei]->child->z-curoct->z)<0.){
		// the neighbor is a periodic mirror
		//printf("wouhou\n");
		vecnei[icur+inei*stride]=curoct->vecpos; // the curoct is its own neighbor
	      }
	    }
#endif

#ifdef TRANSZM
	    if(inei==4){
	      if((curoct->nei[inei]->child->z-curoct->z)>0.5){
		// the neighbor is a periodic mirror
		//printf("wouhou\n");
		vecnei[icur+inei*stride]=curoct->vecpos; // the curoct is its own neighbor
	      }
	    }
#endif

	    if(vecnei[icur+inei*stride]>=stride){
	      printf("error vecnei %d %d\n",stride,curoct->nei[inei]->child->vecpos);
	      abort();
	    }
	  }
	  else{ // we need to interpolate from coarser level
	  
	    coarseoct=cell2oct(curoct->nei[inei]);
	    if(coarseoct->vecpos<0){
	      coarseoct->vecpos=ipos;
	      ipos++;
	    }
	    vecnei[icur+inei*stride]=coarseoct->vecpos; // we assign a neighbor
	  }
	}
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }

  (*nread)=iread;

  return nextoct;
}


//============================================================================

struct OCT *gathervec(struct OCT *octstart, REAL *vec, char var, int *vecl, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos;
  int iread=0;
  int icell;
  
  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      //getting the vector element
     ipos=curoct->vecpos;

      // filling the values
      for(icell=0;icell<8;icell++){
	switch(var){
	case 0:
	  vec[ipos+icell*stride]=curoct->cell[icell].density;
	  break;
	case 1 :
	  vec[ipos+icell*stride]=curoct->cell[icell].pot;
	  break;
	}
      }
      
      vecl[ipos]=curoct->level; // assigning a level
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}

//============================================================================

struct OCT *gathervec2(struct OCT *octstart, REAL *vec, char var, int *vecl, int *vecicoarse, int *veccpu, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos;
  int iread=0;
  int icell;
  
  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      //getting the vector element
     ipos=curoct->vecpos;
     if(ipos<0) continue; // for coarseocts not involved in fine level calculations

      // filling the values
      for(icell=0;icell<8;icell++){
	switch(var){
	case 0:
	  vec[ipos+icell*stride]=curoct->cell[icell].density;
	  break;
	case 1 :
	  vec[ipos+icell*stride]=curoct->cell[icell].pot;
	  break;
#ifdef WHYDRO
	case 10 :
	  vec[ipos+icell*stride]=curoct->cell[icell].d;
	  break;
	case 20 :
	  vec[ipos+icell*stride]=curoct->cell[icell].u;
	  break;
	case 30 :
	  vec[ipos+icell*stride]=curoct->cell[icell].v;
	  break;
	case 40 :
	  vec[ipos+icell*stride]=curoct->cell[icell].w;
	  break;
	case 50 :
	  vec[ipos+icell*stride]=curoct->cell[icell].p;
q	  break;
#endif
	}
      }
      
      vecicoarse[ipos]=curoct->parent->idx; // we store the idx of the parent cell of the current oct
      vecl[ipos]=curoct->level; // assigning a level
      veccpu[ipos]=curoct->cpu; // assigning a cpu
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}



#ifdef WHYDRO2
//=====================================================================================================================
//
//    REFLECHISSONS UN PEU EN FAIT LA RECHERCHE EST RECURSIVE
//

// Structure de base

void recursive_neighbor_gather_oct(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct HGRID *stencil){

  int ix[6]={-1,1,0,0,0,0};
  int iy[6]={0,0,-1,1,0,0};
  int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;
  int face[8]={0,1,2,3,4,5,6,7};
  REAL dist;

  struct Wtype Wi[8];
  struct OCT *oct;
  struct OCT *neioct;
  struct CELL *neicell;

  if(order==1){
    ineiloc=inei;
  }
  else if(order==2){
    ineiloc=inei2;
  }
  else if(order==3){
    ineiloc=inei3;

  }

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


  oct=cell2oct(cell);
  neioct=cell2oct(neicell);

  // ============================ TRANSMISSIVE BOUNDARIES ====================
#ifdef TRANSXP
    if(ineiloc==1){
      dist=neioct->x-oct->x;
      if(dist<0.){
	neicell=cell;
	face[0]=1;
	face[1]=1;
	face[2]=3;
	face[3]=3;
	face[4]=5;
	face[5]=5;
	face[6]=7;
	face[7]=7;
      }
    }
#endif


#ifdef TRANSYP
    if(ineiloc==3){
      dist=neioct->y-oct->y;
      if(dist<0.){
	neicell=cell;
	face[0]=2;
	face[1]=3;
	face[2]=2;
	face[3]=3;
	face[4]=7;
	face[5]=6;
	face[6]=6;
	face[7]=7;
      }
    }
#endif

#ifdef TRANSZP
    if(ineiloc==5){
      dist=neioct->z-oct->z;
      if(dist<0.){
	neicell=cell;
	face[0]=4;
	face[1]=5;
	face[2]=6;
	face[3]=7;
	face[4]=4;
	face[5]=5;
	face[6]=6;
	face[7]=7;
      }
    }
#endif


      
#ifdef TRANSXM
    if(ineiloc==0){
      dist=neioct->x-oct->x;
      if(dist>0.5){
	neicell=cell;
	face[0]=0;
	face[1]=0;
	face[2]=2;
	face[3]=2;
	face[4]=4;
	face[5]=4;
	face[6]=6;
	face[7]=6;
      }
    }
#endif

#ifdef TRANSYM
    if(ineiloc==2){
      dist=neioct->y-oct->y;
      if(dist>0.5){
	neicell=cell;
	face[0]=0;
	face[1]=1;
	face[2]=0;
	face[3]=1;
	face[4]=4;
	face[5]=5;
	face[6]=4;
	face[7]=5;
      }
    }
#endif

#ifdef TRANSZM
    if(ineiloc==4){
      dist=neioct->z-oct->z;
      if(dist>0.5){
	neicell=cell;
	face[0]=0;
	face[1]=1;
	face[2]=2;
	face[3]=3;
	face[4]=0;
	face[5]=1;
	face[6]=2;
	face[7]=3;
      }
    }
#endif


  // ============================ END TRANSMISSIVE BOUNDARIES ====================

  if(neicell->child!=NULL){
    // optimal case
    for(icell=0;icell<8;icell++) memcpy(Wi+icell,&(neicell->child->cell[icell].field),sizeof(struct Wtype));
    }
  else{
    coarse2fine_hydro2(neicell,Wi);

  }

#ifdef WGRAV
    REAL floc[8*3];
    for(icell=0;icell<8;icell++) memcpy(floc+3*icell,neicell->child->cell[icell].f,sizeof(REAL)*3);
#endif



  for(icell=0;icell<8;icell++){
    memcpy(&(stencil->oct[ioct].cell[icell].field),Wi+face[icell],sizeof(struct Wtype)); //
#ifdef WGRAV
    memcpy(stencil->oct[ioct].cell[icell].f,floc+3*face[icell],sizeof(REAL)*3); //
#endif
  }

   
  
  // next order
  if(order==1){
    for(i=0;i<6;i++){
      if((i/2)==(inei/2)) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      recursive_neighbor_gather_oct(ioct2, inei, i, -1, 2, neicell, stencil);
    }
  }
  else if(order==2) {
    for(i=0;i<6;i++){
      if(((i/2)==(inei/2))||((i/2)==(inei2/2))) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      recursive_neighbor_gather_oct(ioct2, inei, inei2, i, 3, neicell, stencil);
    }
  }
}
#endif



#ifdef WHYDRO2
//=====================================================================================================================

struct OCT *gatherstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  
  int ix[6]={-1,1,0,0,0,0};
  int iy[6]={0,0,-1,1,0,0};
  int iz[6]={0,0,0,0,-1,1};
  int ioct;

  //printf("let's gather\n");
  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;


      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(&(stencil[iread].oct[13].cell[icell].field),&(curoct->cell[icell].field),sizeof(struct Wtype)); //
#ifdef WGRAV 
 	memcpy(stencil[iread].oct[13].cell[icell].f,curoct->cell[icell].f,sizeof(REAL)*3); // 
#endif 


      }
      
      //abort();
      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  recursive_neighbor_gather_oct(ioct, inei, -1, -1, 1, cell, stencil+iread);
	}
      
      /* if(curoct==SOCT) { printf(" SOCT=%e %e %e %e ir=%d\n",SOCT->cell[2].field.d,stencil[iread].oct[12].cell[3].field.d,stencil[iread].oct[12].cell[3].field.u,stencil[iread].oct[12].cell[3].field.p,iread);IR=iread;}  */
      /* if(curoct==SOCT2) {printf(" SOCT2=%e %e %e %e ir=%d\n",SOCT2->cell[0].field.d,stencil[iread].oct[12].cell[1].field.d,stencil[iread].oct[12].cell[1].field.u,stencil[iread].oct[12].cell[1].field.p,iread);IR2=iread;}  */


      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}


#endif






//=====================================================================================================================
//=====================================================================================================================

struct OCT *gathervec2_light(struct OCT *octstart, REAL *vec, char var, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos;
  int iread=0;
  int icell;
  
  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      // only the boundaries are relevant
      if(curoct->cpu==cpu->rank) continue;

      //getting the vector element
     ipos=curoct->vecpos;
     if(ipos<0) continue; // for coarseocts not involved in fine level calculations

      // filling the values
      for(icell=0;icell<8;icell++){
	switch(var){
	case 1:
	  vec[ipos+icell*stride]=curoct->cell[icell].pot;
	  break;
	case 0 :
	  vec[ipos+icell*stride]=curoct->cell[icell].density;
	  break;
	}
      }
      
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}

//============================================================================

struct OCT *checknei(struct OCT *octstart, int *vecnei, int stride)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int inei,ipos,iread;
  nextoct=octstart;
  iread=0;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      ipos=curoct->vecpos;
      for(inei=0;inei<6;inei++)
	{
	  if(curoct->nei[inei]->child->vecpos!=vecnei[ipos+inei*stride]){
	    printf("ERROR in checnknei\n");
	    abort();
	  }
	}
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}

//============================================================================

struct OCT *scattervec(struct OCT *octstart, REAL *vec, char var, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  
  nextoct=octstart;
  iread=0;

  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      //getting the vector element
     ipos=curoct->vecpos;

      // filling the values
      for(icell=0;icell<8;icell++){
	switch(var){
	case 1:
	  curoct->cell[icell].pot=vec[ipos+icell*stride];
	  break;
	case 0 :
	  curoct->cell[icell].density=vec[ipos+icell*stride];
	  break;
	case 2 :
	  curoct->cell[icell].temp=vec[ipos+icell*stride];
	  break;
#ifdef WHYDRO
	case 10 :
	  curoct->cell[icell].d=vec[ipos+icell*stride];
	  break;
	case 20 :
	  curoct->cell[icell].u=vec[ipos+icell*stride];
	  break;
	case 30 :
	  curoct->cell[icell].v=vec[ipos+icell*stride];
	  break;
	case 40 :
	  curoct->cell[icell].w=vec[ipos+icell*stride];
	  break;
	case 50 :
	  curoct->cell[icell].p=vec[ipos+icell*stride];
	  break;
#endif
	}
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}

#ifdef WHYDRO
struct OCT *scattervechydro(struct OCT *octstart, struct MULTIVECT *data, int stride, struct CPUINFO *cpu)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  
  nextoct=octstart;
  iread=0;

  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      //getting the vector element
      ipos=curoct->vecpos;

      // filling the values
      for(icell=0;icell<8;icell++){
	curoct->cell[icell].d=data->vec_dnew[ipos+icell*stride];
	curoct->cell[icell].u=data->vec_unew[ipos+icell*stride];
	curoct->cell[icell].v=data->vec_vnew[ipos+icell*stride];
	curoct->cell[icell].w=data->vec_wnew[ipos+icell*stride];
	curoct->cell[icell].p=data->vec_pnew[ipos+icell*stride];
	
	if(data->vec_pnew[ipos+icell*stride]<0){
	  printf("NEgative pressure scatter\n");
	  abort();
	}

/* 	switch(var){ */
/* 	case 1: */
/* 	  curoct->cell[icell].pot=vec[ipos+icell*stride]; */
/* 	  break; */
/* 	case 0 : */
/* 	  curoct->cell[icell].density=vec[ipos+icell*stride]; */
/* 	  break; */
/* 	case 2 : */
/* 	  curoct->cell[icell].temp=vec[ipos+icell*stride]; */
/* 	  break; */
/* #ifdef WHYDRO */
/* 	case 10 : */
/* 	  curoct->cell[icell].d=vec[ipos+icell*stride]; */
/* 	  break; */
/* 	case 20 : */
/* 	  curoct->cell[icell].u=vec[ipos+icell*stride]; */
/* 	  break; */
/* 	case 30 : */
/* 	  curoct->cell[icell].v=vec[ipos+icell*stride]; */
/* 	  break; */
/* 	case 40 : */
/* 	  curoct->cell[icell].w=vec[ipos+icell*stride]; */
/* 	  break; */
/* 	case 50 : */
/* 	  curoct->cell[icell].p=vec[ipos+icell*stride]; */
/* 	  break; */
/* #endif */
/* 	} */
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}
#endif

//=======================================================================================================
//=======================================================================================================

#ifdef WHYDRO2
struct OCT *scatterstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int ioct[7]={12,14,10,16,4,22,13};
  
  nextoct=octstart;
  iread=0;

  //printf("let's scatter\n");
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(&(curoct->cell[icell].flux),&(stencil[iread].new.cell[icell].flux),sizeof(REAL)*NFLUX);
      }
      
      /* if(SOCT==curoct){ */
      /*  	printf(" F SOCT = %e %e\n",stencil[iread].new.cell[2].flux[0],stencil[iread].new.cell[2].flux[0+NVAR]);  */
      /* } */
      /* if(SOCT2==curoct){ */
      /*  	printf(" F SOCT2 = %e %e\n",stencil[iread].new.cell[0].flux[0],stencil[iread].new.cell[0].flux[0+NVAR]);  */
      /* 	//abort(); */
      /* } */

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}
#endif

//=======================================================================================================
//=======================================================================================================

struct OCT *scattervec_light(struct OCT *octstart, REAL *vec, char var, int stride, struct CPUINFO *cpu, int nread, int level)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  
  nextoct=octstart;
  iread=0;

  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if((curoct->border!=1)||(curoct->level!=level)) continue; // only border octs needs to be scattered

      //getting the vector element
     ipos=curoct->vecpos;

      // filling the values
      for(icell=0;icell<8;icell++){
	switch(var){
	case 1:
	  curoct->cell[icell].pot=vec[ipos+icell*stride];
	  break;
	case 0 :
	  curoct->cell[icell].density=vec[ipos+icell*stride];
	  break;
	}
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}





//============================================================================
void remove_valvec(REAL *vec, int nval, int stride, REAL avg, int level, int *vecl)
{
  int icell;
  int i;
  for(icell=0;icell<8;icell++){
    for(i=0;i<stride;i++){
      vec[i+icell*stride]=(vec[i+icell*stride]-avg)*(level==vecl[i])/avg;
    }
  }
}

//============================================================================
REAL square_vec(REAL *vec, int nval, int stride, int level, int curcpu, int *vecl, int *veccpu)
{
  int icell;
  int i;
  REAL sum=0.;
  for(icell=0;icell<8;icell++){
    for(i=0;i<stride;i++){
#ifdef WMPI
      sum+=vec[i+icell*stride]*vec[i+icell*stride]*(level==vecl[i])*(curcpu==veccpu[i]);
#else
      sum+=vec[i+icell*stride]*vec[i+icell*stride]*(level==vecl[i]);
#endif
    }  
  }
  return sum;
}


//============================================================================
REAL laplacian_vec(REAL *vecden,REAL *vecpot,REAL *vecpotnew,int *vecnei,int *vecl, int level, int nread,int stride,REAL dx,REAL omegam,REAL tsim){

  int inei,icell;
  int i;
  REAL temp;
  REAL res,restot=0.;
  int vnei[6],vcell[6];
  int idxnei;
  REAL ominterp=0.2;

  for(i=0;i<stride;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells

      // we skip octs which do not belong to the current level
      if(vecl[i]!=level){
	vecpotnew[i+icell*stride]=vecpot[i+icell*stride]; // nevertheless we copy the (fixed) potential in the new potential
	continue; 
      }

      
      temp=0.;
      getcellnei(icell, vnei, vcell); // we get the neighbors

      // we compute the neighbor part of the laplacian
      for(inei=0;inei<6;inei++){ // we scan the neighbors
	if(vnei[inei]==6){ // the neighbor is in the current oct
	  temp+=vecpot[i+vcell[inei]*stride];
	  }
	else{
	  if(vecl[vecnei[i+vnei[inei]*stride]]==vecl[i]){ // the octs share the same level
	    temp+=vecpot[vecnei[i+vnei[inei]*stride]+vcell[inei]*stride];
	  }
	  else{ // mixing values from two different levels
	    temp+=vecpot[vecnei[i+vnei[inei]*stride]+vcell[inei]*stride]*(1.-ominterp)+vecpot[i+icell*stride]*ominterp;
	  }
	}
      }

      //we setup the residual
      res=temp;

      // we finish the laplacian
      temp=temp/6.0;
#ifndef TESTCOSMO
      temp=temp-dx*dx*vecden[i+icell*stride]/6.*4.*M_PI;
#else
#ifdef SUPERCOMOV
      temp=temp-dx*dx*vecden[i+icell*stride]/6.*6.0*tsim;
#else
      temp=temp-dx*dx*vecden[i+icell*stride]/6.*1.5*omegam/tsim;
#endif
#endif



      // we finsih the residual
      res=res-6.0*vecpot[i+icell*stride];

#ifndef TESTCOSMO
      res=res/(dx*dx)-4.0*M_PI*vecden[i+icell*stride];
#else
#ifdef SUPERCOMOV
      res=res/(dx*dx)-6.0*tsim*vecden[i+icell*stride];
#else
      res=res/(dx*dx)-1.5*omegam/tsim*vecden[i+icell*stride];
#endif
#endif
      restot+=res*res;

      // we store the new value
      vecpotnew[i+icell*stride]=temp;

      // ready for the next cell
    }
    //ready for the next oct
  }

  return restot;
}


//============================================================================
REAL laplacian_vec2(REAL *vecden,REAL *vecpot,REAL *vecpotnew,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,REAL dx,REAL factdens){

  int inei,icell,icellcoarse;
  int i;
  REAL temp;
  REAL res,restot=0.;
  int vnei[6],vcell[6];
  int vneic[6],vcellc[6];
  int idxnei;
  REAL ominterp=0.2;

  for(i=0;i<stride;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells

      // we skip octs which do not belong to the current level
      if(vecl[i]!=level){
	vecpotnew[i+icell*stride]=vecpot[i+icell*stride]; // nevertheless we copy the (fixed) potential in the new potential
	continue; 
      }

      // we skip octs which do not belong to the current cpu
#ifdef WMPI
      if(veccpu[i]!=curcpu){
      	vecpotnew[i+icell*stride]=vecpot[i+icell*stride]; // nevertheless we copy the (fixed) potential in the new potential
      	continue;
      }
#endif


      temp=0.;
      getcellnei(icell, vnei, vcell); // we get the neighbors

      // we compute the neighbor part of the laplacian
      for(inei=0;inei<6;inei++){ // we scan the neighbors
	if(vnei[inei]==6){ // the neighbor is in the current oct
	  temp+=vecpot[i+vcell[inei]*stride];
	  }
	else{
	  if(vecl[vecnei[i+vnei[inei]*stride]]==vecl[i]){ // the octs share the same level
	    temp+=vecpot[vecnei[i+vnei[inei]*stride]+vcell[inei]*stride];
	  }
	  else{ // mixing values from two different levels
	    //
	    icellcoarse=vecicoarse[i];
	    getcellnei(icellcoarse,vneic,vcellc);
	    temp+=vecpot[vecnei[i+vnei[inei]*stride]+vcellc[inei]*stride]*(1.-ominterp)+vecpot[i+icell*stride]*ominterp;
	  }
	}
      }

      //we setup the residual
      res=temp;

      // we finish the laplacian
      temp=temp/6.0;
      temp=temp-dx*dx*vecden[i+icell*stride]/6.*factdens;

      // we finsih the residual
      res=res-6.0*vecpot[i+icell*stride];
      res=res/(dx*dx)-factdens*vecden[i+icell*stride];
      restot+=res*res;

      // we store the new value
      vecpotnew[i+icell*stride]=temp;

      // ready for the next cell
    }
    //ready for the next oct
  }

  return restot;
}



//============================================================================
int residual_vec2(REAL *vecden,REAL *vecpot,REAL *vecres,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,REAL dx,REAL factdens){

  int inei,icell,icellcoarse;
  int i;
  REAL temp;
  REAL res,restot=0.;
  int vnei[6],vcell[6];
  int vneic[6],vcellc[6];
  int idxnei;
  REAL ominterp=0.2;
  int count=0;
  
  for(i=0;i<stride;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells

      // we skip octs which do not belong to the current level
      if(vecl[i]!=level){
	vecres[i+icell*stride]=0.; // nevertheless we copy the (fixed) potential in the new potential
	continue; 
      }

      count++;

      // we skip octs which do not belong to the current cpu
#ifdef WMPI
      if(veccpu[i]!=curcpu){
      	vecres[i+icell*stride]=0.; // nevertheless we copy the (fixed) potential in the new potential
      	continue;
      }
#endif

      
      temp=0.;
      getcellnei(icell, vnei, vcell); // we get the neighbors

      // we compute the neighbor part of the laplacian
      for(inei=0;inei<6;inei++){ // we scan the neighbors
	if(vnei[inei]==6){ // the neighbor is in the current oct
	  temp+=vecpot[i+vcell[inei]*stride];
	  }
	else{
	  if(vecl[vecnei[i+vnei[inei]*stride]]==vecl[i]){ // the octs share the same level
	    temp+=vecpot[vecnei[i+vnei[inei]*stride]+vcell[inei]*stride];
	  }
	  else{ // mixing values from two different levels
	    //
	    icellcoarse=vecicoarse[i];
	    getcellnei(icellcoarse,vneic,vcellc);
	    temp+=vecpot[vecnei[i+vnei[inei]*stride]+vcellc[inei]*stride]*(1.-ominterp)+vecpot[i+icell*stride]*ominterp;
	  }
	}
      }

      //we setup the residual
      res=temp;

      // we finish the residual
      res=res-6.0*vecpot[i+icell*stride];
      res=res/(dx*dx)-factdens*vecden[i+icell*stride];

      // we store the residual value
      vecres[i+icell*stride]=res;

      // ready for the next cell
    }
    //ready for the next oct
  }

  //  printf("count=%d\n",count);
  return 0;
}

//============================================================================
struct OCT* gathercomp(struct OCT *octstart, REAL *vec, char nei, char var, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int icell;
  int vnei[6],vcell[6];
  int ioct=0;
  REAL ominterp=0.2;
  int iread=0;
  int vvnei[8];
  int vvcell[8];
  int j;

  nextoct=octstart;
  if(nei==6){
    if(nextoct!=NULL){
      do{// sweeping level
	ioct++;
	curoct=nextoct;
	nextoct=curoct->next;

	if(cpu->rank!=curoct->cpu) continue; // we consider only current cpu octs
	
	for(icell=0;icell<8;icell++){
	  int kk=0;
	  //getcellnei(icell, vnei, vcell);
	  switch(var){
	    //=================================
	  case 0: //------>  density
	    vec[iread]=curoct->cell[icell].density;
	    break;
	    //=================================
	  case 1: //------>  pot
	    vec[iread]=curoct->cell[icell].pot;
	    break;
	    //=================================
	  case 2: //------>  TEMP
	    vec[iread]=curoct->cell[icell].temp;
	    break;
	  }
	  iread++;
	}
      }while((nextoct!=NULL)&&(iread<stride));
    }
  }
  else{

    /// preparing arrays of neighbor and cells
    for(icell=0;icell<8;icell++){
      getcellnei(icell, vnei, vcell);
      vvnei[icell] =vnei[nei];
      vvcell[icell]=vcell[nei];
    }

    if(nextoct!=NULL){
      do{// sweeping level
	curoct=nextoct;
	nextoct=curoct->next;
	
	if(cpu->rank!=curoct->cpu) continue; // we consider only current cpu octs
       
	for(icell=0;icell<8;icell++){
	  switch(var){
	    //=================================
	  case 1: //------>  pot
	    if(vvnei[icell]==6){
	      vec[iread]=curoct->cell[vvcell[icell]].pot;
	    }
	    else{
	      if(curoct->nei[vvnei[icell]]->child!=NULL){
		vec[iread]=curoct->nei[vvnei[icell]]->child->cell[vvcell[icell]].pot;
	      }
	      else{
		vec[iread]=curoct->nei[vvnei[icell]]->pot*(1.-ominterp)+curoct->cell[icell].pot*(ominterp);
	      }
	    }
	    break;
	    //=================================
	  case 0: //------>  density
	    if(vvnei[icell]==6){
	      vec[iread]=curoct->cell[vvcell[icell]].density;
	    }
	    else{
	      if(curoct->nei[vvnei[icell]]->child!=NULL){
		vec[iread]=curoct->nei[vvnei[icell]]->child->cell[vvcell[icell]].density;
	      }
	      else{
		vec[iread]=curoct->nei[vvnei[icell]]->density;
	      }
	    }
	    break;
	    //=================================
	  case 2: //------>  TEMP
	      if(vvnei[icell]==6){
		vec[iread]=curoct->cell[vvcell[icell]].temp;
	      }
	      else{
		if(curoct->nei[vvnei[icell]]->child!=NULL){
		  vec[iread]=curoct->nei[vvnei[icell]]->child->cell[vvcell[icell]].temp;
		}
		else{
		  vec[iread]=curoct->nei[vvnei[icell]]->temp;
		}
	      }
	    break;
	  }
	  iread++;
	}
      }while((nextoct!=NULL)&&(iread<stride));
    }
  }

  (*nread)=iread; // relevant for subsequent reductions

  return nextoct;
  
}
//============================================================================


struct OCT* scattercomp(struct OCT *octstart, REAL *vec, char nei, char var, int stride, struct CPUINFO *cpu)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int nread=0;
  int icell;
  int vnei[6],vcell[6];
  int vvnei[8];
  int vvcell[8];

  nextoct=octstart;
  
  if(nei==6){
    if(nextoct!=NULL){
      do{// sweeping level
	
	curoct=nextoct;
	nextoct=curoct->next;

	if(cpu->rank!=curoct->cpu) continue;

	for(icell=0;icell<8;icell++){
	  switch(var){
	    //==============================================
	  case 0: //------>  density
	    curoct->cell[icell].density=vec[nread];
	    break;
	    //==============================================
	  case 1: //------>  potential
	    curoct->cell[icell].pot=vec[nread];
	    break;
	    //==============================================
	  case 2: //------>  temp
	    curoct->cell[icell].temp=vec[nread];
	    break;
#ifdef AXLFORCE
	    //==============================================
	  case 3: //------>  temp
	    curoct->cell[icell].fx=vec[nread];
	    break;
	    //==============================================
	  case 4: //------>  temp
	    curoct->cell[icell].fy=vec[nread];
	    break;
	    //==============================================
	  case 5: //------>  temp
	    curoct->cell[icell].fz=vec[nread];
	    break;
#endif
	  }
	  nread++;
	}
      }while((nextoct!=NULL)&&(nread<stride));
    }
  }
  else{

    /// preparing arrays of neighbor and cells
    for(icell=0;icell<8;icell++){
      getcellnei(icell, vnei, vcell);
      vvnei[icell] =vnei[nei];
      vvcell[icell]=vcell[nei];
    }

    if(nextoct!=NULL){
      do{// sweeping level
      
	curoct=nextoct;
	nextoct=curoct->next;

	if(cpu->rank!=curoct->cpu) continue;

	for(icell=0;icell<8;icell++){
	  switch(var){
	    //==============================================
	  case 0: //------>  density
	    if(vvnei[icell]==6){
	      curoct->cell[vvcell[icell]].density=vec[nread];
	    }
	    else{
	      if(curoct->nei[vvnei[icell]]->child!=NULL){
		curoct->nei[vvnei[icell]]->child->cell[vvcell[icell]].density=vec[nread];
	      }
	      else{
		curoct->nei[vvnei[icell]]->density=vec[nread];
	      }
	    }
	    break;
	    //==============================================
	  case 1: //------>  potential
	    if(vvnei[icell]==6){
	      curoct->cell[vvcell[icell]].pot=vec[nread];
	    }
	    else{
	      if(curoct->nei[vvnei[icell]]->child!=NULL){
		curoct->nei[vvnei[icell]]->child->cell[vvcell[icell]].pot=vec[nread];
	      }
	      else{
		curoct->nei[vvnei[icell]]->pot=vec[nread];
	      }
	    }
	    break;
	    //==============================================
	  case 2: //------>  temp
	    if(vvnei[icell]==6){
	      curoct->cell[vvcell[icell]].temp=vec[nread];
	    }
	    else{
	      if(curoct->nei[vvnei[icell]]->child!=NULL){
		curoct->nei[vvnei[icell]]->child->cell[vvcell[icell]].temp=vec[nread];
	      }
	      else{
		curoct->nei[vvnei[icell]]->temp=vec[nread];
	      }
	    }
	    break;
	  }
	  nread++;
	}
      }while((nextoct!=NULL)&&(nread<stride));
    }
  }

  return nextoct;
}
  

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
 

REAL laplacian(REAL **vcomp, int stride, REAL dx, int locres){
  int i;
  REAL temp;
  REAL res,res2=0.;
  REAL tt;
  for(i=0;i<stride;i++){
    temp=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i])/6.-dx*dx*vcomp[7][i]/6.*4*M_PI;
    res=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i]-6.*vcomp[6][i])/(dx*dx)-4.*M_PI*vcomp[7][i];
    res2+=res*res;
    vcomp[locres][i]=temp;
  }
  return res2;
}

REAL laplaciancosmo(REAL **vcomp, int stride, REAL dx, int locres, REAL omegam, REAL a){
  int i;
  REAL temp;
  REAL res,res2=0.;
  REAL tt;

  for(i=0;i<stride;i++){
    temp=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i])/6.-dx*dx*vcomp[7][i]/6.*1.5*omegam/a;
    res=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i]-6.*vcomp[6][i])/(dx*dx)-1.5*omegam*vcomp[7][i]/a;
    res2+=res*res;
    vcomp[locres][i]=temp;
  }
  return res2;
}



void grad(REAL **vcomp, int stride, REAL dx, int dir){

  int i;
  REAL temp;
  REAL om=1.;
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
#ifndef AXLFORCE
    vcomp[6][i]=temp;
#else
    vcomp[dir+6][i]=temp;
#endif
  }
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

void remove_avg(REAL *vcomp, int stride, REAL avg)
{
  int i;
  for(i=0;i<stride;i++){
    vcomp[i]=vcomp[i]-avg;
  }
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

REAL square(REAL *vcomp, int stride)
{
  int i;
  REAL res=0.;
  for(i=0;i<stride;i++){
    res+=vcomp[i]*vcomp[i];
  }
  return res;
}

REAL square_res(REAL **vcomp, int stride, REAL dx)
{
  int i;
  REAL res=0.;
  REAL res2=0.;
  for(i=0;i<stride;i++){
    res2=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i]-6.*vcomp[6][i])/(dx*dx)-vcomp[7][i]*4.*M_PI;
    res+=res2*res2;
  }
  return res;
}
