#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
struct OCT *gathervecnei(struct OCT *octstart, int *vecnei, float *vec, int var, int *vecl, int stride, struct CPUINFO *cpu, int *nread)
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

struct OCT *gathervec(struct OCT *octstart, float *vec, char var, int *vecl, int stride, struct CPUINFO *cpu, int *nread)
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

struct OCT *gathervec2(struct OCT *octstart, float *vec, char var, int *vecl, int *vecicoarse, int *veccpu, int stride, struct CPUINFO *cpu, int *nread)
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
	  break;
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

#ifdef WHYDRO
struct OCT *gathervechydro(struct OCT *octstart, struct MULTIVECT *data, int stride, struct CPUINFO *cpu, int *nread)
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

	data->vec_d[ipos+icell*stride]=curoct->cell[icell].d;
	data->vec_u[ipos+icell*stride]=curoct->cell[icell].u;
	data->vec_v[ipos+icell*stride]=curoct->cell[icell].v;
	data->vec_w[ipos+icell*stride]=curoct->cell[icell].w;
	data->vec_p[ipos+icell*stride]=curoct->cell[icell].p;

#ifdef AXLFORCE
#ifdef SELFGRAV
	// we store the gravitational force in the new fields
	data->vec_unew[ipos+icell*stride]=curoct->cell[icell].fx;
	data->vec_vnew[ipos+icell*stride]=curoct->cell[icell].fy;
	data->vec_wnew[ipos+icell*stride]=curoct->cell[icell].fz;
	//if(curoct->cell[icell].d>0.5) printf("fx=%e\n",data->vec_unew[ipos+icell*stride]);
#endif
#endif
/* 	switch(var){ */
/* 	case 0: */
/* 	  vec[ipos+icell*stride]=curoct->cell[icell].density; */
/* 	  break; */
/* 	case 1 : */
/* 	  vec[ipos+icell*stride]=curoct->cell[icell].pot; */
/* 	  break; */
/* #ifdef WHYDRO */
/* 	case 10 : */
/* 	  vec[ipos+icell*stride]=curoct->cell[icell].d; */
/* 	  break; */
/* 	case 20 : */
/* 	  vec[ipos+icell*stride]=curoct->cell[icell].u; */
/* 	  break; */
/* 	case 30 : */
/* 	  vec[ipos+icell*stride]=curoct->cell[icell].v; */
/* 	  break; */
/* 	case 40 : */
/* 	  vec[ipos+icell*stride]=curoct->cell[icell].w; */
/* 	  break; */
/* 	case 50 : */
/* 	  vec[ipos+icell*stride]=curoct->cell[icell].p; */
/* 	  break; */
/* #endif */
/* 	} */

      }
      
      data->vecicoarse[ipos]=curoct->parent->idx; // we store the idx of the parent cell of the current oct
      data->vecl[ipos]=curoct->level; // assigning a level
      data->veccpu[ipos]=curoct->cpu; // assigning a cpu
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}
#endif

struct OCT *gathervec2_light(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int *nread)
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

struct OCT *scattervec(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int *nread)
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

struct OCT *scattervec_light(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int nread, int level)
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
void remove_valvec(float *vec, int nval, int stride, float avg, int level, int *vecl)
{
  int icell;
  int i;
  for(icell=0;icell<8;icell++){
    for(i=0;i<stride;i++){
      vec[i+icell*stride]=(vec[i+icell*stride]-avg)*(level==vecl[i]);
    }
  }
}

//============================================================================
float square_vec(float *vec, int nval, int stride, int level, int curcpu, int *vecl, int *veccpu)
{
  int icell;
  int i;
  float sum=0.;
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
float laplacian_vec(float *vecden,float *vecpot,float *vecpotnew,int *vecnei,int *vecl, int level, int nread,int stride,float dx,float omegam,float tsim){

  int inei,icell;
  int i;
  float temp;
  float res,restot=0.;
  int vnei[6],vcell[6];
  int idxnei;
  float ominterp=0.2;

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
      temp=temp-dx*dx*vecden[i+icell*stride]/6.*1.5*omegam/tsim;
#endif



      // we finsih the residual
      res=res-6.0*vecpot[i+icell*stride];

#ifndef TESTCOSMO
      res=res/(dx*dx)-4.0*M_PI*vecden[i+icell*stride];
#else
      res=res/(dx*dx)-1.5*omegam/tsim*vecden[i+icell*stride];
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
float laplacian_vec2(float *vecden,float *vecpot,float *vecpotnew,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,float dx,float factdens){

  int inei,icell,icellcoarse;
  int i;
  float temp;
  float res,restot=0.;
  int vnei[6],vcell[6];
  int vneic[6],vcellc[6];
  int idxnei;
  float ominterp=0.2;

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
int residual_vec2(float *vecden,float *vecpot,float *vecres,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,float dx,float factdens){

  int inei,icell,icellcoarse;
  int i;
  float temp;
  float res,restot=0.;
  int vnei[6],vcell[6];
  int vneic[6],vcellc[6];
  int idxnei;
  float ominterp=0.2;
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
struct OCT* gathercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int icell;
  int vnei[6],vcell[6];
  int ioct=0;
  float ominterp=0.2;
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


struct OCT* scattercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu)
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
 

float laplacian(float **vcomp, int stride, float dx, int locres){
  int i;
  float temp;
  float res,res2=0.;
  float tt;
  for(i=0;i<stride;i++){
    temp=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i])/6.-dx*dx*vcomp[7][i]/6.*4*M_PI;
    res=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i]-6.*vcomp[6][i])/(dx*dx)-4.*M_PI*vcomp[7][i];
    res2+=res*res;
    vcomp[locres][i]=temp;
  }
  return res2;
}

float laplaciancosmo(float **vcomp, int stride, float dx, int locres, float omegam, float a){
  int i;
  float temp;
  float res,res2=0.;
  float tt;

  for(i=0;i<stride;i++){
    temp=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i])/6.-dx*dx*vcomp[7][i]/6.*1.5*omegam/a;
    res=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i]-6.*vcomp[6][i])/(dx*dx)-1.5*omegam*vcomp[7][i]/a;
    res2+=res*res;
    vcomp[locres][i]=temp;
  }
  return res2;
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
