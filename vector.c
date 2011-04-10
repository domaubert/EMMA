#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"

//============================================================================
struct OCT *gathervecnei(struct OCT *octstart, int *vecnei, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int iread=0;
  int inei;
  int ipos=0;
  int icur;
  
  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->vecpos==-1){
	// the current oct has not been vectorized yet
	curoct->vecpos=ipos;
	ipos++;
      }
      
      // getting the vector element
      icur=curoct->vecpos;

      // we scan the neighbors
      for(inei=0;inei<6;inei++){
	if(curoct->nei[inei]->child!=NULL){ // the neighbor oct exists (always true for levelcoarse)
	  if(curoct->nei[inei]->child->vecpos==-1){ // the neighbor does not have a vector position
	    curoct->nei[inei]->child->vecpos=ipos;
	    ipos++;
	  }
	  vecnei[icur+inei*stride]=curoct->nei[inei]->child->vecpos; // we assign a neighbor
	}
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }

  (*nread)=iread;
  printf("ipos=%d\n",ipos);

  return nextoct;
}

//============================================================================

struct OCT *gathervec(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int *nread)
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

struct OCT *scattervec(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int nread)
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
	case 0:
	  curoct->cell[icell].density=vec[ipos+icell*stride];
	  break;
	case 1 :
	  curoct->cell[icell].pot=vec[ipos+icell*stride];
	  break;
	}
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}




//============================================================================
void remove_valvec(float *vec, int nval, int stride, float avg)
{
  int icell;
  int i;
  for(icell=0;icell<8;icell++){
    for(i=0;i<nval;i++){
      vec[i+icell*stride]=vec[i+icell*stride]-avg;
    }
  }
}

//============================================================================
float square_vec(float *vec, int nval, int stride)
{
  int icell;
  int i;
  float sum=0.;
  for(icell=0;icell<8;icell++){
    for(i=0;i<nval;i++){
      sum+=vec[i+icell*stride]*vec[i+icell*stride];
    }  
  }
  return sum;
}


//============================================================================
float laplacian_vec(float *vecden,float *vecpot,float *vecpotnew,int *vecnei,int nread,int stride,float dx,float omegam,float tsim){

  int inei,icell;
  int i;
  float temp;
  float res,restot=0.;
  int vnei[6],vcell[6];
  int idxnei;
  
  for(i=0;i<stride;i++){ // we scan the octs
    for(icell=0;icell<8;icell++){ // we scan the cells
      temp=0.;
      getcellnei(icell, vnei, vcell); // we get the neighbors

      // we compute the neighbor part of the laplacian
      for(inei=0;inei<6;inei++){ // we scan the neighbors
	if(vnei[inei]==6){ // the neighbor is in the current oct
	  temp+=vecpot[i+vcell[inei]*stride];
	  }
	else{
	  temp+=vecpot[vecnei[vnei[inei]]+vcell[inei]*stride];
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
