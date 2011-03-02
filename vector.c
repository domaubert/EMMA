#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"

//============================================================================
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


struct OCT* gathercompnew(struct OCT *octstart, float **vec, char *nei, char *var, int stride, struct CPUINFO *cpu, int *nread, int ncomp)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct OCT oct;
  int icell;
  int vnei[6],vcell[6];
  int ioct=0;
  float ominterp=0.2;
  int ncell=0;
  int j;
  double t1,t2,t3;
  int iread=0;
  int idx;

  int vvnei[6*8];
  int vvcell[6*8];
  t1=MPI_Wtime();


  /// preparing arrays of neighbor and cells
  for(icell=0;icell<8;icell++){
    getcellnei(icell, vnei, vcell);
    for(j=0;j<6;j++){
      vvnei[j*8+icell] =vnei[j];
      vvcell[j*8+icell]=vcell[j];
    }
  }
  t2=MPI_Wtime();

  nextoct=octstart;
  //  (*nread)=0; // nread will return the effective number of CELLS read (important for reductions)
  if(nextoct!=NULL){
    do{// sweeping level
	ioct++;
	curoct=nextoct;
	nextoct=curoct->next;

	if(cpu->rank!=curoct->cpu) continue; // we consider only current cpu octs

	for(j=0;j<ncomp;j++){ // we scan all the requested calculations
	  for(icell=0;icell<8;icell++){ // we scan the cells

	    idx=nei[j]*8+icell;
	    switch(var[j]){
	      //=================================
	    case 0: //------>  density
	      if(nei[j]==6){
		vec[j][iread+icell]=curoct->cell[icell].density;
	      }
	      else{
		if(vvnei[idx]==6){
		  vec[j][iread+icell]=curoct->cell[vvcell[idx]].density;
		}
		else{
		  if(curoct->nei[vvnei[idx]]->child!=NULL){
		    vec[j][iread+icell]=curoct->nei[vvnei[idx]]->child->cell[vvcell[idx]].density;
		  }
		  else{
		    vec[j][iread+icell]=curoct->nei[vvnei[idx]]->density;
		  }
		}
	      }
	      break;
	      //=================================
	    case 1: //------>  pot
	      if(nei[j]==6){
		vec[j][iread+icell]=curoct->cell[icell].pot;
	      }
	      else{
		if(vvnei[idx]==6){
		  vec[j][iread+icell]=curoct->cell[vvcell[idx]].pot;
		}
		else{
		  if(curoct->nei[vvnei[idx]]->child!=NULL){
		    vec[j][iread+icell]=curoct->nei[vvnei[idx]]->child->cell[vvcell[idx]].pot;
		  }
		  else{
		    vec[j][iread+icell]=curoct->nei[vvnei[idx]]->pot*(1.-ominterp)+curoct->cell[icell].pot*(ominterp);
		  }
		}
	      }
	      break;
	      //=================================
	    case 2: //------>  TEMP
	      if(nei[j]==6){
		vec[j][iread+icell]=curoct->cell[icell].temp;
	      }
	      else{
		if(vvnei[idx]==6){
		  vec[j][iread+icell]=curoct->cell[vvcell[idx]].temp;
		}
		else{
		  if(curoct->nei[vvnei[idx]]->child!=NULL){
		    vec[j][iread+icell]=curoct->nei[vvnei[idx]]->child->cell[vvcell[idx]].temp;
		  }
		  else{
		    vec[j][iread+icell]=curoct->nei[vvnei[idx]]->temp;
		  }
		}
	      }
	      break;

	    }
	  }
	}
	iread+=8;
    }while((nextoct!=NULL)&&(iread<stride));
  }

  (*nread)=iread;
  t3=MPI_Wtime();
  //if(cpu->rank==0) printf("tgather=%e neical=%e dir=%e\n",t3-t1,t2-t1,t3-t2);
  return nextoct;
  
}




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
