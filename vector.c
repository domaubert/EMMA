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
  int ncell;

  nextoct=octstart;
  (*nread)=0; // nread will return the effective number of CELLS read (important for reductions)

  if(nextoct!=NULL){
    do{// sweeping level
	ioct++;
	curoct=nextoct;
	nextoct=curoct->next;

	if(cpu->rank!=curoct->cpu) continue; // we consider only current cpu octs
	
	for(icell=0;icell<8;icell++){
	  int kk=0;
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
		  kk=1;
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


struct OCT* gathercompnew(struct OCT *octstart, float **vec, char *nei, char *var, int stride, struct CPUINFO *cpu, int *nread, int ncomp)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int icell;
  int vnei[6],vcell[6];
  int ioct=0;
  float ominterp=0.2;
  nextoct=octstart;
  int ncell=0;
  int j;

  (*nread)=0; // nread will return the effective number of CELLS read (important for reductions)

  if(nextoct!=NULL){
    do{// sweeping level
	ioct++;
	curoct=nextoct;
	nextoct=curoct->next;

	if(cpu->rank!=curoct->cpu) continue; // we consider only current cpu octs

	for(icell=0;icell<8;icell++){ // we scan the cells
	  ncell++;
	  getcellnei(icell, vnei, vcell);

	  for(j=0;j<ncomp;j++){ // we scan all the requested calculations
	    switch(var[j]){
	      //=================================
	    case 0: //------>  density
	      if(nei[j]==6){
		vec[j][(*nread)]=curoct->cell[icell].density;
	      }
	      else{
		if(vnei[nei[j]]==6){
		  vec[j][(*nread)]=curoct->cell[vcell[nei[j]]].density;
		}
		else{
		  if(curoct->nei[vnei[nei[j]]]->child!=NULL){
		    vec[j][(*nread)]=curoct->nei[vnei[nei[j]]]->child->cell[vcell[nei[j]]].density;
		  }
		  else{
		    vec[j][(*nread)]=curoct->nei[vnei[nei[j]]]->density;
		  }
		}
	      }
	      break;
	      //=================================
	    case 1: //------>  pot
	      if(nei[j]==6){
		vec[j][(*nread)]=curoct->cell[icell].pot;
	      }
	      else{
		if(vnei[nei[j]]==6){
		  vec[j][(*nread)]=curoct->cell[vcell[nei[j]]].pot;
		}
		else{
		  if(curoct->nei[vnei[nei[j]]]->child!=NULL){
		    vec[j][(*nread)]=curoct->nei[vnei[nei[j]]]->child->cell[vcell[nei[j]]].pot;
		  }
		  else{
		    vec[j][(*nread)]=curoct->nei[vnei[nei[j]]]->pot*(1.-ominterp)+curoct->cell[icell].pot*(ominterp);
		  }
		}
	      }
	      break;
	      //=================================
	    case 2: //------>  TEMP
	      if(nei[j]==6){
		vec[j][(*nread)]=curoct->cell[icell].temp;
	      }
	      else{
		if(vnei[nei[j]]==6){
		  vec[j][(*nread)]=curoct->cell[vcell[nei[j]]].temp;
		}
		else{
		  if(curoct->nei[vnei[nei[j]]]->child!=NULL){
		    vec[j][(*nread)]=curoct->nei[vnei[nei[j]]]->child->cell[vcell[nei[j]]].temp;
		  }
		  else{
		    vec[j][(*nread)]=curoct->nei[vnei[nei[j]]]->temp;
		  }
		}
	      }
	      break;

	    }
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
  for(i=0;i<stride;i++){
    temp=(vcomp[0][i]+vcomp[1][i]+vcomp[2][i]+vcomp[3][i]+vcomp[4][i]+vcomp[5][i])/6.-dx*dx*vcomp[7][i]/6.*4*M_PI;
    vcomp[8][i]=temp;
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
