
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include <string.h>


// returns the cells on a neighbour face given by idx

void getfcell(int idx, int *fcell){
  
  switch(idx){
  case 0:
    fcell[0]=1;
    fcell[1]=3;
    fcell[2]=5;
    fcell[3]=7;
    break;
  case 1:
    fcell[0]=0;
    fcell[1]=2;
    fcell[2]=4;
    fcell[3]=6;
    break;
  case 2:
    fcell[0]=2;
    fcell[1]=3;
    fcell[2]=6;
    fcell[3]=7;
    break;
  case 3:
    fcell[0]=0;
    fcell[1]=1;
    fcell[2]=4;
    fcell[3]=5;
    break;
  case 4:
    fcell[0]=4;
    fcell[1]=5;
    fcell[2]=6;
    fcell[3]=7;
    break;
  case 5:
    fcell[0]=0;
    fcell[1]=1;
    fcell[2]=2;
    fcell[3]=3;
    break;
  }
}


///------------------------------------------------------------------------
// it flips the cell along a certain direction for transmissive boundary conditions
void flipcell(struct OCT *oct, int dir){
  
  int i;
  struct CELL temp;

  switch(dir){
  case 0:
    {
      char c1[4]={0,2,4,6};
      char c2[4]={1,3,5,7};
      for(i=0;i<4;i++){
	temp=oct->cell[c1[i]];
	memcpy(&(oct->cell[c1[i]]),&(oct->cell[c2[i]]),sizeof(struct CELL));
        memcpy(&(oct->cell[c2[i]]),&temp,sizeof(struct CELL));
      }
    }
    break;
  case 1:
    {
    char c1[4]={0,1,4,5};
    char c2[4]={2,3,6,7};
    for(i=0;i<4;i++){
      temp=oct->cell[c1[i]];
      memcpy(&(oct->cell[c1[i]]),&(oct->cell[c2[i]]),sizeof(struct CELL));
      memcpy(&(oct->cell[c2[i]]),&temp,sizeof(struct CELL));
    }
    break;
    }
  case 2:
    {
      char c1[4]={0,1,2,3};
      char c2[4]={4,5,6,7};
      for(i=0;i<4;i++){
	temp=oct->cell[c1[i]];
	memcpy(&(oct->cell[c1[i]]),&(oct->cell[c2[i]]),sizeof(struct CELL));
	memcpy(&(oct->cell[c2[i]]),&temp,sizeof(struct CELL));
      }
      break;
    }
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
//------------------------------------------------------------------------

#ifdef PIC
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

#endif
