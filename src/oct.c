
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include <string.h>
#include "cic.h"

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
  unsigned long int adress;
  struct OCT *poct;
  adress=(unsigned long int) cell;
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

void getneicell_6(struct CELL *cell, struct CELL** neicell ){
/// return a pointer table with the 6 neighbors of a given cell

    int vnei[6],vcell[6];
    getcellnei(cell->idx, vnei, vcell);

    struct OCT* curoct = cell2oct(cell);

    struct OCT* neioct[6];

    int i;
    for(i=0;i<6;i++){
        neicell[i] = NULL;
        if(vnei[i]==6){
            neicell[i] = &curoct->cell[vcell[i]];
        }else{
            struct OCT* neioct = curoct->nei[vnei[i]]->child;
            if(neioct != NULL){
                neicell[i] = &neioct->cell[vcell[i]];
            }
        }
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



////////////////////////////////////////////////////////////////////////////////

void cleanOctList(struct CPUINFO *cpu, struct RUNPARAMS *param, int level){
  int i;
  for(i=0;i<param->ngridmax; i++)  cpu->octList[level-1][i] = NULL;
}


void setOctList(struct OCT *firstoct, struct CPUINFO *cpu, struct RUNPARAMS *param, int level){

 // printf("Building octList\n");

  struct OCT  *nextoct = firstoct;


//  cleanOctList(cpu,param,level);

  int nOct = 0;
  do{ if(nextoct==NULL) 		continue;
        struct OCT  *curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue;

      cpu->octList[level-1][nOct++] = curoct;

  }while(nextoct!=NULL);

  cpu->locNoct[level-1] = nOct;

//  printf("nOct = %d on level %d by cpu %d\n",nOct, level, cpu->rank);
}



////////////////////////////////////////////////////////////////////////////////


void cell2lcell(struct CELL *cell, struct LCELL *lcell){

  lcell->marked=cell->marked;
  lcell->child=(cell->child!=NULL);
#ifdef PIC
  lcell->density=cell->density;
#endif

#ifdef WGRAV
  lcell->den=cell->gdata.d;
  lcell->pot=cell->gdata.p;
  lcell->res=cell->res;

  lcell->f[0]=cell->f[0];
  lcell->f[1]=cell->f[1];
  lcell->f[2]=cell->f[2];
#endif

#ifdef WHYDRO2
  lcell->d=cell->field.d;
  lcell->u=cell->field.u;
  lcell->v=cell->field.v;
  lcell->w=cell->field.w;
  lcell->p=cell->field.p;
#ifdef WRADHYD
  lcell->dX=cell->field.dX;
#ifdef HELIUM
  lcell->dXHE=cell->field.dXHE;
  lcell->dXXHE=cell->field.dXXHE;
#endif
#endif
#endif

#ifdef WRAD
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    lcell->e[igrp]=cell->rfield.e[igrp];
    lcell->fx[igrp]=cell->rfield.fx[igrp];
    lcell->fy[igrp]=cell->rfield.fy[igrp];
    lcell->fz[igrp]=cell->rfield.fz[igrp];
    lcell->src[igrp]=cell->rfield.src[igrp];
  }
#ifdef SUPERNOVAE
  lcell->snfb=cell->rfield.snfb;
#endif
  lcell->xion=cell->rfield.nhplus/cell->rfield.nh;
#ifdef HELIUM
  lcell->xHE=cell->rfield.nheplus/cell->rfield.nh;
  lcell->xxHE=cell->rfield.nhepplus/cell->rfield.nh;
#endif
  lcell->temp=cell->rfield.temp;
#endif

}

void oct2loct(struct OCT *oct, struct LOCT *loct){
  int icell;

  for(icell=0;icell<8;icell++){
    cell2lcell(&oct->cell[icell],&loct->cell[icell]);
    //    memcpy(&loct->cell[icell],&oct->cell[icell],sizeof(struct CELL));
  }

  loct->x=oct->x;
  loct->y=oct->y;
  loct->z=oct->z;

  loct->cpu=oct->cpu;
  loct->level=oct->level;
}

