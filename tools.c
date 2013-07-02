#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hilbert.h"
#include "prototypes.h"

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

//==================================================================


 //------------------------------------------------------------------------
REAL multicheck(struct OCT **firstoct,int npart,int levelcoarse, int levelmax, int rank, struct CPUINFO *cpu){

  int ntot;
  REAL ntotd;
  REAL nlevd;
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dx;
  int nlev,noct;
  int icell;
  struct PART *nexp;
  struct PART *curp;
  REAL mtot;
  
  REAL xc,yc,zc;
  int *vnoct=cpu->noct;
  int *vnpart=cpu->npart;

  //if(rank==0) printf("Check\n");
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
      if(nextoct==NULL){
	vnoct[level-1]=noct;
	continue;
      }
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu!=rank) continue;
	  if(level>=levelcoarse)
	    {
	      for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
	      
	      xc=curoct->x+(icell%2)*dx+dx*0.5;
	      yc=curoct->y+((icell/2)%2)*dx+dx*0.5;
	      zc=curoct->z+(icell/4)*dx+dx*0.5;
#ifdef PIC
	      ntotd+=curoct->cell[icell].density*dx*dx*dx;
	      nlevd+=curoct->cell[icell].density*dx*dx*dx;
	      
	      nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell

	      if((curoct->cell[icell].child!=NULL)&&(curoct->cell[icell].phead!=NULL)){
		printf("check: split cell with particles !\n");
		printf("curoct->cpu = %d curoct->level=%d\n",curoct->cpu,curoct->level);
		abort();
	      }

	      if(nexp==NULL)continue;
	      do{ 
		nlev++;
		ntot++;
		curp=nexp;
		nexp=curp->next;
		
		if((fabs(curp->x-xc)>0.5*dx)*(fabs(curp->y-yc)>0.5*dx)*(fabs(curp->z-zc)>0.5*dx)){
		  printf("proc %d particle %d not in cell %d: abort\n",cpu->rank,curp->idx,icell);
		  printf("xp=%e xc=%e yp=%e yc=%e zp=%e zc=%e \n DX=%e DY=%e DZ=%e dx=%e\n",curp->x,xc,curp->y,yc,curp->z,zc,curp->x-xc,curp->y-yc,curp->z-zc,0.5*dx);
		  abort();
		}
	      }while(nexp!=NULL);
#endif
	      }
	    }
	  noct++;
	}while(nextoct!=NULL);
      //if((noct!=0)&&(level>=levelcoarse)) printf("proc %d level=%d npart=%d npartd=%f noct=%d\n",cpu->rank,level,nlev,nlevd,noct);
      
      if(level==levelcoarse) mtot=nlevd;
      vnoct[level-1]=noct;
      vnpart[level-1]=nlev;
    }
  

  return mtot;
}
 //------------------------------------------------------------------------
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

//============================================================================
//============================================================================


void grid_census(struct RUNPARAMS *param, struct CPUINFO *cpu){

  int level;
  int ltot,gtot=0,nomax,nomin;
  int lpart,ptot=0;

  if(cpu->rank==0){
    printf("===============================================================\n");
  }
  for(level=2;level<=param->lmax;level++){
    ltot=cpu->noct[level-1];
    lpart=cpu->npart[level-1];
    nomax=ltot;
    nomin=ltot;
    gtot+=ltot;
    ptot+=lpart;
#ifdef WMPI
    MPI_Allreduce(&ltot,&nomax,1,MPI_INT,MPI_MAX,cpu->comm);
    MPI_Allreduce(&ltot,&nomin,1,MPI_INT,MPI_MIN,cpu->comm);
    MPI_Allreduce(MPI_IN_PLACE,&ltot,1,MPI_INT,MPI_SUM,cpu->comm);
    MPI_Allreduce(MPI_IN_PLACE,&lpart,1,MPI_INT,MPI_SUM,cpu->comm);
#endif
    if(cpu->rank==0){
      if(ltot!=0) {printf("level=%2d noct=%9d min=%9d max=%9d npart=%9d",level,ltot,nomin,nomax,lpart);
	int I;
	REAL frac=(ltot/(1.0*pow(2,3*(level-1))))*100.;
	printf("[",frac);
	for(I=0;I<12;I++) printf("%c",(I/12.*100<frac?'*':' '));
	printf("]\n");
      }
    }
  }
#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&gtot,1,MPI_INT,MPI_MAX,cpu->comm);
#endif
  if(cpu->rank==0) printf("\n");

  if(cpu->rank==0){
    int I;
    REAL frac=(gtot/(1.0*param->ngridmax))*100.;
    printf("grid occupancy=%6.1f [",frac);
    for(I=0;I<24;I++) printf("%c",(I/24.*100<frac?'*':' '));
    printf("]\n");
  }

#ifdef PIC
  if(cpu->rank==0){
    int I;
    REAL frac=(ptot/(1.0*param->npartmax))*100.;
    printf("part occupancy=%6.1f [",frac);
    for(I=0;I<24;I++) printf("%c",(I/24.*100<frac?'*':' '));
    printf("]\n");
  }
#endif

  if(cpu->rank==0){
    int I;
    REAL frac=(gtot/(1.0*pow(2,(param->lmax-1)*3)))*100.;
    printf("comp factor   =%6.1f [",frac);
    for(I=0;I<24;I++) printf("%c",(I/24.*100<frac?'*':' '));
    printf("]\n\n");
  }


}

