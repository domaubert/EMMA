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
float  multicheck(struct OCT **firstoct,int npart,int levelcoarse, int levelmax, int rank, int *vnoct){

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
  float mtot;
  
  float xc,yc,zc;


  if(rank==0) printf("Check\n");
  ntot=0.;
  ntotd=0.;
  nlevd=0.;

  for(level=levelcoarse;level<=levelmax;level++)
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
	  if(curoct->cpu!=rank) continue;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      ntotd+=curoct->cell[icell].density*dx*dx*dx;
	      nlevd+=curoct->cell[icell].density*dx*dx*dx;
	      
	      xc=curoct->x+(icell%2)*dx+dx*0.5;
	      yc=curoct->y+((icell/2)%2)*dx+dx*0.5;
	      zc=curoct->z+(icell/4)*dx+dx*0.5;

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
		
		/* if(curp->idx==96571){ */
		/*   printf("mcheck %f %f %f lev=%d\n",curp->x,xc,(curp->x-xc)*2./dx,level); */
		/* } */

		if((fabs(curp->x-xc)>0.5*dx)*(fabs(curp->y-yc)>0.5*dx)*(fabs(curp->z-zc)>0.5*dx)){
		  printf("particle not in cell: abort\n");
		  printf("xp=%e xc=%e yp=%e yc=%e zp=%e zc=%e\n",curp->x,xc,curp->y,yc,curp->z,zc);
		  abort();
		}

	      }while(nexp!=NULL);
	    }
	  noct++;
	}while(nextoct!=NULL);
      if(noct!=0) printf("level=%d npart=%d npartd=%f noct=%d\n",level,nlev,nlevd,noct);
      if(level==levelcoarse) mtot=nlevd;
      vnoct[level-1]=noct;
    }
  
  //if(rank==0) printf("CHECK==> RANK # %d total   npart=%d/%d npartd=%f\n",rank,ntot,npart,ntotd);
  if(ntot!=npart) {
    printf("particles number discrepancy ntot=%d npart=%d\n",ntot,npart);
    abort();
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
