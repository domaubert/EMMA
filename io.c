#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "prototypes.h"
#include "friedmann.h"
#include "oct.h"
#include "parameters.h"
#include "restart.h"


//====================================================================================================

void dumpgrid(int levelmax,struct OCT **firstoct, char filename[],REAL tsim, struct RUNPARAMS *param)
{

  int icur,ii,jj,kk;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  struct LOCT loct;
  FILE *fp =NULL;
  int noct=0;

  fp=fopen(filename,"wb");
	if(fp == NULL) printf("Cannot open %s\n", filename);
  //printf("tsim=%f\n",tsim);
	fwrite(&tsim,sizeof(REAL),1,fp);

#ifdef WRAD
  fwrite(&(param->unit.unit_l),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_v),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_t),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_n),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_d),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_N),sizeof(REAL),1,fp);
#endif

  fwrite(&(firstoct[0]),sizeof(struct OCT*),1,fp);


  for(level=param->lcoarse;level<=levelmax;level++) // looping over octs
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];

      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  nextoct=oct.next;

	  oct2loct(&oct,&loct);
	  fwrite(&loct,sizeof(struct LOCT),1,fp);
	  noct++;
/* 	  fwrite(&oct,sizeof(struct OCT),1,fp); */
	}while(nextoct!=NULL);
    }

  //printf("noct=%d\n",noct);
  fwrite(&noct,sizeof(int),1,fp);

  fclose(fp);
}

#ifdef PIC
void dumppart(struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu){

  FILE *fp = NULL;
  float val;
  int vali;
  int ipart=0;
  int level;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  struct PART *nexp;
  struct PART *curp;
  int icell;
  float tsimf=tsim;

  int npart=0;

#ifdef STARS
  int nstar=0;

  char filenamestar[128];							char filenamepart[128];
  FILE *fstar;												FILE *fpart;
#ifdef MULTIFOLDER
  sprintf(filenamestar,"data/%05d/star/star.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
	sprintf(filenamepart,"data/%05d/part/part.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
#else
  sprintf(filenamestar,"data/star.%05d.p%05d",*(cpu->ndumps),cpu->rank);
	sprintf(filenamepart,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER

  fstar=fopen(filenamestar,"wb");						fpart=fopen(filenamepart,"wb");

	if(fstar == NULL) printf("Cannot open %s\n", filenamestar);
	if(fpart == NULL) printf("Cannot open %s\n", filenamepart);

  fwrite(&nstar,1,sizeof(int)  ,fstar);						fwrite(&npart,1,sizeof(int)  ,fpart);
  fwrite(&tsimf,1,sizeof(float),fstar);						fwrite(&tsimf,1,sizeof(float),fpart);



#else
  fp=fopen(filename,"wb");
	if(fp == NULL) printf("Cannot open %s\n", filename);
  fwrite(&npart,1,sizeof(int)  ,fp);
  fwrite(&tsimf,1,sizeof(float),fp);
#endif

  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];
      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
//	  dxcur=1./POW(2,oct.level);
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
		do{
		  curp=nexp;
		  nexp=curp->next;

#ifdef STARS
		  if(curp->isStar) 	{	fp=fstar;	nstar++;	}
		  else 			{	fp=fpart;	npart++;	}
#else
		  npart++;
#endif
		  val=curp->x;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->y;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->z;			fwrite(&val,1,sizeof(float),fp);
#ifndef PARTN
#ifdef PART_EGY
		  val=curp->ekin+curp->epot;	fwrite(&val,1,sizeof(float),fp);
		  val=curp->fx;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;			fwrite(&val,1,sizeof(float),fp);
#else
		  val=curp->vx;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->vy;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->vz;			fwrite(&val,1,sizeof(float),fp);
#endif
#else
		  val=curp->fx;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->fz;			fwrite(&val,1,sizeof(float),fp);
#endif
		  val=(float)(curp->idx);	fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->mass);	fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->epot);	fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->ekin);	fwrite(&val,1,sizeof(float),fp);
#ifdef STARS
		  if(curp->isStar) {
		    val = curp->age;		fwrite(&val,1,sizeof(float),fp);
		  }
#endif
		  ipart++;

		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }

#ifdef STARS
  rewind(fpart);	fwrite(&npart,1,sizeof(int)  ,fpart);	fclose(fpart);
  rewind(fstar);	fwrite(&nstar,1,sizeof(int)  ,fstar);	fclose(fstar);
#else
  rewind(fp);		fwrite(&npart,1,sizeof(int)  ,fp);	fclose(fp);
#endif

  //printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}
#endif // PIC

//=======================================================================================================

void makeFolders(struct CPUINFO *cpu){
/**
  * create output subfolders for the current dump
  */
  char folder_step[128];
  char folder_field[128];

  sprintf(folder_step,"data/%05d/",*(cpu->ndumps));
  mkdir(folder_step, 0755);

  sprintf(folder_field,"%sgrid/",folder_step);
  mkdir(folder_field, 0755);
  copy_param(folder_field);

#ifdef PIC
  sprintf(folder_field,"%spart/",folder_step);
  mkdir(folder_field, 0755);
  copy_param(folder_field);
#endif // PIC

#ifdef STARS
  sprintf(folder_field,"%sstar/",folder_step);
  mkdir(folder_field, 0755);
  copy_param(folder_field);
#endif // STARS

}

void dumpIO(REAL tsim, struct RUNPARAMS *param,struct CPUINFO *cpu, struct OCT **firstoct, REAL *adt, int pdump){

  REAL tdump,adump;
  char filename[128];

  int idir=cpu->rank%8;

#ifndef TESTCOSMO
#ifdef WRAD
	tdump=(tsim)*param->unit.unit_t/MYR;
#else
	tdump=(tsim);
#endif
	adump=tdump;
#else
	tdump=interp_aexp(tsim,(double*)param->cosmo->tab_aexp,(double*)param->cosmo->tab_ttilde);
	adump=tdump;
#endif

#ifdef MULTIFOLDER
  char folder_step[128];
  char folder_field[128];
  sprintf(folder_step,"data/%05d/",*(cpu->ndumps));
  makeFolders(cpu);
#endif

	if(pdump){
	  // === particle dump
#ifdef PIC

#ifdef MULTIFOLDER
    sprintf(folder_field,"%spart/",folder_step);
    sprintf(filename,"%spart.%05d.p%05d",folder_field,*(cpu->ndumps),cpu->rank);
#else
	  sprintf(filename,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif

	  if(cpu->rank==RANK_DISP){
	    printf("Dumping .......");
	    printf("%s %p\n",filename,cpu->part);
	  }

	  dumppart(firstoct,filename,param->lcoarse,param->lmax,adump,cpu);

#endif
	}
	else{
	  // === Hydro dump
#ifdef MULTIFOLDER
    sprintf(folder_field,"%sgrid/",folder_step);
    sprintf(filename,"%sgrid.%05d.p%05d",folder_field,*(cpu->ndumps),cpu->rank);
#else
	  sprintf(filename,"data/grid.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif

	  if(cpu->rank==RANK_DISP){
	    printf("Dumping .......");
	    printf("%s\n",filename);
	  }
	  dumpgrid(param->lmax,firstoct,filename,adump,param);

/// backups for restart
#ifdef BKP
	  if(*(cpu->ndumps)%FBKP==0){

	    if(cpu->rank==RANK_DISP){
	      printf("BACKUP .......#%d\n",*cpu->ndumps%2+1);
	    }
      char folder_bkp[128];
      sprintf(folder_bkp,"data/bkp/");
      mkdir(folder_bkp, 0755);

	    sprintf(filename,"data/bkp/grid.%05d.p%05d",*(cpu->ndumps)%2+1,cpu->rank);
	    save_amr(filename,firstoct,tdump,cpu->tinit,cpu->nsteps,*(cpu->ndumps),param,cpu,cpu->firstpart,adt);

#ifdef PIC
	    // backups for restart
	    sprintf(filename,"data/bkp/part.%05d.p%05d",*(cpu->ndumps)%2+1,cpu->rank);
	    save_part(filename,firstoct,param->lcoarse,param->lmax,tdump,cpu,cpu->firstpart);
#endif // PIC

	  }
#endif // BKP
	}
}
