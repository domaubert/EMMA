#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "prototypes.h"
#include "friedmann.h"

#include "restart.h"

#include "oct.h"
#include "ic.h"
#include "io.h"
#include "parameters.h"

#include "segment.h"
#include "stars.h"
#include "hydro_utils.h"
#include "tools.h"

#ifdef WMPI
#include <mpi.h>
#endif

#ifdef HDF5
#include <hdf5.h>
#endif // HDF5

#define AVGFACT (1.) // set to 0 to get an homogenous cosmo field 1

float assign_grid_field(int field,struct CELL *cell){

/**
  * This function return the appropriate field, depending of a given ID.
  * see also parameters.c for the link between the ID and the input given in param.output.
  */

  float res;

  switch(field){
  case 0:
    res= cell2oct(cell)->x+( (cell->idx)   &1)*POW(0.5,cell2oct(cell)->level);
    break;

  case 1:
    res= cell2oct(cell)->y+( (cell->idx>>1)& 1 )*POW(0.5,cell2oct(cell)->level);
    break;

  case 2:
    res= cell2oct(cell)->z+( (cell->idx)>>2    )*POW(0.5,cell2oct(cell)->level);
    break;

  case 3:
    res= cell2oct(cell)->level;
    break;

  case 4:
    res=cell2oct(cell)->cpu;
    break;

#ifdef WGRAV
  case 5:
    res=cell->gdata.d;
    break;
#ifdef PIC
  case 6:
    res=cell->density;
    break;
#endif // PIC
  case 7:
    res=cell->gdata.p;
    break;
  case 8:
    res=cell->res;
    break;
  case 9:
    res=cell->f[0];
    break;
  case 10:
    res=cell->f[1];
    break;
  case 11:
    res=cell->f[2];
    break;
#endif // WGRAV
  case 12:
    res=cell->marked;
    break;
#ifdef WHYDRO2
  case 13:
    res=cell->field.d;
    break;
  case 14:
    res=cell->field.u;
    break;
  case 15:
    res=cell->field.v;
    break;
  case 16:
    res=cell->field.w;
    break;
  case 17:
    res=cell->field.p;
    break;
#endif // WHYDRO2

#ifdef WRAD
  case 18:
    res=cell->rfield.e[0];
    break;
  case 19:
    res=cell->rfield.fx[0];
    break;
  case 20:
    res=cell->rfield.fy[0];
    break;
  case 21:
    res=cell->rfield.fz[0];
    break;
  case 22:
    res=cell->rfield.e[1];
    break;
  case 23:
    res=cell->rfield.fx[1];
    break;
  case 24:
    res=cell->rfield.fy[1];
    break;
  case 25:
    res=cell->rfield.fz[1];
    break;
  case 26:
#ifdef SUPERNOVAE
//    res=cell->rfield.snfb;
#endif // SUPERNOVAE
    break;


  case 27:
    res=cell->rfield.e[2];
    break;
  case 28:
    res=cell->rfield.fx[2];
    break;
  case 29:
    res=cell->rfield.fy[2];
    break;
  case 30:
    res=cell->rfield.fz[2];
    break;

  case 31:
    res=0;
    int igrp;
    for(igrp=0;igrp<NGRP;igrp++) res+=cell->rfield.src[igrp];
    break;
#ifdef WCHEM
  case 32:
    res=cell->rfield.nhplus/cell->rfield.nh;
    break;
#ifdef WRADHYD
  case 33:
    res=cell->field.dX;
    break;
#ifdef HELIUM
  case 34:
    res=cell->field.dXHE;
    break;
  case 35:
    res=cell->field.dXXHE;
    break;
  case 36:
    res=cell->field.xHE;
    break;
  case 37:
    res=cell->field.xxHE;
    break;
#endif // HELIUM
#endif // WRADHYD
  case 38:
    res=cell->rfield.temp;
    break;
#endif // WCHEM
  case 39:
    res=cell->t_first_xion;
    break;
  case 40:
    res=cell->t_last_xion;
    break;
  case 41:
    res=cell->rho_first_xion;
    break;
#endif // WRAD
  }

  return res;
}

float assign_part_field(int field,struct PART *curp){

/**
  * return the appropriate particle field, depending of a given ID.
  */

  float res=0;
  switch(field){
    case 0:
      res=(float)curp->x;
      break;
    case 1:
      res=(float)curp->y;
      break;
    case 2:
		  res=(float)curp->z;
		  break;

    case 3:
		  res=(float)curp->vx;
		  break;
    case 4:
		  res=(float)curp->vy;
		  break;
    case 5:
		  res=(float)curp->vz;
		  break;

#ifndef TESTCOSMO
    case 6:
		  res=(float)curp->fx;
		  break;
    case 7:
		  res=(float)curp->fy;
		  break;
    case 8:
		  res=(float)curp->fz;
		  break;
#endif // TESTCOSMO

    case 9:
		  res=(float)(curp->idx);
		  break;
#ifdef STARS
    case 10:
		  res=(float)(curp->isStar);
		  break;
#endif // STARS

    case 11:
		  res=(float)(curp->epot);
		  break;
    case 12:
		  res=(float)(curp->ekin);
		  break;

#ifdef STARS
    case 13:
		  res=(float)(curp->mass);
		  break;

    case 14:
      if(curp->isStar) {
		    res=(float)(curp->age);
		  }
    break;
#endif // STARS
  }
  return res;
}

#if defined(MPIIO) || defined(HDF5)
void set_offset(struct RUNPARAMS *param, struct CPUINFO *cpu){
/**
  * This function compute the MPIIO offset
  * It cover the amr tree and count the number of cells in each processor domain.
  * Then it set the offset corresponding to the local processor rank.
  */

  const int debug=0;

  // init n_cell to zero
  int ncells=0;
#ifdef PIC
  int nparts=0;
#endif // PIC
#ifdef STARS
  int nstars=0;
#endif // STARS

  // count the cells, parts and stars
  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){

    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *oct=cpu->octList[level-1][iOct];

      int icell;
      for(icell=0;icell<8;icell++){
        struct CELL * cell= &oct->cell[icell];

        if(((oct->cell[icell].child==0)||(oct->level==param->lmax))){

          ncells++;

#ifdef PIC
          struct PART *nexp=cell->phead;
          if(nexp!=NULL){
            do{
              struct PART *curp=nexp;
              nexp=curp->next;

#ifdef STARS
              if (curp->isStar) nstars++;
              else              nparts++;
#else
                                nparts++;
#endif // STARS
              }while(nexp!=NULL);
          }
#endif // PIC
        }
      }
    }
  }

  //if (debug) printf("br cpu=%d ncells=%d\n",cpu->rank, ncells);
  //if (debug) printf("br cpu=%d nparts=%d\n",cpu->rank, nparts);
  //if (debug) printf("br cpu=%d nstars=%d\n",cpu->rank, nstars);

  // broadcast the result

 // MPI_Barrier(cpu->comm);

  MPI_Allgather(&ncells,1,MPI_INT, cpu->mpiio_ncells,1,MPI_INT, cpu->comm);
#ifdef PIC
  MPI_Allgather(&nparts,1,MPI_INT, cpu->mpiio_nparts,1,MPI_INT, cpu->comm);
#endif // PIC
#ifdef STARS
  MPI_Allgather(&nstars,1,MPI_INT, cpu->mpiio_nstars,1,MPI_INT, cpu->comm);
#endif // STARS

  MPI_Barrier(cpu->comm);

  //if (debug) printf("ar cpu=%d ncells=%d\n",cpu->rank, cpu->mpiio_ncells[cpu->rank]);
  //if (debug) printf("ar cpu=%d nparts=%d\n",cpu->rank, cpu->mpiio_nparts[cpu->rank]);
  //if (debug) printf("ar cpu=%d nstars=%d\n",cpu->rank, cpu->mpiio_nstars[cpu->rank]);

  // compute the offset
  cpu->mpiio_grid_offsets=0;
#ifdef PIC
    cpu->mpiio_part_offsets=0;
#endif // PIC
#ifdef STARS
    cpu->mpiio_star_offsets=0;
#endif // STARS

  int i;
  for (i=0;i<cpu->rank;i++){
    cpu->mpiio_grid_offsets += cpu->mpiio_ncells[i];
#ifdef PIC
    cpu->mpiio_part_offsets += cpu->mpiio_nparts[i];
#endif // PIC
#ifdef STARS
    cpu->mpiio_star_offsets += cpu->mpiio_nstars[i];
#endif // STARS
  }

  //if (debug) printf("cpu=%d grid_offset=%lu\n",cpu->rank, cpu->mpiio_grid_offsets);
  //if (debug) printf("cpu=%d part_offset=%lu\n",cpu->rank, cpu->mpiio_part_offsets);
  //if (debug) printf("cpu=%d star_offset=%lu\n",cpu->rank, cpu->mpiio_star_offsets);

}
#endif // MPIIO

#ifdef PIC
void dumppart_serial(struct RUNPARAMS *param, struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu){

  const int debug =0;

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
  int nstar=0;
  int first=-1;

  float part_xmin=2;
  float part_xmax=-1;
  float part_ymin=2;
  float part_ymax=-1;
  float part_zmin=2;
  float part_zmax=-1;

  float star_xmin=2;
  float star_xmax=-1;
  float star_ymin=2;
  float star_ymax=-1;
  float star_zmin=2;
  float star_zmax=-1;


  FILE **f_part=(FILE **)malloc(param->out_part->n_field*sizeof(FILE *));
#ifdef STARS
  FILE **f_star=(FILE **)malloc(param->out_part->n_field*sizeof(FILE *));
#endif // STARS


  int n_field=0;
  int i;
  for(i=0;i<param->out_part->n_field_tot;i++){
    if(param->out_part->field_id[i]){

      if(first==-1) first=i; // looking for the first non nil field

#ifdef STARS
      char filenamestar[128];
      char filenamepart[128];
#ifdef MULTIFOLDER
      char folder_field_star[128];
      sprintf(folder_field_star,"data/%05d/star_%s/",*(cpu->ndumps),param->out_part->field_name[i]);
      mkdir(folder_field_star, 0755);
      sprintf(filenamestar,"%s%s.%05d.p%05d",folder_field_star,param->out_part->field_name[i],*(cpu->ndumps),cpu->rank);

      char folder_field_part[128];
      sprintf(folder_field_part,"data/%05d/part_%s/",*(cpu->ndumps),param->out_part->field_name[i]);
      mkdir(folder_field_part, 0755);
      sprintf(filenamepart,"%s%s.%05d.p%05d",folder_field_part,param->out_part->field_name[i],*(cpu->ndumps),cpu->rank);
#else
      sprintf(filenamestar,"data/star.%s.%05d.p%05d",param->out_part->field_name[i],*(cpu->ndumps),cpu->rank);
      sprintf(filenamepart,"data/part.%s.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER


      if(debug) printf("openning %s at %p\n",filenamestar, f_star[n_field]);
      f_star[n_field]=fopen(filenamestar,"wb");
      if(f_star[n_field] == NULL) {
        printf("Cannot open %s\n", filenamestar);
        abort();
      }
      fwrite(&nstar,1,sizeof(int)  ,f_star[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_star[n_field]);
      fwrite(&star_xmin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_xmax,sizeof(float),1,f_star[n_field]);
      fwrite(&star_ymin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_ymax,sizeof(float),1,f_star[n_field]);
      fwrite(&star_zmin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_zmax,sizeof(float),1,f_star[n_field]);

      if(debug) printf("openning %s at %p\n",filenamepart,f_part[n_field]);
      f_part[n_field]=fopen(filenamepart,"wb");
      if(f_part[n_field] == NULL){
        printf("Cannot open %s\n", filenamepart);
        abort();
      }

      fwrite(&npart,1,sizeof(int)  ,f_part[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_part[n_field]);
      fwrite(&part_xmin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_xmax,sizeof(float),1,f_part[n_field]);
      fwrite(&part_ymin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_ymax,sizeof(float),1,f_part[n_field]);
      fwrite(&part_zmin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_zmax,sizeof(float),1,f_part[n_field]);

#else
      // NOSTARS CASE
      char filenamepart[128];
#ifdef MULTIFOLDER
      char folder_field_part[128];
      sprintf(folder_field_part,"data/%05d/part_%s/",*(cpu->ndumps),param->out_part->field_name[i]);
      mkdir(folder_field_part, 0755);
      sprintf(filenamepart,"%s%s.%05d.p%05d",folder_field_part,param->out_part->field_name[i],*(cpu->ndumps),cpu->rank);
#else
      sprintf(filenamepart,"data/part.%s.%05d.p%05d",param->out_part->field_name[i],*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER

      if(debug) printf("openning %s at %p\n",filenamepart,f_part[n_field]);

      f_part[n_field]=fopen(filenamepart,"wb");
      if(f_part[n_field] == NULL){
        printf("Cannot open %s\n", filenamepart);
        abort();
      }

      fwrite(&npart,1,sizeof(int)  ,f_part[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_part[n_field]);

#endif // STARS
      n_field++;
    }
  }

  if(debug) printf("opening file OK\n");

  for(level=levelcoarse;level<=levelmax;level++){
    if (debug) printf("entering level %d\n",level);

    nextoct=firstoct[level-1];
    do{
      if(nextoct==NULL) continue; // in case the level is empty
      oct=(*nextoct);
      nextoct=oct.next;

      for(icell=0;icell<8;icell++){ // looping over cells in oct

        nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
        if(nexp!=NULL){
          do{
            curp=nexp;
            nexp=curp->next;

            int ii=0;
            for (i=0;i<param->out_part->n_field_tot; i++){
              if(param->out_part->field_id[i]){

            //    if(debug) printf("field_id=%d\n",param->out_part->field_id[i]);

#ifdef STARS
                if(curp->isStar) 	{	fp=f_star[ii];	if(i==first) nstar++;	}
                else 			        {	fp=f_part[ii];	if(i==first) npart++;	}
#else
                fp=f_part[ii]; if(i==first) npart++;
#endif // STARS

                float dat = (float)assign_part_field(i,curp);
                fwrite(&dat,sizeof(float),1,fp);
                ii++;
              }
            }

      // update file boundaries
#ifdef STARS
      if(curp->isStar) 	{
        if(curp->x<star_xmin) star_xmin=curp->x;
        if(curp->y<star_ymin) star_ymin=curp->y;
        if(curp->z<star_zmin) star_zmin=curp->z;
        if(curp->x>star_xmax) star_xmax=curp->x;
        if(curp->y>star_ymax) star_ymax=curp->y;
        if(curp->z>star_zmax) star_zmax=curp->z;
      }else
#endif // STARS
      {
        if(curp->x<part_xmin) part_xmin=curp->x;
        if(curp->y<part_ymin) part_ymin=curp->y;
        if(curp->z<part_zmin) part_zmin=curp->z;
        if(curp->x>part_xmax) part_xmax=curp->x;
        if(curp->y>part_ymax) part_ymax=curp->y;
        if(curp->z>part_zmax) part_zmax=curp->z;
      }

		  ipart++;

		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }

  if (debug) printf("writing OK\n");

  n_field=0;
  for(i=0;i<param->out_part->n_field_tot;i++){
    if(param->out_part->field_id[i]){

      rewind(f_part[n_field]);
      fwrite(&npart,1,sizeof(int)  ,f_part[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_part[n_field]);
      fwrite(&part_xmin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_xmax,sizeof(float),1,f_part[n_field]);
      fwrite(&part_ymin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_ymax,sizeof(float),1,f_part[n_field]);
      fwrite(&part_zmin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_zmax,sizeof(float),1,f_part[n_field]);
      fclose(f_part[n_field]);
#ifdef STARS
      rewind(f_star[n_field]);
      fwrite(&nstar,1,sizeof(int)  ,f_star[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_star[n_field]);
      fwrite(&star_xmin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_xmax,sizeof(float),1,f_star[n_field]);
      fwrite(&star_ymin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_ymax,sizeof(float),1,f_star[n_field]);
      fwrite(&star_zmin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_zmax,sizeof(float),1,f_star[n_field]);
      fclose(f_star[n_field]);
#endif // STARS

      n_field++;
    }
  }

  if (debug) printf("closing  OK\n");
  //printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);

  free(f_part);
#ifdef STARS
  free(f_star);
#endif // STARS
}

void dumppart_serial_bkp(struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu){

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

  char filenamestar[128];
  char filenamepart[128];

#ifdef MULTIFOLDER
  sprintf(filenamestar,"data/%05d/star/star.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
  sprintf(filenamepart,"data/%05d/part/part.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
#else
  sprintf(filenamestar,"data/star.%05d.p%05d",*(cpu->ndumps),cpu->rank);
  sprintf(filenamepart,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER

  FILE *fstar=fopen(filenamestar,"wb");
  if(fstar == NULL) {
    printf("Cannot open %s\n", filenamestar);
    abort();
  }
  fwrite(&nstar,1,sizeof(int)  ,fstar);
  fwrite(&tsimf,1,sizeof(float),fstar);

  FILE *fpart=fopen(filenamepart,"wb");
  if(fpart == NULL){
    printf("Cannot open %s\n", filenamepart);
    abort();
  }
  fwrite(&npart,1,sizeof(int)  ,fpart);
  fwrite(&tsimf,1,sizeof(float),fpart);

#else

#ifdef MULTIFOLDER
  sprintf(filename,"data/%05d/part/part.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
#else
  sprintf(filename,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER

  fp=fopen(filename,"wb");
  if(fp == NULL){
    printf("Cannot open %s\n", filename);
    abort();
  }
  fwrite(&npart,1,sizeof(int)  ,fp);
  fwrite(&tsimf,1,sizeof(float),fp);
#endif // STARS



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
		  else 			        {	fp=fpart;	npart++;	}
#else
		  npart++;
#endif // STARS
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
#endif // PART_EGY
#else
		  val=curp->fx;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;			fwrite(&val,1,sizeof(float),fp);
		  val=curp->fz;			fwrite(&val,1,sizeof(float),fp);
#endif // PARTN
		  val=(float)(curp->idx);	  fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->mass);	fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->epot);	fwrite(&val,1,sizeof(float),fp);
		  val=(float)(curp->ekin);	fwrite(&val,1,sizeof(float),fp);
#ifdef STARS
		  if(curp->isStar) {
		    val = curp->age;		fwrite(&val,1,sizeof(float),fp);
		  }
#endif // STARS
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
#endif // STARS

  //printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}

#ifdef MPIIO
void dumppart_MPI_bkp(struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct RUNPARAMS *param){

  int ipart=0;
  int istar=0;
  float tsimf=tsim;

//-----------------------------------------------------------------------------------------------

  float *part_tmp = (float*)calloc(cpu->mpiio_nparts[cpu->rank],10*sizeof(float));
  int i_tmp_part = 0;

#ifdef STARS
  float *star_tmp = (float*)calloc(cpu->mpiio_nstars[cpu->rank],11*sizeof(float));
  int i_tmp_star = 0;
#endif // STARS

  int level;
  for(level=levelcoarse;level<=levelmax;level++){
    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *oct=cpu->octList[level-1][iOct];
      int icell;
      for(icell=0;icell<8;icell++){
	      struct PART *nexp=oct->cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
          do{
            struct PART *curp=nexp;
            nexp=curp->next;

      float* tmp=part_tmp;
      int* i_tmp=&i_tmp_part;
#ifdef STARS
		  if(curp->isStar){
        tmp=star_tmp;
        i_tmp=&i_tmp_star;
      }
#endif // STARS

		  tmp[(*i_tmp)++]=(float)curp->x;
      tmp[(*i_tmp)++]=(float)curp->y;
		  tmp[(*i_tmp)++]=(float)curp->z;
#ifndef PARTN
#ifdef PART_EGY
		  tmp[(*i_tmp)++]=(float)curp->ekin+curp->epot;
		  tmp[(*i_tmp)++]=(float)curp->fx;
		  tmp[(*i_tmp)++]=(float)curp->fy;
#else
		  tmp[(*i_tmp)++]=(float)curp->vx;
		  tmp[(*i_tmp)++]=(float)curp->vy;
		  tmp[(*i_tmp)++]=(float)curp->vz;
#endif // PART_EGY
#else
		  tmp[(*i_tmp)++]=(float)curp->fx;
		  tmp[(*i_tmp)++]=(float)curp->fy;
		  tmp[(*i_tmp)++]=(float)curp->fz;
#endif // PARTN
		  tmp[(*i_tmp)++]=(float)(curp->idx);
		  tmp[(*i_tmp)++]=(float)(curp->mass);
		  tmp[(*i_tmp)++]=(float)(curp->epot);
		  tmp[(*i_tmp)++]=(float)(curp->ekin);
#ifdef STARS
		  if(curp->isStar) {
		    tmp[(*i_tmp)++]=(float)(curp->age);
		  }
#endif // STARS
		  ipart++;

          }while(nexp!=NULL);
	      }
	    }
    }
  }

//-----------------------------------------------------------------------------------------------
  int i;
  int npart=0;
  for (i=0;i<cpu->nproc;i++){
   npart+=cpu->mpiio_nparts[i];
  }

  char filenamepart[128];
#ifdef MULTIFOLDER
#ifndef MPIIO
	sprintf(filenamepart,"data/%05d/part/part.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
#else
	sprintf(filenamepart,"data/%05d/part.%05d",*(cpu->ndumps),*(cpu->ndumps));
#endif // MPIIO
#else
	sprintf(filenamepart,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER

  MPI_File fpart;
  MPI_File_open(cpu->comm,filenamepart,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fpart);
	if(fpart==NULL){
    printf("Cannot open %s\n", filenamepart);
    abort();
  }

	MPI_File_write(fpart, &npart,1, MPI_INT, MPI_STATUS_IGNORE);
  MPI_File_write(fpart, &tsimf,1, MPI_FLOAT, MPI_STATUS_IGNORE);
  const size_t part_header_size = sizeof(int)+sizeof(float);

  const size_t part_size = 10;
  MPI_Datatype part_type;
  MPI_Type_contiguous(part_size, MPI_FLOAT, &part_type);
  MPI_Type_commit(&part_type);

  MPI_Offset part_offset = part_size*cpu->mpiio_part_offsets*sizeof(float) + part_header_size ;
  MPI_File_set_view(fpart, part_offset, MPI_FLOAT, part_type, "native", MPI_INFO_NULL);
  MPI_File_write(fpart, part_tmp,i_tmp_part, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_close(&fpart);
  free(part_tmp);

//-----------------------------------------------------------------------------------------------
#ifdef STARS
  int nstar=0;
  for (i=0;i<cpu->nproc;i++){
   nstar+=cpu->mpiio_nstars[i];
  }

  char filenamestar[128];
#ifdef MULTIFOLDER
#ifndef MPIIO
  sprintf(filenamestar,"data/%05d/star/star.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
#else
  sprintf(filenamestar,"data/%05d/star.%05d",*(cpu->ndumps),*(cpu->ndumps));
#endif // MPIIO
#else
  sprintf(filenamestar,"data/star.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER

  MPI_File fstar;
  MPI_File_open(cpu->comm,filenamestar,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fstar);
  if(fstar==NULL){
    printf("Cannot open %s\n", filenamestar);
    abort();
  }
  MPI_File_write(fstar, &nstar,1, MPI_INT, MPI_STATUS_IGNORE);
  MPI_File_write(fstar, &tsimf,1, MPI_FLOAT, MPI_STATUS_IGNORE);
  const size_t star_header_size = sizeof(int)+sizeof(float);

  const size_t star_size = 11;

  MPI_Datatype star_type;
  MPI_Type_contiguous(star_size, MPI_FLOAT, &star_type);
  MPI_Type_commit(&star_type);

  MPI_Offset star_offset = star_size*cpu->mpiio_part_offsets*sizeof(float) + star_header_size ;
  MPI_File_set_view(fstar, star_offset, MPI_FLOAT, star_type, "native", MPI_INFO_NULL);
  MPI_File_write(fstar, star_tmp,i_tmp_star, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_close(&fstar);
  free(star_tmp);
#endif // STARS
//printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}

void dumppart_MPI(struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct RUNPARAMS *param){

  int ipart=0;
  int istar=0;
  float tsimf=tsim;

//-----------------------------------------------------------------------------------------------

  float *part_tmp = (float*)calloc(cpu->mpiio_nparts[cpu->rank],10*sizeof(float));
  int i_tmp_part = 0;

#ifdef STARS
  float *star_tmp = (float*)calloc(cpu->mpiio_nstars[cpu->rank],11*sizeof(float));
  int i_tmp_star = 0;
#endif // STARS

  int level;
  for(level=levelcoarse;level<=levelmax;level++){
    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *oct=cpu->octList[level-1][iOct];
      int icell;
      for(icell=0;icell<8;icell++){
	      struct PART *nexp=oct->cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
          do{
            struct PART *curp=nexp;
            nexp=curp->next;

      float* tmp=part_tmp;
      int* i_tmp=&i_tmp_part;
#ifdef STARS
		  if(curp->isStar){
        tmp=star_tmp;
        i_tmp=&i_tmp_star;
      }
#endif // STARS

		  tmp[(*i_tmp)++]=(float)curp->x;
      tmp[(*i_tmp)++]=(float)curp->y;
		  tmp[(*i_tmp)++]=(float)curp->z;
#ifndef PARTN
#ifdef PART_EGY
		  tmp[(*i_tmp)++]=(float)curp->ekin+curp->epot;
		  tmp[(*i_tmp)++]=(float)curp->fx;
		  tmp[(*i_tmp)++]=(float)curp->fy;
#else
		  tmp[(*i_tmp)++]=(float)curp->vx;
		  tmp[(*i_tmp)++]=(float)curp->vy;
		  tmp[(*i_tmp)++]=(float)curp->vz;
#endif // PART_EGY
#else
		  tmp[(*i_tmp)++]=(float)curp->fx;
		  tmp[(*i_tmp)++]=(float)curp->fy;
		  tmp[(*i_tmp)++]=(float)curp->fz;
#endif // PARTN
		  tmp[(*i_tmp)++]=(float)(curp->idx);
		  tmp[(*i_tmp)++]=(float)(curp->mass);
		  tmp[(*i_tmp)++]=(float)(curp->epot);
		  tmp[(*i_tmp)++]=(float)(curp->ekin);
#ifdef STARS
		  if(curp->isStar) {
		    tmp[(*i_tmp)++]=(float)(curp->age);
		  }
#endif // STARS
		  ipart++;

          }while(nexp!=NULL);
	      }
	    }
    }
  }

//-----------------------------------------------------------------------------------------------
  int i;
  int npart=0;
  for (i=0;i<cpu->nproc;i++){
   npart+=cpu->mpiio_nparts[i];
  }

  char filenamepart[128];
#ifdef MULTIFOLDER
#ifndef MPIIO
	sprintf(filenamepart,"data/%05d/part/part.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
#else
	sprintf(filenamepart,"data/%05d/part.%05d",*(cpu->ndumps),*(cpu->ndumps));
#endif // MPIIO
#else
	sprintf(filenamepart,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER

  MPI_File fpart;
  MPI_File_open(cpu->comm,filenamepart,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fpart);
	if(fpart==NULL){
    printf("Cannot open %s\n", filenamepart);
    abort();
  }

	MPI_File_write(fpart, &npart,1, MPI_INT, MPI_STATUS_IGNORE);
  MPI_File_write(fpart, &tsimf,1, MPI_FLOAT, MPI_STATUS_IGNORE);


  //TODO fix theses values
  int ncell[cpu->nproc];
  float bounds[6*cpu->nproc];

  MPI_File_write(fpart, ncell, cpu->nproc, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_write(fpart, bounds, 6*cpu->nproc, MPI_FLOAT, MPI_STATUS_IGNORE);

  const size_t part_header_size = sizeof(int) +
                                  sizeof(float) +
                                  cpu->nproc *sizeof(int) +
                                  6* cpu->nproc *sizeof(float);


  const size_t part_size = 10;
  MPI_Datatype part_type;
  MPI_Type_contiguous(part_size, MPI_FLOAT, &part_type);
  MPI_Type_commit(&part_type);

  MPI_Offset part_offset = part_size*cpu->mpiio_part_offsets*sizeof(float) + part_header_size ;
  MPI_File_set_view(fpart, part_offset, MPI_FLOAT, part_type, "native", MPI_INFO_NULL);
  MPI_File_write(fpart, part_tmp,i_tmp_part, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_close(&fpart);
  free(part_tmp);

//-----------------------------------------------------------------------------------------------
#ifdef STARS
  int nstar=0;
  for (i=0;i<cpu->nproc;i++){
   nstar+=cpu->mpiio_nstars[i];
  }

  char filenamestar[128];
#ifdef MULTIFOLDER
#ifndef MPIIO
  sprintf(filenamestar,"data/%05d/star/star.%05d.p%05d",*(cpu->ndumps),*(cpu->ndumps),cpu->rank);
#else
  sprintf(filenamestar,"data/%05d/star.%05d",*(cpu->ndumps),*(cpu->ndumps));
#endif // MPIIO
#else
  sprintf(filenamestar,"data/star.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MUTLTIFOLDER

  MPI_File fstar;
  MPI_File_open(cpu->comm,filenamestar,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fstar);
  if(fstar==NULL){
    printf("Cannot open %s\n", filenamestar);
    abort();
  }





  MPI_File_write(fstar, &nstar,1, MPI_INT, MPI_STATUS_IGNORE);
  MPI_File_write(fstar, &tsimf,1, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_write(fpart, ncell, cpu->nproc, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_write(fpart, bounds, 6*cpu->nproc, MPI_FLOAT, MPI_STATUS_IGNORE);

  const size_t star_header_size = sizeof(int) +
                                  sizeof(float) +
                                  cpu->nproc *sizeof(int) +
                                  6* cpu->nproc *sizeof(float);



  const size_t star_size = 11;

  MPI_Datatype star_type;
  MPI_Type_contiguous(star_size, MPI_FLOAT, &star_type);
  MPI_Type_commit(&star_type);

  MPI_Offset star_offset = star_size*cpu->mpiio_part_offsets*sizeof(float) + star_header_size ;
  MPI_File_set_view(fstar, star_offset, MPI_FLOAT, star_type, "native", MPI_INFO_NULL);
  MPI_File_write(fstar, star_tmp,i_tmp_star, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_close(&fstar);
  free(star_tmp);
#endif // STARS
//printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}
#endif // MPIIO
#endif // PIC

void dumpgrid(int levelmax,struct OCT **firstoct, char filename[],REAL tsim, struct RUNPARAMS *param){

  int icur,ii,jj,kk;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  struct LOCT loct;
  FILE *fp =NULL;
  int noct=0;

  fp=fopen(filename,"wb");
  if(fp == NULL) printf("Cannot open %s\n", filename);

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


  //  for(level=param->lcoarse;level<=levelmax;level++) // looping over octs
    for(level=1;level<=levelmax;level++) // looping over octs
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

#ifdef MPIIO
void dumpalloct_MPI(char folder[],REAL tsim, struct RUNPARAMS *param, struct CPUINFO *cpu){

/**
  * This function dump the output data with MPIIO
  * only the most reffined cell are dumped
  *
  * format:
  *     int n : the number of dumped cells
  *     followed by n times:
  *     float field : value of the corresponding field
  *
  *     TODO: for the flux and all the fields depending of NGRP, dump in 1 field instead of NGRP fields
  *
  */
  int i;

  int n_cell_tot=0;
  for (i=0;i<cpu->nproc;i++){
    n_cell_tot+=cpu->mpiio_ncells[i];
  }

  float *tmp = (float*)calloc(cpu->mpiio_ncells[cpu->rank],sizeof(float));

  int n_field=0;
  for (i=0;i<param->out_grid->n_field_tot; i++){
    if(param->out_grid->field_id[i]){

      //reduce data
      int i_tmp=0;
      int level;
      for(level=param->lcoarse;level<=param->lmax;level++){
        int iOct;
        for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
          struct OCT *oct=cpu->octList[level-1][iOct];
          int icell;
          for(icell=0;icell<8;icell++){
            struct CELL * cell= &oct->cell[icell];
            if(((oct->cell[icell].child==0)||(oct->level==param->lmax))){
              tmp[i_tmp++] = (float)assign_grid_field(i,cell);
            }
          }
        }
      }

      //Open the field file
      char dat_name[256];
      sprintf(dat_name,"%s%s.%05d",folder,param->out_grid->field_name[i],*(cpu->ndumps));

      MPI_File f_dat;
      MPI_File_open(cpu->comm,dat_name,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&f_dat);

      if(f_dat == NULL){
       printf("Cannot open %s\n", dat_name);
       abort();
      }

      //TODO fix theses values
      float tsim;
      int ncell[cpu->nproc];
      float bounds[6*cpu->nproc];

      //write header
      MPI_File_write(f_dat, &n_cell_tot,1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(f_dat, &tsim,1, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write(f_dat, ncell, cpu->nproc, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write(f_dat, bounds, 6*cpu->nproc, MPI_FLOAT, MPI_STATUS_IGNORE);

      const size_t grid_header_size = sizeof(int) +
                                      sizeof(float) +
                                      cpu->nproc *sizeof(int) +
                                      6* cpu->nproc *sizeof(float);

      //set view
      MPI_Offset grid_offset = cpu->mpiio_grid_offsets*sizeof(float) + grid_header_size ;
      MPI_File_set_view(f_dat, grid_offset, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);

      //write data
      MPI_File_write(f_dat, tmp, cpu->mpiio_ncells[cpu->rank], MPI_FLOAT, MPI_STATUS_IGNORE);

      // close the field file
      MPI_File_close(&f_dat);
      n_field++;
    }
  }
  free(tmp);
}
#endif // MPIIO

void dumpalloct_serial(char folder[],REAL tsim, struct RUNPARAMS *param, struct CPUINFO *cpu){


/**
  * This function dump the output data with STDIO
  * only the most reffined cell are dumped
  *
  * format:
  *   - N_fields time N_procs files corresponding to the defined field in param.output
  *     int n : the number of dumped cells
  *     followed by n times:
  *     float field : value of the corresponding field
  *
  *     TODO: for the flux and all the fields depending of NGRP, dump in 1 field instead of NGRP field
  *
  */

  const int debug=0;

  int i;
  int n_field=0;
  int n_cell=0;
  float xmin,xmax,ymin,ymax,zmin,zmax;

  xmin=2;
  xmax=-1;
  ymin=2;
  ymax=-1;
  zmin=2;
  zmax=-1;

// Opening all the fields files
  if(debug) printf("Allocating %d file pointers \n", param->out_grid->n_field);

  FILE **f_dat;
  f_dat=(FILE **)malloc(param->out_grid->n_field*sizeof(FILE *));


  for(i=0;i<param->out_grid->n_field_tot;i++){
    if(param->out_grid->field_id[i]){

      char folder_field[128];
      sprintf(folder_field,"%sgrid_%s/",folder, param->out_grid->field_name[i]);
      mkdir(folder_field, 0755);
      char dat_name[256];
      sprintf(dat_name,"%s%s.%05d.p%05d",folder_field,param->out_grid->field_name[i],*(cpu->ndumps),cpu->rank);

      if(debug) printf("Openning : %s",dat_name);

      f_dat[n_field]=fopen(dat_name,"wb");
      if(f_dat[n_field] == NULL){
       printf("Cannot open %s\n", dat_name);
       abort();
      }

      fwrite(&n_cell,sizeof(int),1,f_dat[n_field]);
      fwrite(&tsim,sizeof(float),1,f_dat[n_field]);
      fwrite(&xmin,sizeof(float),1,f_dat[n_field]);
      fwrite(&xmax,sizeof(float),1,f_dat[n_field]);
      fwrite(&ymin,sizeof(float),1,f_dat[n_field]);
      fwrite(&ymax,sizeof(float),1,f_dat[n_field]);
      fwrite(&zmin,sizeof(float),1,f_dat[n_field]);
      fwrite(&zmax,sizeof(float),1,f_dat[n_field]);

      n_field++;
    }
  }

  if(debug) printf("Files open, let's write");

// writing the data
  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){
    double dx=pow(0.5,level-1); // oct size

    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *oct=cpu->octList[level-1][iOct];

      int icell;
      for(icell=0;icell<8;icell++){
        struct CELL * cell= &oct->cell[icell];
        if(((oct->cell[icell].child==0)||(oct->level==param->lmax))){
          n_cell++;

          // get and write the fields
          int ii=0;
          for (i =0;i<param->out_grid->n_field_tot; i++){
            if(param->out_grid->field_id[i]){
              float dat = (float)assign_grid_field(i,cell);
              fwrite(&dat,sizeof(float),1,f_dat[ii]);
              ii++;
            }
          }

	  // update file boundaries
	  if(oct->x<xmin) xmin=oct->x;
	  if(oct->y<ymin) ymin=oct->y;
	  if(oct->z<zmin) zmin=oct->z;

	  if(oct->x+dx>xmax) xmax=oct->x+dx;
	  if(oct->y+dx>ymax) ymax=oct->y+dx;
	  if(oct->z+dx>zmax) zmax=oct->z+dx;

        }
      }
    }
  }

  if(debug) printf("Write OK, header update and close");

  // write n_cells and close the fields files
  n_field=0;
  for(i=0;i<param->out_grid->n_field_tot;i++){
    if(param->out_grid->field_id[i]){
      rewind(f_dat[n_field]);
      fwrite(&n_cell,sizeof(int),1,f_dat[n_field]);
      fwrite(&tsim,sizeof(float),1,f_dat[n_field]);
      fwrite(&xmin,sizeof(float),1,f_dat[n_field]);
      fwrite(&xmax,sizeof(float),1,f_dat[n_field]);
      fwrite(&ymin,sizeof(float),1,f_dat[n_field]);
      fwrite(&ymax,sizeof(float),1,f_dat[n_field]);
      fwrite(&zmin,sizeof(float),1,f_dat[n_field]);
      fwrite(&zmax,sizeof(float),1,f_dat[n_field]);
      fclose(f_dat[n_field]);
      n_field++;
    }
  }

  free(f_dat);
}

void makeFolders(struct CPUINFO *cpu){

  char folder_step[128];
  char folder_field[128];

  sprintf(folder_step,"data/%05d/",*(cpu->ndumps));
  mkdir(folder_step, 0755);

#ifndef ALLOCT
  sprintf(folder_field,"%sgrid/",folder_step);
  mkdir(folder_field, 0755);
//  copy_param(folder_field);
#endif // ALLOCT


#ifndef MPIIO
#ifdef PIC

#ifndef ALLOCT
  sprintf(folder_field,"%spart/",folder_step);
  mkdir(folder_field, 0755);
#endif // ALLOCT

//  copy_param(folder_field);
#endif // PIC

#ifdef STARS
#ifndef ALLOCT
  sprintf(folder_field,"%sstar/",folder_step);
  mkdir(folder_field, 0755);
#endif // ALLOCT

//  copy_param(folder_field);
#endif // STARS
#endif // MPIIO
}

void dump_domain(struct RUNPARAMS *param, struct CPUINFO *cpu){
  /**
    * Dump the processor domain limit information
    * Create the ascii file data/domain.dat
    * called once at the start as domains do not change during runtime
    * this function should be revised when implementing load balancing
    */


  // Reduce domains
  float xmin=2,xmax=-1,ymin=2,ymax=-1,zmin=2,zmax=-1;
  int i_tmp=0;
  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){
    REAL dx = POW(0.5,level-1);
    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *oct=cpu->octList[level-1][iOct];
        if(oct->x<xmin) xmin=oct->x;
        if(oct->y<ymin) ymin=oct->y;
        if(oct->z<zmin) zmin=oct->z;
        if(oct->x+dx>xmax) xmax=oct->x+dx;
        if(oct->y+dx>ymax) ymax=oct->y+dx;
        if(oct->z+dx>zmax) zmax=oct->z+dx;
    }
  }

  float* Xmin= (float*)calloc(cpu->nproc,sizeof(float));
  float* Xmax= (float*)calloc(cpu->nproc,sizeof(float));
  float* Ymin= (float*)calloc(cpu->nproc,sizeof(float));
  float* Ymax= (float*)calloc(cpu->nproc,sizeof(float));
  float* Zmin= (float*)calloc(cpu->nproc,sizeof(float));
  float* Zmax= (float*)calloc(cpu->nproc,sizeof(float));

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI

  MPI_Allgather(&xmin,1,MPI_FLOAT,Xmin,1,MPI_FLOAT, cpu->comm);
  MPI_Allgather(&xmax,1,MPI_FLOAT,Xmax,1,MPI_FLOAT, cpu->comm);
  MPI_Allgather(&ymin,1,MPI_FLOAT,Ymin,1,MPI_FLOAT, cpu->comm);
  MPI_Allgather(&ymax,1,MPI_FLOAT,Ymax,1,MPI_FLOAT, cpu->comm);
  MPI_Allgather(&zmin,1,MPI_FLOAT,Zmin,1,MPI_FLOAT, cpu->comm);
  MPI_Allgather(&zmax,1,MPI_FLOAT,Zmax,1,MPI_FLOAT, cpu->comm);

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI

  if(cpu->rank == RANK_DISP){
    char filename[256];
    sprintf(filename,"data/domains.dat");
    FILE *fp=fopen(filename,"wb");
    if(fp == NULL) printf("Cannot open %s\n", filename);

    fprintf(fp,"ncells       \t");
    fprintf(fp,"xmin         \t");
    fprintf(fp,"xmax         \t");
    fprintf(fp,"ymin         \t");
    fprintf(fp,"ymax         \t");
    fprintf(fp,"zmin         \t");
    fprintf(fp,"zmax         \t");
    fprintf(fp,"\n");

    int i_proc;
    for (i_proc=0; i_proc<cpu->nproc; i_proc++){
      fprintf(fp,"%d\t\t",cpu->mpiio_ncells[i_proc]);
      fprintf(fp,"%f\t",Xmin[i_proc]);
      fprintf(fp,"%f\t",Xmax[i_proc]);
      fprintf(fp,"%f\t",Ymin[i_proc]);
      fprintf(fp,"%f\t",Ymax[i_proc]);
      fprintf(fp,"%f\t",Zmin[i_proc]);
      fprintf(fp,"%f\t",Zmax[i_proc]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  free(Xmin);
  free(Xmax);
  free(Ymin);
  free(Ymax);
  free(Zmin);
  free(Zmax);
}


#ifdef HDF5

#define PHDF5
#define HDF5_METHOD H5FD_MPIO_COLLECTIVE
//#define HDF5_METHOD H5FD_MPIO_INDEPENDENT

void dump_HDF5_grid(char folder[],REAL tsim, struct RUNPARAMS *param, struct CPUINFO *cpu){

/**
  * This function dump the output data with HDF5
  */

  const int debug=0;

  // Create step folder
  char folder_step[128];
  sprintf(folder_step,"data/%05d/",*(cpu->ndumps));
  mkdir(folder_step, 0755);

  float *tmp = (float*)calloc(cpu->mpiio_ncells[cpu->rank],sizeof(float));
  if (tmp == NULL){
    printf("can't allocate tmp vector in dump_HDF5_grid\n");
    abort();
  }

  int ifield;
  for (ifield=0;ifield<param->out_grid->n_field_tot; ifield++){
    if(param->out_grid->field_id[ifield]){

      hid_t plist;

      //Set up file access property list with parallel I/O access
      plist = H5Pcreate(H5P_FILE_ACCESS);
#ifdef PHDF5
      H5Pset_fapl_mpio(plist, cpu->comm, MPI_INFO_NULL);
#endif // PHDF5

      //Create a new file collectively
      char file_name[512];
      sprintf(file_name,"data/%05d/grid_%s_%05d.h5", *cpu->ndumps, param->out_grid->field_name[ifield], *cpu->ndumps);

      hid_t file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
      H5Pclose(plist);

      hsize_t n_cell_tot=0;
      int i;
      for (i=0;i<cpu->nproc;i++){
        n_cell_tot+=cpu->mpiio_ncells[i];
      }

      // Create the data space for the dataset.
      hid_t dataspace = H5Screate_simple(1, &n_cell_tot, NULL);

      //Select hyperslab in the file.
      hsize_t offset = cpu->mpiio_grid_offsets;
      hsize_t n_loc = cpu->mpiio_ncells[cpu->rank];
      if (debug) printf("ncell_tot=%llu nloc=%llu \n",n_cell_tot,n_loc);

      H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&offset,NULL,&n_loc,NULL);

      // Create property list
      plist = H5Pcreate(H5P_DATASET_XFER);

#ifdef PHDF5
      H5Pset_dxpl_mpio(plist, HDF5_METHOD);
#endif // PHDF5

      hid_t	memspace = H5Screate_simple (1, &n_loc, NULL);

      //reduce data
      int i_tmp=0;
      int level;
      for(level=param->lcoarse;level<=param->lmax;level++){
        int iOct;
        for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
          struct OCT *oct=cpu->octList[level-1][iOct];
          int icell;
          for(icell=0;icell<8;icell++){
            struct CELL * cell= &oct->cell[icell];
            if(((oct->cell[icell].child==0)||(oct->level==param->lmax))){
              tmp[i_tmp++] = (float)assign_grid_field(ifield,cell);
            }
          }
        }
      }

      char field_name[256];
      //sprintf(field_name,"%s",param->out_grid->field_name[ifield]);
      sprintf(field_name,"data");
      hid_t  dataset = H5Dcreate(file, field_name, H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, plist, tmp);
      H5Dclose(dataset);

      // tsim attribute
      hsize_t dims = 1;
      hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);
      hid_t attribute_id = H5Acreate2 (file, "a", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      float attr_data=(float)tsim;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(1, &dims, NULL);
      attribute_id = H5Acreate2 (file, "z", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      attr_data= 1./(float)tsim -1;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(1, &dims, NULL);
      attribute_id = H5Acreate2 (file, "t", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      attr_data= (float)param->cosmo->tphy;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      H5Pclose(plist);
      H5Fclose(file);
    }
  }
  free(tmp);
}

#ifdef PIC
void dump_HDF5_part(char filename[],REAL tsim,  struct RUNPARAMS *param, struct CPUINFO *cpu){

  const int debug =0;

  // Create step folder
  char folder_step[128];
  sprintf(folder_step,"data/%05d/",*(cpu->ndumps));
  mkdir(folder_step, 0755);

  if (debug) printf("allocating %d float\n",cpu->mpiio_nparts[cpu->rank]);

  float *tmp = (float*)calloc(cpu->mpiio_nparts[cpu->rank],sizeof(float));
  if (tmp == NULL){
    printf("can't allocate tmp vector in dump_HDF5_part\n");
    abort();
  }

  int ifield;
  for (ifield=0;ifield<param->out_part->n_field_tot-2; ifield++){ // the -2 is to exclude the age which doesnt exist for DM part and the mass which is not relevant
    if(param->out_part->field_id[ifield]){

      hid_t plist;

      //Set up file access property list with parallel I/O access
      plist = H5Pcreate(H5P_FILE_ACCESS);
#ifdef PHDF5
      H5Pset_fapl_mpio(plist, cpu->comm, MPI_INFO_NULL);
#endif // PHDF5

      //Create a new file collectively
      char file_name[256];
      sprintf(file_name,"data/%05d/part_%s_%05d.h5", *cpu->ndumps, param->out_part->field_name[ifield], *cpu->ndumps);
      hid_t file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
      H5Pclose(plist);

      hsize_t n_part_tot=0;
      int i;
      for (i=0;i<cpu->nproc;i++){
        n_part_tot+= (hsize_t)cpu->mpiio_nparts[i];
      }

      // Create the data space for the dataset.
      hid_t dataspace = H5Screate_simple(1, &n_part_tot, NULL);

      //Select hyperslab in the file.
      hsize_t offset = (hsize_t)cpu->mpiio_part_offsets;
      hsize_t n_loc = (hsize_t)cpu->mpiio_nparts[cpu->rank];

      if (debug) printf("npart_tot=%llu nloc=%llu offset=%llu \n",n_part_tot,n_loc,offset );

      H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&offset,NULL,&n_loc,NULL);

      // Create property list
      plist = H5Pcreate(H5P_DATASET_XFER);

#ifdef PHDF5
      H5Pset_dxpl_mpio(plist, HDF5_METHOD);
#endif // PHDF5
      hid_t	memspace = H5Screate_simple (1, &n_loc, NULL);


#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI
    double ti;
    ti =MPI_Wtime();

      //reduce data
      int i_tmp=0;
      int level;
      for(level=param->lcoarse;level<=param->lmax;level++){
        int iOct;
        for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
          struct OCT *oct=cpu->octList[level-1][iOct];
          int icell;
          for(icell=0;icell<8;icell++){ // looping over cells in oct
            struct PART * nexp=oct->cell[icell].phead; //sweeping the particles of the current cell
            if(nexp!=NULL){
              do{
                struct PART *curp=nexp;
                nexp=curp->next;
#ifdef STARS
                if(!(curp->isStar))
#endif
                  tmp[i_tmp++] = (float)assign_part_field(ifield,curp);
              }while(nexp!=NULL);
            }
          }
	      }
	    }

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI
      //if(cpu->rank == RANK_DISP) printf("reduce time %e\n",MPI_Wtime() - ti);

      char field_name[256];
      //sprintf(field_name,"%s",param->out_part->field_name[ifield]);
      sprintf(field_name,"data");


#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI
    ti =MPI_Wtime();

      hid_t dataset = H5Dcreate(file, field_name , H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, plist, tmp);
      H5Dclose(dataset);
#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI
//          if(cpu->rank == RANK_DISP) printf("write time %e\n",MPI_Wtime() - ti);


#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif // WMPI
    ti =MPI_Wtime();
      // tsim attribute
      hsize_t dims = 1;
      hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);
      hid_t attribute_id = H5Acreate2 (file, "a", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      float attr_data=(float)tsim;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(1, &dims, NULL);
      attribute_id = H5Acreate2 (file, "z", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      attr_data= 1./(float)tsim -1;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(1, &dims, NULL);
      attribute_id = H5Acreate2 (file, "t", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      attr_data= (float)param->cosmo->tphy;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      H5Pclose(plist);
      H5Fclose(file);
    }
  }

  free(tmp);
}

#ifdef STARS
void dump_HDF5_star(char filename[],REAL tsim,  struct RUNPARAMS *param, struct CPUINFO *cpu){

  const int debug =0;

  // Create step folder
  char folder_step[128];
  sprintf(folder_step,"data/%05d/",*(cpu->ndumps));
  mkdir(folder_step, 0755);

  float *tmp = (float*)calloc(cpu->mpiio_nstars[cpu->rank],sizeof(float));
  if (tmp == NULL){
    printf("can't allocate tmp vector in dump_HDF5_star\n");
    abort();
  }

  if (debug) printf("Starting main reducing/writting loop\n");
  int ifield;
  for (ifield=0;ifield<param->out_part->n_field_tot; ifield++){
    if(param->out_part->field_id[ifield]){

      hid_t plist;

      //Set up file access property list with parallel I/O access
      plist = H5Pcreate(H5P_FILE_ACCESS);
#ifdef PHDF5
      H5Pset_fapl_mpio(plist, cpu->comm, MPI_INFO_NULL);
#endif // PHDF5

      //Create a new file collectively
      char file_name[256];
      sprintf(file_name,"data/%05d/star_%s_%05d.h5", *cpu->ndumps, param->out_part->field_name[ifield], *cpu->ndumps);
      hid_t file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
      H5Pclose(plist);

      if (debug) printf("getting n stars tot \n");
      hsize_t n_star_tot=0;
      int i;
      for (i=0;i<cpu->nproc;i++){
        n_star_tot+=cpu->mpiio_nstars[i];
      }

      // Create the data space for the dataset.
      hid_t dataspace = H5Screate_simple(1, &n_star_tot, NULL);

      if (debug) printf("hyperslab \n");
      //Select hyperslab in the file.
      hsize_t offset = cpu->mpiio_star_offsets;
      hsize_t n_loc = cpu->mpiio_nstars[cpu->rank];
      H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&offset,NULL,&n_loc,NULL);

      // Create property list
      plist = H5Pcreate(H5P_DATASET_XFER);

#ifdef PHDF5
      H5Pset_dxpl_mpio(plist, HDF5_METHOD);
#endif // PHDF5

      hid_t	memspace = H5Screate_simple (1, &n_loc, NULL);

      if (debug>1) printf("begin reduce for field %s \n", param->out_part->field_name[ifield]);

      //reduce data
      int i_tmp=0;
      int level;
      for(level=param->lcoarse;level<=param->lmax;level++){
        int iOct;
        for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
          struct OCT *oct=cpu->octList[level-1][iOct];
          int icell;
          for(icell=0;icell<8;icell++){ // looping over cells in oct
            struct PART * nexp=oct->cell[icell].phead; //sweeping the particles of the current cell
            if(nexp!=NULL){
              do{
                struct PART *curp=nexp;
                nexp=curp->next;
                if(curp->isStar)
                  tmp[i_tmp++] = (float)assign_part_field(ifield,curp);
              }while(nexp!=NULL);
            }
          }
	      }
	    }

      char field_name[256];
      //sprintf(field_name,"%s",param->out_part->field_name[ifield]);
      sprintf(field_name,"data");

      if (debug>1) printf("creating dataset for field %s \n", param->out_part->field_name[ifield]);
      hid_t dataset = H5Dcreate(file, field_name , H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (debug>1) printf("writing dataset for field %s \n", param->out_part->field_name[ifield]);
      H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, plist, tmp);
      if (debug>1) printf("closing dataset for field %s \n", param->out_part->field_name[ifield]);
      H5Dclose(dataset);

      // tsim attribute
      hsize_t dims = 1;
      hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);
      hid_t attribute_id = H5Acreate2 (file, "a", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      float attr_data=(float)tsim;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(1, &dims, NULL);
      attribute_id = H5Acreate2 (file, "z", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      attr_data= 1./(float)tsim -1;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(1, &dims, NULL);
      attribute_id = H5Acreate2 (file, "t", H5T_NATIVE_FLOAT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      attr_data= (float)param->cosmo->tphy;
      H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      H5Pclose(plist);
      H5Fclose(file);
    }
  }

  free(tmp);
}
#endif // STARS

#endif // PIC
#endif // HDF5


void dumpIO(REAL tsim, struct RUNPARAMS *param,struct CPUINFO *cpu, struct OCT **firstoct, REAL *adt, int pdump){

  const int debug=0;
#ifdef WMPI
  MPI_Barrier(cpu->comm);
  double  tstart=MPI_Wtime();
#endif // WMPI

  REAL tdump,adump;
  char filename[128];

  int idir=cpu->rank%8;


#if defined(MPIIO) || defined(HDF5)
  set_offset(param,cpu);
#endif

#ifndef TESTCOSMO
#ifdef WRAD
	tdump=(tsim)*param->unit.unit_t/MYR;
#else
	tdump=(tsim);
#endif // WRAD
	adump=tdump;
#else
	tdump=interp_aexp(tsim,(double*)param->cosmo->tab_aexp,(double*)param->cosmo->tab_ttilde);
	adump=tdump;
#endif // TESTCOSMO

#ifdef MULTIFOLDER
  char folder_step[128];
  char folder_field[128];
  sprintf(folder_step,"data/%05d/",*(cpu->ndumps));
  makeFolders(cpu);
#endif // MULTIFOLDER

	if(pdump){
	  // === particle dump
#ifdef PIC

#ifdef MULTIFOLDER
    sprintf(folder_field,"%spart/",folder_step);
    sprintf(filename,"%spart.%05d.p%05d",folder_field,*(cpu->ndumps),cpu->rank);
#else
	  sprintf(filename,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MULTIFOLDER

	  if(cpu->rank==RANK_DISP){
	    printf("Dumping .......");
	    printf("%s %p\n",filename,cpu->part);
	 //   printf("%spart.%05d\n",folder_field,*(cpu->ndumps));
	  }

#ifdef WMPI
    MPI_Barrier(cpu->comm);
    REAL tt1=MPI_Wtime();
#endif // WMPI

#ifdef HDF5
  dump_HDF5_part(folder_field,adump,param, cpu);
#ifdef STARS
  if (param->stars->n) dump_HDF5_star(folder_field,adump,param, cpu);
#endif // STARS
#else
#ifndef MPIIO
    dumppart_serial(param,firstoct,filename,param->lcoarse,param->lmax,adump,cpu);
#else
    dumppart_MPI(firstoct,filename,param->lcoarse,param->lmax,adump,cpu, param);
#endif // MPIIO
#endif // HDF5

#ifdef WMPI
  MPI_Barrier(cpu->comm);
  REAL tt2=MPI_Wtime();
#endif // WMPI

//  printf("Total time for output writing %f \n", tt2-tt1);

#endif // PIC
	}
	else{
	  // === Hydro dump

#ifdef ALLOCT
#ifdef MULTIFOLDER
    sprintf(folder_field,"%s",folder_step);
    sprintf(filename,"%sgrid",folder_field);
#else
	  sprintf(filename,"data/grid.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MULTIFOLDER

#else
#ifdef MULTIFOLDER
    sprintf(folder_field,"%sgrid/",folder_step);
    sprintf(filename,"%sgrid.%05d.p%05d",folder_field,*(cpu->ndumps),cpu->rank);
#else
	  sprintf(filename,"data/grid.%05d.p%05d",*(cpu->ndumps),cpu->rank);
#endif // MULTIFOLDER
#endif // ALLOCT


	  if(cpu->rank==RANK_DISP){
	    printf("Dumping .......");
	    printf("%s\n",filename);
	  }
#ifdef WMPI
MPI_Barrier(cpu->comm);
REAL tt1=MPI_Wtime();
#endif // WMPI


#ifdef HDF5
  dump_HDF5_grid(folder_field,adump,param, cpu);
#else
#ifdef ALLOCT
#ifdef MPIIO
	  dumpalloct_MPI(folder_field,adump,param, cpu);
#else
    dumpalloct_serial(folder_field,adump,param, cpu);
#endif // MPIIO
#else
    dumpgrid(param->lmax,firstoct,filename,adump,param);
#endif // ALLOCT
#endif // HDF5

//#ifdef WMPI
//  MPI_Barrier(cpu->comm);
//  REAL tt2=MPI_Wtime();
//#endif // WMPI
//  printf("Total time for output writing %f \n", tt2-tt1);

#ifdef WMPI
  MPI_Barrier(cpu->comm);
  double tend=MPI_Wtime();
  if(cpu->rank==RANK_DISP) printf("===CPU DUMP TOTAL TIME =%fs \n", tend-tstart);
#endif // WMPI


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

#ifdef WMPI
  MPI_Barrier(cpu->comm);
  if(cpu->rank==RANK_DISP) printf("===CPU BKP TOTAL TIME =%fs \n", MPI_Wtime() - tend);
#endif // WMPI

#endif // BKP

	}



}

