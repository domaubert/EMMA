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
#include "parameters.h"

#include "segment.h"
#include "stars.h"
#include "hydro_utils.h"
#include "atomic_data/Atomic.h"
#include "tools.h"

#ifdef WMPI
#include <mpi.h>
#endif

#define AVGFACT (1.) // set to 0 to get an homogenous cosmo field 1

float assign_field(int field,struct CELL *cell){

/**
  * This function return the appropriate field, depending of a given ID.
  * see also parameters.c for the link between the ID and the input given in param.output.
  */

  float res;

  switch(field){

  case 0:
    res= (float)(cell2oct(cell)->x+( cell->idx   %2)*POW(0.5,cell2oct(cell)->level));
    break;

  case 1:
    res= (float)(cell2oct(cell)->y+( cell->idx/2 %2)*POW(0.5,cell2oct(cell)->level));
    break;

  case 2:
    res= (float)(cell2oct(cell)->z+( cell->idx/4   )*POW(0.5,cell2oct(cell)->level));
    break;

  case 3:
    res= (float)(cell2oct(cell)->level);
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
    res=cell->rfield.snfb;
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
    res=cell->rfield.src;
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
#endif // WRAD

  }

  return res;
}

void set_offset(struct RUNPARAMS *param, struct CPUINFO *cpu){
/**
  * This function compute the MPIIO offset
  * It cover the amr tree and count the number of cells in each processor domain.
  * Then it set the offset corresponding to the local processor rank.
  */

  // init n_cell to zero
  int i;
  for (i=0;i<cpu->nproc ;i++){
    cpu->mpiio_ncells[i]=0;
#ifdef PIC
    cpu->mpiio_nparts[i]=0;
#endif // PIC
#ifdef STARS
    cpu->mpiio_nstars[i]=0;
#endif // STARS
  }

  // count the cells
  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){

    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *oct=cpu->octList[level-1][iOct];

      int icell;
      for(icell=0;icell<8;icell++){
        struct CELL * cell= &oct->cell[icell];

        if(((oct->cell[icell].child==0)||(oct->level==param->lmax))){
          cpu->mpiio_ncells[cpu->rank]++;
        }

#ifdef PIC
        struct PART *nexp=cell->phead;
	      if(nexp!=NULL){
          do{
            struct PART *curp=nexp;
            nexp=curp->next;

#ifdef STARS
            if (curp->isStar) cpu->mpiio_nstars[cpu->rank]++;
            else              cpu->mpiio_nparts[cpu->rank]++;
#else
                              cpu->mpiio_nparts[cpu->rank]++;
#endif // STARS

        		}while(nexp!=NULL);
        }
#endif // PIC
      }
    }
  }

  // braodcast the result
  MPI_Allreduce(MPI_IN_PLACE,cpu->mpiio_ncells,cpu->nproc,MPI_INT,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,cpu->mpiio_nparts,cpu->nproc,MPI_INT,MPI_SUM,cpu->comm);

  // compute the offset
  cpu->mpiio_grid_offsets=0;

#ifdef PIC
    cpu->mpiio_part_offsets=0;
#endif // PIC
#ifdef STARS
    cpu->mpiio_star_offsets=0;
#endif // STARS
  for (i=1;i<=cpu->rank;i++){
    cpu->mpiio_grid_offsets += cpu->mpiio_ncells[i-1];
#ifdef PIC
    cpu->mpiio_part_offsets += cpu->mpiio_nparts[i-1];
#endif // PIC
#ifdef STARS
    cpu->mpiio_part_offsets += cpu->mpiio_nstars[i-1];
#endif // STARS
  }


  //printf("cpu=%d grid_offset=%d\n",cpu->rank, cpu->mpiio_grid_offsets);
  //printf("cpu=%d part_offset=%d\n",cpu->rank, cpu->mpiio_part_offsets);
  //printf("cpu=%d star_offset=%d\n",cpu->rank, cpu->mpiio_star_offsets);

}

#ifdef PIC
void dumppart_serial(struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu){

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

  FILE *fstar=NULL;
  fstar=fopen(filenamestar,"wb");
  if(fstar == NULL) {
    printf("Cannot open %s\n", filenamestar);
    abort();
  }
  fwrite(&nstar,1,sizeof(int)  ,fstar);
  fwrite(&tsimf,1,sizeof(float),fstar);

  FILE *fpart=NULL;
  fpart=fopen(filenamepart,"wb");
	if(fpart == NULL){
    printf("Cannot open %s\n", filenamepart);
    abort();
	}
  fwrite(&npart,1,sizeof(int)  ,fpart);
  fwrite(&tsimf,1,sizeof(float),fpart);

#else
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
		  else 			{	fp=fpart;	npart++;	}
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

void dumppart_MPIIO(struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct RUNPARAMS *param){

  float val;
  int vali;
  int ipart=0;
  int istar=0;
  int level;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  struct PART *nexp;
  struct PART *curp;
  int icell;
  float tsimf=tsim;
  int i;

  set_offset(param,cpu);

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

  MPI_File fpart=NULL;
  MPI_File_open(cpu->comm,filenamepart,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fpart);
	if(fpart == NULL){
    printf("Cannot open %s\n", filenamepart);
    abort();
	}
	MPI_File_write(fpart, &npart,1, MPI_INT, MPI_STATUS_IGNORE);
  MPI_File_write(fpart, &tsimf,1, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_seek(fpart, (MPI_Offset)(16*cpu->mpiio_part_offsets*sizeof(float)), MPI_SEEK_CUR);

#ifdef STARS
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

  int nstar=0;
  for (i=0;i<cpu->nproc;i++){
   nstar+=cpu->mpiio_nstars[i];
  }
  MPI_File fstar=NULL;
  MPI_File_open(cpu->comm,filenamestar,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fstar);
  if(fstar == NULL) {
    printf("Cannot open %s\n", filenamestar);
    abort();
  }
  MPI_File_write(fstar, &nstar,1, MPI_INT, MPI_STATUS_IGNORE);
  MPI_File_write(fstar, &tsimf,1, MPI_FLOAT, MPI_STATUS_IGNORE);
  MPI_File_seek(fstar, (MPI_Offset)(17*cpu->mpiio_star_offsets*sizeof(float)), MPI_SEEK_CUR);
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

      MPI_File *fp = NULL;
#ifdef STARS
		  if(curp->isStar) 		fp=&fstar;
		  else 			        	fp=&fpart;
#else
                          fp=&fpart;
#endif // STARS

		  val=(float)curp->x;			  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
      val=(float)curp->y;			  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)curp->z;			  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
#ifndef PARTN
#ifdef PART_EGY
		  val=(float)curp->ekin+curp->epot;
                                MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)curp->fx;		  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)curp->fy;		  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
#else
		  val=(float)curp->vx;		  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)curp->vy;		  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)curp->vz;		  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
#endif // PART_EGY
#else
		  val=(float)curp->fx;			MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)curp->fy;			MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)curp->fz;			MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
#endif // PARTN
		  val=(float)(curp->idx);	  MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)(curp->mass);	MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)(curp->epot);	MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  val=(float)(curp->ekin);	MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
#ifdef STARS
		  if(curp->isStar) {
		    val = curp->age;		    MPI_File_write(*fp, &val,1, MPI_FLOAT, MPI_STATUS_IGNORE);
		  }
#endif // STARS
//		  ipart++;

		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }


  MPI_File_close(&fpart);
#ifdef STARS
  MPI_File_close(&fstar);
#endif // STARS

  //printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);

}
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

  set_offset(param,cpu);

  int n_cell_tot=0;
  for (i=0;i<cpu->nproc;i++){
    n_cell_tot+=cpu->mpiio_ncells[i];
  }

// Open all the fields files
  MPI_File *f_dat;
  f_dat=(MPI_File *)malloc(param->out->n_field*sizeof(MPI_File));

  int n_field=0;
  for(i=0;i<param->out->n_field_tot;i++){
    if(param->out->field_id[i]){

      f_dat[n_field]=NULL;
      char dat_name[256];
      sprintf(dat_name,"%s%s.%05d",folder,param->out->field_name[i],*(cpu->ndumps));
      MPI_File_open(cpu->comm,dat_name,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&f_dat[n_field]);

      if(f_dat[n_field] == NULL){
       printf("Cannot open %s\n", dat_name);
       abort();
      }

      MPI_File_write(f_dat[n_field], &n_cell_tot,1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_seek(f_dat[n_field], cpu->mpiio_grid_offsets*sizeof(float), MPI_SEEK_CUR);

      n_field++;
    }
  }

// write the data
  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){

    int iOct;
    for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      struct OCT *oct=cpu->octList[level-1][iOct];

      int icell;
      for(icell=0;icell<8;icell++){
        struct CELL * cell= &oct->cell[icell];
        if(((oct->cell[icell].child==0)||(oct->level==param->lmax))){

          // get and write the fields
          int n_field=0;
          for (i =0;i<param->out->n_field_tot; i++){
            if(param->out->field_id[i]){
              float dat = (float)assign_field(i,cell);
              MPI_File_write(f_dat[n_field], &dat,1, MPI_FLOAT, MPI_STATUS_IGNORE);
              n_field++;
            }
          }
        }
      }
    }
  }

  // close all the fields files
  n_field=0;
  for(i=0; i<param->out->n_field; i++){
    if(param->out->field_id[i]){

      //MPI_File_sync( f_dat[n_field] ) ;
      MPI_File_close(&f_dat[n_field]);
      n_field++;
    }
  }

  free(f_dat);
}

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

  int i;
  int n_field=0;
  int n_cell=0;

// Opening all the fields files
  static FILE **f_dat;
  f_dat=(FILE **)malloc(param->out->n_field*sizeof(FILE *));

  for(i=0;i<param->out->n_field_tot;i++){
    if(param->out->field_id[i]){

      f_dat[i]=NULL;

      char folder_field[128];
      sprintf(folder_field,"%s%s/",folder, param->out->field_name[i]);
      mkdir(folder_field, 0755);
      char dat_name[256];
      sprintf(dat_name,"%s%s.%05d.p%05d",folder_field,param->out->field_name[i],*(cpu->ndumps),cpu->rank);
      f_dat[n_field]=fopen(dat_name,"wb");

      if(f_dat[n_field] == NULL){
       printf("Cannot open %s\n", dat_name);
       abort();
      }

      fwrite(&n_cell,1,sizeof(int),f_dat[n_field]);
      n_field++;
    }
  }

// writing the data
  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){

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
          for (i =0;i<param->out->n_field_tot; i++){
            if(param->out->field_id[i]){
              float dat = (float)assign_field(ii,cell);
              fwrite(&dat ,1,sizeof(float),f_dat[ii]);
              ii++;
            }
          }
        }
      }
    }
  }

  // write n_cells and close the fields files
  for (i=0; i<param->out->n_field; i++){
    rewind(f_dat[i]);
    fwrite(&n_cell,1,sizeof(int),f_dat[i]);
    fclose(f_dat[i]);
  }
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
  sprintf(folder_field,"%spart/",folder_step);
  mkdir(folder_field, 0755);
//  copy_param(folder_field);
#endif // PIC

#ifdef STARS
  sprintf(folder_field,"%sstar/",folder_step);
  mkdir(folder_field, 0755);
//  copy_param(folder_field);
#endif // STARS
#endif // MPIIO
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
	  }

#ifndef MPIIO
	  dumppart_serial(firstoct,filename,param->lcoarse,param->lmax,adump,cpu);
#else
	  dumppart_MPIIO(firstoct,filename,param->lcoarse,param->lmax,adump,cpu, param);
#endif // MPIIO

#endif // PIC
	}
	else{
	  // === Hydro dump

#ifdef ALLOCT
#ifdef MULTIFOLDER
    sprintf(folder_field,"%s",folder_step);
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
#ifdef ALLOCT
#ifdef MPIIO
	  dumpalloct_MPI(folder_field,adump,param, cpu);
#else
    dumpalloct_serial(folder_field,adump,param, cpu);
#endif // MPIIO
#else
    dumpgrid(param->lmax,firstoct,filename,adump,param);
#endif // ALLOCT

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
