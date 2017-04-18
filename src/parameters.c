#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <sys/stat.h>


#include "prototypes.h"
#include "io.h" // assign_grid_field

#if defined(UVBKG) || defined(STARS_TO_UVBKG)
#include "src_utils.h" //setUVBKG
#endif

int copy_file(char const * const source, char const * const destination) {
/**
  * copy a file name 'source' into 'destination'
  */
  FILE* fSrc;
  FILE* fDest;
  char buffer[512];
  int NbLus;

  if ((fSrc=fopen(source,"rb"))==NULL) {return 1;}

  if ((fDest=fopen(destination,"wb"))==NULL){
      fclose(fSrc);
      return 2;
  }

  while ((NbLus = fread(buffer, 1, 512, fSrc)) != 0)
      fwrite(buffer, 1, NbLus, fDest);

  fclose(fDest);
  fclose(fSrc);

  return 0;
}

void copy_param(const char *folder){
/**
  * copy all parameters files into each sudfolders
  */

  char param[128];
  char param_src[128];
  char param_dest[128];

  sprintf(param,"param.run");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);

  sprintf(param,"src/param.mk");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);

  sprintf(param,"param.info");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);

#ifdef GPUAXL
  sprintf(param,"param.avg.gpu");
#else
  sprintf(param,"param.avg.cpu");
#endif

  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);

  sprintf(param,"param.h");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);
#ifdef ALLOCT
  sprintf(param,"param.output");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);
  #endif // ALLOCT
}

void printFileOnScreen(char *filename){
/**
  * print file on screen
  */

  FILE* fp=fopen(filename,"r");
  if(fp == NULL){
    printf("Cannot open %s : you may want to link it in the current directory\n", filename);
  }
  char ch;
  while((ch=fgetc(fp))!=EOF) printf("%c",ch);
}

void dumpFile(char *filename_in, char *filename_out){
/**
  * copy filename_in into filename_out and print it on screen
  */
  int fileok=1;
  FILE *fps[2] = {stdout, NULL};
  fps[1]=fopen(filename_out,"w");
  if(fps[1] == NULL) {
    printf("Cannot open %s\n", filename_out);
    fileok=0;
  }

  FILE* buf=NULL;
  buf=fopen(filename_in,"r");
  if(buf == NULL){
    printf("Cannot open %s : you may want to link it in the current directory\n", filename_in);
    fileok=0;
  }

  int i;
  if(fileok){
    for(i=0;i<2;i++){
      FILE *fp = fps[i];
      char ch;
      fseek(buf,0,SEEK_SET);
      while((ch=fgetc(buf))!=EOF) fprintf(fp,"%c",ch);
    }

    fclose(fps[1]);
    fclose(buf);
  }
}

void readOutputParam_grid(char *fparam, struct RUNPARAMS *param){
  int debug=0;

  char *field_name [] ={
    // The field order has to be the same as in param.output for consistency
    "x",
    "y",
    "z",
    "l",
    "cpu",
    "gdata_d",
    "density",
    "gdata.p",
    "res",
    "f0",
    "f1",
    "f2",
    "marked",
    "field_d",
    "field_u",
    "field_v",
    "field_w",
    "field_p",
    "rfield_e0",
    "rfield_fx0",
    "rfield_fy0",
    "rfield_fz0",
    "rfield_e1",
    "rfield_fx1",
    "rfield_fy1",
    "rfield_fz1",
    "rfield_snfb",
    "rfield_e2",
    "rfield_fx2",
    "rfield_fy2",
    "rfield_fz2",
    "rfield_src",
    "xion",
    "field_dX",
    "field_dXHE",
    "field_dXXHE",
    "field_xHE",
    "field_xxHE",
    "rfield_temp",
    "z_first_xion",
    "z_last_xion"
  };

  int n_field=0;
  int n_field_tot=0;
  param->out_grid->n_field_movie=0;

  FILE *f=fopen(fparam,"r");
  if(f==NULL){
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
  }

  if (debug) printf("opening OK\n");

  char stream[256];

  while (fscanf(f,"%s",stream)!= EOF) {
    param->out_grid->field_name[n_field_tot]= field_name[n_field_tot];
    if( !strncmp(stream, "#" ,1)){
      // skip line if first char == #
      char* line=NULL;
      size_t len=0;
      size_t status= getline(&line,&len,f);
      if (debug) printf("%s\t",line);
      free(line);
      continue;
    }

#ifdef SINGLEPRECISION
    char* type= "%d %d %d %e %e\n";
#else
    char* type= "%d %d %d %le %le\n";
#endif // SINGLEPRECISION

    size_t status=fscanf(f, type,
          &(param->out_grid->field_state_grid[n_field_tot]),
          &(param->out_grid->field_state_movie[n_field_tot]),
          &(param->out_grid->field_state_stat[n_field_tot]),
          &(param->physical_state->field[n_field_tot].bin_min),
          &(param->physical_state->field[n_field_tot].bin_max)
          );

    if (debug) printf("%s\t",stream);
    if (debug) printf("grid=%d stat=%d\t",param->out_grid->field_state_grid[n_field_tot],  param->out_grid->field_state_stat[n_field_tot]);
    if (debug) printf("bin_min=%e bin_max=%e\n", param->physical_state->field[n_field_tot].bin_min,  param->physical_state->field[n_field_tot].bin_max);


    if (param->out_grid->field_state_movie[n_field_tot]) param->out_grid->n_field_movie++;

    if (param->out_grid->field_state_grid[n_field_tot]){
      param->out_grid->field_id[n_field_tot]=1;

      n_field++;
    }else{
      param->out_grid->field_id[n_field_tot]=0;
    }
    n_field_tot++;
  }

  if (debug) printf("n_field_tot=%d\n",n_field_tot);

  fclose(f);

  if (debug) printf("read OK\n");

  param->out_grid->n_field=n_field;
  param->out_grid->n_field_tot=n_field_tot;

  if(debug){
    int i;
    for (i=0;i<n_field_tot; i++){
      if (param->out_grid->field_state_movie[i])
      printf("%d\t%s\n",param->out_grid->field_id[i], param->out_grid->field_name[i]);
    }
  abort();
  }

}

void readOutputParam_part(char *fparam, struct RUNPARAMS *param){

  int debug=0;

  char *field_name [] ={
    // The field order has to be the same as in param.output for consistency
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "fx",
    "fy",
    "fz",
    "idx",
    "isStar",
    "epot",
    "ekin",
    "mass",
    "age"
  };

  int n_field=0;
  int n_field_tot=0;

  FILE *f=NULL;
  f=fopen(fparam,"r");
  if(f==NULL){
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
  }
  char stream[256];
  int state;

  while (fscanf(f,"%s",stream)!= EOF) {
  param->out_part->field_name[n_field_tot]= field_name[n_field_tot];
  if( !strncmp(stream, "#" ,1)){
    // skip line if first char == #
    char* line=NULL;
    size_t len=0;
    size_t status= getline(&line,&len,f);
    if (debug) printf("%s\t",line);
    free(line);
    continue;
  }

  size_t status=fscanf(f, "%d\n",&state);
    if (state){
      param->out_part->field_id[n_field_tot] = 1;
      param->out_part->field_name[n_field_tot] = field_name[n_field_tot];
      n_field++;
    }else{
      param->out_part->field_id[n_field_tot] = 0;
    }
    n_field_tot++;
  }
  fclose(f);


  param->out_part->n_field=n_field;
  param->out_part->n_field_tot=n_field_tot;

  if(debug){
    int i;
    for (i=0;i<n_field_tot; i++){
      if (param->out_part->field_id[i])
      printf("%d\t%s\n",param->out_part->field_id[i], param->out_part->field_name[i]);
    }
    abort();
  }

}

#ifdef SUPERNOVAE
void read_egy_loss(struct RUNPARAMS *param){
/**
  * Read parameter file
  * egy_loss.dat is a rewriting of the Starburst99 figure of the integrated energy loss
  * fig 115 : http://www.stsci.edu/science/starburst99/figs/energy_inst_e.html
  * input are in joules / kilogram
  */

  int size = 10000;

  param->sn->egy_loss_t =	(REAL *)calloc(size,sizeof(REAL ));
  param->sn->egy_loss_egy =	(REAL *)calloc(size,sizeof(REAL ));

  char* fname = "src/phys_data/egy_loss.dat";
  FILE *f=fopen(fname,"r");
  if(f==NULL){
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fname);
      abort();
  }

  int i=0;
  double t,egy;
  while (fscanf(f, "%lf %lf\n", &t,&egy ) != EOF) {

    param->sn->egy_loss_t[i]=(REAL)t;
    param->sn->egy_loss_egy[i]=(REAL)egy;
    i++;
  }
  fclose(f);

  int debug=0;
  if(debug){
    for(i=0;i<size;i++){
      t=param->sn->egy_loss_t[i];
      egy=param->sn->egy_loss_egy[i];
      printf("%e %e\n", t,egy);
    }
    abort();
  }
}

void read_mass_loss(struct RUNPARAMS *param){
/**
  * Read parameter file
  * mass_loss.dat is a rewriting of the Starburst99 figure of the integrated mass loss
  * fig 109 : http://www.stsci.edu/science/starburst99/figs/mass_inst_e.html
  * input are in percent of the particle mass.
  */

  int size = 10000;

  param->sn->mass_loss_t =	(REAL *)calloc(size,sizeof(REAL ));
  param->sn->mass_loss_mass=	(REAL *)calloc(size,sizeof(REAL ));

  char* fname = "src/phys_data/mass_loss.dat";
  FILE *f=fopen(fname,"r");
  if(f==NULL){
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fname);
      abort();
  }

  int i=0;
  double t,mass;
  while (fscanf(f, "%lf %lf\n", &t,&mass) != EOF) {

    param->sn->mass_loss_t[i]= (REAL)t;
    param->sn->mass_loss_mass[i]= (REAL)mass;
    i++;
  }
  fclose(f);

  int debug=0;
  if(debug){
    for(i=0;i<size;i++){
      t=param->sn->mass_loss_t[i];
      mass=param->sn->mass_loss_mass[i];
      printf("%e %e\n", t,mass);
    }
    abort();
  }
}

#endif // SUPERNOVAE

#ifdef WRAD
void readAtomic(struct RUNPARAMS *param){
/**
  * Read the atomic data from file defined in param.run
  */
  int debug=0; //print what was readed and abort if !=0

#ifdef SINGLEPRECISION
    char* type= "%e\n";
#else
    char* type= "%le\n";
#endif // SINGLEPRECISION

  //openfile
  FILE *buf=fopen(param->atomic.path,"r");
  if(buf==NULL){
    printf("ERROR : cannot open the parameter file (%s given), please check\n",param->atomic.path);
    abort();
  }
  if(debug) printf("reading %s\n",param->atomic.path);
  size_t rstat;
  char stream[256];

  //alloc space
  rstat=fscanf(buf,"%s %d",stream,&param->atomic.ngrp_space);
  if(debug) printf("ngrp_space=%d\n", param->atomic.ngrp_space);
  param->atomic.space_bound= (REAL*)calloc(param->atomic.ngrp_space+1,sizeof(REAL));
  param->atomic.space_bound[param->atomic.ngrp_space]=INFINITY;

  //read space bound
  rstat=fscanf(buf,"%s",stream); // read SPACE_BOUND(eV)
  int i_space;
  for (i_space=0; i_space<param->atomic.ngrp_space; i_space++){
    rstat=fscanf(buf,type,&param->atomic.space_bound[i_space]);
    if(debug) printf("space_bound[%d]%e\n",i_space, param->atomic.space_bound[i_space]);
  }
  rstat=fscanf(buf,"%s",stream); //read inf


  //alloc time
  rstat=fscanf(buf,"%s %d",stream,&param->atomic.ngrp_time);
  if(debug) printf("ngrp_time=%d\n", param->atomic.ngrp_time );
  param->atomic.time_bound= (REAL*)calloc(param->atomic.ngrp_time+1,sizeof(REAL));
  param->atomic.time_bound[param->atomic.ngrp_time]=INFINITY;

  //read time bound
  rstat=fscanf(buf,"%s",stream); //read TIME_BOUND(MYR)
  int i_time;
  for (i_time=0; i_time<param->atomic.ngrp_time; i_time++){
    rstat=fscanf(buf,type,&param->atomic.time_bound[i_time]);
    if(debug) printf("time_bound[%d]%e\n",i_time, param->atomic.time_bound[i_time]);
    param->atomic.time_bound[i_time]*=1e6; //Myrs in yrs
  }
  rstat=fscanf(buf,"%s",stream); //read inf

  // alloc grp
  param->atomic.n = param->atomic.ngrp_space * param->atomic.ngrp_time;
  param->atomic.hnu= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphae= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphai= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.factgrp= (REAL*)calloc(param->atomic.n,sizeof(REAL));
#ifdef HESIMPLE
  param->atomic.alphaeHE= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphaiHE= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphaeHE2= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphaiHE2= (REAL*)calloc(param->atomic.n,sizeof(REAL));
#endif
  //skip header
  int i;
  for(i=0;i<4;i++)  rstat=fscanf(buf,"%s",stream);

  //read grp
  double dummy; // to read recombination coeff
  for(i=0;i<param->atomic.n;i++){
    rstat=fscanf(buf,"%lf",&dummy); param->atomic.hnu[i]=dummy*ELECTRONVOLT;
    rstat=fscanf(buf,"%lf",&dummy); param->atomic.alphae[i]=param->clightorg*(dummy*LIGHT_SPEED_IN_M_PER_S);
    rstat=fscanf(buf,"%lf",&dummy); param->atomic.alphai[i]=param->clightorg*(dummy*LIGHT_SPEED_IN_M_PER_S);
#ifdef HESIMPLE
    rstat=fscanf(buf,"%lf",&dummy); param->atomic.alphaeHE[i]=param->clightorg*(dummy*LIGHT_SPEED_IN_M_PER_S);
    rstat=fscanf(buf,"%lf",&dummy); param->atomic.alphaiHE[i]=param->clightorg*(dummy*LIGHT_SPEED_IN_M_PER_S);
    rstat=fscanf(buf,"%lf",&dummy); param->atomic.alphaeHE2[i]=param->clightorg*(dummy*LIGHT_SPEED_IN_M_PER_S);
    rstat=fscanf(buf,"%lf",&dummy); param->atomic.alphaiHE2[i]=param->clightorg*(dummy*LIGHT_SPEED_IN_M_PER_S);
#endif
    rstat=fscanf(buf,"%lf",&dummy);param->atomic.factgrp[i]=dummy;
    if(debug) printf("%e %e %e %e \n", param->atomic.hnu[i], param->atomic.alphae[i], param->atomic.alphai[i], param->atomic.factgrp[i]);
#ifdef HESIMPLE
    if(debug) printf("%e %e \n",param->atomic.alphaeHE[i], param->atomic.alphaiHE[i]);
#endif

  }
  fclose(buf);
  if(debug) abort();
}
#endif // WRAD

void ReadParameters(char *fparam, struct RUNPARAMS *param){
  int debug=0;
  int i;
  FILE *buf=NULL;
  char stream[256];
  size_t rstat;
  double dummyf;
  char RF[]="%s %lf";

  buf=fopen(fparam,"r");
  if(buf==NULL)
    {
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
    }
  else
    {

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ngridmax);  if (debug) printf("param->ngridmax=%d\n", param->ngridmax);
      rstat=fscanf(buf,"%s %d",stream,&param->npartmax);  if (debug) printf("param->npartmax=%d\n", param->npartmax);
      rstat=fscanf(buf,"%s %d",stream,&param->nbuff);     if (debug) printf("param->nbuff=%d\n", param->nbuff);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ndumps);    if (debug) printf("param->ndumps=%d\n", param->ndumps);
      rstat=fscanf(buf,RF,stream,&param->dt_dump);        if (debug) printf("param->dt_dump=%e\n", param->dt_dump);
      rstat=fscanf(buf,"%s %d",stream,&param->nsteps);    if (debug) printf("param->nsteps=%d\n", param->nsteps);
      rstat=fscanf(buf,RF,stream,&dummyf);param->dt=(REAL)dummyf;   if (debug) printf("param->dt=%e\n", param->dt);
      rstat=fscanf(buf,RF,stream,&dummyf);param->tmax=(REAL)dummyf; if (debug) printf("param->tmax=%e\n", param->tmax);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->lcoarse); if (debug) printf("param->lcoarse=%d\n", param->lcoarse);
      rstat=fscanf(buf,"%s %d",stream,&param->lmax);    if (debug) printf("param->lmax=%d\n", param->lmax);
      rstat=fscanf(buf,RF,stream,&dummyf);param->amrthresh0=(REAL)dummyf; if (debug) printf("param->amrthresh=%e\n", param->amrthresh0);
      rstat=fscanf(buf,"%s %d",stream,&param->nsmooth); if (debug) printf("param->nsmooth=%d\n", param->nsmooth);


      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->DM_res);  if (debug) printf("param->DM_res=%d\n", param->DM_res);
      rstat=fscanf(buf,RF,stream,&dummyf);param->dx_res=(REAL)dummyf; if (debug) printf("param->dx_res=%e\n", param->dx_res);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->niter);   if (debug) printf("param->niter=%d\n", param->niter);
      rstat=fscanf(buf,RF,stream,&dummyf);param->poissonacc=(REAL)dummyf; if (debug) printf("param->poissonacc=%e\n", param->poissonacc);
      rstat=fscanf(buf,"%s %d",stream,&param->mgridlmin); if (debug) printf("param->mgridlmin=%d\n", param->mgridlmin);
      if(param->mgridlmin<0){
        param->mgridlmin=param->lcoarse-param->lcoarse;
      }

      rstat=fscanf(buf,"%s %d",stream,&param->nvcycles);  if (debug) printf("param->nvcycles=%d\n", param->nvcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nrelax);    if (debug) printf("param->nrelax=%d\n", param->nrelax);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->nrestart);  if (debug) printf("param->nrestart=%d\n", param->nrestart);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->gstride);    if (debug) printf("param->gstride=%d\n", param->gstride);
      rstat=fscanf(buf,"%s %d",stream,&param->hstride);    if (debug) printf("param->hstride=%d\n", param->hstride);
      rstat=fscanf(buf,"%s %d",stream,&param->nsubcycles); if (debug) printf("param->nsubcycles=%d\n", param->nsubcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nthread);    if (debug) printf("param->nthread=%d\n", param->nthread);
      rstat=fscanf(buf,"%s %d",stream,&param->nstream);    if (debug) printf("param->nstream=%d\n", param->nstream);
      rstat=fscanf(buf,"%s %d",stream,&param->ompthread);  if (debug) printf("param->ompthread=%d\n", param->ompthread);

      rstat=fscanf(buf,"%s",stream);
#ifdef WRAD
      rstat=fscanf(buf,RF,stream,&dummyf);param->clight=(REAL)dummyf;param->clightorg=(REAL)dummyf; if (debug) printf("param->clight=%e\n", param->clight);
      rstat=fscanf(buf,RF,stream,&dummyf);param->denthresh=(REAL)dummyf;  if (debug) printf("param->denthresh=%e\n", param->denthresh);
      rstat=fscanf(buf,RF,stream,&dummyf);param->tmpthresh=(REAL)dummyf;  if (debug) printf("param->tmpthresh=%e\n", param->tmpthresh);
      rstat=fscanf(buf,RF,stream,&dummyf);param->srcint=(REAL)dummyf;     if (debug) printf("param->srcint=%e\n", param->srcint);
      rstat=fscanf(buf,RF,stream,&dummyf);param->fesc=(REAL)dummyf;       if (debug) printf("param->fesc=%e\n", param->fesc);

      param->srcint*=param->fesc;

      char filename[256];
      rstat=fscanf(buf,"%s %s",stream, filename);
      sprintf(param->atomic.path,"./SRC/src/atomic_data/%s",filename);   if (debug) printf("param->atomic.path=%s\n", param->atomic.path);
      param->fudgecool=1.0;
      param->ncvgcool=0;
#else
      int i;
	for (i=0; i<5; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
                      rstat=fscanf(buf,"%s %s",stream, stream);
#endif

      rstat=fscanf(buf,"%s",stream);
#ifdef STARS
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->overdensity_cond=(REAL)dummyf;   if (debug) printf("param->stars->overdensity_cond=%e\n", param->stars->overdensity_cond);
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->density_cond=(REAL)dummyf;       if (debug) printf("param->stars->density_cond=%e\n", param->stars->density_cond);
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->efficiency=(REAL)dummyf;         if (debug) printf("param->stars->efficiency=%e\n", param->stars->efficiency);
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->tlife=(REAL)dummyf;              if (debug) printf("param->stars->tlife=%e\n", param->stars->tlife);
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->mass_res=(REAL)dummyf;           if (debug) printf("param->stars->mass_res=%e\n", param->stars->mass_res);

#else
	for (i=0; i<5; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
#endif

      rstat=fscanf(buf,"%s",stream);
#ifdef SUPERNOVAE
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->feedback_eff	=(REAL)dummyf;        if (debug) printf("param->sn->feedback_eff=%e\n", param->sn->feedback_eff);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->feedback_frac	=(REAL)dummyf;      if (debug) printf("param->sn->feedback_frac=%e\n", param->sn->feedback_frac);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->ejecta_proportion	=(REAL)dummyf;  if (debug) printf("param->sn->ejecta_proportion=%e\n",param->sn->ejecta_proportion);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->sn_egy	=(REAL)dummyf;              if (debug) printf("param->sn->sn_egy=%e\n", param->sn->sn_egy);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->tlife	=(REAL)dummyf;              if (debug) printf("param->sn->tlife=%e\n", param->sn->tlife);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->load_factor =(REAL)dummyf;         if (debug) printf("param->sn->load_factor=%e\n", param->sn->load_factor);

#else
	for (i=0; i<6; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
#endif

#ifdef MOVIE
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->movie->lmap);                    if (debug) printf("param->movie->lmap=%d\n", param->movie->lmap);
      if (param->movie->lmap>param->lmax) param->movie->lmap=param->lmax;
      rstat=fscanf(buf,"%s %s",stream,param->movie->mode_str);                if (debug) printf("param->movie->mode_str=%s\n", param->movie->mode_str);

      rstat=fscanf(buf,RF,stream,&param->movie->xmin);    if (debug) printf("param->movie->xmin=%e\n", param->movie->xmin);
      rstat=fscanf(buf,RF,stream,&param->movie->xmax);    if (debug) printf("param->movie->xmax=%e\n", param->movie->xmax);
      rstat=fscanf(buf,RF,stream,&param->movie->ymin);    if (debug) printf("param->movie->ymin=%e\n", param->movie->ymin);
      rstat=fscanf(buf,RF,stream,&param->movie->ymax);    if (debug) printf("param->movie->ymax=%e\n", param->movie->ymax);
      rstat=fscanf(buf,RF,stream,&param->movie->zmin);    if (debug) printf("param->movie->zmin=%e\n", param->movie->zmin);
      rstat=fscanf(buf,RF,stream,&param->movie->zmax);    if (debug) printf("param->movie->zmax=%e\n", param->movie->zmax);
#endif
      fclose(buf);
    }

 if (debug) abort();
}

void GetParameters(char *fparam, struct RUNPARAMS *param){

  ReadParameters(fparam, param);

  // computing the maxhash
  int val=(POW(2,param->lmax-1)<=512?POW(2,param->lmax-1):512); // limit to 2097152 octs in hash table i.e. 16e6 cells
  param->maxhash=POW(val,3);
  //printf("maxhash=%d\n",param->maxhash);

  // ====================== some checks

  // stencil/streams conformity
#ifdef GPUAXL
  if(param->hstride<(param->nthread*param->nstream)){
    printf(" Stream Thread granulosity too high : nt=%d ns=%d stencil=%d\n",param->hstride,param->nthread,param->nstream);
    abort();
  }
#endif

#ifdef STARS
    param->stars->n		= 0;
#ifdef AGN
    param->stars->nagn		= 0;
#endif
#endif

#ifdef SRCINT
  param->srcint*=SRCINT;
#endif

#if defined(UVBKG) || defined(STARS_TO_UVBKG)
  setUVBKG(param, "SRC/src/phys_data/uvbkg.dat");
#endif // UVBKG


#ifdef WRAD
  readAtomic(param);
#endif // WRAD

#ifdef SUPERNOVAE
  //read_egy_loss(param);
  //read_mass_loss(param);
#endif // SUPERNOVAE

#if defined(WRADTEST) || defined(SNTEST)
/*
  param->unitary_stars_test->lifetime = 3.673e6;
  param->unitary_stars_test->mass=2e2;
  param->unitary_stars_test->src_pos_x=0.;
  param->unitary_stars_test->src_pos_y=0.;
  param->unitary_stars_test->src_pos_z=0.;
*/

  param->unitary_stars_test->lifetime = 8e6;
  param->unitary_stars_test->mass=2e2;
  param->unitary_stars_test->src_pos_x=0.5;
  param->unitary_stars_test->src_pos_y=0.5;
  param->unitary_stars_test->src_pos_z=0.5;
#endif // defined

}

void dumpInfo(char *filename_info, struct RUNPARAMS *param, struct CPUINFO *cpu){
/**
  * Write a file containing information about the simulation like the cosmology or the resolution.
  */
  FILE *fps[2] = {stdout, NULL};

  fps[1]=fopen(filename_info,"w");
  if(fps[1] == NULL) printf("Cannot open %s\n", filename_info);

  char* int_format = "%-24s%d\n";
  char* float_format = "%-24s%.2f\n";
  char* real_format = "%-24s%e\n";

  int i;
  for(i=0;i<2;i++){
    FILE *fp = fps[i];

    fprintf(fp, int_format,"nproc",(cpu->nproc)); 		// number of processor
#ifdef TESTCOSMO
    fprintf(fp, float_format,"box_size_hm1_Mpc",(param->unit.unit_l/PARSEC/1e6*param->cosmo->H0/100));
#else
    fprintf(fp, float_format,"box_size_Kpc",(param->unit.unit_l/PARSEC/1e3));
#endif // TESTCOSMO
    fprintf(fp, int_format,"level_min",(param->lcoarse) );
    fprintf(fp, int_format,"level_max",(param->lmax) );

  fprintf(fp,"##=Unit_code->SI====================\n" );

  fprintf(fp, real_format,"unit_l",(param->unit.unit_l) );		// comoving length size of the box [meters]
  fprintf(fp, real_format,"unit_v",(param->unit.unit_v) );		// unit velocity
  fprintf(fp, real_format,"unit_t",(param->unit.unit_t) );		// unit time [seconds]
  fprintf(fp, real_format,"unit_N",(param->unit.unit_N) );		// unit number [moles typically]
  fprintf(fp, real_format,"unit_mass",(param->unit.unit_mass) );	// unit mass [in kg, total mass is equal to one in unit codes]
  //  fprintf(fp,"\n");

  #ifdef TESTCOSMO
  fprintf(fp,"##=Cosmology========================\n" );
    fprintf(fp, real_format,"om",(param->cosmo->om) );			// Omega matter
    fprintf(fp, real_format,"ov",(param->cosmo->ov) );			// Omega vacuum
    fprintf(fp, real_format,"ob",(param->cosmo->ob) );			// Omega baryon
    fprintf(fp, real_format,"H0",(param->cosmo->H0) );			// Hubble constant
  //  fprintf(fp,"\n");
  #endif

#ifdef PIC
    fprintf(fp,"##=Mass_resolution_(Mo)=============\n" );
#ifdef TESTCOSMO
    //REAL mass_res_DM =  (1.- param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse+param->DM_res))*param->unit.unit_mass/SOLAR_MASS;

    // below hack to get the correct result even in the case of SP calculation
    double munpercell=param->cosmo->om*(3.*pow(param->cosmo->H0*1e3/PARSEC/1e6,2)/(8.*M_PI*NEWTON_G))*pow(param->unit.unit_l/pow(2.0,param->lcoarse),3);
    REAL mass_res_DM=(REAL)((1.- param->cosmo->ob/param->cosmo->om)*munpercell/pow(2.0,3.*param->DM_res)/SOLAR_MASS);

    fprintf(fp, real_format,"mass_res_DM",mass_res_DM );
#endif // TESTCOSMO
#ifdef STARS
    REAL res = param->stars->mass_res;
    if(res>100){
        fprintf(fp, real_format,"mass_res_star",param->stars->mass_res);
    }else{

      if(res>=0){
        /* REAL mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse+param->stars->mass_res)); */
        /* REAL mass_res_star = mstars_level * param->unit.unit_mass /SOLAR_MASS; */

#ifdef TESTCOSMO
	double mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->stars->mass_res));
	REAL mass_res_star =(REAL)( mstars_level * munpercell/SOLAR_MASS);
#else
  //TODO considere non cosmo case
  double mstars_level=0;
  REAL mass_res_star=0;
#endif // TESTCOSMO

        fprintf(fp, real_format,"mass_res_star",mass_res_star);
      }else{
        int level;
        for(level=param->lcoarse;level<=param->lmax;level++){
          REAL mlevel=level-1;
          REAL restmp=-param->stars->mass_res;
#ifdef TESTCOSMO
          REAL mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(mlevel+restmp));
          REAL mass_res_star = mstars_level * param->unit.unit_mass /SOLAR_MASS;
#else
  //TODO considere non cosmo case
  double mstars_level=0;
  REAL mass_res_star=0;
#endif // TESTCOSMO
          char mlev[128];
          sprintf(mlev,"mass_star_L%d",level);
          fprintf(fp, real_format,mlev,mass_res_star );
        }
      }
    }
#endif // STARS
#endif // PIC

    fprintf(fp,"##=Spatial_resolution_(Kpc)=========\n" );
    int level;
    for(level=param->lcoarse;level<=param->lmax;level++){
      REAL dx = param->unit.unit_l * POW(2,-level) /1e3/PARSEC;
      char dxlev[128];
      sprintf(dxlev,"dx_L%d",level);
      fprintf(fp, real_format,dxlev,dx);
    }

  }
  fclose(fps[1]);
}

void dumpHeader(struct RUNPARAMS *param, struct CPUINFO *cpu,char *fparam){
/**
  * Dump on screen and into files, a set of parameters
  */

  printf("\n");
  printf("--------------------------------------------------------------\n");
  dumpInfo("data/param.info", param, cpu);
  printf("\n");
  printf("--------------------------------------------------------------\n");
  printFileOnScreen("SRC/param.mk");
  printf("\n");
  printf("--------------------------------------------------------------\n");
  printFileOnScreen("SRC/src/param.h");
  printf("\n");
  printf("--------------------------------------------------------------\n");
  printFileOnScreen(param->paramrunfile);
  printf("\n");
  printf("--------------------------------------------------------------\n");
#ifdef ALLOCT
  char partoutput[512];
  strcpy(partoutput,param->paramrunfile);
  strcat(partoutput,".part_output");
  printFileOnScreen(partoutput);
  printf("\n");
  printf("--------------------------------------------------------------\n");
  char gridoutput[512];
  strcpy(gridoutput,param->paramrunfile);
  strcat(gridoutput,".grid_output");
  printf("\n");
  printf("--------------------------------------------------------------\n");
#endif // ALLOCT


#ifdef TESTCOSMO
  REAL threshold=param->amrthresh0;
#ifndef ZOOM
  threshold*=POW(2.0,-3.0*param->lcoarse);
#else
  threshold*=POW(2.0,-3.0*param->lmaxzoom);
#endif
  if(cpu->rank==RANK_DISP)
    printf("amrthresh : maximum number of part in a cell before refinement : %d -> compute density thresold of %e \n ", (int)param->amrthresh0, threshold);
#endif

#ifdef WRAD
    printf("SRCINT set to %e\n",param->srcint);

#ifndef SRCINT
  if(param->srcint!=0) // usefull to unactivate sources
  if(param->srcint<2.){ // it is likely to be an error
    if(cpu->rank==RANK_DISP){
      printf("ERROR FLAG SRCINT NOT DEFINED: param->srcint = %e Photons/s/kg\n. ",param->srcint);
    }
    abort();
  }
#endif
#endif

#ifdef UVBKG
    printf("UVBKG is defined -> clight set to 1e-4 \n");
#endif // UVBKG
  //abort();

  printf("\n");
}

void printFieldInfo(struct FIELD_INFO *field){
  printf("************\n");
  printf("min=%e\n",field->min);
  printf("max=%e\n",field->max);
  printf("mean=%e\n",field->mean);
  printf("sigma=%e\n",field->sigma);
}

void initFieldInfo(struct FIELD_INFO *field, REAL pdf_min, REAL pdf_max){

  int debug=0;

  field->mean=0.;
  field->min=INFINITY;
  field->max=0.;
  field->sigma=0.;

  //int type = 0; //lin bins
  int type = 1; //log bins

  if (type){
    double bin= (log10(pdf_max)-log10(pdf_min))/N_BIN_PDF;
    if(debug)printf("bin=%e\n",bin);
    int i;
    for(i=0;i<N_BIN_PDF+1;i++){
      field->bins_edges[i]=POW(10, log10(pdf_min) + i*bin  );
      if(debug) printf("field->bins[i]=%e\n",field->bins_edges[i]);
    }
  }else{
    double bin= (pdf_max-pdf_min)/N_BIN_PDF;
    int i;
    for(i=0;i<N_BIN_PDF+1;i++){
      field->bins_edges[i]=i*bin;
      //printf("field->bins[i]=%e\n",field->bins[i]);
    }
  }

  int i;
  for(i=0;i<N_BIN_PDF;i++){
    field->pdf[i]=0;
  }
  if(debug) abort();
}

void getFieldInfo(struct FIELD_INFO *field, double value, double vweight){
  field->mean+=value*vweight;
  field->sigma+=value*value*vweight;
  field->min=fmin(field->min, value);
  field->max=fmax(field->max, value);
/*
  int i;
  for(i=0;i<N_BIN_PDF;i++){
    if( (value>field->bins_edges[i]) && (value<=field->bins_edges[i+1]) ){
      field->pdf[i]+=value*vweight;
      break;
    }
  }
  */
}

#ifdef WMPI
void comFieldInfo(struct CPUINFO *cpu, struct FIELD_INFO *field){
  MPI_Allreduce(MPI_IN_PLACE,&(field->mean  ),1,MPI_DOUBLE,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&(field->min   ),1,MPI_DOUBLE,MPI_MIN,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&(field->max   ),1,MPI_DOUBLE,MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&(field->sigma ),1,MPI_DOUBLE,MPI_SUM,cpu->comm);

  MPI_Allreduce(MPI_IN_PLACE,field->pdf ,N_BIN_PDF,MPI_DOUBLE,MPI_SUM,cpu->comm);
}
#endif // WMPI

void setSigmaFieldInfo(struct FIELD_INFO *field){
  field->sigma = sqrt(field->sigma - field->mean*field->mean);
}

void writeFieldInfoHeader(FILE* fp,struct FIELD_INFO *field){

  fprintf(fp,"Nbin %d\n",N_BIN_PDF);
  fprintf(fp,"bins_edges \n");

  int i;
  for(i=0;i<N_BIN_PDF+1;i++){
    fprintf(fp, "%e\t",field->bins_edges[i]);
  }
  fprintf(fp,"\n\n");

  fprintf(fp,"mean        \t");
  fprintf(fp,"sigma       \t");
  fprintf(fp,"min         \t");
  fprintf(fp,"max         \t");
  fprintf(fp,"pdf         \t");
  fprintf(fp,"\n");
}

void writeFieldInfo(struct FIELD_INFO *field, FILE* fp, REAL t){
  char* format = "%e\t";
  fprintf(fp, format,t);
  fprintf(fp, format,field->mean);
  fprintf(fp, format,field->sigma);
  fprintf(fp, format,field->min);
  fprintf(fp, format,field->max);
  fprintf(fp,"\n");
}

void dumpStepInfoField(struct RUNPARAMS *param, char* field_name, struct FIELD_INFO *field, int nsteps, REAL t){

    mkdir("data/avg/", 0755);
    char filename[256];
    sprintf(filename,"data/avg/%s",field_name);

    FILE* fp;
    if (nsteps==0){
      fp=fopen(filename,"w");
      if(fp == NULL) printf("Cannot open %s\n", filename);
      //writeFieldInfoHeader(fp,field);
    }else{
      fp=fopen(filename,"a+");
      if(fp == NULL) printf("Cannot open %s\n", filename);
    }
    writeFieldInfo(field,fp,t);
    fclose(fp);
}

void getStepInfo(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu){

  const int debug =0;

  param->physical_state->sfr=0;
  param->physical_state->Nsn=0;
  param->physical_state->src=0;


#ifdef STARS
  double mass_star;

  if(param->stars->mass_res>100){
    mass_star=param->stars->mass_res;
  }
  else if(param->stars->mass_res>=0){
  #ifdef TESTCOSMO
    double munpercell=param->cosmo->om*(3.*pow(param->cosmo->H0*1e3/PARSEC/1e6,2)/(8.*M_PI*NEWTON_G))*pow(param->unit.unit_l/pow(2.0,param->lcoarse),3);
    double mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->stars->mass_res));
    mass_star =(REAL)( mstars_level * munpercell /SOLAR_MASS);
#else
  //TODO considere non cosmo case
  double mstars_level=0;
  mass_star =0;
#endif // TESTCOSMO
  }
  else{
    mass_star=1.;
  }
#endif

  double pre_mstar=(cpu->nsteps>0)? param->physical_state->mstar:0;
  param->physical_state->mstar=0;

  double pre_mstar_sfr=(cpu->nsteps>0)? param->physical_state->mstar_sfr:0;
  param->physical_state->mstar_sfr=0;

#ifdef TESTCOSMO
  double prev_t =(cpu->nsteps>0)? param->physical_state->t:0;
  param->physical_state->t = param->cosmo->tphy;
  double dt_yr = param->cosmo->tphy - prev_t;
#endif // TESTCOSMO
  int i;
  for (i=0;i<param->out_grid->n_field_tot; i++){
     if (param->out_grid->field_state_stat[i]){
        initFieldInfo(&(param->physical_state->field[i]), param->physical_state->field[i].bin_min, param->physical_state->field[i].bin_max);
    }
  }

  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){
    struct OCT *nextoct=firstoct[level-1];
    double vweight=POW(0.5,3*level); // volume d'une cellule de niveau level
    do{
      if(nextoct==NULL) continue;
      struct OCT *curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue;// ne pas compter les quantites des cellules de bord
      param->physical_state->max_level=level;
      int icell;
      for(icell=0;icell<8;icell++) {
        struct CELL *curcell = &curoct->cell[icell];
        if( curcell->child==NULL){

//------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------//

  for (i=0;i<param->out_grid->n_field_tot; i++){
    if (param->out_grid->field_state_stat[i]){
      getFieldInfo(&(param->physical_state->field[i]),assign_grid_field(i,curcell),vweight);
    }
  }

#ifdef WRAD
    int igrp;
	  for(igrp=0;igrp<NGRP;igrp++){
      param->physical_state->src+=curcell->rfield.src[igrp]*vweight;
    }
#endif // WRAD


#ifdef STARS
    struct PART *curp;
    struct PART *nexp=curcell->phead;
    do{
      if(nexp==NULL) continue;
      curp=nexp;
      nexp=curp->next;
      if((curp->isStar)&&(curp->isStar!=100)){
        param->physical_state->mstar += curp->mass;

#ifdef SUPERNOVAE
        if (param->sn->feedback_eff){
          if (param->cosmo->tphy - curp->age < param->sn->tlife){
            param->physical_state->mstar_sfr += curp->mass;
          }else{
            param->physical_state->mstar_sfr += curp->mass/(1.-param->sn->ejecta_proportion);
          }

        }else{
          param->physical_state->mstar_sfr += curp->mass;
        }
#else
          param->physical_state->mstar_sfr += curp->mass;
#endif // SUPERNOVAE

        if(curp->isStar==5||(curp->isStar==7||curp->isStar==8)){
          param->physical_state->Nsn++;
        }
      }
    }while(nexp!=NULL);
#endif // STARS

//------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------//
        }
      }
    }while(nextoct!=NULL);
  }

#ifdef WMPI
  for (i=0;i<param->out_grid->n_field_tot; i++){
    if (param->out_grid->field_state_stat[i]){
      comFieldInfo(cpu,&(param->physical_state->field[i]));
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,&param->physical_state->src,      1,MPI_DOUBLE,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&param->physical_state->max_level,1,MPI_INT,   MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&param->physical_state->Nsn,      1,MPI_INT,   MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&param->physical_state->mstar,    1,MPI_DOUBLE,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&param->physical_state->mstar_sfr,1,MPI_DOUBLE,MPI_SUM,cpu->comm);

#endif

#ifdef TESTCOSMO
  double dm_M0 = (param->physical_state->mstar_sfr - pre_mstar_sfr)*param->unit.unit_mass/SOLAR_MASS;
#endif // TESTCOSMO

  for (i=0;i<param->out_grid->n_field_tot; i++){
    if (param->out_grid->field_state_stat[i]){
      setSigmaFieldInfo(&(param->physical_state->field[i]));
    }
  }

  if (debug){
    int i;
    for (i=0;i<param->out_grid->n_field_tot; i++){
      if (param->out_grid->field_state_stat[i]){
        printFieldInfo(&(param->physical_state->field[i]));
      }
    }
  }

#ifdef TESTCOSMO
  REAL h=param->cosmo->H0/100.;
  REAL l= param->unit.unit_l/(1e6*PARSEC);
  REAL V_Mpc = POW(l,3);

  if (debug) printf("dm_M0%lf dt_yr=%lf V_Mpc=%lf\n",dm_M0,dt_yr,V_Mpc);

  param->physical_state->sfr = dm_M0/dt_yr/V_Mpc;

  if (debug) printf("param->physical_state->sfr = %lf\n",param->physical_state->sfr);
  if (debug) abort();
#endif // TESTCOSMO
}

void dumpStepInfo(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, int nsteps,REAL dt, REAL t){
/**
  * At each timestep, write information about the average state of several physical quantities
  */

  getStepInfo(firstoct, param, cpu);


  if(cpu->rank==RANK_DISP) printf("Dumping step info\n");

//  char* int_format = "%-8d\t";
  char* real_format = "%e\t";

  if(cpu->rank==RANK_DISP){

#ifdef GPUAXL
    char* filename = "data/param.avg.gpu";
#else
    char* filename = "data/param.avg.cpu";
#endif

    char* write_type;
    if (nsteps==0){
      write_type="w";
    }else{
      write_type="a+";
    }

    FILE* fp=fopen(filename,write_type);
    if(fp == NULL) printf("Cannot open %s\n", filename);

    if (nsteps==0){
      fprintf(fp,"step\t");

      fprintf(fp,"aexp\t\t");
      fprintf(fp,"z\t\t");
      fprintf(fp,"t_yrs\t\t");

      fprintf(fp,"dt\t\t");
      fprintf(fp,"max_level\t");

      fprintf(fp,"Nstars\t\t");
      fprintf(fp,"Mstars\t\t");
      fprintf(fp,"SFR\t\t");
      fprintf(fp,"SN\t\t");
      fprintf(fp,"src\t\t");
      fprintf(fp,"xion\t\t");
      fprintf(fp,"temp\t\t");
      fprintf(fp,"densb**2\t");
      fprintf(fp,"pressb\t\t");
      fprintf(fp,"Nagn\t\t");
      fprintf(fp,"dtff\t\t");
      fprintf(fp,"dthydro\t\t");
      fprintf(fp,"dtcosmo\t\t");
      fprintf(fp,"dtpic\t\t");
      fprintf(fp,"dtrad\t\t");

      fprintf(fp,"\n");
    }

    fprintf(fp, "%d\t",nsteps);
#ifdef TESTCOSMO
    fprintf(fp, real_format,(float)param->cosmo->aexp);
    fprintf(fp, real_format,(float)(1./param->cosmo->aexp-1.));
    fprintf(fp, real_format,(float)param->cosmo->tphy);
#else
    fprintf(fp, real_format,(float)0.);
    fprintf(fp, real_format,(float)0.);
    fprintf(fp, real_format,(float)0.);
#endif // TESTCOSMO

    fprintf(fp, real_format,(float)dt);
    fprintf(fp, real_format ,(float)param->physical_state->max_level);

#ifdef STARS

    int nstar;
#ifndef AGN
    nstar=param->stars->n;
#else
    nstar=param->stars->n-param->stars->nagn;
#endif
    fprintf(fp, real_format ,(float)nstar);
    fprintf(fp, real_format ,(float)param->physical_state->mstar);
    fprintf(fp, real_format ,(float)param->physical_state->sfr);
#else
    fprintf(fp, real_format ,(float)0.);
    fprintf(fp, real_format ,(float)0.);
    fprintf(fp, real_format ,(float)0.);
#endif // STARS

    fprintf(fp, real_format ,(float)param->physical_state->Nsn);
    fprintf(fp, real_format ,(float)param->physical_state->src);

    fprintf(fp, real_format ,(float)param->physical_state->field[32].mean); // mass weighted ionization fraction
    fprintf(fp, real_format ,(float)param->physical_state->field[38].mean); // temperature

    fprintf(fp, real_format ,(float)param->physical_state->field[13].sigma); // density
    fprintf(fp, real_format ,(float)param->physical_state->field[17].mean); // pressure

#ifdef AGN
    fprintf(fp, real_format ,(float)param->stars->nagn);
#else
    fprintf(fp, real_format ,(float)0);
#endif

    fprintf(fp, real_format ,(float)param->physical_state->dt_ff);
    fprintf(fp, real_format ,(float)param->physical_state->dt_hydro);
    fprintf(fp, real_format ,(float)param->physical_state->dt_cosmo);
    fprintf(fp, real_format ,(float)param->physical_state->dt_pic);
    fprintf(fp, real_format ,(float)param->physical_state->dt_rad);


    fprintf(fp,"\n");
    fclose(fp);
  }

  if(cpu->rank==RANK_DISP){
    int i;
    REAL tloc;
#ifdef TESTCOSMO
    tloc=(REAL)param->cosmo->aexp;
#else
    tloc=t;
#endif // TESTCOSMO
    printf("tloc=%e\n",tloc);

    for (i=0;i<param->out_grid->n_field_tot; i++){
      if (param->out_grid->field_state_stat[i]){
	dumpStepInfoField(param,param->out_grid->field_name[i], &(param->physical_state->field[i]), nsteps,tloc);

      }
    }

  }
}
