#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "prototypes.h"

#ifdef UVBKG
#include "src_utils.h"
#endif // UVBKG

#define KPLANCK (6.62617e-34)	//J/s
#define ELECTRONVOLT (1.602176565e-19) //J


char *field_name [] ={
// The field order has to be the same as in param.output for consistency
"x",
"y",
"z",
"l",
"cpu",
"gdata.d",
"density",
"gdata.p",
"res",
"f[0]",
"f[1]",
"f[2]",
"marked",
"field.d",
"field.u",
"field.v",
"field.w",
"field.p",
"rfield.e[0]",
"rfield.fx[0]",
"rfield.fy[0]",
"rfield.fz[0]",
"rfield.e[1]",
"rfield.fx[1]",
"rfield.fy[1]",
"rfield.fz[1]",
"rfield.snfb",
"rfield.e[2]",
"rfield.fx[2]",
"rfield.fy[2]",
"rfield.fz[2]",
"rfield.src",
"xion",
"field.dX",
"field.dXHE",
"field.dXXHE",
"field.xHE",
"field.xxHE",
"rfield.temp"
};

int copy(char const * const source, char const * const destination) {
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
  copy(param_src, param_dest);

  sprintf(param,"param.mk");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy(param_src, param_dest);

  sprintf(param,"param.info");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy(param_src, param_dest);

  sprintf(param,"param.avg");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy(param_src, param_dest);

  sprintf(param,"param.h");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy(param_src, param_dest);
#ifdef ALLOCT
  sprintf(param,"param.output");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy(param_src, param_dest);
  #endif // ALLOCT
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

void readOutputParam(char *fparam, struct RUNPARAMS *param){

  int n_field=0;
  int n_field_tot=0;

  FILE *f=NULL;
  f=fopen(fparam,"r");
  if(f==NULL){
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
  }else{
    char stream[256];
    while (fscanf(f, "%s\n", stream) != EOF) {

      if ( strncmp( stream,"#", 1) ){
        param->out->field_id[n_field_tot] = 1;
        param->out->field_name[n_field_tot] = field_name[n_field_tot];
        n_field++;
      }else{
        param->out->field_id[n_field_tot] = 0;
      }
      n_field_tot++;
    }
    fclose(f);
  }

  param->out->n_field=n_field;
  param->out->n_field_tot=n_field_tot;

/*
  int i;
  for (i=0;i<n_field_tot; i++){
    if (param->out->field_id[i])
    printf("%d\t%s\n",param->out->field_id[i], param->out->field_name[i]);
  }
*/

//abort();

}

void readAtomic(struct RUNPARAMS *param){
/**
  * Read the atomic data from file
  */

  //openfile
  FILE *buf=NULL;
  buf=fopen(param->atomic.path,"r");
  if(buf==NULL){
    printf("ERROR : cannot open the parameter file (%s given), please check\n",param->atomic.path);
    abort();
  }

  size_t rstat;
  char stream[256];

  //alloc space
  rstat=fscanf(buf,"%s %d",stream,&param->atomic.ngrp_space);
 // printf("ngrp_space=%d\n", param->atomic.ngrp_space );
  param->atomic.space_bound= (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));

  //read space bound
  rstat=fscanf(buf,"%s",stream);
  int i_space;
  for (i_space=0; i_space<param->atomic.ngrp_space; i_space++){
    rstat=fscanf(buf,"%lf",&param->atomic.space_bound[i_space]);
    //printf("time space[%d]%e\n",i_space, param->atomic.space_bound[i_space]);
  }
  rstat=fscanf(buf,"%s",stream);

  //alloc time
  rstat=fscanf(buf,"%s %d",stream,&param->atomic.ngrp_time);
  //printf("ngrp_time=%d\n", param->atomic.ngrp_time );
  param->atomic.time_bound= (REAL*)calloc(param->atomic.ngrp_time,sizeof(REAL));

  //read time bound
  rstat=fscanf(buf,"%s",stream);
  int i_time;
  for (i_time=0; i_time<param->atomic.ngrp_time; i_time++){
    rstat=fscanf(buf,"%lf",&param->atomic.time_bound[i_time]);
    //printf("time bound[%d]%e\n",i_time, param->atomic.time_bound[i_time]);
  }

  // alloc grp
  param->atomic.n = param->atomic.ngrp_space * param->atomic.ngrp_time;
  param->atomic.hnu= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphae= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphai= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.factgrp= (REAL*)calloc(param->atomic.n,sizeof(REAL));

  //skip header
  rstat=fscanf(buf,"%s",stream);
  int i;
  for(i=0;i<4;i++){
    rstat=fscanf(buf,"%s",stream);
  }

  //read grp
  for(i=0;i<param->atomic.n;i++){
    rstat=fscanf(buf,"%lf",&param->atomic.hnu[i]);    param->atomic.hnu[i]*=1.6022e-19;
    rstat=fscanf(buf,"%lf",&param->atomic.alphae[i]); param->atomic.alphae[i]=param->clightorg*LIGHT_SPEED_IN_M_PER_S;
    rstat=fscanf(buf,"%lf",&param->atomic.alphai[i]); param->atomic.alphai[i]=param->clightorg*LIGHT_SPEED_IN_M_PER_S;
    rstat=fscanf(buf,"%lf",&param->atomic.factgrp[i]);
    //printf("%e %e %e %e \n", param->atomic.hnu[i], param->atomic.alphae[i], param->atomic.alphai[i], param->atomic.factgrp[i]);
  }
  fclose(buf);
}

REAL cross_section(REAL egy){
  if(egy >= 13.6){
    REAL P=2.963;
    REAL x=egy/4.298e-1;

    REAL a=POW(x-1.,2);
    REAL b=POW(x,P/2.-5.5);
    REAL c=POW(1+SQRT(x/3.288e1),-P);

    return 5.475e4*a*b*c*1e-22;

  }else {
    return 0;
  }
}

REAL eV2lambda(REAL eV){
/**
  * Convert electron Volt to wavelenght in meter
  */
	REAL h = KPLANCK ;
	REAL c = LIGHT_SPEED_IN_M_PER_S;
	REAL Ej = eV*ELECTRONVOLT;
	return h*c/Ej;
}

REAL lambda2E(REAL l){
  REAL h = KPLANCK ;
	REAL c = LIGHT_SPEED_IN_M_PER_S;
	return h*c/l;
}

void integ_dwave(struct RUNPARAMS *param, REAL* y, REAL space_bound_min, REAL space_bound_max, REAL *out){

  *out = 0;
  int ispace;
  for(ispace=0;ispace<param->spectrum.nwave;ispace++){
    if ((param->spectrum.wavelength[ispace]<=space_bound_min)
    &&  (param->spectrum.wavelength[ispace]>space_bound_max)) {

      REAL dx = param->spectrum.wavelength[ispace+1]-param->spectrum.wavelength[ispace];
      REAL y_mean= (y[ispace+1]+y[ispace])/2.;

      *out+=dx*y_mean;
    }
  }
}

void integ_dtime(struct RUNPARAMS *param, REAL* y, REAL time_bound_min, REAL time_bound_max, REAL *out){

  *out=0;
  int i_time;
  for(i_time=0;i_time<param->spectrum.ntime;i_time++){
    if ((param->spectrum.time[i_time]>time_bound_min)
    &&  (param->spectrum.time[i_time]<=time_bound_max)){

      REAL dx = param->spectrum.time[i_time+1]-param->spectrum.time[i_time];
      REAL y_mean = (y[i_time+1]+y[i_time])/2.;

      printf("dx=%e dy=%e out=%e\n", dx, y_mean, *out);

      *out+=dx*y_mean;
    }
  }
}

struct GRP{
  REAL *s_N;
  REAL *s_E;
  REAL *Epho;
  REAL *nphot;
  REAL *factgrp;
};

void get_space_GRP(struct RUNPARAMS *param, int i_time, struct GRP *grp){

  grp->s_N =  (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  grp->s_E =  (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  grp->Epho = (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  grp->nphot= (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  grp->factgrp= (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));

  REAL *s_lambda = (REAL*)calloc(param->spectrum.nwave,sizeof(REAL));
  REAL *N_lambda = (REAL*)calloc(param->spectrum.nwave,sizeof(REAL));
  REAL *tmp_space= (REAL*)calloc(param->spectrum.nwave,sizeof(REAL));

  REAL nphot_tot=0;

  int i_grp_space;
  for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){

    REAL space_bound_min=eV2lambda(param->atomic.space_bound[i_grp_space]);
    REAL space_bound_max=eV2lambda(param->atomic.space_bound[i_grp_space+1]);

    //Eion
    REAL Eion = 0;
    integ_dwave(param, param->spectrum.flux[i_time], space_bound_min, space_bound_max, &Eion);

    //Nphot_per_sec
    REAL Nphot_per_sec=0;
    int ispace;
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      REAL E = lambda2E(param->spectrum.wavelength[ispace]);
      tmp_space[ispace]=param->spectrum.flux[i_time][ispace]/E;
    }
    integ_dwave(param, tmp_space, space_bound_min, space_bound_max, &Nphot_per_sec);

    //s_lambda
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      REAL E = lambda2E(param->spectrum.wavelength[ispace]);
      s_lambda[ispace] = cross_section(E/ELECTRONVOLT);
    }

    //N_lambda
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      REAL E = lambda2E(param->spectrum.wavelength[ispace]);
      N_lambda[ispace] = param->spectrum.flux[i_time][ispace]/E;
    }

    //s_N
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      tmp_space[ispace]=s_lambda[ispace]*N_lambda[ispace];
    }
    integ_dwave(param, tmp_space, space_bound_min, space_bound_max, &grp->s_N[i_grp_space]);
    grp->s_N[i_grp_space]/=Nphot_per_sec;

    //s_E
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      tmp_space[ispace]=s_lambda[ispace]*param->spectrum.flux[i_time][ispace];
    }
    integ_dwave(param, tmp_space, space_bound_min, space_bound_max, &grp->s_E[i_grp_space]);
    grp->s_E[i_grp_space]/=Eion;

    grp->Epho[i_grp_space]= Eion/Nphot_per_sec;
    grp->nphot[i_grp_space] = Nphot_per_sec/(1e6*SOLAR_MASS);

    nphot_tot+=grp->nphot[i_grp_space];
  }

  for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){
    grp->factgrp[i_grp_space]=grp->nphot[i_grp_space]/nphot_tot;
  }

  free(s_lambda);
  free(N_lambda);
  free(tmp_space);
}

void get_time_GRP(struct RUNPARAMS *param, struct GRP *grp){

  REAL *tmp_time= (REAL*)calloc(param->spectrum.ntime,sizeof(REAL));


  int i_grp_time;
  for(i_grp_time=0;i_grp_time<param->atomic.ngrp_time;i_grp_time++){

    REAL time_bound_min=param->atomic.time_bound[i_grp_time];
    REAL time_bound_max=param->atomic.time_bound[i_grp_time+1];
    REAL nphot_tot=0;

    int i_grp_space;
    for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){

      int id= i_grp_time*param->atomic.ngrp_space+i_grp_space;
      //printf("id=%d \n", id);

      //Nphot
      int itmp;
      for(itmp=0;itmp<param->spectrum.ntime;itmp++){
        tmp_time[itmp]= grp[itmp].nphot[i_grp_space];
  //      printf("grp[%d].nphot[%d]=%e \n", itmp,i_grp_space,tmp_time[itmp]);
      }
      REAL Nphot=0;
      integ_dtime(param, tmp_time, time_bound_min, time_bound_max, &Nphot);
      printf("Nphot=%e \n", Nphot);

//      nphot_tot+=Nphot;

      //E_tot
      for(itmp=0;itmp<param->spectrum.ntime;itmp++){
        tmp_time[itmp]=grp[itmp].Epho[i_grp_space]*grp[itmp].nphot[i_grp_space];
      }
      integ_dtime(param, tmp_time, time_bound_min, time_bound_max, &param->atomic.hnu[id]);
      param->atomic.hnu[id]/=Nphot;

//    printf("hnu=%e \n", param->atomic.hnu[id]/ELECTRONVOLT);

      //s_N_tot
      for(itmp=0;itmp<param->spectrum.ntime;itmp++){
        tmp_time[itmp]=grp[itmp].s_N[i_grp_space]*grp[itmp].nphot[i_grp_space];
      }
      integ_dtime(param, tmp_time, time_bound_min, time_bound_max, &param->atomic.alphai[id]);
      param->atomic.alphai[id]/=Nphot;
  //   printf("ai=%e \n", param->atomic.alphai[id]);

      //s_E_tot
      for(itmp=0;itmp<param->spectrum.ntime;itmp++){
        tmp_time[itmp]=grp[itmp].s_E[i_grp_space]*grp[itmp].nphot[i_grp_space];
      }
      integ_dtime(param, tmp_time, time_bound_min, time_bound_max, &param->atomic.alphae[id]);
      param->atomic.alphae[id]/=Nphot;
    //  printf("ae=%e \n", param->atomic.alphae[id]);

    }


//			factgrp = np.sum(nphot[igrp])/np.sum(nphot)



    int i_grp_time;
    for(i_grp_time=0;i_grp_time<param->atomic.ngrp_time;i_grp_time++){
      REAL time_bound_min=param->atomic.time_bound[i_grp_time];
      REAL time_bound_max=param->atomic.time_bound[i_grp_time+1];

      REAL nphot_tot=0;

      for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){

        int id= i_grp_time*param->atomic.ngrp_space+i_grp_space;

        int i_time;
        for(i_time=0;i_time<param->spectrum.ntime;i_time++){
          nphot_tot+=grp[i_time].nphot[i_grp_space];
        }
        for(i_time=0;i_time<param->spectrum.ntime;i_time++){
          if ((param->spectrum.time[i_time]>time_bound_min)
          &&  (param->spectrum.time[i_time]<time_bound_max)){
            param->atomic.factgrp[id] += grp[i_time].nphot[i_grp_space]/nphot_tot;
          }
        }
      }
    }
  }
}


void setAtomic(struct RUNPARAMS *param){
/**
  * reduce spectrum to mean energy and cross section, depending of NGRP
  */

  //  param->atomic.ngrp_space=1;   REAL space_bound[] = {13.6,INFINITY};
  param->atomic.ngrp_space=3;   REAL space_bound[] = {13.6, 24.587387, 54.417760,INFINITY};
  //param->atomic.ngrp_time=1;    REAL time_bound[] = {0.,INFINITY};
  param->atomic.ngrp_time=3;    REAL time_bound[] = {0.,10,100,INFINITY};

  int i_time;
  for(i_time=0;i_time<param->atomic.ngrp_time;i_time++){
    time_bound[i_time]*=1e6;
  }

  param->atomic.n = param->atomic.ngrp_space * param->atomic.ngrp_time;
  param->atomic.space_bound= (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  param->atomic.time_bound= (REAL*)calloc(param->atomic.ngrp_time,sizeof(REAL));

  param->atomic.space_bound= space_bound;
  param->atomic.time_bound= time_bound;

  struct GRP *grp = (struct GRP*)calloc(param->spectrum.ntime,sizeof(struct GRP));

  for(i_time=0;i_time<param->spectrum.ntime;i_time++){
    get_space_GRP(param, i_time,&grp[i_time]);
    //printf("\n");
  }
/*
  for(i_time=0;i_time<param->spectrum.ntime;i_time++){
    int i_grp_space;
    for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){

      printf("hnu=%e \t",grp[i_time].Epho[i_grp_space] /ELECTRONVOLT);
      printf("sN=%e \t", grp[i_time].s_N[i_grp_space]);
      printf("sE=%e \t", grp[i_time].s_E[i_grp_space]);
      printf("fact=%e \n", grp[i_time].factgrp[i_grp_space]);

    }
  }
*/
  //alloc atomic
  param->atomic.hnu=    (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphae= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphai= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.factgrp=(REAL*)calloc(param->atomic.n,sizeof(REAL));

  get_time_GRP(param, grp);

  int i_grp_time;
  for(i_grp_time=0;i_grp_time<param->atomic.ngrp_time;i_grp_time++){
      int i_grp_space;
      for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){
        int id= i_grp_time*param->atomic.ngrp_space+i_grp_space;
        printf("hnu=%e \t", param->atomic.hnu[id]/ELECTRONVOLT);
        printf("sE=%e \t", param->atomic.alphae[id]);
        printf("sN=%e \t", param->atomic.alphai[id]);
        printf("fact=%e \n", param->atomic.factgrp[id]);
      }
      printf("\n");
    }
}




void readStarburst99(struct RUNPARAMS *param){
/**
  * Read the Starburst99 spectrum database
  * http://www.stsci.edu/science/starburst99/data/fig1e.dat
  */

  sprintf(param->spectrum.path,"./src/atomic_data/Starburst99_spectrum.dat");

  FILE *buf=fopen(param->spectrum.path,"r");
  if(buf==NULL){
    printf("ERROR : cannot open the parameter file (%s given), please check\n",param->spectrum.path);
    abort();
  }

  //skip human readable header
  int nmax = 2048;
  char stream[nmax];
  fgets(stream, nmax, buf);
  fgets(stream, nmax, buf);
  fgets(stream, nmax, buf);

  //read header
  size_t rstat;
  rstat=fscanf(buf,"%s %d",stream,&param->spectrum.ntime);
  rstat=fscanf(buf,"%s %d",stream,&param->spectrum.nwave);
  printf("%d \n", param->spectrum.ntime);printf("%d \n", param->spectrum.nwave);

  //alloc
  param->spectrum.time=(REAL*)calloc(param->spectrum.ntime,sizeof(REAL));
  param->spectrum.wavelength=(REAL*)calloc(param->spectrum.nwave,sizeof(REAL));

  param->spectrum.flux=(REAL**)calloc(param->spectrum.ntime,sizeof(REAL*));
  int iwave;
  for(iwave=0;iwave<param->spectrum.nwave;iwave++)
    param->spectrum.flux[iwave]=(REAL*)calloc(param->spectrum.nwave,sizeof(REAL));

  //read time
  int itime;
  for(itime=0;itime<param->spectrum.ntime;itime++)
    rstat=fscanf(buf,"%lf",&param->spectrum.time[itime]);

  //read wavelenght and flux
  for(iwave=0;iwave<param->spectrum.nwave; iwave++){
    rstat=fscanf(buf,"%lf",&param->spectrum.wavelength[iwave]);
    for(itime=0;itime<param->spectrum.ntime;itime++){
      rstat=fscanf(buf,"%lf",&param->spectrum.flux[itime][iwave]);
    }
  }

  //set time (Myr in yr)
  for(itime=0;itime<param->spectrum.ntime;itime++){
    param->spectrum.time[itime]*=1e6;
  }

  //set Wavelenght (Angstrom in meter)
  for(iwave=0;iwave<param->spectrum.nwave; iwave++){
    param->spectrum.wavelength[iwave]*=1e-10;
  }

  //set flux (log(Erg/A) in J/m/s)
  for(iwave=0;iwave<param->spectrum.nwave; iwave++){
    for(itime=0;itime<param->spectrum.ntime;itime++){
      //printf("%e\n",param->spectrum.flux[itime][iwave] );
      param->spectrum.flux[itime][iwave] = POW(10.,param->spectrum.flux[itime][iwave]) * 1e-7/1e-10;
      //printf("%e\n",param->spectrum.flux[itime][iwave] );
    }
  }
}

void GetParameters(char *fparam, struct RUNPARAMS *param){
  FILE *buf=NULL;
  char stream[256];
  size_t rstat;
  double dummyf;
  char RF[]="%s %lf";
  int i;
  buf=fopen(fparam,"r");
  if(buf==NULL)
    {
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
    }
  else
    {
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ngridmax);
      rstat=fscanf(buf,"%s %d",stream,&param->npartmax);
      rstat=fscanf(buf,"%s %d",stream,&param->nbuff);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ndumps);
      rstat=fscanf(buf,RF,stream,&param->dt_dump);
      rstat=fscanf(buf,"%s %d",stream,&param->nsteps);
      rstat=fscanf(buf,RF,stream,&dummyf);param->dt=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->tmax=(REAL)dummyf;

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->lcoarse);
      rstat=fscanf(buf,"%s %d",stream,&param->lmax);
      rstat=fscanf(buf,RF,stream,&dummyf);param->amrthresh0=(REAL)dummyf;
      rstat=fscanf(buf,"%s %d",stream,&param->nsmooth);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->DM_res);
      rstat=fscanf(buf,RF,stream,&param->dx_res);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->niter);
      rstat=fscanf(buf,RF,stream,&dummyf);param->poissonacc=(REAL)dummyf;
      rstat=fscanf(buf,"%s %d",stream,&param->mgridlmin);
      if(param->mgridlmin<0){
        param->mgridlmin=param->lcoarse-param->lcoarse;
      }

      rstat=fscanf(buf,"%s %d",stream,&param->nvcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nrelax);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->nrestart);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->gstride);
      rstat=fscanf(buf,"%s %d",stream,&param->hstride);
      rstat=fscanf(buf,"%s %d",stream,&param->nsubcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nthread);
      rstat=fscanf(buf,"%s %d",stream,&param->nstream);
      rstat=fscanf(buf,"%s %d",stream,&param->ompthread);

      rstat=fscanf(buf,"%s",stream);
#ifdef WRAD
      rstat=fscanf(buf,RF,stream,&dummyf);param->clight=(REAL)dummyf;param->clightorg=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->denthresh=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->tmpthresh=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->srcint=(REAL)dummyf;
      rstat=fscanf(buf,"%s %s",stream, param->atomic.path);
      param->fudgecool=1.0;
      param->ncvgcool=0;
#else
	for (i=0; i<4; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
#endif

      rstat=fscanf(buf,"%s",stream);
#ifdef STARS
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->overdensity_cond=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->density_cond=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->efficiency=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->tlife=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->mass_res=(REAL)dummyf;

#else
	for (i=0; i<5; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
#endif

      rstat=fscanf(buf,"%s",stream);
#ifdef SUPERNOVAE
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->feedback_eff	=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->feedback_frac	=(REAL)dummyf;
#else
	for (i=0; i<2; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
#endif

#ifdef MOVIE
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->movie->lmap);
      rstat=fscanf(buf,RF,stream,&param->movie->xmin);
      rstat=fscanf(buf,RF,stream,&param->movie->xmax);
      rstat=fscanf(buf,RF,stream,&param->movie->ymin);
      rstat=fscanf(buf,RF,stream,&param->movie->ymax);
      rstat=fscanf(buf,RF,stream,&param->movie->zmin);
      rstat=fscanf(buf,RF,stream,&param->movie->zmax);
#endif
      fclose(buf);
    }




  // computing the maxhash
  int val=(POW(2,param->lmax-1)<256?POW(2,param->lmax-1):256); // limit to 2097152 octs in hash table i.e. 16e6 cells
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
#endif

#ifdef SRCINT
  param->srcint*=SRCINT;
#endif

#ifdef UVBKG
  setUVBKG(param, "src/phys_data/uvbkg.dat");
#endif // UVBKG

#ifdef ALLOCT
  readOutputParam("param.output", param);
#endif // ALLOCT

  //readAtomic(param);
  readStarburst99(param);
  setAtomic(param);

  abort();
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
    fprintf(fp, float_format,"box_size_Mpc/h",(param->unit.unit_l/PARSEC/1e6*param->cosmo->H0/100));
    fprintf(fp, int_format,"level_min",(param->lcoarse) );
    fprintf(fp, int_format,"level_max",(param->lmax) );

  fprintf(fp,"##=Unit_code->SI====================\n" );
  #ifdef WRAD
    fprintf(fp, real_format,"unit_l",(param->unit.unit_l) );		// comoving length size of the box [meters]
    fprintf(fp, real_format,"unit_v",(param->unit.unit_v) );		// unit velocity
    fprintf(fp, real_format,"unit_t",(param->unit.unit_t) );		// unit time [seconds]
    fprintf(fp, real_format,"unit_n",(param->unit.unit_n) );		// unit number [moles typically]
    fprintf(fp, real_format,"unit_mass",(param->unit.unit_mass) );	// unit mass [in kg, total mass is equal to one in unit codes]
  //  fprintf(fp,"\n");
  #endif

  fprintf(fp,"##=Cosmology========================\n" );
  #ifdef TESTCOSMO
    fprintf(fp, real_format,"om",(param->cosmo->om) );			// Omega matter
    fprintf(fp, real_format,"ov",(param->cosmo->ov) );			// Omega vacuum
    fprintf(fp, real_format,"ob",(param->cosmo->ob) );			// Omega baryon
    fprintf(fp, real_format,"H0",(param->cosmo->H0) );			// Hubble constant
  //  fprintf(fp,"\n");
  #endif



    fprintf(fp,"##=Mass_resolution_(Mo)=============\n" );
    REAL mass_res_DM =  (1.- param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse+param->DM_res))*param->unit.unit_mass/SOLAR_MASS;
    fprintf(fp, real_format,"mass_res_DM",mass_res_DM );
#ifdef STARS
    REAL res = param->stars->mass_res;
    if(res>100){
        fprintf(fp, real_format,"mass_res_star",param->stars->mass_res);
    }else{

      if(res>=0){
        REAL mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse+param->stars->mass_res));
        REAL mass_res_star = mstars_level * param->unit.unit_mass /SOLAR_MASS;
        fprintf(fp, real_format,"mass_res_star",mass_res_star);
      }else{
        int level;
        for(level=param->lcoarse;level<=param->lmax;level++){
          REAL mlevel=level-1;
          REAL restmp=-param->stars->mass_res;
          REAL mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(mlevel+restmp));
          REAL mass_res_star = mstars_level * param->unit.unit_mass /SOLAR_MASS;
          char mlev[128];
          sprintf(mlev,"mass_star_L%d",level);
          fprintf(fp, real_format,mlev,mass_res_star );
        }
      }
    }
#endif // STARS

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

  printf("Dumping parameter files \n\n");

  dumpInfo("data/param.info", param, cpu);
  printf("\n");
  dumpFile("param.mk", "data/param.mk");
  printf("\n");
  dumpFile("param.h", "data/param.h");
  printf("\n");
  dumpFile(fparam, "data/param.run");
  printf("\n");
#ifdef ALLOCT
  dumpFile("param.output", "data/param.output");
  printf("\n");
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

void dumpStepInfo(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, int nsteps,REAL dt){
/**
  * At each timestep, compute and write information about the average state of several physical quantities
  */

  if(cpu->rank==RANK_DISP) printf("Dumping step info\n");

  int max_level=0;
  int Nsn=0;
  REAL mean_xion=0;
  REAL mean_T=0;
  REAL max_T=0;
  REAL max_rho=0;
  REAL src=0;

  int level;
  REAL vweight;
  for(level=param->lcoarse;level<=param->lmax;level++){
    struct OCT *nextoct=firstoct[level-1];
    vweight=POW(0.5,3*level); // volume d'une cellule de niveau level
    do{
      if(nextoct==NULL) continue;
      struct OCT *curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue;// ne pas compter les quantites des cellules de bord
      max_level=level;
      int icell;
      for(icell=0;icell<8;icell++) {
	struct CELL *curcell = &curoct->cell[icell];
	if( curcell->child==NULL){
	  // ajout d'une ponderation en volume, plus conforme a l'idee de moyenne
	  // note somme(vweight)=1.
#ifdef WRADHYD
	  mean_xion+=curcell->field.dX/curcell->field.d*vweight;
#endif // WRADHYD

#ifdef WRAD
	  mean_T+=curcell->rfield.temp*vweight;

    int igrp;
	  for(igrp=0;igrp<NGRP;igrp++){
      src+=curcell->rfield.src[igrp]*vweight;
    }

	  max_T=FMAX(max_T,curcell->rfield.temp);
#endif // WRAD
	  max_rho=FMAX(max_rho,curcell->field.d);


#ifdef STARS
    struct PART *curp;
    struct PART *nexp=curcell->phead;
    do{
      if(nexp==NULL) continue;
      curp=nexp;
      nexp=curp->next;
      if(curp->isStar==2||curp->isStar==3){
        Nsn++;
      }
        //------------------------------------------------//

    }while(nexp!=NULL);
#endif // STARS
	}
      }
    }while(nextoct!=NULL);
  }

#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&mean_xion,1,MPI_REEL,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&mean_T,1,MPI_REEL,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&max_T,1,MPI_REEL,MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&max_rho,1,MPI_REEL,MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&src,1,MPI_REEL,MPI_SUM,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&max_level,1,MPI_INT,MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&Nsn,1,MPI_INT,MPI_MAX,cpu->comm);
#endif

//  char* int_format = "%-8d\t";
  char* real_format = "%e\t";

  if(cpu->rank==RANK_DISP){

    char* filename = "data/param.avg";

    FILE* fp=NULL;

    if (nsteps==0){
      fp=fopen(filename,"w");
      if(fp == NULL) printf("Cannot open %s\n", filename);
      fprintf(fp,"step\t");
      fprintf(fp,"aexp\t\t");
      fprintf(fp,"z\t\t");
      fprintf(fp,"dt\t\t");
      fprintf(fp,"max_level\t");
      fprintf(fp,"max_rho\t\t");
      fprintf(fp,"mean_xion\t");
      fprintf(fp,"mean_T\t\t");
      fprintf(fp,"max_T\t\t");
      fprintf(fp,"stars\t\t");
      fprintf(fp,"SN\t\t");
      fprintf(fp,"src");
      fprintf(fp,"\n");

    }else{
      fp=fopen(filename,"a+");
      if(fp == NULL) printf("Cannot open %s\n", filename);
    }

    fprintf(fp, "%d\t",nsteps);
    fprintf(fp, real_format,param->cosmo->aexp);
    fprintf(fp, real_format,1./param->cosmo->aexp-1.);
    fprintf(fp, real_format,dt);
    fprintf(fp, real_format ,(float)max_level);
    fprintf(fp, real_format,max_rho);

    fprintf(fp, real_format,mean_xion);
    fprintf(fp, real_format,mean_T);
    fprintf(fp, real_format,max_T);

#ifdef STARS
    fprintf(fp, real_format ,(float)param->stars->n);
#else
    fprintf(fp, real_format ,0.);
#endif // STARS

    fprintf(fp, real_format ,(float)Nsn);
    fprintf(fp, real_format , src);

    fprintf(fp,"\n");
    fclose(fp);
  }
}
