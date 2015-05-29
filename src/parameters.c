#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"

#ifdef UVBKG
#include "src_utils.h"
#endif // UVBKG



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

  sprintf(param,"param.output");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy(param_src, param_dest);
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

void GetParameters(char *fparam, struct RUNPARAMS *param)
{
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


  readOutputParam("param.output", param);


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
  dumpFile("param.output", "data/param.output");
  printf("\n");

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
    src+=curcell->rfield.src*vweight;
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


