#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include <sys/types.h>
#include <sys/stat.h>

#include "friedmann.h"
#include "segment.h"
#include "stars.h"

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
#endif
#endif

#ifdef WRAD
  int igrp;
  for(igrp=0;igrp<NGRP;igrp++){
    lcell->e[igrp]=cell->rfield.e[igrp];
    lcell->fx[igrp]=cell->rfield.fx[igrp];
    lcell->fy[igrp]=cell->rfield.fy[igrp];
    lcell->fz[igrp]=cell->rfield.fz[igrp];
  }
  lcell->src=cell->rfield.src;
#ifdef STARS
  lcell->snfb=cell->rfield.snfb;
#endif
  lcell->xion=cell->rfield.nhplus/cell->rfield.nh;
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


void dumpFile(char *filename_in, char *filename_out){

  FILE *fps[2] = {stdout, NULL};
  fps[1]=fopen(filename_out,"w");
  if(fps[1] == NULL) printf("Cannot open %s\n", filename_out);

  FILE* buf=NULL;
  buf=fopen(filename_in,"r");
  if(buf == NULL) printf("Cannot open %s\n", filename_in);

  int i;
  for(i=0;i<2;i++){
    FILE *fp = fps[i];
    char ch;
    fseek(buf,0,SEEK_SET);
    while((ch=fgetc(buf))!=EOF) fprintf(fp,"%c",ch);
  }

  fclose(fps[1]);
  fclose(buf);
}

void dumpInfo(char *filename_info, struct RUNPARAMS *param, struct CPUINFO *cpu){

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
    REAL mass_res_DM =  ( 1.- param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse))*param->unit.unit_mass /1.98e30;
    fprintf(fp, real_format,"mass_res_DM",mass_res_DM );

    REAL res = param->stars->mass_res;
    if(res>=0){
      REAL mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse+param->stars->mass_res));
      fprintf(fp, real_format,"mass_res_star",mstars_level);
    }else{
      int level;
      for(level=param->lcoarse;level<=param->lmax;level++){
        REAL mlevel=0;
        if(res>=0){mlevel=param->lcoarse;}
        else{mlevel=level+1;}
        REAL mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(mlevel+res));

        REAL mass_res_star = mstars_level * param->unit.unit_mass /1.98e30;
        char mlev[128];
        sprintf(mlev,"mass_star_L%d",level);
        fprintf(fp, real_format,mlev,mass_res_star );
      }
    }

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
  printf("Dumping parameter files \n\n");

  dumpInfo("data/param.info", param, cpu);
  dumpFile("param.mk", "data/param.mk");
  dumpFile(fparam, "data/param.run");

  printf("\n");
  //abort();
}

void dumpStepInfo(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, int nsteps,REAL dt){

  if(cpu->rank==RANK_DISP) printf("Dumping step info\n");

  int ncell=0;
  int max_level=0;

  REAL mean_xion=0;
  REAL mean_T=0;
  REAL max_T=0;
  REAL max_rho=0;

  REAL sfr=0;


  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){
	struct OCT *nextoct=firstoct[level-1];
    do{if(nextoct==NULL) continue;
      max_level=level;
      struct OCT *curoct=nextoct;
      nextoct=curoct->next;
      int icell;
      for(icell=0;icell<8;icell++) {
        struct CELL *curcell = &curoct->cell[icell];
        if( curcell->child==NULL){
          mean_xion+=curcell->field.dX/curcell->field.d;
          mean_T+=curcell->rfield.temp;

          max_T=FMAX(max_T,curcell->rfield.temp);
          max_rho=FMAX(max_rho,curcell->field.d);
        }
      }
      ncell+=8;
    }while(nextoct!=NULL);
  }

#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&mean_xion,1,MPI_REEL,MPI_SUM,cpu->comm); mean_xion/=cpu->nproc*ncell;
  MPI_Allreduce(MPI_IN_PLACE,&mean_T,1,MPI_REEL,MPI_SUM,cpu->comm); mean_T/=cpu->nproc*ncell;

  MPI_Allreduce(MPI_IN_PLACE,&max_T,1,MPI_REEL,MPI_MAX,cpu->comm);
  MPI_Allreduce(MPI_IN_PLACE,&max_rho,1,MPI_REEL,MPI_MAX,cpu->comm);
//  MPI_Allreduce(MPI_IN_PLACE,&sfr,1,MPI_REEL,MPI_SUM,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&ncell,1,MPI_INT,MPI_SUM,cpu->comm);
	MPI_Allreduce(MPI_IN_PLACE,&max_level,1,MPI_INT,MPI_MAX,cpu->comm);
#endif



  char* int_format = "%-8d\t";
  char* real_format = "%e\t";

  if(cpu->rank==RANK_DISP){

    char* filename = "data/param.avg";

    FILE* fp=NULL;

    if (nsteps==0){
      fp=fopen(filename,"w");
      if(fp == NULL) printf("Cannot open %s\n", filename);
      fprintf(fp,"step\taexp\t\tz\t\tmax_level\tmax_rho\t\tmean_xion\tmean_T\t\tmax_T\t\tstars\n");

    }else{
      fp=fopen(filename,"a+");
      if(fp == NULL) printf("Cannot open %s\n", filename);
    }

    fprintf(fp, "%d\t",nsteps);
    fprintf(fp, real_format,param->cosmo->aexp);
    fprintf(fp, real_format,1./param->cosmo->aexp-1.);
    fprintf(fp, int_format ,max_level);
    fprintf(fp, real_format,max_rho);

    fprintf(fp, real_format,mean_xion);
    fprintf(fp, real_format,mean_T);
    fprintf(fp, real_format,max_T);

    fprintf(fp, int_format ,param->stars->n);

//   fprintf(fp, real_format,sfr);


    fprintf(fp,"\n");
    fclose(fp);

  }
}

//====================================================================================================
void save_amr(char filename[], struct OCT **firstoct,REAL tsim, REAL tinit,int nsteps, int ndumps, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt){

  FILE *fp = NULL;
  int level;
  struct OCT * nextoct;
  struct OCT * root;
  struct CELL * rootcell;
  struct PART * rootpart;
  struct OCT *oct;
  int ic=0;


  root=firstoct[0];
  rootcell=&(root->cell[0]);
  rootpart=proot;

  fp=fopen(filename,"wb");
  if(fp == NULL) printf("Cannot open %s\n", filename);

  fwrite(&tsim,sizeof(REAL),1,fp);
  fwrite(&tinit,sizeof(REAL),1,fp);
  fwrite(&nsteps,sizeof(int),1,fp);
  fwrite(&ndumps,sizeof(int),1,fp);
  fwrite(&(param->lcoarse),sizeof(int),1,fp);
  fwrite(&(param->lmax),sizeof(int),1,fp);

  // dumping timesteps
  for(level=1;level<=param->lmax;level++){
    fwrite(adt+level-1,sizeof(REAL),1,fp);
  }

#ifdef WRAD
  fwrite(&(param->unit.unit_l),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_v),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_t),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_n),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_mass),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_d),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_N),sizeof(REAL),1,fp);
#endif


#ifdef TESTCOSMO
  fwrite(&(param->cosmo->om),sizeof(REAL),1,fp);
  fwrite(&(param->cosmo->ov),sizeof(REAL),1,fp);
  fwrite(&(param->cosmo->ob),sizeof(REAL),1,fp);
  fwrite(&(param->cosmo->H0),sizeof(REAL),1,fp);
#endif



  // writing pointer informations
  fwrite(&root,sizeof(struct OCT*),1,fp);
  fwrite(&rootcell,sizeof(struct CELL*),1,fp);
  fwrite(&rootpart,sizeof(struct PART*),1,fp);


  for(level=1;level<=param->lmax;level++) // looping over octs
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];

      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=nextoct;
	  nextoct=oct->next;

	  fwrite(&oct,sizeof(struct OCT*),1,fp);
	  fwrite(oct,sizeof(struct OCT),1,fp);
	  ic++;
	}while(nextoct!=NULL);
    }


  //printf("%d octs dumped by proc %d\n",ic,cpu->rank);


  fclose(fp);

}



//====================================================================================================
struct OCT * restore_amr(char filename[], struct OCT **firstoct,struct OCT **lastoct, REAL *tsim, REAL *tinit, int *nsteps, int *ndumps,struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt, struct CELL *root){

  FILE *fp = NULL;
  int level,lcoarse,lmax;
  struct OCT * curoct;
  struct OCT * nextoct;

  struct OCT * root_sna;
  struct CELL * rootcell_sna;
  struct PART * rootpart_sna;

  struct OCT * root_mem;
  struct CELL * rootcell_mem;
  struct PART * rootpart_mem;


  struct OCT oct;
  struct OCT * oct_ad;
  int ic,ioct=0;

  struct OCT *freeoct=NULL;



  root_mem=firstoct[0]; // the root oct of the grid in memory
  rootcell_mem=&(root_mem->cell[0]); // the root cell of the grid in memory
  rootpart_mem=proot; // the root cell of the grid in memory

  // reset global pointers

  for(level=1;level<=param->lmax;level++){
    firstoct[level-1]=NULL;
    lastoct[level-1]=NULL;
  }


  // opening the file
  fp=fopen(filename,"rb");
  if(fp == NULL) printf("Cannot open %s\n", filename);

  size_t outf;

  // reading snapshot time
  outf=fread(tsim,sizeof(REAL),1,fp);
  outf=fread(tinit,sizeof(REAL),1,fp);
  outf=fread(nsteps,sizeof(int),1,fp);
  outf=fread(ndumps,sizeof(int),1,fp);
  outf=fread(&lcoarse,sizeof(int),1,fp);
  outf=fread(&lmax,sizeof(int),1,fp);

  for(level=1;level<=lmax;level++){
    outf=fread(adt+level-1,sizeof(REAL),1,fp);
  }

#ifdef WRAD
  // reading units
  outf=fread(&(param->unit.unit_l),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_v),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_t),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_n),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_mass),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_d),sizeof(REAL),1,fp);
  outf=fread(&(param->unit.unit_N),sizeof(REAL),1,fp);
#endif
  //printf("UNIT L=%e\n",param->unit.unit_l);

#ifdef TESTCOSMO
  outf=fread(&(param->cosmo->om),sizeof(REAL),1,fp);
  outf=fread(&(param->cosmo->ov),sizeof(REAL),1,fp);
  outf=fread(&(param->cosmo->ob),sizeof(REAL),1,fp);
  outf=fread(&(param->cosmo->H0),sizeof(REAL),1,fp);
#endif

  // reading pointer informations in the snapshot
  outf=fread(&root_sna,sizeof(struct OCT*),1,fp);
  outf=fread(&rootcell_sna,sizeof(struct CELL*),1,fp);
  outf=fread(&rootpart_sna,sizeof(struct PART*),1,fp);


  if(cpu->rank==RANK_DISP) printf(" STARTING OCT READ%p %p %p\n",root_sna,rootcell_sna,rootpart_sna);
  // reading the octs sequence

  outf=fread(&oct_ad,sizeof(struct OCT *),1,fp);
  outf=fread(&oct,sizeof(struct OCT),1,fp);

  while(!feof(fp)){

    // do stuff
    ioct++;
    // 1 copy the content of the oct at the right location
    curoct=root_mem+(oct_ad-root_sna);
    //printf("cpu=%d ioct=%d curoct=%p  oct-ad=%p root_sna=%p dif=%ld lev=%d chi=%p\n",cpu->rank,ioct,curoct,oct_ad,root_sna,(unsigned long int) (oct_ad-root_sna),oct.level,oct.cell[0].child);
    //if(curoct-root_mem>param->ngridmax) printf("ERROR BKP\n");
    memcpy(curoct,&oct,sizeof(struct OCT));

    //if(cpu->rank==RANK_DISP) printf("le=%d ic=%d\n",curoct->level,ic);
    // 2.a modify the oct pointers within curoct

    curoct->next=(curoct->next==NULL?NULL:(curoct->next-root_sna)+root_mem);
    curoct->prev=(curoct->prev==NULL?NULL:(curoct->prev-root_sna)+root_mem);
    curoct->nexthash=(curoct->nexthash==NULL?NULL:(curoct->nexthash-root_sna)+root_mem);

    // 2.c modifity the oct pointers within cells
    for(ic=0;ic<8;ic++){
      curoct->cell[ic].child=(curoct->cell[ic].child==NULL?NULL:(curoct->cell[ic].child-root_sna)+root_mem);
    }


    // 2.b modify the cell pointers within curoct

    for(ic=0;ic<6;ic++){

      if(curoct->nei[ic]!=NULL){

	curoct->nei[ic]=(struct CELL *)(((unsigned long long int)(curoct->nei[ic])-(unsigned long long int)rootcell_sna)+(unsigned long long int)rootcell_mem);
      }


    }

    if(curoct->parent!=NULL){

	struct CELL *co;
	if(curoct->level>1){
	  co=(struct CELL *)((unsigned long long int)curoct->parent-(unsigned long long int) rootcell_sna+(unsigned long long int) rootcell_mem);
	  curoct->parent=co;
	}
	else{
	  curoct->parent=root;
	}

    }


#ifdef PIC
    // 2.c modifity the particle pointers within cells

    for(ic=0;ic<8;ic++){
      curoct->cell[ic].phead=(curoct->cell[ic].phead==NULL?NULL:(curoct->cell[ic].phead-rootpart_sna)+rootpart_mem);
    }
#endif

    // Overall pointer management
    if(firstoct[curoct->level-1]==NULL){
      if(curoct->prev==NULL){
	firstoct[curoct->level-1]=curoct;
      }
    }

    if(lastoct[curoct->level-1]==NULL){
      if(curoct->next==NULL){
	lastoct[curoct->level-1]=curoct;
      }
    }

    // read next oct

    outf=fread(&oct_ad,sizeof(struct OCT *),1,fp);
    outf=fread(&oct,sizeof(struct OCT),1,fp);

  }

  //printf("%d octs recovered by proc %d with root=%p firstoct=%p\n",ioct,cpu->rank,root_mem,firstoct[0]);

  // done
  fclose(fp);


  // Building back the freeoct chained list

  struct OCT * loct;
  for(curoct=root_mem;curoct<root_mem+param->ngridmax;curoct++){
    if(curoct->level==0){
      if(freeoct==NULL){
	freeoct=curoct;
	curoct->prev=NULL;
	curoct->next=NULL;
	loct=curoct;
      }
      else{
	curoct->next=NULL;
	curoct->prev=loct;
	loct->next=curoct;
	loct=curoct;
      }
    }
  }

  //abort();

  /* // check and count the number of freeocts */

  /* nextoct=freeoct; */
  /* int ifree=0; */
  /* while(nextoct!=NULL){ */
  /*   curoct=nextoct; */
  /*   ifree++; */
  /*   nextoct=curoct->next; */
  /* } */

  /* printf("ifree=%d\n",ifree); */



  return freeoct;

}


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


//====================================================================================================
//=================================================================================================

  //------------------------------------------------------------------------

  //------------------------------------------------------------------------



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

//=======================================================================================================

void save_part(char filename[],struct OCT **firstoct, int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct PART* proot){

  FILE *fp = NULL;
  int vali;
  int ipart=0;
  int level;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  struct PART *nexp;
  struct PART *curp;
  int icell;


  int npart=0;
  for(level=levelcoarse;level<=levelmax;level++) npart+=cpu->npart[level-1];

  fp=fopen(filename,"wb");
	if(fp == NULL) printf("Cannot open %s\n", filename);

  fwrite(&npart,1,sizeof(int),fp);
  fwrite(&tsim,1,sizeof(REAL),fp);
  fwrite(&proot,1,sizeof(struct PART *),fp);

//	printf("cpu %d \t%d\t%e\t%p \n", cpu->rank, npart, tsim,proot );

  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];
      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./POW(2,oct.level);
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
		do{
		  curp=nexp;
		  nexp=curp->next;

		  fwrite(&curp,1,sizeof(struct PART *),fp);
		  fwrite(curp,1,sizeof(struct PART),fp);

		  ipart++;
		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }

  fclose(fp);
  //printf("wrote %d particles (%d expected) in %s on proc %d\n",ipart,npart,filename,cpu->rank);

}




// ===================================================================================================
// ===================================================================================================
// ===================================================================================================

struct PART * restore_part(char filename[], struct OCT **firstoct, REAL *tsim, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot){

  FILE *fp = NULL;
  int level;

  struct PART * rootpart_sna;
  struct PART * rootpart_mem;

  rootpart_mem=proot;

  int ic,ipart=0;
  int npart;
  int sp;

  struct PART part;
  struct PART *part_ad;
  struct PART *curp;

#ifdef STARS
  int nstar=0;
#endif

  rootpart_mem=proot; // the root cell of the grid in memory

  // opening the file
  fp=fopen(filename,"rb");
	if(fp == NULL) printf("Cannot open %s\n", filename);

  // reading snapshot time
  size_t outf;
  outf=fread(&npart,1,sizeof(int),fp);
  outf=fread(tsim,1,sizeof(REAL),fp);
  outf=fread(&rootpart_sna,1,sizeof(struct PART *),fp);

//	printf("cpu %d \t%d\t%e\t%p \n", cpu->rank, npart, *tsim,rootpart_sna);

  // reading the particle sequence

  outf=fread(&part_ad,sizeof(struct PART*),1,fp);
  outf=fread(&part,sizeof(struct PART),1,fp);


 // printf("debut de lecture %ld\t%p\t%p\t%p\n ", (unsigned long long int)(part_ad-rootpart_sna), part_ad,rootpart_sna,rootpart_mem );

  while(!feof(fp)){

    /* if(cpu->rank==RANK_DISP){ */
    /*   printf("ipart=%d\n",ipart); */
    /* } */
    // do stuff
    ipart++;

    // 1 copy the content of the particle at the right location

    curp=(part_ad-rootpart_sna)+rootpart_mem;

    memcpy(curp,&part,sizeof(struct PART));

#ifdef STARS
    if(curp->isStar)  nstar++;
#endif

 //   printf("memcpy OK \n");

    // 2.a modify the particle pointers
    curp->next=(curp->next==NULL?NULL:(curp->next-rootpart_sna)+rootpart_mem);
    curp->prev=(curp->prev==NULL?NULL:(curp->prev-rootpart_sna)+rootpart_mem);

  //  printf("*part ok \n");

    // read next particle
    outf=fread(&part_ad,sizeof(struct PART*),1,fp);
    outf=fread(&part,sizeof(struct PART),1,fp);

//    printf("%d\t ", ipart);
  }

//  printf("READ OK \n ");
  // Building back the freepart chained list

  struct PART * lpart;
  struct PART *freepart;
  freepart=NULL;

  for(curp=rootpart_mem;curp<rootpart_mem+param->npartmax;curp++){
    if(curp->mass==-1.){ // flag empty particles
      if(freepart==NULL){
	freepart=curp;
	curp->prev=NULL;
	curp->next=NULL;
	lpart=curp;
      }
      else{
	curp->next=NULL;
	curp->prev=lpart;
	lpart->next=curp;
	lpart=curp;
      }
    }
    /* else if(cpu->rank==RANK_DISP){ */
    /*   printf("%e\n",curp->mass); */
    /* } */
  }

  //printf("%d/%d part recovered by proc %d with freepart=%p\n",ipart,param->npartmax,cpu->rank,freepart);

  // done
  fclose(fp);

#ifdef STARS
#ifdef WMPI
	MPI_Allreduce(MPI_IN_PLACE,&nstar,1,MPI_INT,   MPI_SUM,cpu->comm);
#endif
	//printf("nstar=%d\n",nstar);
  param->stars->n=nstar;
#endif

  return freepart;
}

#endif

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

//------------------------------------------------------------------------

void GetParameters(char *fparam, struct RUNPARAMS *param)
{
  FILE *buf;
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
      rstat=fscanf(buf,"%s %d",stream,&param->niter);
      rstat=fscanf(buf,RF,stream,&dummyf);param->poissonacc=(REAL)dummyf;
      rstat=fscanf(buf,"%s %d",stream,&param->mgridlmin);
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
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->tcar=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->tlife=(REAL)dummyf;
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->mass_res=(REAL)dummyf;

#else
	for (i=0; i<4; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
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

}

//==================================================================================
//==================================================================================
#ifdef PIC
#ifdef GRAFIC
struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param,int level)
{
  FILE *fx = NULL;
  FILE *fy = NULL;
  FILE *fz = NULL;

  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
  int dummy;
  struct PART *lastpart;
  int ip;
  char filename[256];

  if(cpu->rank==0){
    sprintf(filename,"./level_%03d/ic_velbx",level);
    fx=fopen(filename,"rb"); 	if(fx == NULL) printf("Cannot open %s\n", filename);
    sprintf(filename,"./level_%03d/ic_velby",level);
    fy=fopen(filename,"rb"); 	if(fy == NULL) printf("Cannot open %s\n", filename);
    sprintf(filename,"./level_%03d/ic_velbz",level);
    fz=fopen(filename,"rb");	if(fz == NULL) printf("Cannot open %s\n", filename);


    // reading the headers
    size_t outf;

    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&np1,1,4,fx);
    outf=fread(&np2,1,4,fx);
    outf=fread(&np3,1,4,fx);
    outf=fread(&dx,1,4,fx);
    printf("DX=%e\n",dx);
    outf=fread(&x1o,1,4,fx);
    outf=fread(&x2o,1,4,fx);
    outf=fread(&x3o,1,4,fx);
    outf=fread(&astart,1,4,fx);
    outf=fread(&om,1,4,fx);
    outf=fread(&ov,1,4,fx);
    outf=fread(&h0,1,4,fx);
    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&dummy,1,sizeof(dummy),fy);
    outf=fread(&np1,1,4,fy);
    outf=fread(&np2,1,4,fy);
    outf=fread(&np3,1,4,fy);
    outf=fread(&dx,1,4,fy);
    outf=fread(&x1o,1,4,fy);
    outf=fread(&x2o,1,4,fy);
    outf=fread(&x3o,1,4,fy);
    outf=fread(&astart,1,4,fy);
    outf=fread(&om,1,4,fy);
    outf=fread(&ov,1,4,fy);
    outf=fread(&h0,1,4,fy);
    outf=fread(&dummy,1,sizeof(dummy),fy);

    outf=fread(&dummy,1,sizeof(dummy),fz);
    outf=fread(&np1,1,4,fz);
    outf=fread(&np2,1,4,fz);
    outf=fread(&np3,1,4,fz);
    outf=fread(&dx,1,4,fz);
    outf=fread(&x1o,1,4,fz);
    outf=fread(&x2o,1,4,fz);
    outf=fread(&x3o,1,4,fz);
    outf=fread(&astart,1,4,fz);
    outf=fread(&om,1,4,fz);
    outf=fread(&ov,1,4,fz);
    outf=fread(&h0,1,4,fz);
    outf=fread(&dummy,1,sizeof(dummy),fz);
  }

#ifdef WMPI
  // Massive broadcast of grafic header
  MPI_Bcast(&np1,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np2,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np3,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&dx,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x1o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x2o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x3o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&astart,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&om,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&ov,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&h0,1,MPI_FLOAT,0,cpu->comm);
  MPI_Barrier(cpu->comm);
#endif

  if(cpu->rank==0){
    printf("============================================\n");
    printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
    printf("om=%f ov=%f h0=%f\n",om,ov,h0);
    printf("dx=%f np1*dx=%f\n",dx,np1*dx);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }

  if(level==param->lcoarse){
  if((np1*np2*np3)/(cpu->nproc)*1.>param->npartmax){
    printf("Error : Number of particles greater than npartmax Np(grafic, est. )=%d Npartmax=%d!\n",(np1*np2*np3)/(cpu->nproc),param->npartmax);
    abort();
  }
  }
  //setting omegab

  ob=OMEGAB;


  // computing Zeldovich displacement quantities

  double vfact;
  vfact=fomega(astart,om,ov)*h0*dladt(astart,om,ov)/astart;
  if(cpu->rank==0) printf("vfact=%f\n",vfact);
  // reading the grafic planes

  float *velx;
  float *vely;
  float *velz;

  double x;
  double y;
  double z;

  double vx;
  double vy;
  double vz;

  double x0,y0,z0;

  int i1,i2,i3;
  int offidx=0;
  int keep;
  double mass;
  double veloff=0.;

#ifdef WHYDRO2
  mass=(1.-ob/om)/(np1*np2*np3);
#else
  mass=1./(np1*np2*np3);
#endif

#ifdef ZOOM
  if(level>param->lcoarse){
    offidx=POW(2,3*(level-param->lcoarse)); // for ids of particles
  }
#endif



  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);

  //REAL rmin=2.;

  ip=0;
  size_t outf;
  for(i3=1;i3<=np3;i3++){

    if(cpu->rank==0){
      printf("\r %f percent done",(i3*1.0)/np3*100.);
      outf=fread(&dummy,1,sizeof(dummy),fx);
      outf=fread(velx,np1*np2,sizeof(float),fx);
      outf=fread(&dummy,1,sizeof(dummy),fx);

      outf=fread(&dummy,1,sizeof(dummy),fy);
      outf=fread(vely,np1*np2,sizeof(float),fy);
      outf=fread(&dummy,1,sizeof(dummy),fy);

      outf=fread(&dummy,1,sizeof(dummy),fz);
      outf=fread(velz,np1*np2,sizeof(float),fz);
      outf=fread(&dummy,1,sizeof(dummy),fz);
    }


#ifdef WMPI
    MPI_Barrier(cpu->comm);
    MPI_Bcast(velx,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(vely,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velz,np1*np2,MPI_FLOAT,0,cpu->comm);
#endif


    z0=(i3-0.5)*dx;
    for(i2=1;i2<=np2;i2++){
      y0=(i2-0.5)*dx;
      for(i1=1;i1<=np1;i1++){
	x0=(i1-0.5)*dx;
	// computing the displacements
	x=(x0+velx[(i1-1)+(i2-1)*np1]/vfact)/(np1*dx);
	y=(y0+vely[(i1-1)+(i2-1)*np1]/vfact)/(np2*dx);
	z=(z0+velz[(i1-1)+(i2-1)*np1]/vfact)/(np3*dx);

	// periodic boundary conditions

	x+=((x<0.)-(x>1.))*1.; 
	y+=((y<0.)-(y>1.))*1.; 
	z+=((z<0.)-(z>1.))*1.; 

	// ugly fix for huge config in SINGLE FLOAT precision
	// generally affects a tiny fraction of particle (like 1 over 1e7)
 	if(x>0.99999) x=0.;
	if(y>0.99999) y=0.;
	if(z>0.99999) z=0.;

	// computing the velocities
	vx=(velx[(i1-1)+(i2-1)*np1]+veloff/astart*1e-3)*astart/(np1*dx*h0)/(sqrt(om)*0.5);
	vy=(vely[(i1-1)+(i2-1)*np1])*astart/(np2*dx*h0)/(sqrt(om)*0.5);
	vz=(velz[(i1-1)+(i2-1)*np1])*astart/(np3*dx*h0)/(sqrt(om)*0.5);

	// if it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,cpu,cpu->levelcoarse)){

	  keep=1;
#ifdef ZOOM
	  // is the current particle at the correct level?
	  int lzoom;
	  lzoom=pos2levelzoom(x,y,z,param);
	  REAL rloc;
	  rloc=sqrt(((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)));

	  if(lzoom!=level){
	    keep=0;
	  }
 #endif

	  if(keep) {
 	    part[ip].x=x;
	    part[ip].y=y;
	    part[ip].z=z;

 	    //rmin=(rloc<rmin?rloc:rmin);

	    part[ip].vx=vx;
	    part[ip].vy=vy;
	    part[ip].vz=vz;
	    part[ip].level=level;

	    part[ip].mass=mass;
	    part[ip].idx=(i1-1)+(i2-1)*np1+(i3-1)*np1*np2+offidx;
	    lastpart=part+ip;
	    ip++;
	  }
	}
      }
    }
  }

  free(velx);
  free(vely);
  free(velz);

  if(cpu->rank==0){
    fclose(fx);
    fclose(fy);
    fclose(fz);
  }

  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=ip;
#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif

  if(cpu->rank==RANK_DISP){
    printf("Grafic Particle Read ok\n");
  }
  return lastpart;
}
#endif

// ========================== ZELDOVICH
// ====================================
// ====================================

#ifdef ZELDOVICH
struct PART * read_zeldovich_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param, struct OCT **firstoct)
{
  FILE *fd;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
  int dummy;
  struct PART *lastpart;
  int ip;

  int nploc;
  float munit_z;
  float lbox;
  float ainit_z;
  REAL mass;
  size_t outf;

  fd=fopen("utils/grafic_src/ZEL.PM.0","rb");

  outf=fread(&dummy,sizeof(dummy),1,fd);
  outf=fread(&nploc,sizeof(int),1,fd);
  outf=fread(&munit_z,sizeof(float),1,fd);
  outf=fread(&ainit_z,sizeof(float),1,fd);
  outf=fread(&lbox,sizeof(float),1,fd);
  outf=fread(&om,sizeof(float),1,fd);
  outf=fread(&ov,sizeof(float),1,fd);
  outf=fread(&h0,sizeof(float),1,fd);
  outf=fread(&dummy,sizeof(dummy),1,fd);

  astart=ainit_z;


  if(cpu->rank==0){
    printf("============================================\n");
    printf("ntot%d\n",nploc);
    printf("om=%f ov=%f h0=%f\n",om,ov,h0);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }

  if((nploc)/(cpu->nproc)*1.>param->npartmax){
    printf("Error : Number of particles greater than npartmax Np(grafic, est. )=%d Npartmax=%d!\n",(nploc)/(cpu->nproc),param->npartmax);
    abort();
  }

  //setting omegab

  ob=OMEGAB;


#ifdef WHYDRO2
  mass=(1.-ob/om)/(nploc);
#else
  mass=1./(nploc);
#endif


  // reading the grafic planes
  float *pos;
  float *vel;
  int nread=nploc; // quick fixes
  int npatch=1.; // quick fixes
  int ipatch;
  int i;
  pos=(float *)malloc(sizeof(REAL)*3*nread);
  vel=(float *)malloc(sizeof(REAL)*3*nread);

  double x;
  double y;
  double z;

  double vx;
  double vy;
  double vz;
  size_t rstat;

  int pstart=ftell(fd);

  ip=0.;
  for(ipatch=0;ipatch<npatch;ipatch++) {
    //    rstat=outf=fread(&dummy,sizeof(dummy),1,fd);
    //    fseek(fd,pstart,SEEK_SET);
    fseek(fd,pstart+(0*nploc+ipatch*nread)*sizeof(float)+1*sizeof(dummy),SEEK_SET);
    outf=fread(pos,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(1*nploc+ipatch*nread)*sizeof(float)+3*sizeof(dummy),SEEK_SET);
    outf=fread(pos+nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(2*nploc+ipatch*nread)*sizeof(float)+5*sizeof(dummy),SEEK_SET);
    outf=fread(pos+2*nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(3*nploc+ipatch*nread)*sizeof(float)+7*sizeof(dummy),SEEK_SET);
    outf=fread(vel,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(4*nploc+ipatch*nread)*sizeof(float)+9*sizeof(dummy),SEEK_SET);
    outf=fread(vel+nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(5*nploc+ipatch*nread)*sizeof(float)+11*sizeof(dummy),SEEK_SET);
    outf=fread(vel+2*nread,sizeof(float),nread,fd);
    //outf=fread(&dummy,sizeof(dummy),1,fd);


    for(i=0;i<nread;i++)
      {
	x=pos[i];
	y=pos[i+nread];
	z=pos[i+2*nread];
	vx=vel[i];
	vy=vel[i+nread];
	vz=vel[i+2*nread];
	// periodic boundary conditions

	x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x);
	y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y);
	z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z);
	// it it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,cpu,cpu->levelcoarse)){

	  part[ip].x=x;
	  part[ip].y=y;
	  part[ip].z=z;

	  part[ip].vx=vx;
	  part[ip].vy=vy;
	  part[ip].vz=vz;

	  part[ip].mass=mass;
	  part[ip].idx=i;
	  lastpart=part+ip;
	  ip++;
	}
      }
  }


  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=ip;

  if(cpu->rank==RANK_DISP){
    printf("Zeldovich Particle Read ok\n");
  }

#ifdef WHYDRO2

  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;
  REAL ZC=10; // hard coded collapse of redshift
  REAL ZI=1./astart-1.;

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
 		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		W.d=(1.+(1.+ZC)/(1.+ZI)*cos(2.*M_PI*(xc-0.5)))*ob/om;
		W.p=PMIN;
		W.u=-(1.+ZC)/POW(1.+ZI,1.5)*sin(2.*M_PI*(xc-0.5))/(M_PI); // for omegam=1. only
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);
		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

#endif

  return lastpart;
}
#endif



// =====================================================================================
// =====================================================================================

#ifdef EDBERT
struct PART * read_edbert_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param, struct OCT **firstoct)
{
  float astart,om,ov,h0,ob;
  int dummy;

  struct PART *lastpart;
  int ip,iploc;
  int i,j,k;

  float lbox;
  float ainit_z;
  REAL mass;
  REAL x,y,z,r;
  REAL delta;
  REAL lsphere;
  printf("Start EDBERT\n");

  lsphere=0.05;
  om=0.99999;
  ov=0.00001;
  ob=OMEGAB;
  delta=0.2;
  h0=70.;
  lbox=1.;
  astart=1e-3;

  int n1d=POW(2,param->lcoarse);
  REAL dx=1./n1d;
  iploc=0;
  int nin=0;
  int nout=0;
  REAL m;

#ifdef PIC

  REAL mout,mint;
  REAL vsphere=0.;
  for(k=0;k<n1d;k++)
    {
      for(j=0;j<n1d;j++)
	{
	  for(i=0;i<n1d;i++)
	    {
	      x=(i+0.5)*dx;
	      y=(j+0.5)*dx;
	      z=(k+0.5)*dx;

	      r=sqrt(POW(x-0.5,2)+POW(y-0.5,2)+POW(z-0.5,2));
	      if(r<lsphere){
		nin++;
		mout=-1;
		vsphere+=POW(dx,3);
	      }
	      else{
		nout++;
		mout=1;
	      }

	      if(segment_part(x,y,z,cpu,cpu->levelcoarse)){
		part[iploc].x=x;
		part[iploc].y=y;
		part[iploc].z=z;

		part[iploc].vx=0;
		part[iploc].vy=0;
		part[iploc].vz=0;

		part[iploc].mass=mout;
		part[iploc].idx=-nin;
		lastpart=part+iploc;
		iploc++;
	      }
	    }
	}
    }

  mint=(om-ob)*(1.+delta)*vsphere/nin;
  mout=((om-ob)-mint*nin)/nout;
  printf("mint=%e mout=%e\n",mint,mout);

  for(i=0;i<iploc;i++){
    part[i].mass=(part[i].mass<0?mint:mout);
  }

#endif

  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=iploc;

  if(cpu->rank==RANK_DISP){
    printf("Edbert Particle Read ok\n");
  }

#ifdef WHYDRO2

  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		rad=sqrt((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5));
		if(rad<lsphere){
		  W.d=(1.+delta)*ob;
		}
		else{
		  W.d=ob*(1.-(1.+delta)*vsphere)/(1.-vsphere);
		}
		W.p=PMIN;
		W.u=0.; // vstar is expressed in m/s and grafic vel are in km/s
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);

		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

#endif

  return lastpart;
}
#endif


#endif
//==================================================================================
//==================================================================================

#ifdef WHYDRO2

#ifdef TUBE
// =====================================================================================================
// =====================================================================================================

void read_shocktube(struct CPUINFO *cpu, REAL *ainit, struct RUNPARAMS *param, struct OCT **firstoct)
{
  FILE *fd = NULL;

  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;


  struct Wtype WL, WR;
  if(cpu->rank==RANK_DISP) printf("Init Hydro\n");

  /* /\*  /\\* // TEST 1 *\\/ *\/ */

  WL.d=1.;
  WL.u=0.;
  WL.v=0.;
  WL.w=0.;
  WL.p=1.0;
  WL.a=sqrt(GAMMA*WL.p/WL.d);
  getE(&WL);

  WR.d=0.125;
  WR.u=0.;
  WR.v=0.;
  WR.w=0.;
  WR.p=0.1;
  WR.a=sqrt(GAMMA*WR.p/WR.d);
  getE(&WR);

  REAL X0=0.3125;

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		if(xc<X0){
		  memcpy(&(curoct->cell[icell].field),&WL,sizeof(struct Wtype));
		}
		else{
		  memcpy(&(curoct->cell[icell].field),&WR,sizeof(struct Wtype));
		}
	      }
	  }while(nextoct!=NULL);
      }
}
#endif

#ifdef EVRARD
int read_evrard_hydro(struct CPUINFO *cpu,struct OCT **firstoct, struct RUNPARAMS *param){

  int level;
  REAL dxcur;
  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  struct Wtype W;
  REAL rad;
  REAL xc,yc,zc;

  //==== parameters of evrard sphere
  REAL R=0.35;
  REAL M=1.;
  REAL rhostar=M/(4./3.*M_PI*R*R*R);
  REAL estar=M/R; //assuming G=1
  REAL pstar=rhostar*estar;
  REAL tstar=sqrt(M_PI*M_PI/8.)*POW(R,1.5)/POW(M,0.5);
  if(cpu->rank==RANK_DISP) printf("Generating Evrard Test Case ts=%e, rhostar=%e\n",tstar,rhostar);

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		rad=sqrt((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5))/R;
		if(rad<1.){
		  W.d=rhostar/rad;
		  W.p=pstar/rad*0.05;
		}
		else{
		  W.d=1e-3;
		  W.p=1e-5;
		}

		W.u=0.; // vstar is expressed in m/s and grafic vel are in km/s
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);

		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

  }
#endif


#ifdef TESTCOSMO
#ifdef GRAFIC
 int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, struct RUNPARAMS *param,int level){

  FILE *fx = NULL;
  FILE *fy = NULL;
  FILE *fz = NULL;
  FILE *fdx = NULL;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,ob,h0;
  int dummy;
  int ip;
  struct Wtype W;
  size_t outf;
  char filename[256];

  // Note only the rank 0 reads the file.

  if(cpu->rank==0){

    sprintf(filename,"./level_%03d/ic_deltab",level);
    fdx=fopen(filename,"rb");		if(fdx == NULL) printf("Cannot open %s\n", filename);
    sprintf(filename,"./level_%03d/ic_velbx",level);
    fx=fopen(filename,"rb");		if(fx == NULL) printf("Cannot open %s\n", filename);
    sprintf(filename,"./level_%03d/ic_velby",level);
    fy=fopen(filename,"rb");		if(fy == NULL) printf("Cannot open %s\n", filename);
    sprintf(filename,"./level_%03d/ic_velbz",level);
    fz=fopen(filename,"rb");		if(fz == NULL) printf("Cannot open %s\n", filename);

    /* fdx=fopen("./ic_deltab","rb"); */
    /* fx=fopen("./ic_velcx","rb"); */
    /* fy=fopen("./ic_velcy","rb"); */
    /* fz=fopen("./ic_velcz","rb"); */

    // reading the headers

    outf=fread(&dummy,1,sizeof(dummy),fdx);
    outf=fread(&np1,1,4,fdx);
    outf=fread(&np2,1,4,fdx);
    outf=fread(&np3,1,4,fdx);
    outf=fread(&dx,1,4,fdx);
    outf=fread(&x1o,1,4,fdx);
    outf=fread(&x2o,1,4,fdx);
    outf=fread(&x3o,1,4,fdx);
    outf=fread(&astart,1,4,fdx);
    outf=fread(&om,1,4,fdx);
    outf=fread(&ov,1,4,fdx);
    outf=fread(&h0,1,4,fdx);
    outf=fread(&dummy,1,sizeof(dummy),fdx);

    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&np1,1,4,fx);
    outf=fread(&np2,1,4,fx);
    outf=fread(&np3,1,4,fx);
    outf=fread(&dx,1,4,fx);
    outf=fread(&x1o,1,4,fx);
    outf=fread(&x2o,1,4,fx);
    outf=fread(&x3o,1,4,fx);
    outf=fread(&astart,1,4,fx);
    outf=fread(&om,1,4,fx);
    outf=fread(&ov,1,4,fx);
    outf=fread(&h0,1,4,fx);
    outf=fread(&dummy,1,sizeof(dummy),fx);

    outf=fread(&dummy,1,sizeof(dummy),fy);
    outf=fread(&np1,1,4,fy);
    outf=fread(&np2,1,4,fy);
    outf=fread(&np3,1,4,fy);
    outf=fread(&dx,1,4,fy);
    outf=fread(&x1o,1,4,fy);
    outf=fread(&x2o,1,4,fy);
    outf=fread(&x3o,1,4,fy);
    outf=fread(&astart,1,4,fy);
    outf=fread(&om,1,4,fy);
    outf=fread(&ov,1,4,fy);
    outf=fread(&h0,1,4,fy);
    outf=fread(&dummy,1,sizeof(dummy),fy);

    outf=fread(&dummy,1,sizeof(dummy),fz);
    outf=fread(&np1,1,4,fz);
    outf=fread(&np2,1,4,fz);
    outf=fread(&np3,1,4,fz);
    outf=fread(&dx,1,4,fz);
    outf=fread(&x1o,1,4,fz);
    outf=fread(&x2o,1,4,fz);
    outf=fread(&x3o,1,4,fz);
    outf=fread(&astart,1,4,fz);
    outf=fread(&om,1,4,fz);
    outf=fread(&ov,1,4,fz);
    outf=fread(&h0,1,4,fz);
    outf=fread(&dummy,1,sizeof(dummy),fz);
  }

  // setting baryon density parameter
  ob=OMEGAB;

#ifdef WMPI
  // Massive broadcast of grafic header
  MPI_Bcast(&np1,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np2,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np3,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&dx,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x1o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x2o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x3o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&astart,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&om,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&ov,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&h0,1,MPI_FLOAT,0,cpu->comm);
  MPI_Barrier(cpu->comm);

#endif

  if(cpu->rank==0){
    printf("============================================\n");
    printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
    printf("om=%f ov=%f ob=%f h0=%f\n",om,ov,ob,h0);
    printf("dx=%f np1*dx=%f\n",dx,np1*dx);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }



  if(np1!=(int)POW(2,level)){
    printf("ERROR !ABORT! Grafic hydro  file not compliant with parameter file : ngrafic=%d nquartz=%d\n",np1,(int)POW(2,level));
    abort();
  }


  // reading the grafic planes

  float *deltab;
  float *velz;
  float *vely;
  float *velx;

  int i1,i2,i3;
  int icx,icy,icz,icell;
  unsigned long long key;
  struct OCT *curoct;
  struct OCT *nextoct;
  unsigned long hidx;
  int found;
  float z0,y0,x0;
  int ifound=0;


  deltab=(float*)malloc(sizeof(float)*np1*np2);
  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);

  double rhob,pressure;
  double H0=h0*1e3/3.08568025e22; // Hubble constant (s-1)
  double rhoc=3.*H0*H0/(8.*M_PI*NEWTON_G); // comoving critical density (kg/m3)
  double zstart=1./astart-1.;

  // ---------- Setting the initial temperature ---- //
  double temp;
  //  double temp=550.*((1.0+zstart)*(1.0+zstart)); // baryon temperature (to check) in K
  //double temp=2.7*(1+zstart);
  //double temp=1e4;
  //double temp=170.*(1.+zstart)*(1.+zstart)/10000.;

  //double temp=0.0874545+0.0302621*zstart+0.00675076*zstart*zstart; // recfast ob fit
#ifdef COOLING
  if(om==1.) {
    temp=33.64/POW(41.,2)*POW(1.+zstart,2);
    if(cpu->rank==RANK_DISP) printf("WARNING: YOU ARE USING SCDM COSMOLOGY\n");
  }
  else{
    if(cpu->rank==RANK_DISP) printf("No temperature law for cosmologies other than SCDM -> F** it\n");
    temp=33.64/POW(41.,2)*POW(1.+zstart,2);
    //    abort();
  }
#else
  temp=1e4;
#endif

  // supercomoving unit values
  double rhostar;
  double rstar;
  double vstar;
  double tstar;
  double tstar2;
  double pstar;
  double mpc=3.08568025e22; // Mpc in m
  double veloff=0.;
#ifdef BULKFLOW
  veloff=VBC; // relative motion in km/s at z=999
#endif

  rstar= np1*dx*mpc; // box size in m
  rhostar=rhoc*om;
  tstar=2./H0/sqrt(om); // sec
  tstar2=2./h0/sqrt(om); // Mpc sec / km
  vstar=rstar/tstar; //m/s
  pstar=rhostar*vstar*vstar;

  if(cpu->rank==RANK_DISP) printf("rhoc=%e temperature=%lf rstar=%e(%e) pstar=%e tstar=%e vstar=%e rhostar=%e\n",rhoc,temp,rstar,np1*dx,pstar,tstar,vstar,rhostar);


  for(i3=0;i3<np3;i3++){

    if(cpu->rank==0){
      printf("\r %f percent done",(i3*1.0)/np3*100.);

      outf=fread(&dummy,1,sizeof(dummy),fdx);
      outf=fread(deltab,np1*np2,sizeof(float),fdx);
      outf=fread(&dummy,1,sizeof(dummy),fdx);

      outf=fread(&dummy,1,sizeof(dummy),fx);
      outf=fread(velx,np1*np2,sizeof(float),fx);
      outf=fread(&dummy,1,sizeof(dummy),fx);

      outf=fread(&dummy,1,sizeof(dummy),fy);
      outf=fread(vely,np1*np2,sizeof(float),fy);
      outf=fread(&dummy,1,sizeof(dummy),fy);

      outf=fread(&dummy,1,sizeof(dummy),fz);
      outf=fread(velz,np1*np2,sizeof(float),fz);
      outf=fread(&dummy,1,sizeof(dummy),fz);
    }

#ifdef WMPI
    MPI_Barrier(cpu->comm);
    MPI_Bcast(deltab,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velx,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(vely,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velz,np1*np2,MPI_FLOAT,0,cpu->comm);
#endif

    z0=(i3*1.0)/(np3);
    for(i2=0;i2<np2;i2++){
      y0=(i2*1.0)/(np2);
      for(i1=0;i1<np1;i1++){
	x0=(i1*1.0)/(np1);

	key=pos2key(x0,y0,z0,level);

	// first we compute the adress from the hashfunction
	hidx=hfun(key,cpu->maxhash);
	nextoct=cpu->htable[hidx];

	// looking for the oct
	found=0;
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,level)==key)&&(curoct->level==level));
	  }while((nextoct!=NULL)&&(!found));
	}

	// filling the cell

	if(found){
	  icx=i1%2;
	  icy=i2%2;
	  icz=i3%2;

	  icell=icx+icy*2+icz*4;

	  rhob=(deltab[i1+i2*np1]+1.0)*ob*rhoc/POW(astart,3); // physical baryon density in kg/m3
	  pressure=(GAMMA-1.0)*1.5*(rhob/(PROTON_MASS*MOLECULAR_MU))*KBOLTZ*temp; // physical pressure

	  //printf("pres=%e\n",pressure);
	  // filling the cells using supercomoving values

	  //abort();

	  W.d=(deltab[i1+i2*np1]+1.0)*ob/om;
	  W.u=((velx[i1+i2*np1]+veloff)*1e3)*astart/vstar; // vstar is expressed in m/s and grafic vel are in km/s
	  W.v=(vely[i1+i2*np1]*1e3)*astart/vstar;
	  W.w=(velz[i1+i2*np1]*1e3)*astart/vstar;
	  W.p=pressure/pstar*POW(astart,5);
	  W.a=sqrt(GAMMA*W.p/W.d);
	  getE(&W);

	  if(W.p<PMIN) {printf(" YOU SHOULD RECONSIDER PMIN %e %e %e\n",W.p,pressure,pstar);abort();}

#ifdef WRADHYD
	  // Testing ADVECTION
	  //W.X=(i1/6)%2+((i2+1)/6)%2;
	  W.dX=0.2e-3*W.d;
#endif
	  memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	  ifound++;
	}
	else{

	  // this branch corresponds to cell out of the domain
	  //	  printf("euh pas trouve! hidx=%d %p",hidx,cpu->htable[hidx]);
	  //	  abort();
	}
      }
    }
  }

  if(cpu->rank==0){
    fclose(fdx);
    fclose(fx);
    fclose(fy);
    fclose(fz);
  }


  free(deltab);
  free(velx);
  free(vely);
  free(velz);

  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  param->cosmo->H0=h0;
  param->cosmo->unit_l=rstar;

#ifdef WRAD
  param->unit.unit_l=rstar;
  param->unit.unit_v=vstar;
  param->unit.unit_t=param->unit.unit_l/param->unit.unit_v;
  param->unit.unit_n=1.;//(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS)*rhostar; //
  param->unit.unit_mass=rhostar*POW(param->unit.unit_l,3);
  param->unit.unit_d=rhostar; // kg/m3
  param->unit.unit_N=rhostar/PROTON_MASS; // atom/m3
#endif

//  REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS*MOLECULAR_MU/param->unit.unit_mass)/POW(param->unit.unit_l,3);
 /*   REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS*MOLECULAR_MU)*rhostar; */

/* if(cpu->rank==RANK_DISP) printf("navg=%e \n",navg); */

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif
  if(cpu->rank==RANK_DISP) printf("Grafic hydro read ok\n");
  return ifound;
}
#endif
#endif
#endif

// ======================================================================
// ======================================================================
int copy(char const * const source, char const * const destination) {
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
  mkdir(folder_step, 0755);
#endif

	if(pdump){
	  // === particle dump
#ifdef PIC

#ifdef MULTIFOLDER
#ifdef STARS
    sprintf(folder_field,"%sstar/",folder_step);
    mkdir(folder_field, 0755);
    copy_param(folder_field);
#endif // STARS
    sprintf(folder_field,"%spart/",folder_step);
    mkdir(folder_field, 0755);
    copy_param(folder_field);
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
    mkdir(folder_field, 0755);
    copy_param(folder_field);
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

	    sprintf(filename,"bkp/grid.%05d.p%05d",*(cpu->ndumps)%2+1,cpu->rank);
	    save_amr(filename,firstoct,tdump,cpu->tinit,cpu->nsteps,*(cpu->ndumps),param,cpu,cpu->firstpart,adt);

#ifdef PIC
	    // backups for restart
	    sprintf(filename,"bkp/part.%05d.p%05d",*(cpu->ndumps)%2+1,cpu->rank);
	    save_part(filename,firstoct,param->lcoarse,param->lmax,tdump,cpu,cpu->firstpart);
#endif

	  }
#endif
	}
}
