#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "friedmann.h"
#include "segment.h"
#include <string.h>





//====================================================================================================
void save_amr(char filename[], struct OCT **firstoct,REAL tsim, REAL tinit,int nsteps, int ndumps, struct RUNPARAMS *param, struct CPUINFO *cpu,struct PART *proot, REAL *adt){
  
  FILE *fp;
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
struct OCT * restore_amr(char filename[], struct OCT **firstoct,struct OCT **lastoct, REAL *tsim, REAL *tinit, int *nsteps, int *ndumps,struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt){
  
  FILE *fp;
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
#endif

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
  
  
  //if(cpu->rank==0) printf(" %p %p %p\n",root_sna,rootcell_sna,rootpart_sna);
  // reading the octs sequence
   
  outf=fread(&oct_ad,sizeof(struct OCT *),1,fp);
  outf=fread(&oct,sizeof(struct OCT),1,fp);

  while(!feof(fp)){

    // do stuff
    ioct++;
    // 1 copy the content of the oct at the right location
    curoct=root_mem+(oct_ad-root_sna);
    memcpy(curoct,&oct,sizeof(struct OCT));

    //if(cpu->rank==0) printf("le=%d ic=%d\n",curoct->level,ic);
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
	curoct->nei[ic]=(struct CELL *)(((unsigned long int)(curoct->nei[ic])-(unsigned long int)rootcell_sna)+(unsigned long int)rootcell_mem);
      }
      
    }
    
    if(curoct->parent!=NULL){
      curoct->parent=(struct CELL *)(((unsigned long int)(curoct->parent)-(unsigned long int)rootcell_sna)+(unsigned long int)rootcell_mem);
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
  FILE *fp;


   fp=fopen(filename,"wb");

  //printf("tsim=%f\n",tsim);
  fwrite(&tsim,sizeof(REAL),1,fp); 

#ifdef WRAD
  fwrite(&(param->unit.unit_l),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_v),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_t),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_n),sizeof(REAL),1,fp);
  fwrite(&(param->unit.unit_mass),sizeof(REAL),1,fp);
#endif

  fwrite(&(firstoct[0]),sizeof(struct OCT*),1,fp);


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
	  
	  fwrite(&oct,sizeof(struct OCT),1,fp);
	}while(nextoct!=NULL);
    }

  
  fclose(fp);
}


//====================================================================================================
//=================================================================================================

  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
#ifdef PIC
void dumppart(struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu){

  FILE *fp;
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
  for(level=levelcoarse;level<=levelmax;level++) npart+=cpu->npart[level-1];

  fp=fopen(filename,"wb");
  fwrite(&npart,1,sizeof(int),fp);
  fwrite(&tsimf,1,sizeof(float),fp);
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      //printf("level=%d\n",level);
      // setting the first oct
      
      nextoct=firstoct[level-1];
      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
		do{ 
		  curp=nexp; 
		  nexp=curp->next; 
		  val=curp->x;fwrite(&val,1,sizeof(float),fp);
		  val=curp->y;fwrite(&val,1,sizeof(float),fp);
		  val=curp->z;fwrite(&val,1,sizeof(float),fp);
#ifndef PARTN
#ifdef PART_EGY
		  val=curp->ekin+curp->epot;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;fwrite(&val,1,sizeof(float),fp);
#else
		  val=curp->vx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vy;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vz;fwrite(&val,1,sizeof(float),fp);
#endif
#else
		  val=curp->fx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fz;fwrite(&val,1,sizeof(float),fp);
#endif
		  val=(REAL)(curp->idx);fwrite(&val,1,sizeof(float),fp);
		  val=(REAL)(curp->mass);fwrite(&val,1,sizeof(float),fp);
		  val=(REAL)(curp->epot);fwrite(&val,1,sizeof(float),fp);
		  val=(REAL)(curp->ekin);fwrite(&val,1,sizeof(float),fp);
		  ipart++;
		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }
  fclose(fp);
  //printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}

//=======================================================================================================

void save_part(char filename[],struct OCT **firstoct, int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct PART* proot){

  FILE *fp;
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
  fwrite(&npart,1,sizeof(int),fp);
  fwrite(&tsim,1,sizeof(REAL),fp);
  fwrite(&proot,1,sizeof(struct PART *),fp);

  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      //printf("level=%d\n",level);
      // setting the first oct
      
      nextoct=firstoct[level-1];
      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
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
  
  FILE *fp;
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

  rootpart_mem=proot; // the root cell of the grid in memory

  // opening the file
  fp=fopen(filename,"rb");
  

  // reading snapshot time
  size_t outf;
  outf=fread(&npart,1,sizeof(int),fp);
  outf=fread(tsim,1,sizeof(REAL),fp);
  outf=fread(&rootpart_sna,1,sizeof(struct PART *),fp);

  // reading the particle sequence
   
  outf=fread(&part_ad,sizeof(struct PART*),1,fp);
  outf=fread(&part,sizeof(struct PART),1,fp);

  while(!feof(fp)){

    // do stuff
    ipart++;
    // 1 copy the content of the particle at the right location
    curp=(part_ad-rootpart_sna)+rootpart_mem;
    memcpy(curp,&part,sizeof(struct PART));

    // 2.a modify the particle pointers 
    curp->next=(curp->next==NULL?NULL:(curp->next-rootpart_sna)+rootpart_mem);
    curp->prev=(curp->prev==NULL?NULL:(curp->prev-rootpart_sna)+rootpart_mem);
    


    // read next particle

    outf=fread(&part_ad,sizeof(struct PART*),1,fp);
    outf=fread(&part,sizeof(struct PART),1,fp);

  }


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
    /* else if(cpu->rank==0){ */
    /*   printf("%e\n",curp->mass); */
    /* } */
  }

  printf("%d/%d part recovered by proc %d with freepart=%p\n",ipart,param->npartmax,cpu->rank,freepart);

  // done
  fclose(fp);

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
  REAL dummyf;

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
      rstat=fscanf(buf,"%s %d",stream,&param->nsteps);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->dt=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->tmax=dummyf;

       rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->lcoarse);
      rstat=fscanf(buf,"%s %d",stream,&param->lmax);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->amrthresh=dummyf;
      rstat=fscanf(buf,"%s %d",stream,&param->nsmooth);

       rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->niter);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->poissonacc=dummyf;
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

#ifdef WRAD
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->clight=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->denthresh=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->tmpthresh=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->srcint=dummyf;
      param->fudgecool=1.0;
      param->ncvgcool=0;

#endif
      fclose(buf);
    }

  // computing the maxhash
  int val=(pow(2,param->lmax-1)<256?pow(2,param->lmax-1):256); // limit to 2097152 octs in hash table i.e. 16e6 cells
  param->maxhash=pow(val,3);
  //printf("maxhash=%d\n",param->maxhash);

  // ====================== some checks

  // stencil/streams conformity
#ifdef GPUAXL
  if(param->hstride<(param->nthread*param->nstream)){
    printf(" Stream Thread granulosity too high : nt=%d ns=%d stencil=%d\n",param->hstride,param->nthread,param->nstream);
    abort();
  }
#endif
    
}

//==================================================================================
//==================================================================================
#ifdef PIC
#ifdef GRAFIC
struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param)
{
  FILE *fx, *fy, *fz;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
  int dummy;
  struct PART *lastpart;
  int ip;


  if(cpu->rank==0){
    fx=fopen("./ic_velcx","rb");
    fy=fopen("./ic_velcy","rb");
    fz=fopen("./ic_velcz","rb");

    // reading the headers
    size_t outf;

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

#ifdef WMPI
  // Massive broadcast of grafic header
  MPI_Bcast(&np1,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np2,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np3,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&dx,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x1o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x2o,1,MPI_FLOAT,0,cpu->comm);
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

  if((np1*np2*np3)/(cpu->nproc)*1.>param->npartmax){
    printf("Error : Number of particles greater than npartmax Np(grafic, est. )=%d Npartmax=%d!\n",(np1*np2*np3)/(cpu->nproc),param->npartmax);
    abort();
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

  double mass;

#ifdef WHYDRO2
  mass=(1.-ob/om)/(np1*np2*np3);
#else
  mass=1./(np1*np2*np3);
#endif

  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);


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

	x+=(x<=0.)*((int)(-x)+1.)-(x>1.)*((int)x); 
	y+=(y<=0.)*((int)(-y)+1.)-(y>1.)*((int)y); 
	z+=(z<=0.)*((int)(-z)+1.)-(z>1.)*((int)z); 

	// computing the velocities
	vx=velx[(i1-1)+(i2-1)*np1]*astart/(np1*dx*h0)/(sqrt(om)*0.5);
	vy=vely[(i1-1)+(i2-1)*np1]*astart/(np2*dx*h0)/(sqrt(om)*0.5);
	vz=velz[(i1-1)+(i2-1)*np1]*astart/(np3*dx*h0)/(sqrt(om)*0.5);

	// if it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,cpu,cpu->levelcoarse)){
	  
	  part[ip].x=x;
	  part[ip].y=y;
	  part[ip].z=z;
	  
	  part[ip].vx=vx;
	  part[ip].vy=vy;
	  part[ip].vz=vz;
	  
	  part[ip].mass=mass;
	  part[ip].idx=(i1-1)+(i2-1)*np1+(i3-1)*np1*np2;
	  lastpart=part+ip;
	  ip++;
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

  if(cpu->rank==0){
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

  if(cpu->rank==0){
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
	dxcur=pow(0.5,level);
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
		W.u=-(1.+ZC)/pow(1.+ZI,1.5)*sin(2.*M_PI*(xc-0.5))/(M_PI); // for omegam=1. only
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

  int n1d=pow(2,param->lcoarse);
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
	      
	      r=sqrt(pow(x-0.5,2)+pow(y-0.5,2)+pow(z-0.5,2));
	      if(r<lsphere){
		nin++;
		mout=-1;
		vsphere+=pow(dx,3);
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

  if(cpu->rank==0){
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
	dxcur=pow(0.5,level);
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
  FILE *fd;
  
  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;


  struct Wtype WL, WR;
  if(cpu->rank==0) printf("Init Hydro\n");
  
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
	dxcur=pow(0.5,level);
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
  REAL tstar=sqrt(M_PI*M_PI/8.)*pow(R,1.5)/pow(M,0.5);
  if(cpu->rank==0) printf("Generating Evrard Test Case ts=%e, rhostar=%e\n",tstar,rhostar);

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=pow(0.5,level);
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
int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, struct RUNPARAMS *param){
  
  FILE *fx;
  FILE *fy;
  FILE *fz;
  FILE *fdx;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,ob,h0;
  int dummy;
  int ip;
  struct Wtype W;
  size_t outf;

  // Note only the rank 0 reads the file.
 
  if(cpu->rank==0){
    fdx=fopen("./ic_deltab","rb");
    fx=fopen("./ic_velcx","rb");
    fy=fopen("./ic_velcy","rb");
    fz=fopen("./ic_velcz","rb");
  
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



  if(np1!=(int)pow(2,cpu->levelcoarse)){
    printf("ERROR !ABORT! Grafic file not compliant with parameter file : ngrafic=%d nquartz=%d\n",np1,(int)pow(2,cpu->levelcoarse));
    abort();
  }


  // reading the grafic planes

  float *deltab;
  float *velz;
  float *vely;
  float *velx;

  int i1,i2,i3;
  int icx,icy,icz,icell;
  unsigned long key;
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
  if(om==1.) {
    temp=33.64/pow(41.,2)*pow(1.+zstart,2);
    if(cpu->rank==0) printf("WARNING: YOU ARE USING SCDM COSMOLOGY\n");
  }
  else{
    printf("No temperature law for cosmologies other than SCDM -> F** it\n");
    temp=33.64/pow(41.,2)*pow(1.+zstart,2);
    //    abort();
  }


  // supercomoving unit values
  double rhostar;
  double rstar;
  double vstar;
  double tstar;
  double tstar2;
  double pstar;
  double mpc=3.08568025e22; // Mpc in m

  rstar= np1*dx*mpc; // box size in m
  rhostar=rhoc*om;
  tstar=2./H0/sqrt(om); // sec
  tstar2=2./h0/sqrt(om); // Mpc sec / km
  vstar=rstar/tstar; //m/s
  pstar=rhostar*vstar*vstar;
  
  if(cpu->rank==0) printf("rhoc=%e temperature=%lf rstar=%e(%e) pstar=%e tstar=%e vstar=%e rhostar=%e\n",rhoc,temp,rstar,np1*dx,pstar,tstar,vstar,rhostar);

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
	
	key=pos2key(x0,y0,z0,cpu->levelcoarse);

	// first we compute the adress from the hashfunction
	hidx=hfun(key,cpu->maxhash);
	nextoct=cpu->htable[hidx];

	// looking for the oct
	found=0;
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,cpu->levelcoarse)==key)&&(curoct->level==cpu->levelcoarse));
	  }while((nextoct!=NULL)&&(!found));
	}

	// filling the cell

	if(found){
	  icx=i1%2;
	  icy=i2%2;
	  icz=i3%2;
	  
	  icell=icx+icy*2+icz*4;
	
	  rhob=(deltab[i1+i2*np1]+1.0)*ob*rhoc/pow(astart,3); // physical baryon density in kg/m3
	  pressure=(GAMMA-1.0)*1.5*(rhob/(MOLECULAR_MU*PROTON_MASS))*KBOLTZ*temp; // physical pressure
	  
	  //printf("pres=%e\n",pressure);
	  // filling the cells using supercomoving values
	  
	  //abort();

	  W.d=(deltab[i1+i2*np1]+1.0)*ob/om;
	  W.u=(velx[i1+i2*np1]*1e3)*astart/vstar; // vstar is expressed in m/s and grafic vel are in km/s
	  W.v=(vely[i1+i2*np1]*1e3)*astart/vstar;
	  W.w=(velz[i1+i2*np1]*1e3)*astart/vstar;
	  W.p=pressure/pstar*pow(astart,5);
	  W.a=sqrt(GAMMA*W.p/W.d);
	  getE(&W);

#ifdef WRADHYD
	  // Testing ADVECTION
	  //W.X=(i1/6)%2+((i2+1)/6)%2;
	  W.X=0.2e-3;
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

#ifdef WRAD
  param->unit.unit_l=rstar;
  param->unit.unit_v=vstar;
  param->unit.unit_t=param->unit.unit_l/param->unit.unit_v;
  param->unit.unit_n=1.;
  param->unit.unit_mass=rhostar*pow(param->unit.unit_l,3);
#endif

  
#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif
  if(cpu->rank==0) printf("Grafic hydro read ok\n");
  return ifound;
}
#endif
#endif
#endif

// ======================================================================
// ======================================================================

void dumpIO(REAL tsim, struct RUNPARAMS *param,struct CPUINFO *cpu, struct OCT **firstoct, REAL *adt, int pdump){

  REAL tdump,adump;
  char filename[128]; 
  
#ifndef TESTCOSMO
#ifdef WRAD
	tdump=(tsim)*param->unit.unit_t/MYR;
#else
	tdump=(tsim);
#endif
	adump=tdump;
#else
	tdump=interp_aexp(tsim,param->cosmo->tab_aexp,param->cosmo->tab_ttilde);
	adump=tdump;
#endif



	if(pdump){
	  // === particle dump
#ifdef PIC
	  sprintf(filename,"data/part.%05d.p%05d",*(cpu->ndumps),cpu->rank);
	  if(cpu->rank==0){
	    printf("Dumping .......");
	    printf("%s\n",filename);
	  }
	  dumppart(firstoct,filename,param->lcoarse,param->lmax,tdump,cpu);
	  
	  // backups for restart
	  sprintf(filename,"bkp/part.%05d.p%05d",*(cpu->ndumps),cpu->rank); 
	  save_part(filename,firstoct,param->lcoarse,param->lmax,tdump,cpu,cpu->part);
#endif
	}
	else{
	  // === Hydro dump
    
	  sprintf(filename,"data/grid.%05d.p%05d",*(cpu->ndumps),cpu->rank); 
	  if(cpu->rank==0){
	    printf("Dumping .......");
	    printf("%s\n",filename);
	  }
	  dumpgrid(param->lmax,firstoct,filename,adump,param); 

	  // backups for restart
	  sprintf(filename,"bkp/grid.%05d.p%05d",*(cpu->ndumps),cpu->rank); 
	  save_amr(filename,firstoct,tdump,cpu->tinit,cpu->nsteps,*(cpu->ndumps),param,cpu,cpu->part,adt);

	  *(cpu->ndumps)=*(cpu->ndumps)+1;
	}
}
