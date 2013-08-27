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

  
  printf("%d octs dumped by proc %d\n",ic,cpu->rank);


  fclose(fp);

}



//====================================================================================================
struct OCT * restore_amr(char filename[], struct OCT **firstoct,struct OCT **lastoct, REAL *tsim, REAL *tinit, int *nsteps, int *ndumps,struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt){
  
  FILE *fp;
  int level,lcoarse,lmax;
  struct OCT * curoct;

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
  

  // reading snapshot time
  fread(tsim,sizeof(REAL),1,fp); 
  fread(tinit,sizeof(REAL),1,fp); 
  fread(nsteps,sizeof(int),1,fp); 
  fread(ndumps,sizeof(int),1,fp); 
  fread(&lcoarse,sizeof(int),1,fp); 
  fread(&lmax,sizeof(int),1,fp); 

  for(level=1;level<=lmax;level++){
    fread(adt+level-1,sizeof(REAL),1,fp); 
  }
  
#ifdef WRAD
  // reading units
  fread(&(param->unit.unit_l),sizeof(REAL),1,fp);
  fread(&(param->unit.unit_v),sizeof(REAL),1,fp);
  fread(&(param->unit.unit_t),sizeof(REAL),1,fp);
  fread(&(param->unit.unit_n),sizeof(REAL),1,fp);
  fread(&(param->unit.unit_mass),sizeof(REAL),1,fp);
#endif

#ifdef TESTCOSMO
  fread(&(param->cosmo->om),sizeof(REAL),1,fp);
  fread(&(param->cosmo->ov),sizeof(REAL),1,fp);
  fread(&(param->cosmo->ob),sizeof(REAL),1,fp);
  fread(&(param->cosmo->H0),sizeof(REAL),1,fp);
#endif
  // reading pointer informations in the snapshot
  fread(&root_sna,sizeof(struct OCT*),1,fp);
  fread(&rootcell_sna,sizeof(struct CELL*),1,fp);
  fread(&rootpart_sna,sizeof(struct PART*),1,fp);
  
  
  if(cpu->rank==0) printf(" %p %p %p\n",root_sna,rootcell_sna,rootpart_sna);
  // reading the octs sequence
   
  fread(&oct_ad,sizeof(struct OCT *),1,fp);
  fread(&oct,sizeof(struct OCT),1,fp);

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
    

    // 2.c modifity the particle pointers within cells
    
    for(ic=0;ic<8;ic++){
      curoct->cell[ic].phead=(curoct->cell[ic].phead==NULL?NULL:(curoct->cell[ic].phead-rootpart_sna)+rootpart_mem);
      if(curoct->cell[ic].phead!=NULL)
	if(curoct->cell[ic].phead->level!=7)
	  abort();
      
    }


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

    fread(&oct_ad,sizeof(struct OCT *),1,fp);
    fread(&oct,sizeof(struct OCT),1,fp);

  }

  printf("%d octs recovered by proc %d with root=%p firstoct=%p\n",ioct,cpu->rank,root_mem,firstoct[0]);

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
	loct=curoct;
      }
    }
  }
  
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
  printf("wrote %d particles (%d expected) in %s on proc %d\n",ipart,npart,filename,cpu->rank);

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
  struct PART *lastp;

  rootpart_mem=proot; // the root cell of the grid in memory

  // opening the file
  fp=fopen(filename,"rb");
  

  // reading snapshot time
  fread(&npart,1,sizeof(int),fp);
  fread(tsim,1,sizeof(REAL),fp);
  fread(&rootpart_sna,1,sizeof(struct PART *),fp);

  // reading the particle sequence
   
  fread(&part_ad,sizeof(struct PART*),1,fp);
  fread(&part,sizeof(struct PART),1,fp);

  lastp=proot;

  while(!feof(fp)){

    // do stuff
    ipart++;
    // 1 copy the content of the particle at the right location
    curp=(part_ad-rootpart_sna)+rootpart_mem;
    memcpy(curp,&part,sizeof(struct PART));

    // 2.a modify the particle pointers 
    curp->next=(curp->next==NULL?NULL:(curp->next-rootpart_sna)+rootpart_mem);
    curp->prev=(curp->prev==NULL?NULL:(curp->prev-rootpart_sna)+rootpart_mem);
    
    // 3. setting the last particle
    if(curp>lastp) lastp=curp;


    // read next particle

    fread(&part_ad,sizeof(struct PART*),1,fp);
    fread(&part,sizeof(struct PART),1,fp);

  }

  printf("%d part recovered by proc %d\n",ipart,cpu->rank);

  // done
  fclose(fp);

  return lastp;
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
      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ngridmax);
      rstat=fscanf(buf,"%s %d",stream,&param->npartmax);
      rstat=fscanf(buf,"%s %d",stream,&param->nbuff);
      
      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ndumps);
      rstat=fscanf(buf,"%s %d",stream,&param->nsteps);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->dt=dummyf;

      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->lcoarse);
      rstat=fscanf(buf,"%s %d",stream,&param->lmax);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->amrthresh=dummyf;

      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->niter);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->poissonacc=dummyf;
      rstat=fscanf(buf,"%s %d",stream,&param->mgridlmin);
      rstat=fscanf(buf,"%s %d",stream,&param->nvcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nrelax);

      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->nrestart);

      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->gstride);
      rstat=fscanf(buf,"%s %d",stream,&param->hstride);
      rstat=fscanf(buf,"%s %d",stream,&param->nsubcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nthread);
      rstat=fscanf(buf,"%s %d",stream,&param->nstream);

#ifdef WRAD
      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->clight=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->denthresh=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->tmpthresh=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->srcint=dummyf;
#endif
      fclose(buf);
    }

  // computing the maxhash
  int val=(pow(2,param->lmax-1)<256?pow(2,param->lmax-1):256); // limit to 2097152 octs in hash table i.e. 16e6 cells
  param->maxhash=pow(val,3);
  //printf("maxhash=%d\n",param->maxhash);
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
    fx=fopen("utils/grafic_src/ic_velcx","rb");
    fy=fopen("utils/grafic_src/ic_velcy","rb");
    fz=fopen("utils/grafic_src/ic_velcz","rb");

    // reading the headers

    fread(&dummy,1,sizeof(dummy),fx);
    fread(&np1,1,4,fx);
    fread(&np2,1,4,fx);
    fread(&np3,1,4,fx);
    fread(&dx,1,4,fx);
    fread(&x1o,1,4,fx);
    fread(&x2o,1,4,fx);
    fread(&x3o,1,4,fx);
    fread(&astart,1,4,fx);
    fread(&om,1,4,fx);
    fread(&ov,1,4,fx);
    fread(&h0,1,4,fx);
    fread(&dummy,1,sizeof(dummy),fx);

    fread(&dummy,1,sizeof(dummy),fy);
    fread(&np1,1,4,fy);
    fread(&np2,1,4,fy);
    fread(&np3,1,4,fy);
    fread(&dx,1,4,fy);
    fread(&x1o,1,4,fy);
    fread(&x2o,1,4,fy);
    fread(&x3o,1,4,fy);
    fread(&astart,1,4,fy);
    fread(&om,1,4,fy);
    fread(&ov,1,4,fy);
    fread(&h0,1,4,fy);
    fread(&dummy,1,sizeof(dummy),fy);

    fread(&dummy,1,sizeof(dummy),fz);
    fread(&np1,1,4,fz);
    fread(&np2,1,4,fz);
    fread(&np3,1,4,fz);
    fread(&dx,1,4,fz);
    fread(&x1o,1,4,fz);
    fread(&x2o,1,4,fz);
    fread(&x3o,1,4,fz);
    fread(&astart,1,4,fz);
    fread(&om,1,4,fz);
    fread(&ov,1,4,fz);
    fread(&h0,1,4,fz);
    fread(&dummy,1,sizeof(dummy),fz);
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

  if((np1*np2*np3)/(cpu->nproc)*1.1>param->npartmax){
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
  for(i3=1;i3<=np3;i3++){

    if(cpu->rank==0){
      printf("\r %f percent done",(i3*1.0)/np3*100.);
      fread(&dummy,1,sizeof(dummy),fx);
      fread(velx,np1*np2,sizeof(float),fx);
      fread(&dummy,1,sizeof(dummy),fx);
      
      fread(&dummy,1,sizeof(dummy),fy);
      fread(vely,np1*np2,sizeof(float),fy);
      fread(&dummy,1,sizeof(dummy),fy);
      
      fread(&dummy,1,sizeof(dummy),fz);
      fread(velz,np1*np2,sizeof(float),fz);
      fread(&dummy,1,sizeof(dummy),fz);
      

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
#endif
//==================================================================================
//==================================================================================

#ifdef WHYDRO2
#ifdef TESTCOSMO
int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, struct RUNPARAMS *param){
  
  FILE *fx;
  FILE *fy;
  FILE *fz;
  FILE *fdx;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,ob,h0;
  int dummy;
  struct PART *lastpart;
  int ip;
  struct Wtype W;

  // Note only the rank 0 reads the file.
 
  if(cpu->rank==0){
    fdx=fopen("utils/grafic_src/ic_deltab","rb");
    fx=fopen("utils/grafic_src/ic_velcx","rb");
    fy=fopen("utils/grafic_src/ic_velcy","rb");
    fz=fopen("utils/grafic_src/ic_velcz","rb");
  
    // reading the headers

    fread(&dummy,1,sizeof(dummy),fdx);
    fread(&np1,1,4,fdx);
    fread(&np2,1,4,fdx);
    fread(&np3,1,4,fdx);
    fread(&dx,1,4,fdx);
    fread(&x1o,1,4,fdx);
    fread(&x2o,1,4,fdx);
    fread(&x3o,1,4,fdx);
    fread(&astart,1,4,fdx);
    fread(&om,1,4,fdx);
    fread(&ov,1,4,fdx);
    fread(&h0,1,4,fdx);
    fread(&dummy,1,sizeof(dummy),fdx);

    fread(&dummy,1,sizeof(dummy),fx);
    fread(&np1,1,4,fx);
    fread(&np2,1,4,fx);
    fread(&np3,1,4,fx);
    fread(&dx,1,4,fx);
    fread(&x1o,1,4,fx);
    fread(&x2o,1,4,fx);
    fread(&x3o,1,4,fx);
    fread(&astart,1,4,fx);
    fread(&om,1,4,fx);
    fread(&ov,1,4,fx);
    fread(&h0,1,4,fx);
    fread(&dummy,1,sizeof(dummy),fx);

    fread(&dummy,1,sizeof(dummy),fy);
    fread(&np1,1,4,fy);
    fread(&np2,1,4,fy);
    fread(&np3,1,4,fy);
    fread(&dx,1,4,fy);
    fread(&x1o,1,4,fy);
    fread(&x2o,1,4,fy);
    fread(&x3o,1,4,fy);
    fread(&astart,1,4,fy);
    fread(&om,1,4,fy);
    fread(&ov,1,4,fy);
    fread(&h0,1,4,fy);
    fread(&dummy,1,sizeof(dummy),fy);

    fread(&dummy,1,sizeof(dummy),fz);
    fread(&np1,1,4,fz);
    fread(&np2,1,4,fz);
    fread(&np3,1,4,fz);
    fread(&dx,1,4,fz);
    fread(&x1o,1,4,fz);
    fread(&x2o,1,4,fz);
    fread(&x3o,1,4,fz);
    fread(&astart,1,4,fz);
    fread(&om,1,4,fz);
    fread(&ov,1,4,fz);
    fread(&h0,1,4,fz);
    fread(&dummy,1,sizeof(dummy),fz);
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
    printf("No temperature law for cosmologies other than SCDM\n");
    abort();
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

      fread(&dummy,1,sizeof(dummy),fdx);
      fread(deltab,np1*np2,sizeof(float),fdx);
      fread(&dummy,1,sizeof(dummy),fdx);

      fread(&dummy,1,sizeof(dummy),fx);
      fread(velx,np1*np2,sizeof(float),fx);
      fread(&dummy,1,sizeof(dummy),fx);

      fread(&dummy,1,sizeof(dummy),fy);
      fread(vely,np1*np2,sizeof(float),fy);
      fread(&dummy,1,sizeof(dummy),fy);

      fread(&dummy,1,sizeof(dummy),fz);
      fread(velz,np1*np2,sizeof(float),fz);
      fread(&dummy,1,sizeof(dummy),fz);
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

