#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prototypes.h"


//====================================================================================================
//            AMR OUT
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
//            AMR IN
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
//            PART OUT
//====================================================================================================

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

//====================================================================================================
//            PART IN
//====================================================================================================

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


