
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include "oct.h"
#include "cic.h"
#include "segment.h"
#include "particle.h"

#ifdef WMPI
#include <mpi.h>
#endif

void setup_mpi(struct CPUINFO *cpu, struct OCT **firstoct, int levelmax, int levelcoarse, int ngridmax,int loadb){

  int nnei=0;
  /*int *mpinei; // the COMPLETE list of neighbors cpu (can be redundant) */
   int *neicpu; // the reduced list of neighbors CPU (not redundant); 
  int hidx;
  int inidx; // counts the number of inner boundaries octs
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *desoct;
  struct OCT *newoct;
  unsigned long key;
  int inei;
  int i,j;
  int nbnd;
  int icell;
  int *flagcpu;

  /* mpinei =(int*)calloc(ngridmax,sizeof(int)); */
  /* neicpu =(int*)calloc(ngridmax,sizeof(int)); */
  flagcpu=(int*) calloc(cpu->nproc,sizeof(REAL));

  if(cpu->bndoct!=NULL) free(cpu->bndoct);
  cpu->bndoct=(struct OCT**)calloc(cpu->nbuff,sizeof(struct OCT*));

  // we clean the hash table
  memset(cpu->htable,0,cpu->maxhash*sizeof(struct OCT*));

  // looking for neighbors

  for(level=1;level<=levelmax;level++)
    {
      nextoct=firstoct[level-1];
      int level_in_cpu=0;

      if(nextoct!=NULL){
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    //if(level>=levelcoarse){
	      if((level<=levelcoarse)&&(loadb)){
		assigncpu2coarseoct(curoct, cpu, levelcoarse);
	      }

	      if(cpu->rank==curoct->cpu) level_in_cpu=1;
	      key=oct2key(curoct,level); // getting the key of the current oct

	      //filling the hash table for levelcoarse octs
	      hidx=hfun(key,cpu->maxhash); // getting the hash index from the key
	      
	      if(cpu->htable[hidx]==NULL){ //this bucket is empty
		cpu->htable[hidx]=curoct;
		curoct->nexthash=NULL;
	      }
	      else{ // we solve for collisions by looking for the last oct the chained list
		newoct=cpu->htable[hidx];
		do{
		  desoct=newoct;
		  newoct=desoct->nexthash;
		}while(newoct!=NULL);
		desoct->nexthash=curoct;
		curoct->nexthash=NULL;
	      }
	  }while(nextoct!=NULL);

	if((level_in_cpu)||(level>=levelcoarse)){
	  nextoct=firstoct[level-1];
	  do{
	    curoct=nextoct;
	    nextoct=curoct->next;
	    if((curoct->cpu!=cpu->rank)&&(curoct->cpu!=-1)){
	      // a neighbor has been found
	      //mpinei[nnei]=curoct->cpu;
	      flagcpu[curoct->cpu]=1;
	      cpu->bndoct[nnei]=curoct; // contains the oct adresses for emission
	      //if(cpu->rank==0) printf("levl=%d nei=%d\n",curoct->level,curoct->cpu);
	      nnei++;
	    }
	  }while(nextoct!=NULL);
	}
      }
    }


  // =========================================== SETTING UP THE COMMUNICATIONS BETWEEN NEIGHBORS


  // computing the mpi neighbor list
  neicpu=(int*) calloc(cpu->nproc,sizeof(REAL));
  for(i=0;i<cpu->nproc;i++) neicpu[i]=500000+i; // large enough to be at the end
  j=0;
  for(i=0;i<cpu->nproc;i++){ // we scan the list
    if(flagcpu[i]==1){
      neicpu[j]=i;
      j++;
    }
  }

  free(flagcpu);
  myradixsort(neicpu,cpu->nproc); // we sort the neighbors list

  /* myradixsort(mpinei,nnei); // we sort the neighbors list */
  /* neicpu[0]=mpinei[0]; */
  /* j=0; */
  /* for(i=1;i<nnei;i++){ // we scan the list */
  /*   if(mpinei[i]!=neicpu[j]){ */
  /*     j++; */
  /*     neicpu[j]=mpinei[i]; */
  /*   } */
  /* } */

  nbnd=nnei;
  nnei=j;



  cpu->nebnd=nbnd;
  cpu->nnei=nnei;
  if(cpu->mpinei!=NULL) free(cpu->mpinei);
  cpu->mpinei=(int*)calloc(nnei,sizeof(int)); // we reallocate the array to free some memory
  for(i=0;i<cpu->nnei;i++) cpu->mpinei[i]=neicpu[i];
  free(neicpu);



  // AT THIS STAGE: 
  // nbnd contains the number of boundary octs
  // nnei contains the number of neighbor procs
  // mpinei contains the rank of the neighbor procs and has a size nnei


  // some displays
  printf("Found %d neighbors and %d bnd octs  on rank %d :",cpu->nnei,cpu->nebnd,cpu->rank);
  for(i=0;i<cpu->nnei;i++) printf("%d ",cpu->mpinei[i]);
  printf("\n");

#ifdef WMPI
  int nvois_loc=cpu->nebnd;
  int nvois_max;
  MPI_Allreduce(&nvois_loc,&nvois_max,1,MPI_INT,MPI_MAX,cpu->comm);
  if(cpu->rank==0) printf("Max number of neighbors=%d\n",nvois_max);
#endif
  
#ifdef WMPI
  if(cpu->nebnd>cpu->nbuff){
    int err=777;
    printf("ERROR: mpi buffer size should be increased nbnd=%d nbuff=%d\n",cpu->nebnd,cpu->nbuff);
    MPI_Abort(cpu->comm,err);
  }
#endif
  
  // creating a cpu dictionnary to translate from cpu number to inei
  
  if(cpu->dict!=NULL) free(cpu->dict);
  cpu->dict=(int*)calloc(cpu->nproc,sizeof(int));
  for(i=0;i<cpu->nproc;i++) cpu->dict[i]=-1;
  for(i=0;i<cpu->nnei;i++){
    cpu->dict[cpu->mpinei[i]]=i;
  }
}
#ifdef WMPI

//======================================================================================
void gather_ex(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field){
  
  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
   */
  int *countpacket;
  int icpu;
  struct PACKET *pack;
  int i,ii;

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer
    
    // assigning the values
    pack->level=cpu->bndoct[i]->level;
    pack->key=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    for(ii=0;ii<8;ii++){
      switch(field){
      case 0:
	pack->data[ii]=cpu->bndoct[i]->cell[ii].density;
	break;
      case 1:
	pack->data[ii]=cpu->bndoct[i]->cell[ii].marked;
	break;
      }
    }
    
    // counting the number of packets for icpu
    countpacket[icpu]++;
  }
  
  //for(i=0;i<cpu->nnei;i++) printf("rank=%d cpu %d nbnd=%d\n",cpu->rank,cpu->mpinei[i],countpacket[i]);
  free(countpacket);

}

//======================================================================================
void gather_ex_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer){
  
  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
    USed for coarse cell update of hydro at boundaries
   */
  int *countpacket;
  int icpu;
  struct HYDRO_MPI *pack;
  int i,ii;

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer
    
    // assigning the values
    pack->level=cpu->bndoct[i]->level;
    pack->key=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    for(ii=0;ii<8;ii++) pack->data[ii]=cpu->bndoct[i]->cell[ii].fieldnew;
    // counting the number of packets for icpu
    countpacket[icpu]++;
  }
  
  //for(i=0;i<cpu->nnei;i++) printf("rank=%d cpu %d nbnd=%d\n",cpu->rank,cpu->mpinei[i],countpacket[i]);
  free(countpacket);

}

// ====================================================================================
// ====================================================================================

void gather_ex_rad(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer){
  
  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
    USed for coarse cell update of hydro at boundaries
   */
  int *countpacket;
  int icpu;
  struct RAD_MPI *pack;
  int i,ii;

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer
    
    // assigning the values
    pack->level=cpu->bndoct[i]->level;
    pack->key=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    for(ii=0;ii<8;ii++) pack->data[ii]=cpu->bndoct[i]->cell[ii].rfieldnew;
    // counting the number of packets for icpu
    countpacket[icpu]++;
  }
  
  //for(i=0;i<cpu->nnei;i++) printf("rank=%d cpu %d nbnd=%d\n",cpu->rank,cpu->mpinei[i],countpacket[i]);
  free(countpacket);
}


//======================================================================================
int gather_ex_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer,struct PART **lastp){
  
  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
   */
  int *countpacket;
  int icpu;
  struct PART_MPI *part;
  struct PART *curp;
  struct PART *nexp;
  int i,ii;
  int nrem=0;

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    for(ii=0;ii<8;ii++){
      nexp=cpu->bndoct[i]->cell[ii].phead;
      // sweeping the particles of the boundary cells
      if(nexp!=NULL){
	do{ 
	  curp=nexp; 
	  nexp=curp->next; 
	  
	  part=psendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer

	  // assigning the values
	  part->level=cpu->bndoct[i]->level;
	  part->key=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct (is the eventual destination)
	  part->icell=ii;

	  // we assign the data to the communicator
	  part->x=curp->x;
	  part->y=curp->y;
	  part->z=curp->z;

	  part->vx=curp->vx;
	  part->vy=curp->vy;
	  part->vz=curp->vz;

	  part->mass=curp->mass;
	  part->idx=curp->idx;
	  
	  // counting the number of packets for icpu
	  countpacket[icpu]++;

	  // switching the mass to -1 to flag exiting particles
	  curp->mass=-1.;
	  nrem++;

	  // disconnecting the particle from its list

	  curp->prev=NULL;
	  curp->next=NULL;
	  
	  if(nexp==NULL){ // reached the last particle from the cell
	    cpu->bndoct[i]->cell[ii].phead=NULL; // we "clear" the cell
	  }

	  // is it the global last particle ?
	  if(curp==(*lastp)){
	    (*lastp)=curp-1;
	    while((*lastp)->mass<0){ // if true this particle does not exist anymore
	      (*lastp)=(*lastp)-1;
	    }
	  }

	}while(nexp!=NULL);
      }
    }
  }
  free(countpacket);

  return nrem; // we return the number of exiting particles
}



//======================================================================================

void gather_mpi(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field){

  int i,j;
  int found=0;
  int hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  int *countpacket;

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

    for(j=0;j<cpu->nnei;j++){
      for(i=0;i<cpu->nbuff;i++){
	pack=sendbuffer[j]+i; // we assume that the sendbuffer already contains the keys
	if(pack->level!=0){ // we do something
	  //t[0]=MPI_Wtime();
	  countpacket[j]++;

	  // first we compute the adress from the hashfunction
	  hidx=hfun(pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      curoct->border=1;

	      for(icell=0;icell<8;icell++){
		switch(field){
		case 0:
		  pack->data[icell]=curoct->cell[icell].density; // density
		  break;
		case 1:
		  pack->data[icell]=curoct->cell[icell].density; // density again we reproduce the case 1 in order to be consistent with scatter_mpi
		  break;
		case 2:
		  pack->data[icell]=curoct->cell[icell].gdata.p; // potential
		  break;
		case 3:
		  pack->data[icell]=curoct->cell[icell].marked*(curoct->level>=cpu->levelcoarse); //refinment mark
		  break;
		case 4:
		  pack->data[icell]=curoct->cell[icell].temp; //temp field for force calculation
		  break;
#ifdef WGRAV
		case 5:
		  pack->data[icell]=curoct->cell[icell].f[0]; //temp field for force calculation
		  break;
		case 6:
		  pack->data[icell]=curoct->cell[icell].f[1]; //temp field for force calculation
		  break;
		case 7:
		  pack->data[icell]=curoct->cell[icell].f[2]; //temp field for force calculation
		  break;
#endif
		}
	      }
	    }
	    else{
	      printf("error no reception oct found !");
	      abort();
	    }
	  }else{
	    printf("error no hash key obtained !!\n");
	    abort();
	  }
	}
      }
    }
    

    free(countpacket);
}


//====================================================================

void gather_mpi_level(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field, int level){

  int i,j;
  int found=0;
  int hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  int *countpacket;

/*   double t[10]; */
/*   double th=0,tg=0,tt=0,tc=0; */
/*   int nt=0; */

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

    for(j=0;j<cpu->nnei;j++){
      for(i=0;i<cpu->nbuff;i++){
	pack=sendbuffer[j]+i; // we assume that the sendbuffer already contains the keys
	if(pack->level==level){ // we do something
	  //t[0]=MPI_Wtime();
	  countpacket[j]++;

	  // first we compute the adress from the hashfunction
	  hidx=hfun(pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      curoct->border=1;

	      for(icell=0;icell<8;icell++){
		switch(field){
		case 0:
		  pack->data[icell]=curoct->cell[icell].density; // density
		  break;
		case 1:
		  pack->data[icell]=curoct->cell[icell].density; // density again we reproduce the case 1 in order to be consistent with scatter_mpi
		  break;
		case 2:
		  pack->data[icell]=curoct->cell[icell].gdata.p; // potential
		  break;
		case 3:
		  pack->data[icell]=curoct->cell[icell].marked; //refinment mark
		  break;
		case 4:
		  pack->data[icell]=curoct->cell[icell].temp; //temp field for force calculation
		  break;
#ifdef WGRAV
		case 5:
		  pack->data[icell]=curoct->cell[icell].f[0]; //temp field for force calculation
		  break;
		case 6:
		  pack->data[icell]=curoct->cell[icell].f[1]; //temp field for force calculation
		  break;
		case 7:
		  pack->data[icell]=curoct->cell[icell].f[2]; //temp field for force calculation
		  break;
#endif
		}
	      }
	    }
	    else{
	      printf("error no reception oct found !");
	      abort();
	    }
	    //t[2]=MPI_Wtime();

	    
	  }else{
	    printf("error no hash key obtained !!\n");
	    abort();
	  }
	}
      }
    }
    

    free(countpacket);
}

 //------------------------------------------------------------------------
void scatter_mpi(struct CPUINFO *cpu, struct PACKET **recvbuffer,  int field){

  int i,j;
  int found=0;
  int hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun(pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      switch(field){
	      case 0:
		curoct->cell[icell].density+=pack->data[icell]; // density += for CIC correction
		break;
	      case 1:
		curoct->cell[icell].density =pack->data[icell]; // density
		break;
	      case 2:
		curoct->cell[icell].gdata.p=pack->data[icell]; // potential
		break;
	      case 3:
		curoct->cell[icell].marked=fmax(pack->data[icell],(REAL)curoct->cell[icell].marked); // refinement mark
		break;
	      case 4:
		curoct->cell[icell].temp=pack->data[icell]; // temp field for force calculation
		break;
#ifdef WGRAV
	      case 5:
		curoct->cell[icell].f[0]=pack->data[icell]; // temp field for force calculation
		break;
	      case 6:
		curoct->cell[icell].f[1]=pack->data[icell]; // temp field for force calculation
		break;
	      case 7:
		curoct->cell[icell].f[2]=pack->data[icell]; // temp field for force calculation
		break;
#endif
	      }
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%ld\n",i,pack->level,pack->key);
	    abort();
	  }
	    
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
      
}


 //------------------------------------------------------------------------
void scatter_mpi_level(struct CPUINFO *cpu, struct PACKET **recvbuffer,  int field, int level){

  int i,j;
  int found=0;
  int hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level==level){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun(pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      switch(field){
	      case 0:
		curoct->cell[icell].density+=pack->data[icell]; // density += for CIC correction
		break;
	      case 1:
		curoct->cell[icell].density =pack->data[icell]; // density
		break;
	      case 2:
		curoct->cell[icell].gdata.p=pack->data[icell]; // potential
		break;
	      case 3:
		curoct->cell[icell].marked=fmax(pack->data[icell],(REAL)curoct->cell[icell].marked); // refinement mark
		break;
	      case 4:
		curoct->cell[icell].temp=pack->data[icell]; // temp field for force calculation
		break;
#ifdef WGRAV
	      case 5:
		curoct->cell[icell].f[0]=pack->data[icell]; // temp field for force calculation
		break;
	      case 6:
		curoct->cell[icell].f[1]=pack->data[icell]; // temp field for force calculation
		break;
	      case 7:
		curoct->cell[icell].f[2]=pack->data[icell]; // temp field for force calculation
		break;
#endif
	      }
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%ld\n",i,pack->level,pack->key);
	    abort();
	  }
	    
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
      
}


 //------------------------------------------------------------------------
int scatter_mpi_part(struct CPUINFO *cpu, struct PART_MPI **precvbuffer, struct PART **lastp){

  int i,j;
  int found=0;
  int hidx;
  struct PART_MPI *part;
  struct PART *curp;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  int nadd=0;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      part=precvbuffer[j]+i;
      //      if(part==NULL) continue;
      if(part->level!=0){ // we do something
	// first we compute the adress from the hashfunction
	hidx=hfun(part->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==part->key)&&(part->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    // we assign the new particle to the global particle list of the client
	    *lastp=*lastp+1; // getting the next slot in the global particle list
	    

	    // copying the data
	    (*lastp)->x =part->x;
	    (*lastp)->y =part->y;
	    (*lastp)->z =part->z;

	    (*lastp)->vx=part->vx;
	    (*lastp)->vy=part->vy;
	    (*lastp)->vz=part->vz;
	      
	    (*lastp)->mass=part->mass;
	    (*lastp)->idx=part->idx;
	    (*lastp)->level=part->level;

	    nadd++;

	    // if recepction cell is refined
	    struct CELL *newcell=&(curoct->cell[part->icell]);

	    if(newcell->child!=NULL){
	      struct OCT *newoct=newcell->child;
	      REAL dxcur2=1./pow(2.,newoct->level);
				
	      int xp=(int)(2*(part->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
	      int yp=(int)(2*(part->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
	      int zp=(int)(2*(part->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
	      int ip=xp+yp*2+zp*4;
	      

	      newcell=&(newoct->cell[ip]);
	      part->icell=ip;
	      (*lastp)->level=part->level+1;
	    }

	    // looking for the tail particle in the destination cell
	    curp=findlastpart(newcell->phead);
	    if(curp!=NULL){
	      curp->next=(*lastp);
	      (*lastp)->next=NULL;
	      (*lastp)->prev=curp;
	    }
	    else{
	      newcell->phead=(*lastp);
	      (*lastp)->next=NULL;
	      (*lastp)->prev=NULL;
	    }
	  }
	  else{
	    printf("error no reception oct found !");
	    abort();
	  }
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
  return nadd; // returning the new number of particles
}

//------------------------------------------------------------------------
void compute_bndkeys(struct CPUINFO *cpu, struct PACKET **recvbuffer){

  // we create a counter of values for each neighbor
  int *countpacket;
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct PACKET *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    int keyloc;
    int cpuloc;
    int inei;
    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }  

  free(countpacket);
}


void compute_bndkeys_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **recvbuffer){

  // we create a counter of values for each neighbor
  int *countpacket;
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct HYDRO_MPI *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    int keyloc;
    int cpuloc;
    int inei;
    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }  

  free(countpacket);
}

#ifdef WRAD
void compute_bndkeys_rad(struct CPUINFO *cpu, struct RAD_MPI **recvbuffer){

  // we create a counter of values for each neighbor
  int *countpacket;
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct RAD_MPI *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    int keyloc;
    int cpuloc;
    int inei;
    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }  

  free(countpacket);
}
#endif


void compute_bndkeys_flux(struct CPUINFO *cpu, struct FLUX_MPI **recvbuffer){

  // we create a counter of values for each neighbor
  int *countpacket;
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct FLUX_MPI *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    int keyloc;
    int cpuloc;
    int inei;
    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }  

  free(countpacket);
}

//------------------------------------------------------------------------
void compute_bndkeys_level(struct CPUINFO *cpu, struct PACKET **recvbuffer, int level, int *countpacket){

  /* // we create a counter of values for each neighbor */

  /* int *countpacket; */
  /* countpacket=(int*)calloc(cpu->nnei,sizeof(int)); */

  memset(countpacket,0,cpu->nnei*sizeof(int));


  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct PACKET *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    int keyloc;
    int cpuloc;
    int inei;

    if(cpu->bndoct[i]->level != level) continue; // we skip the borders which do not belong to current level

    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }  

  //free(countpacket);
}

//------------------------------------------------------------------------

void  clean_mpibuff(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer){
  int i;
  for(i=0;i<cpu->nnei;i++) {
    memset(sendbuffer[i],0,cpu->nbuff*sizeof(struct PACKET));
    memset(recvbuffer[i],0,cpu->nbuff*sizeof(struct PACKET));
  }
}

void  clean_mpibuff_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer, struct HYDRO_MPI **recvbuffer){
  int i;
  
  for(i=0;i<cpu->nnei;i++) {
    memset(sendbuffer[i],0,cpu->nbuff*sizeof(struct HYDRO_MPI));
    memset(recvbuffer[i],0,cpu->nbuff*sizeof(struct HYDRO_MPI));
  }

}

#ifdef WRAD
void  clean_mpibuff_rad(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer, struct RAD_MPI **recvbuffer){
  int i;
  
  for(i=0;i<cpu->nnei;i++) {
    memset(sendbuffer[i],0,cpu->nbuff*sizeof(struct RAD_MPI));
    memset(recvbuffer[i],0,cpu->nbuff*sizeof(struct RAD_MPI));
  }

}
#endif

void  clean_mpibuff_flux(struct CPUINFO *cpu, struct FLUX_MPI **sendbuffer, struct FLUX_MPI **recvbuffer){
  int i;
  
  for(i=0;i<cpu->nnei;i++) {
    memset(sendbuffer[i],0,cpu->nbuff*sizeof(struct FLUX_MPI));
    memset(recvbuffer[i],0,cpu->nbuff*sizeof(struct FLUX_MPI));
  }

}


void  clean_mpibuff_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer){
  int i;
  for(i=0;i<cpu->nnei;i++) {
    memset(psendbuffer[i],0,cpu->nbuff*sizeof(struct PART_MPI));
    memset(precvbuffer[i],0,cpu->nbuff*sizeof(struct PART_MPI));
  }
}
 //------------------------------------------------------------------------

#ifdef WMPI
void mpi_exchange(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field, int cmp_keys)
{
  int i;
  int icpu;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_PACKET);
  int mpitag=1;
  
  double t[10];
  double tot;
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));
  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));


  // ---------- The key calculation may already been computed (cmp_key=0) or must be recomputed (cmp_key=1)

  if(cmp_keys){
    // ----------- 0  / we clean the mpi buffers
    clean_mpibuff(cpu,sendbuffer,recvbuffer);

    // ----------- I  / we compute the boundary keys and store them in recvbuffer
    compute_bndkeys(cpu,recvbuffer);

    // ----------- II / we send the keys to the server

    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Isend(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
    }
    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Irecv(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
    }
    MPI_Waitall(2*cpu->nnei,req,stat);
    MPI_Barrier(cpu->comm);
  }

  // ----------- III/ the server gather the data
  MPI_Barrier(cpu->comm);
  gather_mpi(cpu, sendbuffer, field);
  MPI_Barrier(cpu->comm);

  //if(cpu->rank==0) printf("--- X ---\n");
  memset(req,0,2*cpu->nnei*sizeof(MPI_Request));
  memset(stat,0,2*cpu->nnei*sizeof(MPI_Status));
  MPI_Barrier(cpu->comm);

  // ----------- IV / the server send the data back to the client


  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank         ,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);



  t[5]=MPI_Wtime();
  // ----------- V  / the client scatter the data back in the oct tree
  scatter_mpi(cpu,recvbuffer,field);
  
  t[6]=MPI_Wtime();

  //
  free(req);
  free(stat);
  t[7]=MPI_Wtime();
  tot=t[7]-t[0];
  //if(cpu->rank==0) printf("clean=%e keys=%e sendkeys=%e gather=%e senddata=%e scatter=%e free=%e\n",(t[1]-t[0])/tot,(t[2]-t[1])/tot,(t[3]-t[2])/tot,(t[4]-t[3])/tot,(t[5]-t[4])/tot,(t[6]-t[5])/tot,(t[7]-t[6])/tot);

}

//===========================================================================================
//===========================================================================================
//===========================================================================================

#ifdef WHYDRO2
void gather_mpi_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer){

  int i,j;
  int found=0;
  int hidx;
  struct HYDRO_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  int *countpacket;

/*   double t[10]; */
/*   double th=0,tg=0,tt=0,tc=0; */
/*   int nt=0; */

  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

    for(j=0;j<cpu->nnei;j++){
      for(i=0;i<cpu->nbuff;i++){
	pack=sendbuffer[j]+i; // we assume that the sendbuffer already contains the keys
	if(pack->level!=0){ // we do something
	  //t[0]=MPI_Wtime();
	  countpacket[j]++;

	  // first we compute the adress from the hashfunction
	  hidx=hfun(pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      curoct->border=1;

	      for(icell=0;icell<8;icell++){
		memcpy(&(pack->data[icell]),&(curoct->cell[icell].field),sizeof(struct Wtype));
	      }
	    }
	    else{
	      printf("error no reception oct found !");
	      abort();
	    }
	    //t[2]=MPI_Wtime();

	    
	  }else{
	    printf("error no hash key obtained !!\n");
	    abort();
	  }
	}
      }
    }
    

    free(countpacket);
}


// ================================================================================================================================
/* void gather_mpi_flux(struct CPUINFO *cpu, struct FLUX_MPI **sendbuffer){ */

/*   int i,j; */
/*   int found=0; */
/*   int hidx; */
/*   struct FLUX_MPI *pack; */
/*   struct OCT *curoct; */
/*   struct OCT *nextoct; */
/*   int icell; */
/*   int *countpacket; */

/* /\*   double t[10]; *\/ */
/* /\*   double th=0,tg=0,tt=0,tc=0; *\/ */
/* /\*   int nt=0; *\/ */

/*   // we create a counter of values for each neighbor */
/*   countpacket=(int*)calloc(cpu->nnei,sizeof(int)); */

/*     for(j=0;j<cpu->nnei;j++){ */
/*       for(i=0;i<cpu->nbuff;i++){ */
/* 	pack=sendbuffer[j]+i; // we assume that the sendbuffer already contains the keys */
/* 	if(pack->level!=0){ // we do something */
/* 	  //t[0]=MPI_Wtime(); */
/* 	  countpacket[j]++; */

/* 	  // first we compute the adress from the hashfunction */
/* 	  hidx=hfun(pack->key,cpu->maxhash); */
/* 	  nextoct=cpu->htable[hidx]; */
/* 	  //t[3]=MPI_Wtime(); */
/* 	  if(nextoct!=NULL){ */
/* 	    do{ // resolving collisions */
/* 	      curoct=nextoct; */
/* 	      nextoct=curoct->nexthash; */
/* 	      found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level)); */
/* 	    }while((nextoct!=NULL)&&(!found)); */

/* 	    //t[1]=MPI_Wtime(); */
/* 	    if(found){ // the reception oct has been found */

/* 	      // we set the current oct as a border one (for vector based communications) */
/* 	      curoct->border=1; */

/* 	      for(icell=0;icell<8;icell++){ */
/* 		memcpy(&(pack->data[icell*NFLUX]),curoct->cell[icell].flux,sizeof(REAL)*NFLUX); */
/* 	      } */
/* 	    } */
/* 	    else{ */
/* 	      printf("error no reception oct found !"); */
/* 	      abort(); */
/* 	    } */
/* 	    //t[2]=MPI_Wtime(); */

	    
/* 	  }else{ */
/* 	    printf("error no hash key obtained !!\n"); */
/* 	    abort(); */
/* 	  } */
/* 	} */
/*       } */
/*     } */
    

/*     free(countpacket); */
/* } */


 //------------------------------------------------------------------------
void scatter_mpi_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **recvbuffer){

  int i,j;
  int found=0;
  int hidx;
  struct HYDRO_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun(pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      memcpy(&(curoct->cell[icell].field),&(pack->data[icell]),sizeof(struct Wtype));
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%ld\n",i,pack->level,pack->key);
	    abort();
	  }
	    
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
      
}

// ============================================================
// ============================================================
// ============================================================
void scatter_mpi_hydro_ext(struct CPUINFO *cpu, struct HYDRO_MPI **recvbuffer,int level){

  int i,j;
  int found=0;
  int hidx;
  struct HYDRO_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  struct Utype U;
  struct Utype Ue;
  struct Wtype W;
  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun(pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    if(curoct->level!=(level-1)) continue; // we update only coarse neighbours relative to the current level
	    for(icell=0;icell<8;icell++){
	      
	      W2U(&(curoct->cell[icell].fieldnew),&U);
	      W2U(&(pack->data[icell]),&Ue);

	      U.d+=Ue.d;
	      U.du+=Ue.du;
	      U.dv+=Ue.dv;
	      U.dw+=Ue.dw;
	      U.E+=Ue.E;
#ifdef DUAL_E
	      U.eint+=Ue.eint;
#endif	      

#ifdef WRADHYD
	      U.dX+=Ue.dX;
#endif
	      U2W(&U,&W);
	      getE(&W);

	      memcpy(&(curoct->cell[icell].fieldnew),&W,sizeof(struct Wtype));
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%ld\n",i,pack->level,pack->key);
	    abort();
	  }
	    
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
      
}


//=====================================================================

void scatter_mpi_rad_ext(struct CPUINFO *cpu, struct HYDRO_MPI **recvbuffer,int level){

  int i,j;
  int found=0;
  int hidx;
  struct HYDRO_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  struct Rtype *R;
  struct Rtype *Re;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun(pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    if(curoct->level!=(level-1)) continue; // we update only coarse neighbours relative to the current level
	    for(icell=0;icell<8;icell++){
	      
	      R=&(curoct->cell[icell].rfieldnew);
	      Re=&(pack->data[icell]);

	      for(igrp=0;igrp<NGRP;igrp++){
		// update
		R.e[igrp] += Re.e[igrp];
		R.fx[igrp]+= Re.fx[igrp];
		R.fy[igrp]+= Re.fy[igrp];
		R.fz[igrp]+= Re.fz[igrp];
	      }
 	      
	      memcpy(&(curoct->cell[icell].rfieldnew),&R,sizeof(struct Rtype));
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%ld\n",i,pack->level,pack->key);
	    abort();
	  }
	    
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
      
}


//=============================================================================================

/* void scatter_mpi_flux(struct CPUINFO *cpu, struct FLUX_MPI **recvbuffer){ */

/*   int i,j; */
/*   int found=0; */
/*   int hidx; */
/*   struct FLUX_MPI *pack; */
/*   struct OCT *curoct; */
/*   struct OCT *nextoct; */
/*   int icell; */

/*   for(j=0;j<cpu->nnei;j++){ */
/*     for(i=0;i<cpu->nbuff;i++){ */
/*       pack=recvbuffer[j]+i; */
/*       if(pack->level!=0){ // we do something */

/* 	// first we compute the adress from the hashfunction */
/* 	hidx=hfun(pack->key,cpu->maxhash); */
/* 	nextoct=cpu->htable[hidx]; */
/* 	if(nextoct!=NULL){ */
/* 	  do{ // resolving collisions */
/* 	    curoct=nextoct; */
/* 	    nextoct=curoct->nexthash; */
/* 	    found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level)); */
/* 	  }while((nextoct!=NULL)&&(!found)); */

/* 	  if(found){ // the reception oct has been found */
/* 	    for(icell=0;icell<8;icell++){ */
/* 	      memcpy(curoct->cell[icell].flux,&(pack->data[icell*NFLUX]),NFLUX*sizeof(REAL)); */
/* 	    } */
/* 	  } */
/* 	  else{ */
/* 	    printf("error no reception oct found ! for buff #%d lev=%d key=%ld\n",i,pack->level,pack->key); */
/* 	    abort(); */
/* 	  } */
	    
/* 	}else{ */
/* 	  printf("error no hash key obtained !!\n"); */
/* 	  abort(); */
/* 	} */
/*       } */
/*     } */
/*   } */
    
      
/* } */


//=======================================================================================================
//=======================================================================================================

void mpi_exchange_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer, struct HYDRO_MPI **recvbuffer, int cmp_keys)
{
  int i;
  int icpu;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_HYDRO);
  int mpitag=1;
  
  double t[10];
  double tot;
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));
  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));


  // ---------- The key calculation may already been computed (cmp_key=0) or must be recomputed (cmp_key=1)

  if(cmp_keys){
    // ----------- 0  / we clean the mpi buffers
    clean_mpibuff_hydro(cpu,sendbuffer,recvbuffer);

    // ----------- I  / we compute the boundary keys and store them in recvbuffer
    compute_bndkeys_hydro(cpu,recvbuffer);

    // ----------- II / we send the keys to the server

    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Isend(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
    }
    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Irecv(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
    }
    MPI_Waitall(2*cpu->nnei,req,stat);
    MPI_Barrier(cpu->comm);
  }

  // ----------- III/ the server gather the data
  gather_mpi_hydro(cpu, sendbuffer);

  //if(cpu->rank==0) printf("--- X ---\n");
  memset(req,0,2*cpu->nnei*sizeof(MPI_Request));
  memset(stat,0,2*cpu->nnei*sizeof(MPI_Status));
  MPI_Barrier(cpu->comm);

  // ----------- IV / the server send the data back to the client


  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank         ,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);



  t[5]=MPI_Wtime();
  // ----------- V  / the client scatter the data back in the oct tree
  scatter_mpi_hydro(cpu,recvbuffer);
  MPI_Barrier(cpu->comm);
  t[6]=MPI_Wtime();

  //
  free(req);
  free(stat);
  t[7]=MPI_Wtime();
  tot=t[7]-t[0];
  //if(cpu->rank==0) printf("clean=%e keys=%e sendkeys=%e gather=%e senddata=%e scatter=%e free=%e\n",(t[1]-t[0])/tot,(t[2]-t[1])/tot,(t[3]-t[2])/tot,(t[4]-t[3])/tot,(t[5]-t[4])/tot,(t[6]-t[5])/tot,(t[7]-t[6])/tot);

}

// =======================================================================
// =======================================================================


//=======================================================================================================

/* void mpi_exchange_flux(struct CPUINFO *cpu, struct FLUX_MPI **sendbuffer, struct FLUX_MPI **recvbuffer, int cmp_keys) */
/* { */
/*   int i; */
/*   int icpu; */
/*   MPI_Status *stat; */
/*   MPI_Request *req; */
/*   MPI_Datatype MPI_PACKET=*(cpu->MPI_FLUX); */
/*   int mpitag=1; */
  
/*   double t[10]; */
/*   double tot; */
/*   req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request)); */
/*   stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status)); */


/*   // ---------- The key calculation may already been computed (cmp_key=0) or must be recomputed (cmp_key=1) */

/*   if(cmp_keys){ */
/*     // ----------- 0  / we clean the mpi buffers */
/*     clean_mpibuff_flux(cpu,sendbuffer,recvbuffer); */

/*     // ----------- I  / we compute the boundary keys and store them in recvbuffer */
/*     compute_bndkeys_flux(cpu,recvbuffer); */

/*     // ----------- II / we send the keys to the server */

/*     for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys */
/*       MPI_Isend(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  ); */
/*     } */
/*     for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys */
/*       MPI_Irecv(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]); */
/*     } */
/*     MPI_Waitall(2*cpu->nnei,req,stat); */
/*     MPI_Barrier(cpu->comm); */
/*   } */

/*   // ----------- III/ the server gather the data */
/*   gather_mpi_flux(cpu, sendbuffer); */

/*   //if(cpu->rank==0) printf("--- X ---\n"); */
/*   memset(req,0,2*cpu->nnei*sizeof(MPI_Request)); */
/*   memset(stat,0,2*cpu->nnei*sizeof(MPI_Status)); */
/*   MPI_Barrier(cpu->comm); */

/*   // ----------- IV / the server send the data back to the client */


/*   for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys */
/*     MPI_Isend(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank         ,cpu->comm,&req[i]  ); */
/*   } */

/*   for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys */
/*     MPI_Irecv(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]); */
/*   } */
/*   MPI_Waitall(2*cpu->nnei,req,stat); */
/*   MPI_Barrier(cpu->comm); */



/*   t[5]=MPI_Wtime(); */
/*   // ----------- V  / the client scatter the data back in the oct tree */
/*   scatter_mpi_flux(cpu,recvbuffer); */
  
/*   t[6]=MPI_Wtime(); */

/*   // */
/*   free(req); */
/*   free(stat); */
/*   t[7]=MPI_Wtime(); */
/*   tot=t[7]-t[0]; */
/*   //if(cpu->rank==0) printf("clean=%e keys=%e sendkeys=%e gather=%e senddata=%e scatter=%e free=%e\n",(t[1]-t[0])/tot,(t[2]-t[1])/tot,(t[3]-t[2])/tot,(t[4]-t[3])/tot,(t[5]-t[4])/tot,(t[6]-t[5])/tot,(t[7]-t[6])/tot); */

/* } */

#endif


#ifdef WRAD

// ============================================================================================
// ============================================================================================

void gather_mpi_rad(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer){

  int i,j;
  int found=0;
  int hidx;
  struct RAD_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  int *countpacket;


  // we create a counter of values for each neighbor
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

    for(j=0;j<cpu->nnei;j++){
      for(i=0;i<cpu->nbuff;i++){
	pack=sendbuffer[j]+i; // we assume that the sendbuffer already contains the keys
	if(pack->level!=0){ // we do something
	  //t[0]=MPI_Wtime();
	  countpacket[j]++;

	  // first we compute the adress from the hashfunction
	  hidx=hfun(pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      curoct->border=1;

	      for(icell=0;icell<8;icell++){
		memcpy(&(pack->data[icell]),&(curoct->cell[icell].rfield),sizeof(struct Rtype));
	      }
	    }
	    else{
	      printf("error no reception oct found !");
	      abort();
	    }
	    //t[2]=MPI_Wtime();

	    
	  }else{
	    printf("error no hash key obtained !!\n");
	    abort();
	  }
	}
      }
    }
    

    free(countpacket);
}

// =======================================================================
void scatter_mpi_rad(struct CPUINFO *cpu, struct RAD_MPI **recvbuffer){

  int i,j;
  int found=0;
  int hidx;
  struct RAD_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun(pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      memcpy(&(curoct->cell[icell].rfield),&(pack->data[icell]),sizeof(struct Rtype));
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%ld\n",i,pack->level,pack->key);
	    abort();
	  }
	    
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }
    
      
}
// =======================================================================

void mpi_exchange_rad(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer, struct RAD_MPI **recvbuffer, int cmp_keys)
{
  int i;
  int icpu;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_RAD);
  int mpitag=1;
  
  double t[10];
  double tot;
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));
  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));


  // ---------- The key calculation may already been computed (cmp_key=0) or must be recomputed (cmp_key=1)
  if(cmp_keys){
    // ----------- 0  / we clean the mpi buffers
    clean_mpibuff_rad(cpu,sendbuffer,recvbuffer);

    // ----------- I  / we compute the boundary keys and store them in recvbuffer
    compute_bndkeys_rad(cpu,recvbuffer);

    // ----------- II / we send the keys to the server

    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Isend(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
    }
    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Irecv(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
    }
    MPI_Waitall(2*cpu->nnei,req,stat);
    MPI_Barrier(cpu->comm);
  }

  // ----------- III/ the server gather the data
  gather_mpi_rad(cpu, sendbuffer);

  //if(cpu->rank==0) printf("--- X ---\n");
  memset(req,0,2*cpu->nnei*sizeof(MPI_Request));
  memset(stat,0,2*cpu->nnei*sizeof(MPI_Status));
  MPI_Barrier(cpu->comm);

  // ----------- IV / the server send the data back to the client


  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank         ,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);



  t[5]=MPI_Wtime();
  // ----------- V  / the client scatter the data back in the oct tree
  scatter_mpi_rad(cpu,recvbuffer);
  MPI_Barrier(cpu->comm);
  t[6]=MPI_Wtime();

  //
  free(req);
  free(stat);
  t[7]=MPI_Wtime();
  tot=t[7]-t[0];

}

#endif

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------

void mpi_exchange_level(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field, int cmp_keys, int level)
{
  int i;
  int icpu;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_PACKET);
  int mpitag=1;

  double t[10];
  double tot;
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));
  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));

  // ---------- The key calculation may already been computed (cmp_key=0) or must be recomputed (cmp_key=1)

  if(cmp_keys){
    // ----------- 0  / we clean the mpi buffers
    clean_mpibuff(cpu,sendbuffer,recvbuffer);

    // ----------- 0.5  / we allocate the number of octs to transmit
    if(cpu->nsend!=NULL){
      free(cpu->nsend);
      free(cpu->nrecv);
    }

    cpu->nsend=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to send to each neighbor
    cpu->nrecv=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to receive from each neighbor


    // ----------- I  / we compute the boundary keys and store them in recvbuffer
    compute_bndkeys_level(cpu,recvbuffer,level,cpu->nrecv);

    // ----------- I bis/ we send the number of requested packets to each server

    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Isend(cpu->nrecv+i,1,MPI_INT,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
    }
    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Irecv(cpu->nsend+i,1,MPI_INT,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
    }
    MPI_Waitall(2*cpu->nnei,req,stat);
    MPI_Barrier(cpu->comm);
    

    // ----------- II / we send the keys to the server

    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Isend(recvbuffer[i],*(cpu->nrecv+i),MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
    }
    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Irecv(sendbuffer[i],*(cpu->nsend+i),MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
    }
    MPI_Waitall(2*cpu->nnei,req,stat);
    MPI_Barrier(cpu->comm);
  }

  // ----------- III/ the server gather the data
  gather_mpi_level(cpu, sendbuffer, field,level);

  //if(cpu->rank==0) printf("--- X ---\n");
  memset(req,0,2*cpu->nnei*sizeof(MPI_Request));
  memset(stat,0,2*cpu->nnei*sizeof(MPI_Status));
  MPI_Barrier(cpu->comm);

  // ----------- IV / the server send the data back to the client


  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],*(cpu->nsend+i),MPI_PACKET,cpu->mpinei[i],cpu->rank         ,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],*(cpu->nrecv+i),MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);



  t[5]=MPI_Wtime();
  // ----------- V  / the client scatter the data back in the oct tree
  scatter_mpi_level(cpu,recvbuffer,field,level);
  
  t[6]=MPI_Wtime();

  //
  free(req);
  free(stat);

  t[7]=MPI_Wtime();
  tot=t[7]-t[0];
}


//------------------------------------------------------------------------
void mpi_cic_correct(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field)
{
  int i,icpu;
  //MPI_Status stat;
  int mpitag=1;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_PACKET);

  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));

  printf("correcting CIC on rank %d\n",cpu->rank);

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff(cpu,sendbuffer,recvbuffer);

  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  MPI_Barrier(cpu->comm);
  gather_ex(cpu,sendbuffer,field);
  MPI_Barrier(cpu->comm);

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // ---------  third we scatter the data back to the INTERNAL boundary octs
  if(field==1) field=3; // fix for an offset of marked scatter
  scatter_mpi(cpu,recvbuffer,field);
  MPI_Barrier(cpu->comm);

  //
  free(stat);
  free(req);

  
}


// ==============================================================================
// ==============================================================================

void mpi_hydro_correct(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer, struct HYDRO_MPI **recvbuffer,int level)
{
  int i,icpu;
  //MPI_Status stat;
  int mpitag=1;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_HYDRO);

  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));

  printf("correcting hydro CIC on rank %d\n",cpu->rank);

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff_hydro(cpu,sendbuffer,recvbuffer);

  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  MPI_Barrier(cpu->comm);
  gather_ex_hydro(cpu,sendbuffer);
  MPI_Barrier(cpu->comm);

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // ---------  third we scatter the data back to the INTERNAL boundary octs
  scatter_mpi_hydro_ext(cpu,recvbuffer,level);
  MPI_Barrier(cpu->comm);

  //
  free(stat);
  free(req);

  
}

// ==================================================================
#ifdef WRAD
void mpi_rad_correct(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer, struct RAD_MPI **recvbuffer,int level)
{
  int i,icpu;
  //MPI_Status stat;
  int mpitag=1;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_RAD);

  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));

  printf("correcting Rad CIC on rank %d\n",cpu->rank);

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff_rad(cpu,sendbuffer,recvbuffer);

  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  MPI_Barrier(cpu->comm);
  gather_ex_rad(cpu,sendbuffer);
  MPI_Barrier(cpu->comm);

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nbuff,MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // ---------  third we scatter the data back to the INTERNAL boundary octs
  scatter_mpi_rad_ext(cpu,recvbuffer,level);
  MPI_Barrier(cpu->comm);

  // ============================
  free(stat);
  free(req);
}

#endif

#ifdef PIC
int mpi_exchange_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer, struct PART **lastpart){

  int nrem,nadd;
  int mpitag=1;
  int i;

  clean_mpibuff_part(cpu,psendbuffer,precvbuffer);
  nrem=gather_ex_part(cpu,psendbuffer,lastpart);
  MPI_Barrier(cpu->comm);

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors
    mpitag=cpu->rank+cpu->mpinei[i];
    MPI_Sendrecv(psendbuffer[i],cpu->nbuff,*cpu->MPI_PART,cpu->mpinei[i],mpitag,precvbuffer[i],cpu->nbuff,*cpu->MPI_PART,cpu->mpinei[i],mpitag,cpu->comm,MPI_STATUS_IGNORE);
  }
  
  nadd=scatter_mpi_part(cpu,precvbuffer,lastpart);
  MPI_Barrier(cpu->comm);

  // Return delta part
  
  return nadd-nrem;
}
#endif

#endif
#endif 
