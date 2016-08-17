#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include "oct.h"
#include "cic.h"
#include "segment.h"
#include "particle.h"
#include "hydro_utils.h"
#include "tools.h"

#ifdef WMPI
#include <mpi.h>
#endif

void setup_mpi(struct CPUINFO *cpu, struct OCT **firstoct, int levelmax, int levelcoarse, int ngridmax,int loadb){

  const int debug=0;

  int nnei=0;
  int *neicpu; // the reduced list of neighbors CPU (not redundant);
  unsigned long long hidx;
  int inidx; // counts the number of inner boundaries octs
  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *desoct;
  struct OCT *newoct;
  unsigned long long key;
  int inei;
  int i,j;
  int nbnd;
  int icell;
  int *flagcpu;
  int *nfromcpu;


  flagcpu=(int*) calloc(cpu->nproc,sizeof(int));
  nfromcpu=(int*) calloc(cpu->nproc,sizeof(int));
  neicpu=(int*) calloc(cpu->nproc,sizeof(int));

#ifdef PIC
  int *npartfromcpu;
  npartfromcpu=(int*) calloc(cpu->nproc,sizeof(int));
#endif


  if (debug) printf("Cleaning the hash table\n");

  memset(cpu->htable,0,cpu->maxhash*sizeof(struct OCT*));
  memset(cpu->bndoct,0,cpu->nbufforg*sizeof(struct OCT*));

#ifdef WMPI
  MPI_Barrier(cpu->comm);
  if (debug) if(cpu->rank==RANK_DISP) printf("Cleaning MPI Buffers\n");

  // we clean the comm buffers

  if(cpu->sendbuffer!=NULL) {
    for(i=0;i<cpu->nnei;i++){
      free(cpu->sendbuffer[i]);
      cpu->sendbuffer[i]=NULL;
    }
    for(i=0;i<cpu->nnei;i++){
      free(cpu->recvbuffer[i]);
      cpu->recvbuffer[i]=NULL;
    }
    free(cpu->sendbuffer);
    free(cpu->recvbuffer);

    cpu->sendbuffer=NULL;
    cpu->recvbuffer=NULL;
  }

#ifdef PIC

  if(cpu->psendbuffer!=NULL) {
    for(i=0;i<cpu->nnei;i++) free(cpu->psendbuffer[i]);
    for(i=0;i<cpu->nnei;i++) free(cpu->precvbuffer[i]);
    free(cpu->psendbuffer);
    free(cpu->precvbuffer);

    cpu->psendbuffer=NULL;
    cpu->precvbuffer=NULL;
  }

  MPI_Barrier(cpu->comm);
#endif

#ifdef WHYDRO2
  if(cpu->hsendbuffer!=NULL) {
    for(i=0;i<cpu->nnei;i++) free(cpu->hsendbuffer[i]);
    for(i=0;i<cpu->nnei;i++) free(cpu->hrecvbuffer[i]);

    free(cpu->hsendbuffer);
    free(cpu->hrecvbuffer);
    cpu->hsendbuffer=NULL;
    cpu->hrecvbuffer=NULL;

  }

#endif

#ifdef WRAD

  if(cpu->Rsendbuffer!=NULL) {
    for(i=0;i<cpu->nnei;i++) free(cpu->Rsendbuffer[i]);
    for(i=0;i<cpu->nnei;i++) free(cpu->Rrecvbuffer[i]);

    free(cpu->Rsendbuffer);
    free(cpu->Rrecvbuffer);
    cpu->Rsendbuffer=NULL;
    cpu->Rrecvbuffer=NULL;

  }


#endif

  MPI_Barrier(cpu->comm);
#endif // WMPI

  // looking for neighbors
  if (debug) if(cpu->rank==RANK_DISP) printf("Getting MPI Neigbourhs\n");

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
	      nfromcpu[curoct->cpu]++;
	      cpu->bndoct[nnei]=curoct; // contains the oct adresses for emission

	      //if(cpu->rank==RANK_DISP) printf("levl=%d nei=%d\n",curoct->level,curoct->cpu);
	      nnei++;
 	      if(nfromcpu[curoct->cpu]==cpu->nbufforg) {
		printf("ERROR on nbufforg being too small for octs on level =%d!\n",level);
		abort();
	      }
 	      if(nnei==cpu->nbufforg) {
		printf("ERROR on nbufforg being too small for octs on level =%d!\n",level);
		abort();
	      }
#ifdef PIC
	      // counting the particles
	      for(icell=0;icell<8;icell++) npartfromcpu[curoct->cpu]+=countpart(curoct->cell[icell].phead);
	      if(npartfromcpu[curoct->cpu]==cpu->nbufforg) {
		printf("ERROR on nbufforg being too small for particles on level =%d!\n",level);
		abort();
	      }
#endif // PIC

	    }
	  }while(nextoct!=NULL);
	}
      }

    }


  // =========================================== SETTING UP THE COMMUNICATIONS BETWEEN NEIGHBORS

#ifdef WMPI
  MPI_Barrier(cpu->comm);
  //if(cpu->rank==RANK_DISP) printf("Set up MPI Comms\n");
#endif


  // computing the mpi neighbor list
  for(i=0;i<cpu->nproc;i++) neicpu[i]=500000+i; // large enough to be at the end
  j=0;
  int maxn=0;
  int maxnpart=0;
  for(i=0;i<cpu->nproc;i++){ // we scan the list
    if(flagcpu[i]==1){
      neicpu[j]=i;
      j++;
      if(nfromcpu[i]>maxn) maxn=nfromcpu[i];
#ifdef PIC
      if(npartfromcpu[i]>maxnpart) maxnpart=npartfromcpu[i];
#endif
    }
  }


  free(flagcpu);
  free(nfromcpu);
#ifdef PIC
  free(npartfromcpu);
#endif
  myradixsort(neicpu,cpu->nproc); // we sort the neighbors list


  nbnd=nnei;
  nnei=j;

  cpu->nebnd=nbnd;
  cpu->nnei=nnei;
  if(cpu->mpinei!=NULL){
    free(cpu->mpinei);
    cpu->mpinei=NULL;
  }

  cpu->mpinei=(int*)calloc(nnei,sizeof(int)); // we reallocate the array to free some memory
  for(i=0;i<cpu->nnei;i++) cpu->mpinei[i]=neicpu[i];
  free(neicpu);

  // AT THIS STAGE:
  // nbnd contains the number of boundary octs
  // nnei contains the number of neighbor procs
  // mpinei contains the rank of the neighbor procs and has a size nnei


  // some displays
  /* printf("Found %d neighbors and %d bnd octs  on rank %d :",cpu->nnei,cpu->nebnd,cpu->rank); */
  /* for(i=0;i<cpu->nnei;i++) printf("%d ",cpu->mpinei[i]); */
  /* printf("\n"); */

#ifdef WMPI
  MPI_Barrier(cpu->comm);
  int nvois_loc=cpu->nebnd;
  int nvois_max;
  int nvois_maxpart;
  MPI_Allreduce(&nvois_loc,&nvois_max,1,MPI_INT,MPI_MAX,cpu->comm);
  //if(cpu->rank==RANK_DISP) printf("Max total number of neighbors octs=%d\n",nvois_max);
  MPI_Allreduce(&maxn,&nvois_max,1,MPI_INT,MPI_MAX,cpu->comm);
  //if(cpu->rank==RANK_DISP) printf("Max number of neighbors octs from 1 nei=%d\n",nvois_max);
  //if(cpu->rank==RANK_DISP) printf("Overriding param nbuff=%d with maxn =%d\n",cpu->nbuff,nvois_max);
  cpu->nbuff=nvois_max;

#ifdef PIC
  MPI_Allreduce(&maxnpart,&nvois_maxpart,1,MPI_INT,MPI_MAX,cpu->comm);
  //if(cpu->rank==RANK_DISP) printf("Max number of neighbors particles from 1 nei=%d\n",nvois_maxpart);
  //if(cpu->rank==RANK_DISP) printf("Overriding param nbuffpart=%d with maxnpart =%d\n",cpu->nbuffpart,nvois_maxpart);
  cpu->nbuffpart=nvois_maxpart;
#endif

  MPI_Barrier(cpu->comm);
  // ----------- allocating the communication buffers
  //  printf("%p %p\n", cpu->sendbuffer, cpu->recvbuffer);

#if 1
  cpu->sendbuffer=(struct PACKET **)(calloc(cpu->nnei,sizeof(struct PACKET*)));
  cpu->recvbuffer=(struct PACKET **)(calloc(cpu->nnei,sizeof(struct PACKET*)));
  for(i=0;i<cpu->nnei;i++) {
    *(cpu->sendbuffer+i)=(struct PACKET *) (calloc(cpu->nbuff,sizeof(struct PACKET)));
    *(cpu->recvbuffer+i)=(struct PACKET *) (calloc(cpu->nbuff,sizeof(struct PACKET)));
  }

#ifdef PIC
  cpu->psendbuffer=(struct PART_MPI **)(calloc(cpu->nnei,sizeof(struct PART_MPI*)));
  cpu->precvbuffer=(struct PART_MPI **)(calloc(cpu->nnei,sizeof(struct PART_MPI*)));
  for(i=0;i<cpu->nnei;i++) {
    *(cpu->psendbuffer+i)=(struct PART_MPI *) (calloc(cpu->nbuffpart,sizeof(struct PART_MPI)));
    *(cpu->precvbuffer+i)=(struct PART_MPI *) (calloc(cpu->nbuffpart,sizeof(struct PART_MPI)));
  }
#endif

#ifdef WHYDRO2
  cpu->hsendbuffer=(struct HYDRO_MPI **)(calloc(cpu->nnei,sizeof(struct HYDRO_MPI*)));
  cpu->hrecvbuffer=(struct HYDRO_MPI **)(calloc(cpu->nnei,sizeof(struct HYDRO_MPI*)));
  for(i=0;i<cpu->nnei;i++) {
    cpu->hsendbuffer[i]=(struct HYDRO_MPI *) (calloc(cpu->nbuff,sizeof(struct HYDRO_MPI)));
    cpu->hrecvbuffer[i]=(struct HYDRO_MPI *) (calloc(cpu->nbuff,sizeof(struct HYDRO_MPI)));
  }
#endif

#ifdef WRAD

  cpu->Rsendbuffer=(struct RAD_MPI **)(calloc(cpu->nnei,sizeof(struct RAD_MPI*)));
  cpu->Rrecvbuffer=(struct RAD_MPI **)(calloc(cpu->nnei,sizeof(struct RAD_MPI*)));
  for(i=0;i<cpu->nnei;i++) {
    cpu->Rsendbuffer[i]=(struct RAD_MPI *) (calloc(cpu->nbuff,sizeof(struct RAD_MPI)));
    cpu->Rrecvbuffer[i]=(struct RAD_MPI *) (calloc(cpu->nbuff,sizeof(struct RAD_MPI)));
  }

#endif // WRAD
#endif // 1

#endif // WMPI


  // creating a cpu dictionnary to translate from cpu number to inei

  if(cpu->dict!=NULL){
    free(cpu->dict);
    cpu->dict=NULL;
  }
  cpu->dict=(int*)calloc(cpu->nproc,sizeof(int));
  for(i=0;i<cpu->nproc;i++) cpu->dict[i]=-1;
  for(i=0;i<cpu->nnei;i++){
    cpu->dict[cpu->mpinei[i]]=i;
  }

#ifdef WMPI
MPI_Barrier(cpu->comm);
#endif
 if(cpu->rank==RANK_DISP) printf("setup mpi Done\n");

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
    pack->key=(double) oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    for(ii=0;ii<8;ii++){
      switch(field){
#ifdef PIC
      case 0:
	pack->data[ii]=cpu->bndoct[i]->cell[ii].density;
	break;
#endif
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


//========================================================================
//  ========================================================================

void gather_ex_level(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field,int level){
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
    if(cpu->bndoct[i]->level!=level) continue; // only l-1 octs must be corrected from flux diffusion
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer

    // assigning the values
    pack->level=cpu->bndoct[i]->level;
    pack->key=(double) oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    for(ii=0;ii<8;ii++){
      switch(field){
#ifdef PIC
      case 0:
	pack->data[ii]=cpu->bndoct[i]->cell[ii].density;
	break;
#endif
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


#ifdef WGRAV
void gather_ex_dens(struct CPUINFO *cpu, struct PACKET **sendbuffer, int level){

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
    if(cpu->bndoct[i]->level!=(level)) continue; // only l-1 octs must be corrected from flux diffusion
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer

    // assigning the values
    pack->level=cpu->bndoct[i]->level;
    pack->key=(double) oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    for(ii=0;ii<8;ii++){
#ifdef PIC
	pack->data[ii]=cpu->bndoct[i]->cell[ii].density;
#endif // PIC
    }

    // counting the number of packets for icpu
    countpacket[icpu]++;
  }

  //for(i=0;i<cpu->nnei;i++) printf("rank=%d cpu %d nbnd=%d\n",cpu->rank,cpu->mpinei[i],countpacket[i]);
  free(countpacket);

}
#endif

void gather_ex_mark(struct CPUINFO *cpu, struct PACKET **sendbuffer, int level){

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
    if(cpu->bndoct[i]->level!=(level)) continue; // only l-1 octs must be corrected from flux diffusion
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer

    // assigning the values
    pack->level=cpu->bndoct[i]->level;
    pack->key=(double) oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    for(ii=0;ii<8;ii++){
	pack->data[ii]=cpu->bndoct[i]->cell[ii].marked;
    }

    // counting the number of packets for icpu
    countpacket[icpu]++;
  }

  //for(i=0;i<cpu->nnei;i++) printf("rank=%d cpu %d nbnd=%d\n",cpu->rank,cpu->mpinei[i],countpacket[i]);
  free(countpacket);

}




#ifdef WHYDRO2
void gather_ex_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer,int level, int *countpacket){
  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
    USed for coarse cell update of hydro at boundaries
   */
  int icpu;
  struct HYDRO_MPI *pack;
  int i,icell;

  memset(countpacket,0,cpu->nnei*sizeof(int));

  if(cpu->noct[level-1]>0){ // correction must be done only if the current level exists
    for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
      if(cpu->bndoct[i]->level!=(level-1)) continue; // only l-1 octs must be corrected from flux diffusion
      icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
      pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer


      // assigning the values
      pack->level=cpu->bndoct[i]->level;
      pack->key=(double)oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
      //for(ii=0;ii<8;ii++) pack->data[ii]=cpu->bndoct[i]->cell[ii].fieldnew;

      for(icell=0;icell<8;icell++){
	memcpy(&(pack->data[icell]),&(cpu->bndoct[i]->cell[icell].fieldnew),sizeof(struct Wtype));
      }
      // counting the number of packets for icpu
      countpacket[icpu]++;
    }
  }
}
#endif
// ====================================================================================
// ====================================================================================

#ifdef WRAD
void gather_ex_rad(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer,int level, int *countpacket){

  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
    USed for coarse cell update of hydro at boundaries
   */

  int icpu;
  struct RAD_MPI *pack;
  int i,icell;

  memset(countpacket,0,cpu->nnei*sizeof(int));

  if(cpu->noct[level-1]>0){ // correction must be done only if the current level exists
  for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
    if(cpu->bndoct[i]->level!=(level-1)) continue; // only l-1 octs must be corrected from flux diffusion
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    pack=sendbuffer[icpu]+countpacket[icpu]; // taking a free slot in the send buffer

    // assigning the values
    pack->level=cpu->bndoct[i]->level;
    pack->key=(double)oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current oct
    /* for(ii=0;ii<8;ii++){ */
    /*   pack->data[ii]=cpu->bndoct[i]->cell[ii].rfieldnew; */
    /* } */
    for(icell=0;icell<8;icell++){
      memcpy(&(pack->data[icell]),&(cpu->bndoct[i]->cell[icell].rfieldnew),sizeof(struct Rtype));
    }

    // counting the number of packets for icpu
    countpacket[icpu]++;
  }
  }

}
#endif

//======================================================================================
#ifdef PIC
void gather_ex_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, int *nrem){

  /*
    NOTE: this one is peculiar since the keys are directly computed from the source of data
   */
  int *countpacket;
  int icpu;
  struct PART_MPI *part;
  struct PART *curp;
  struct PART *nexp;
  int i,ii;
  nrem[0]=0;
  nrem[1]=0;

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
	  part->key=(double) oct2key(cpu->bndoct[i],cpu->bndoct[i]->level); // getting the key of the current o509ct (is the eventual destination)
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
	  part->is=curp->is;

#ifdef STARS
	  part->rhocell=curp->rhocell;
	  part->age=curp->age;
	  part->isStar=curp->isStar;
    part->radiative_state=curp->radiative_state;
#endif
	  //if(cpu->rank==1) printf(" id send =%d\n",part->idx);

	  // counting the number of packets for icpu
	  countpacket[icpu]++;
	  if(countpacket[icpu]>cpu->nbuffpart){
	    printf("cpu->nbuffpart %d blasted by particles number %d",cpu->nbuffpart,countpacket[icpu]);
	    abort();
	  }
	  // switching the mass to -1 to flag exiting particles
	  curp->mass=-1.;

#ifdef STARS
	if(part->isStar){
		nrem[1] ++;
	//	printf("cpu %d nremstar %d\n",cpu->rank,nrem[1]);
	}
#endif
	  nrem[0]++;

	  // disconnecting the particle from its list

	  curp->prev=NULL;
	  curp->next=NULL;

	  if(nexp==NULL){ // reached the last particle from the cell
	    cpu->bndoct[i]->cell[ii].phead=NULL; // we "clear" the cell
	  }

	  // fixing the freepart list
	  curp->prev=NULL;
	  cpu->freepart->prev=curp;
	  curp->next=cpu->freepart;
	  cpu->freepart=curp;


	}while(nexp!=NULL);
      }
    }
  }
  free(countpacket);

  //printf("nrem=%d\n",nrem[0]);
  //return nrem; // we return the number of exiting particles
}
#endif


//======================================================================================

void gather_mpi(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field){

  int i,j;
  int found=0;
  unsigned long long hidx;
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
	  hidx=hfun((unsigned long long) pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==(unsigned long long) pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      //curoct->border=1;

	      for(icell=0;icell<8;icell++){
		switch(field){
#ifdef WGRAV
#ifdef PIC
		case 0:
		  pack->data[icell]=curoct->cell[icell].density; // density
		  break;
		case 1:
		  pack->data[icell]=curoct->cell[icell].density; // density again we reproduce the case 1 in order to be consistent with scatter_mpi
		  break;
		case 4:
		  //pack->data[icell]=curoct->cell[icell].temp; //temp field for force calculation
		  break;
#endif
		case 2:
		  pack->data[icell]=curoct->cell[icell].gdata.p; // potential
		  break;
#endif
		case 3:
		  pack->data[icell]=curoct->cell[icell].marked*(curoct->level>=cpu->levelcoarse); //refinment mark
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
  unsigned long long hidx;
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
	  hidx=hfun((unsigned long long) pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==(unsigned long long) pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      //curoct->border=1;

	      for(icell=0;icell<8;icell++){
		switch(field){
#ifdef WGRAV
#ifdef PIC
		case 0:
		  pack->data[icell]=curoct->cell[icell].density; // density
		  break;
		case 1:
		  pack->data[icell]=curoct->cell[icell].density; // density again we reproduce the case 1 in order to be consistent with scatter_mpi
		  break;
		case 4:
		  //pack->data[icell]=curoct->cell[icell].temp; //temp field for force calculation
		  break;
#endif
		case 2:
		  pack->data[icell]=curoct->cell[icell].gdata.p; // potential
		  break;
#endif
		case 3:
		  pack->data[icell]=curoct->cell[icell].marked; //refinment mark
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

// ======================================================
// ======================================================
#ifdef WGRAV
void gather_mpi_pot_level(struct CPUINFO *cpu, struct PACKET **sendbuffer, int level){
  int i,j;
  int found=0;
  unsigned long long hidx;
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
	  hidx=hfun((unsigned long long) pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==(unsigned long long) pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      //curoct->border=1;

	      for(icell=0;icell<8;icell++){
#ifdef WGRAV
		pack->data[icell]=curoct->cell[icell].gdata.p; // potential
#endif // WGRAV
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
#endif


 //------------------------------------------------------------------------
void scatter_mpi(struct CPUINFO *cpu, struct PACKET **recvbuffer,  int field){

  int i,j;
  int found=0;
  unsigned long long hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long)pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long) pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      switch(field){
#ifdef WGRAV
#ifdef PIC
	      case 0:
		curoct->cell[icell].density+=pack->data[icell]; // density += for CIC correction
		break;
	      case 1:
		curoct->cell[icell].density =pack->data[icell]; // density
		break;
	      case 4:
		//curoct->cell[icell].temp=pack->data[icell]; // temp field for force calculation
		break;
#endif
	      case 2:
		curoct->cell[icell].gdata.p=pack->data[icell]; // potential
		break;
#endif
	      case 3:
		curoct->cell[icell].marked=FMAX(pack->data[icell],(REAL)curoct->cell[icell].marked); // refinement mark
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
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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
  unsigned long long hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level==level){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long)pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      switch(field){
#ifdef WGRAV
#ifdef PIC
	      case 0:
		curoct->cell[icell].density+=pack->data[icell]; // density += for CIC correction
		break;
	      case 1:
		curoct->cell[icell].density =pack->data[icell]; // density
		break;
	      case 4:
		//curoct->cell[icell].temp=pack->data[icell]; // temp field for force calculation
		break;
#endif
	      case 2:
		curoct->cell[icell].gdata.p=pack->data[icell]; // potential
		break;
#endif
	      case 3:
		curoct->cell[icell].marked=FMAX(pack->data[icell],(REAL)curoct->cell[icell].marked); // refinement mark
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
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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

// ==================================================================================
#ifdef WGRAV
void scatter_mpi_pot_level(struct CPUINFO *cpu, struct PACKET **recvbuffer,  int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level==level){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long)pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
#ifdef WGRAV
		curoct->cell[icell].gdata.p=pack->data[icell]; // potential
#endif // WGRAV
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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
#endif

 //------------------------------------------------------------------------
#ifdef PIC
void scatter_mpi_part(struct CPUINFO *cpu, struct PART_MPI **precvbuffer, int *nadd, int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
  struct PART_MPI *part;
  struct PART *curp;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  nadd[0]=0;
  nadd[1]=0;

  struct PART *lastp;


  /* if(SOCT!=NULL){ */
  /*   printf("test SOCT level=%d key=%llu on rank %d\n",level,oct2key(SOCT,12),cpu->rank); */
  /* } */
  /* else{ */
  /*   printf("SOCT NULL at level=%d on rank=%d\n",level,cpu->rank); */
  /* } */

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuffpart;i++){
      part=precvbuffer[j]+i;
      //      if(part==NULL) continue;

      //if(cpu->rank==0) if(i<50) printf("i=%d l=%d idx=%d on rank%d\n",i,part->level,part->idx,cpu->rank);

      if(part->level!=0){ // we do something

	/* if(part->level!=level){ */
	/*   printf("level=%d plevel=%d rg part=%d cpu=%d key=%e\n",level,part->level,i,cpu->rank,part->key); */
	/*   abort(); */
	/* } */

	/* if(part->level==12){ */
	/*   printf("RECEPTION level 12 particle found with key=%e with pos= %e %e %e\n",part->key,part->x,part->y,part->z); */
	/* } */

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long) part->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)part->key)&&(part->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    // we assign the new particle to the global particle list of the client

	    lastp=cpu->freepart;// getting the next slot in the global particle list
	    if(cpu->freepart->next==NULL){
	      // particles full
	      printf("Particle array full, increase npartmax\n");
	      abort();
	    }
	    cpu->freepart=cpu->freepart->next;
	    cpu->freepart->prev=NULL;

	    //if(cpu->firstoct[0]->next!=NULL) abort();

	    // copying the data
	    (lastp)->x =part->x;
	    (lastp)->y =part->y;
	    (lastp)->z =part->z;

	    (lastp)->vx=part->vx;
	    (lastp)->vy=part->vy;
	    (lastp)->vz=part->vz;


	    (lastp)->mass=part->mass;

	    (lastp)->idx=part->idx;
	    (lastp)->level=part->level;
	    (lastp)->is=part->is;

#ifdef STARS
	    (lastp)->rhocell=part->rhocell;
	    (lastp)->age=part->age;
	    (lastp)->isStar=part->isStar;
	    (lastp)->radiative_state=part->radiative_state;

	if ((lastp)->isStar) {
		nadd[1] ++;
	//	printf("cpu %d naddstar %d\n",cpu->rank, nadd[1]);
	}
#endif
	    nadd[0]++;

	    //if(cpu->firstoct[0]->next!=NULL) abort();

	    // if recepction cell is refined
	    struct CELL *newcell=&(curoct->cell[part->icell]);

	    if(newcell->child!=NULL){
	      struct OCT *newoct=newcell->child;
	      REAL dxcur2=1./POW(2.,newoct->level);

	      int xp=(int)(DFACT*(part->x-newoct->x)/dxcur2);xp=(xp>1?1:xp);xp=(xp<0?0:xp);
	      int yp=(int)(DFACT*(part->y-newoct->y)/dxcur2);yp=(yp>1?1:yp);yp=(yp<0?0:yp);
	      int zp=(int)(DFACT*(part->z-newoct->z)/dxcur2);zp=(zp>1?1:zp);zp=(zp<0?0:zp);
	      int ip=xp+yp*2+zp*4;

	      newcell=&(newoct->cell[ip]);
	      part->icell=ip;
	      (lastp)->level=part->level+1;
	    }

	    //if(cpu->firstoct[0]->next!=NULL) abort();

	    // looking for the tail particle in the destination cell
	    curp=findlastpart(newcell->phead);
	    if(curp!=NULL){
	      curp->next=(lastp);
	      (lastp)->next=NULL;
	      (lastp)->prev=curp;
	    }
	    else{
	      newcell->phead=(lastp);
	      (lastp)->next=NULL;
	      (lastp)->prev=NULL;
	    }
	    //if(cpu->firstoct[0]->next!=NULL) abort();

	  }
	  else{
	    printf("error no reception oct found !");

	    hidx=hfun((unsigned long long)part->key,cpu->maxhash);
	    nextoct=cpu->htable[hidx];
	    found=0;
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==(unsigned long long) part->key)&&(part->level==curoct->level));
	      printf("o2key=%llu part.key=%llu part.level=%d curoct.level=%d hidx=%llu hidx ex=%llu maxhash=%d, curoct cpu=%d\n",oct2key(curoct,curoct->level),(unsigned long long) part->key,part->level,curoct->level,hidx,hfun(oct2key(curoct,curoct->level),cpu->maxhash),cpu->maxhash,curoct->cpu);
#ifdef STARS
	      printf("is=%d is star=%d age=%e\n",part->is,part->isStar,part->age);
#endif
	    }while((nextoct!=NULL)&&(!found));

	    abort();
	  }
	}else{
	  printf("error no hash key obtained !!\n");
	  abort();
	}
      }
    }
  }

  // --- check
  int icpu;
  struct PART *nexp;
  int ii;
  for(i=0;i<cpu->nebnd;i++){ // scanning all the external boundary octs
    icpu=cpu->dict[cpu->bndoct[i]->cpu]; // getting the local destination cpu through the dictionnary
    for(ii=0;ii<8;ii++){
      nexp=cpu->bndoct[i]->cell[ii].phead;
      // sweeping the particles of the boundary cells
      if(nexp!=NULL){
	printf("hey!\n");
	abort();
      }
    }
  }

  //printf("nadd=%d\n",nadd[0]);
  //return nadd; // returning the new number of particles
}
#endif

//------------------------------------------------------------------------
void compute_bndkeys(struct CPUINFO *cpu, struct PACKET **recvbuffer){

  // we create a counter of values for each neighbor
  int *countpacket;
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct PACKET *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    unsigned long long keyloc;
    int cpuloc;
    int inei;
    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=(double) keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }

  free(countpacket);
}

// ===========================================================================
// ===========================================================================

void compute_bndkeys_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **recvbuffer){

  // we create a counter of values for each neighbor
  int *countpacket;
  countpacket=(int*)calloc(cpu->nnei,sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct HYDRO_MPI *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    unsigned long long keyloc;
    int cpuloc;
    int inei;
    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=(double)keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }

  free(countpacket);
}


void compute_bndkeys_hydro_level(struct CPUINFO *cpu, int level, struct HYDRO_MPI **recvbuffer, int *countpacket){

  /* // we create a counter of values for each neighbor */
  /* int *countpacket; */
  /* countpacket=(int*)calloc(cpu->nnei,sizeof(int)); */

  memset(countpacket,0,cpu->nnei*sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct HYDRO_MPI *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    unsigned long long keyloc;
    int cpuloc;
    int inei;

    if(cpu->bndoct[i]->level != level) continue; // we skip the borders which do not belong to current level

    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=(double)keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }
}


#ifdef WRAD

// ======================================================================================
// ======================================================================================

void compute_bndkeys_rad_level(struct CPUINFO *cpu, int level, struct RAD_MPI **recvbuffer, int *countpacket){

  // we create a counter of values for each neighbor
  /* int *countpacket; */
  /* countpacket=(int*)calloc(cpu->nnei,sizeof(int)); */

  memset(countpacket,0,cpu->nnei*sizeof(int));

  // filling the keys in the reception buffer (which will be processed by the remote cpus)
  struct RAD_MPI *pack;

  int i;
  for(i=0;i<cpu->nebnd;i++) {
    unsigned long long keyloc;
    int cpuloc;
    int inei;

    if(cpu->bndoct[i]->level != level) continue; // we skip the borders which do not belong to current level

    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=(double)keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }

}

#endif

// ================================================================================
// ================================================================================


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
    unsigned long long keyloc;
    int cpuloc;
    int inei;

    if(cpu->bndoct[i]->level != level) continue; // we skip the borders which do not belong to current level

    keyloc=oct2key(cpu->bndoct[i],cpu->bndoct[i]->level);
    cpuloc=cpu->bndoct[i]->cpu;
    inei=cpu->dict[cpuloc]; // we recover the local neighbor index by using the dictionnary

    pack=recvbuffer[inei]+countpacket[inei]; // we get the pack
    pack->key=(double)keyloc;
    pack->level=cpu->bndoct[i]->level;

    countpacket[inei]++;
  }

  //free(countpacket);
}

//------------------------------------------------------------------------
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



void  clean_mpibuff_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer){
  int i;
  for(i=0;i<cpu->nnei;i++) {
    memset(psendbuffer[i],0,cpu->nbuffpart*sizeof(struct PART_MPI));
    memset(precvbuffer[i],0,cpu->nbuffpart*sizeof(struct PART_MPI));
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

  //if(cpu->rank==RANK_DISP) printf("--- X ---\n");
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
  //if(cpu->rank==RANK_DISP) printf("clean=%e keys=%e sendkeys=%e gather=%e senddata=%e scatter=%e free=%e\n",(t[1]-t[0])/tot,(t[2]-t[1])/tot,(t[3]-t[2])/tot,(t[4]-t[3])/tot,(t[5]-t[4])/tot,(t[6]-t[5])/tot,(t[7]-t[6])/tot);

}

//===========================================================================================
//===========================================================================================
//===========================================================================================
// ============================================================
void scatter_mpi_mark_ext(struct CPUINFO *cpu, struct PACKET **recvbuffer,int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long)pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    if(curoct->level!=(level)) continue; // we update only coarse neighbours relative to the current level
	    for(icell=0;icell<8;icell++){
	      curoct->cell[icell].marked=FMAX(pack->data[icell],curoct->cell[icell].marked);
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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

#ifdef WGRAV
void scatter_mpi_dens_ext(struct CPUINFO *cpu, struct PACKET **recvbuffer,int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
  struct PACKET *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long)pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    if(curoct->level!=(level)) continue; // we update only coarse neighbours relative to the current level
	    for(icell=0;icell<8;icell++){
#ifdef PIC
	      curoct->cell[icell].density+=(pack->data[icell],curoct->cell[icell].density);
#endif // PIC
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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
#endif

#ifdef WHYDRO2
void gather_mpi_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer){

  int i,j;
  int found=0;
  unsigned long long hidx;
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
	  hidx=hfun((unsigned long long)pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      //curoct->border=1;

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

// ===============================================================
// ===============================================================

void gather_mpi_hydro_level(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer, int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
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
	if(pack->level==level){ // we do something
	  //t[0]=MPI_Wtime();
	  countpacket[j]++;

	  // first we compute the adress from the hashfunction
	  hidx=hfun((unsigned long long) pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==(unsigned long long) pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      //curoct->border=1;

	      for(icell=0;icell<8;icell++){
		memcpy(&(pack->data[icell]),&(curoct->cell[icell].field),sizeof(struct Wtype));
		if(pack->data[icell].d==0.){
		  printf("NULL dens packet %e %e %e %d %d\n",pack->data[icell].d,pack->data[icell].d,pack->data[icell].p, curoct->level,cpu->rank);
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
//------------------------------------------------------------------------
void scatter_mpi_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **recvbuffer){

  int i,j;
  int found=0;
  unsigned long long hidx;
  struct HYDRO_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level!=0){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long) pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      memcpy(&(curoct->cell[icell].field),&(pack->data[icell]),sizeof(struct Wtype));
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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


// =========================================================================
// =========================================================================

void scatter_mpi_hydro_level(struct CPUINFO *cpu, struct HYDRO_MPI **recvbuffer, int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
  struct HYDRO_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level==level){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long) pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      memcpy(&(curoct->cell[icell].field),&(pack->data[icell]),sizeof(struct Wtype));
	      if(pack->data[icell].d==0.){
		printf("SCATTER NULL dens packet %e %e %e %d %d\n",pack->data[icell].d,pack->data[icell].d,pack->data[icell].p, curoct->level,cpu->rank);
	      }
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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
void scatter_mpi_hydro_ext(struct CPUINFO *cpu, struct HYDRO_MPI **recvbuffer,int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
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
	hidx=hfun((unsigned long long)pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    if(curoct->level!=(level-1)) {
	      continue; // we update only coarse neighbours relative to the current level
	    }
	    for(icell=0;icell<8;icell++){

	      W2U(&(curoct->cell[icell].fieldnew),&U);
	      W2U(&(pack->data[icell]),&Ue);

 	      /* if(level==7) if((Ue.d!=0)&&(Ue.E==0.)) {  */
	      /* 	  printf("rank=%d d=%e E=%e level=%d\n",cpu->rank,Ue.d,Ue.E,pack->level); */
	      /* 	  abort(); */
	      /* 	} */

	      U.d +=Ue.d;
	      U.du+=Ue.du;
	      U.dv+=Ue.dv;
	      U.dw+=Ue.dw;
	      U.E +=Ue.E;
	      U.eint+=Ue.eint;

#ifdef WRADHYD
	      U.dX+=Ue.dX;
#ifdef HELIUM
	      U.dXHE+=Ue.dXHE;
	      U.dXXHE+=Ue.dXXHE;

#endif
#endif
	      U2W(&U,&W);
	      memcpy(&(curoct->cell[icell].fieldnew),&W,sizeof(struct Wtype));

 	if(isnan(curoct->cell[icell].fieldnew.u)){
	  printf("HH %e %e %e %e %e %e | %e %e %e %e %e| %d %d %d\n",curoct->cell[icell].fieldnew.d,curoct->cell[icell].fieldnew.u,curoct->cell[icell].fieldnew.v,curoct->cell[icell].fieldnew.w,curoct->cell[icell].fieldnew.p,curoct->cell[icell].fieldnew.E,Ue.d,Ue.du,Ue.dv,Ue.dw,Ue.E,curoct->cpu,cpu->rank,pack->level);
	  abort();
	}

	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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

    // ----------- 0.5  / we allocate the number of octs to transmit
    if(cpu->nsend!=NULL){
      free(cpu->nsend);
      free(cpu->nrecv);

      free(cpu->nsend_coarse);
      free(cpu->nrecv_coarse);
    }

    cpu->nsend=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to send to each neighbor
    cpu->nrecv=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to receive from each neighbor

    cpu->nsend_coarse=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to send to each neighbor
    cpu->nrecv_coarse=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to receive from each neighbor
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

  //if(cpu->rank==RANK_DISP) printf("--- X ---\n");
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
  //if(cpu->rank==RANK_DISP) printf("clean=%e keys=%e sendkeys=%e gather=%e senddata=%e scatter=%e free=%e\n",(t[1]-t[0])/tot,(t[2]-t[1])/tot,(t[3]-t[2])/tot,(t[4]-t[3])/tot,(t[5]-t[4])/tot,(t[6]-t[5])/tot,(t[7]-t[6])/tot);

}

// ====================================================================================================
// ====================================================================================================

void mpi_exchange_hydro_level(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer, struct HYDRO_MPI **recvbuffer, int cmp_keys, int level)
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

    // ----------- 0.5  / we allocate the number of octs to transmit
    if(cpu->nsend!=NULL){
      free(cpu->nsend);
      free(cpu->nrecv);

      free(cpu->nsend_coarse);
      free(cpu->nrecv_coarse);
    }

    cpu->nsend=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to send to each neighbor
    cpu->nrecv=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to receive from each neighbor

    cpu->nsend_coarse=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to send to each neighbor
    cpu->nrecv_coarse=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to receive from each neighbor


    // ----------- I  / we compute the boundary keys and store them in recvbuffer
    compute_bndkeys_hydro_level(cpu,level,recvbuffer,cpu->nrecv);

    // ----------- II / we send the number of octs required from each neighbors

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
      MPI_Isend(recvbuffer[i],cpu->nrecv[i],MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
    }
    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Irecv(sendbuffer[i],cpu->nsend[i],MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
    }
    MPI_Waitall(2*cpu->nnei,req,stat);
    MPI_Barrier(cpu->comm);

  }

  // ----------- III/ the server gather the data
  gather_mpi_hydro_level(cpu, sendbuffer,level);

  //if(cpu->rank==RANK_DISP) printf("--- X ---\n");
  memset(req,0,2*cpu->nnei*sizeof(MPI_Request));
  memset(stat,0,2*cpu->nnei*sizeof(MPI_Status));
  MPI_Barrier(cpu->comm);

  // ----------- IV / the server send the data back to the client

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nsend[i],MPI_PACKET,cpu->mpinei[i],cpu->rank         ,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nrecv[i],MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }

  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);



  t[5]=MPI_Wtime();
  // ----------- V  / the client scatter the data back in the oct tree
  scatter_mpi_hydro_level(cpu,recvbuffer,level);
  MPI_Barrier(cpu->comm);
  t[6]=MPI_Wtime();

  //
  free(req);
  free(stat);
  t[7]=MPI_Wtime();
  tot=t[7]-t[0];
  //if(cpu->rank==RANK_DISP) printf("clean=%e keys=%e sendkeys=%e gather=%e senddata=%e scatter=%e free=%e\n",(t[1]-t[0])/tot,(t[2]-t[1])/tot,(t[3]-t[2])/tot,(t[4]-t[3])/tot,(t[5]-t[4])/tot,(t[6]-t[5])/tot,(t[7]-t[6])/tot);

}


#endif

//=====================================================================
#ifdef WRAD
void scatter_mpi_rad_ext(struct CPUINFO *cpu, struct RAD_MPI **recvbuffer,int level){

  int i,j,igrp;
  int found=0;
  unsigned long long hidx;
  struct RAD_MPI *pack;
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
	hidx=hfun((unsigned long long)pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long) pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    if(curoct->level!=(level-1)) continue; // we update only coarse neighbours relative to the current level
	    for(icell=0;icell<8;icell++){

	      R=&(curoct->cell[icell].rfieldnew);
	      Re=&(pack->data[icell]);

	      REAL vali=R->e[0];
	      for(igrp=0;igrp<NGRP;igrp++){
		// update
		R->e[igrp] += Re->e[igrp];
		R->fx[igrp]+= Re->fx[igrp];
		R->fy[igrp]+= Re->fy[igrp];
		R->fz[igrp]+= Re->fz[igrp];
	      }

	      /* R->nh=curoct->cell[icell].rfield.nh; */
	      /* R->eint=curoct->cell[icell].rfield.eint; */
	      /* R->src=curoct->cell[icell].rfield.src; */

	      if(R->e[0]<0){ printf(" WARNING new = %e delta=%e old=%e vali=%e\n",R->e[0],Re->e[0],curoct->cell[icell].rfield.e[0],vali);
		printf("current cpu=%d curoct->cpu=%d\n",cpu->rank,curoct->cpu);
		printf("level=%d x=%e y=%e z=%e\n",curoct->level,curoct->x,curoct->y,curoct->z);
		printf("key=%e\n",pack->key);
	      }
	      //memcpy(&(curoct->cell[icell].rfieldnew),R,sizeof(struct Rtype)); // inutile ?
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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

void gather_mpi_rad_level(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer, int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
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
	if(pack->level==level){ // we do something
	  //t[0]=MPI_Wtime();
	  countpacket[j]++;

	  // first we compute the adress from the hashfunction
	  hidx=hfun((unsigned long long) pack->key,cpu->maxhash);
	  nextoct=cpu->htable[hidx];
	  //t[3]=MPI_Wtime();
	  if(nextoct!=NULL){
	    do{ // resolving collisions
	      curoct=nextoct;
	      nextoct=curoct->nexthash;
	      found=((oct2key(curoct,curoct->level)==(unsigned long long) pack->key)&&(pack->level==curoct->level));
	    }while((nextoct!=NULL)&&(!found));

	    //t[1]=MPI_Wtime();
	    if(found){ // the reception oct has been found

	      // we set the current oct as a border one (for vector based communications)
	      //curoct->border=1;

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


// =======================================================================
// =======================================================================
void scatter_mpi_rad_level(struct CPUINFO *cpu, struct RAD_MPI **recvbuffer,int level){

  int i,j;
  int found=0;
  unsigned long long hidx;
  struct RAD_MPI *pack;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  for(j=0;j<cpu->nnei;j++){
    for(i=0;i<cpu->nbuff;i++){
      pack=recvbuffer[j]+i;
      if(pack->level==level){ // we do something

	// first we compute the adress from the hashfunction
	hidx=hfun((unsigned long long) pack->key,cpu->maxhash);
	nextoct=cpu->htable[hidx];
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,curoct->level)==(unsigned long long)pack->key)&&(pack->level==curoct->level));
	  }while((nextoct!=NULL)&&(!found));

	  if(found){ // the reception oct has been found
	    for(icell=0;icell<8;icell++){
	      memcpy(&(curoct->cell[icell].rfield),&(pack->data[icell]),sizeof(struct Rtype));
	    }
	  }
	  else{
	    printf("error no reception oct found ! for buff #%d lev=%d key=%e\n",i,pack->level,pack->key);
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



// =======================================================================================================

void mpi_exchange_rad_level(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer, struct RAD_MPI **recvbuffer, int cmp_keys,int level)
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


    // ----------- 0.5  / we allocate the number of octs to transmit
    if(cpu->nsend!=NULL){
      free(cpu->nsend);
      free(cpu->nrecv);
      free(cpu->nsend_coarse);
      free(cpu->nrecv_coarse);
    }

    cpu->nsend=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to send to each neighbor
    cpu->nrecv=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to receive from each neighbor
    cpu->nsend_coarse=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to send to each neighbor
    cpu->nrecv_coarse=(int *)calloc(cpu->nnei,sizeof(int)); // the number of packets to receive from each neighbor

    // ----------- I  / we compute the boundary keys and store them in recvbuffer
    compute_bndkeys_rad_level(cpu,level,recvbuffer,cpu->nrecv);

    // ----------- II / we send the number of octs required from each neighbors

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
      MPI_Isend(recvbuffer[i],cpu->nrecv[i],MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
    }
    for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
      MPI_Irecv(sendbuffer[i],cpu->nsend[i],MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
    }
    MPI_Waitall(2*cpu->nnei,req,stat);
    MPI_Barrier(cpu->comm);


  }

  // ----------- III/ the server gather the data
  gather_mpi_rad_level(cpu, sendbuffer,level);

  //if(cpu->rank==RANK_DISP) printf("--- X ---\n");
  memset(req,0,2*cpu->nnei*sizeof(MPI_Request));
  memset(stat,0,2*cpu->nnei*sizeof(MPI_Status));
  MPI_Barrier(cpu->comm);

  // ----------- IV / the server send the data back to the client


  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nsend[i],MPI_PACKET,cpu->mpinei[i],cpu->rank         ,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nrecv[i],MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);



  t[5]=MPI_Wtime();
  // ----------- V  / the client scatter the data back in the oct tree
  scatter_mpi_rad_level(cpu,recvbuffer,level);
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

  //if(cpu->rank==RANK_DISP) printf("--- X ---\n");
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

// ===================================================
// ===================================================
#ifdef WGRAV
void mpi_exchange_pot_level(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int cmp_keys, int level)
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
  gather_mpi_pot_level(cpu, sendbuffer,level);

  //if(cpu->rank==RANK_DISP) printf("--- X ---\n");
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
  scatter_mpi_pot_level(cpu,recvbuffer,level);

  t[6]=MPI_Wtime();

  //
  free(req);
  free(stat);

  t[7]=MPI_Wtime();
  tot=t[7]-t[0];
}


#endif


//------------------------------------------------------------------------
void mpi_cic_correct(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field)
{
  int i,icpu;
  //MPI_Status stat;
  int mpitag=1;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_PACKET);

  MPI_Barrier(cpu->comm);
  //if(cpu->rank==RANK_DISP) printf("correcting CIC on rank %d with %d\n",cpu->rank,cpu->nnei);

  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));


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


void mpi_cic_correct_level(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field, int level)
{
  int i,icpu;
  //MPI_Status stat;
  int mpitag=1;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_PACKET);

  MPI_Barrier(cpu->comm);
  //if(cpu->rank==RANK_DISP) printf("correcting CIC on rank %d with %d\n",cpu->rank,cpu->nnei);

  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));


  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff(cpu,sendbuffer,recvbuffer);

  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  MPI_Barrier(cpu->comm);
  gather_ex_level(cpu,sendbuffer,field,level);
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

void mpi_mark_correct(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer,int level)
{
  int i,icpu;
  //MPI_Status stat;
  int mpitag=1;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_PACKET);

  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));

//  if(cpu->rank==RANK_DISP) printf("correcting hydro CIC on rank %d\n",cpu->rank);

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff(cpu,sendbuffer,recvbuffer);

  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  MPI_Barrier(cpu->comm);
  gather_ex_mark(cpu,sendbuffer,level);
  MPI_Barrier(cpu->comm);

  // --------------------------------------------
  // --------------------------------------------

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(cpu->nsend_coarse+i,1,MPI_INT,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }
  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(cpu->nrecv_coarse+i,1,MPI_INT,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // --------------------------------------------
  // --------------------------------------------

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nsend_coarse[i],MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nrecv_coarse[i],MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // ---------  third we scatter the data back to the INTERNAL boundary octs
  scatter_mpi_mark_ext(cpu,recvbuffer,level);
  MPI_Barrier(cpu->comm);

  //
  free(stat);
  free(req);


}


#ifdef WGRAV
void mpi_dens_correct(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer,int level)
{
  int i,icpu;
  //MPI_Status stat;
  int mpitag=1;
  MPI_Status *stat;
  MPI_Request *req;
  MPI_Datatype MPI_PACKET=*(cpu->MPI_PACKET);

  stat=(MPI_Status*)calloc(cpu->nnei*2,sizeof(MPI_Status));
  req=(MPI_Request*)calloc(cpu->nnei*2,sizeof(MPI_Request));

//  if(cpu->rank==RANK_DISP) printf("correcting hydro CIC on rank %d\n",cpu->rank);

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff(cpu,sendbuffer,recvbuffer);
  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  MPI_Barrier(cpu->comm);

  gather_ex_dens(cpu,sendbuffer,level);
  MPI_Barrier(cpu->comm);

  // --------------------------------------------
  // --------------------------------------------
#if 1

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(cpu->nsend_coarse+i,1,MPI_INT,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }
  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(cpu->nrecv_coarse+i,1,MPI_INT,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // --------------------------------------------
  // --------------------------------------------

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nsend_coarse[i],MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nrecv_coarse[i],MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // ---------  third we scatter the data back to the INTERNAL boundary octs
  scatter_mpi_dens_ext(cpu,recvbuffer,level);
  MPI_Barrier(cpu->comm);
#endif
  //
  free(stat);
  free(req);


}
#endif


#ifdef WHYDRO2
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

  if(cpu->rank==RANK_DISP) printf("correcting hydro CIC on rank %d\n",cpu->rank);

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff_hydro(cpu,sendbuffer,recvbuffer);

  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  MPI_Barrier(cpu->comm);
  gather_ex_hydro(cpu,sendbuffer,level,cpu->nsend_coarse);
  MPI_Barrier(cpu->comm);

  // --------------------------------------------
  // --------------------------------------------

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(cpu->nsend_coarse+i,1,MPI_INT,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }
  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(cpu->nrecv_coarse+i,1,MPI_INT,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // --------------------------------------------
  // --------------------------------------------

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nsend_coarse[i],MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nrecv_coarse[i],MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
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
#endif

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

  //if(cpu->rank==RANK_DISP) printf("correcting Rad CIC on rank %d\n",cpu->rank);

  // ----------- 0  / we clean the mpi buffers
  clean_mpibuff_rad(cpu,sendbuffer,recvbuffer);

  // ---------  first we collect the data from EXTERNAL boundaries (keys are computed by remote)
  MPI_Barrier(cpu->comm);
  gather_ex_rad(cpu,sendbuffer,level,cpu->nsend_coarse);
  MPI_Barrier(cpu->comm);


  // -----------

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(cpu->nsend_coarse+i,1,MPI_INT,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }
  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(cpu->nrecv_coarse+i,1,MPI_INT,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
  }
  MPI_Waitall(2*cpu->nnei,req,stat);
  MPI_Barrier(cpu->comm);

  // ----------------------- second we send the data

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Isend(sendbuffer[i],cpu->nsend_coarse[i],MPI_PACKET,cpu->mpinei[i],cpu->rank,cpu->comm,&req[i]  );
  }

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors to send the keys
    MPI_Irecv(recvbuffer[i],cpu->nrecv_coarse[i],MPI_PACKET,cpu->mpinei[i],cpu->mpinei[i],cpu->comm,&req[i+cpu->nnei]);
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
void mpi_exchange_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer, int *delta,int level){

  //int* delta = (int*)calloc(2,sizeof(int));
  int val;
  int mpitag=1;
  int i;
  static int nrem[2];
  static int nadd[2];
  MPI_Datatype MPI_PACKET=*(cpu->MPI_PART);

  delta[0]=0;
  delta[1]=0;

  clean_mpibuff_part(cpu,psendbuffer,precvbuffer);
  gather_ex_part(cpu,psendbuffer,nrem);
  MPI_Barrier(cpu->comm);

  for(i=0;i<cpu->nnei;i++){ // we scan all the neighbors
    mpitag=cpu->rank+cpu->mpinei[i];

    MPI_Sendrecv(psendbuffer[i],cpu->nbuffpart,*cpu->MPI_PART,cpu->mpinei[i],mpitag,precvbuffer[i],cpu->nbuffpart,*cpu->MPI_PART,cpu->mpinei[i],mpitag,cpu->comm,MPI_STATUS_IGNORE);
  }

  MPI_Barrier(cpu->comm);
  scatter_mpi_part(cpu,precvbuffer,nadd,level);
  MPI_Barrier(cpu->comm);

  //printf("nadd=%d nrem=%d on rank %d\n",nadd[0],nrem[0],cpu->rank);
  delta[0] = nadd[0]-nrem[0];
  delta[1] = nadd[1]-nrem[1];

  //return delta;		// Return delta part
}
#endif // PIC

#endif // WMPI


void init_MPI(struct CPUINFO *cpu,
              MPI_Datatype *MPI_PACKET,
              MPI_Datatype *MPI_PART,
              MPI_Datatype *MPI_WTYPE,
              MPI_Datatype *MPI_HYDRO,
              MPI_Datatype *MPI_RTYPE,
              MPI_Datatype *MPI_RAD ){

  MPI_Status stat;

  MPI_Comm_size(MPI_COMM_WORLD,&(cpu->nproc));
  MPI_Comm_rank(MPI_COMM_WORLD,&(cpu->rank));

#ifdef WDBG
  //if(cpu.rank==4) gdb_debug();
  MPI_Barrier(MPI_COMM_WORLD);
#endif

 //========= creating a PACKET MPI type =======
  MPI_Datatype oldtypes[16];
  int          blockcounts[16];
  struct PACKET _info_pack;


  /* MPI_Aint type used to be consistent with syntax of */
  /* MPI_Type_extent routine */
  MPI_Aint    offsets[5], extent;
  MPI_Aint base;


  MPI_Address(&_info_pack.data,offsets);
  base=offsets[0];
  MPI_Address(&_info_pack.key,offsets+1);
  MPI_Address(&_info_pack.level,offsets+2);

  offsets[0]=offsets[0]-base;
  offsets[1]=offsets[1]-base;
  offsets[2]=offsets[2]-base;

  oldtypes[0]=MPI_REEL;
  oldtypes[1]=MPI_DOUBLE;
  oldtypes[2]=MPI_INT;

  blockcounts[0]=8;
  blockcounts[1]=1;
  blockcounts[2]=1;

  MPI_Type_struct(3, blockcounts, offsets, oldtypes, MPI_PACKET); // MODKEY WARNING TO 2 AND 3
  MPI_Type_commit(MPI_PACKET);


#ifdef PIC
  //========= creating a PART MPI type =======
  struct PART_MPI _info;


  MPI_Address(&_info.x,offsets);
  base=offsets[0];
  MPI_Address(&_info.key,offsets+1);
  MPI_Address(&_info.idx,offsets+2);

  offsets[0]=offsets[0]-base;
  offsets[1]=offsets[1]-base;
  offsets[2]=offsets[2]-base;

  oldtypes[0]=MPI_REEL;
  oldtypes[1]=MPI_DOUBLE;
  oldtypes[2]=MPI_INT;

#ifdef STARS
  blockcounts[0]=9;
#else
  blockcounts[0]=7;
#endif
  blockcounts[1]=1;

#ifdef STARS
  blockcounts[2]=6;
#else
  blockcounts[2]=4;
#endif

  MPI_Type_struct(3, blockcounts, offsets, oldtypes, MPI_PART); // MODKEY WARNING TO 2 AND 3
  MPI_Type_commit(MPI_PART);
#endif // PIC

#ifdef WHYDRO2
  //========= creating a WTYPE MPI type =======

  struct Wtype _info_hyd;

  MPI_Address(&_info_hyd.d,offsets);
  base=offsets[0];

  offsets[0]=offsets[0]-base;

  oldtypes[0]=MPI_REEL;

#ifdef WRADHYD
#ifdef HELIUM
  blockcounts[0]=10;
#else
  blockcounts[0]=8;
#endif // HELIUM
#else
  blockcounts[0]=7;
#endif // WRADHYD

  /* Now define structured type and commit it */
  MPI_Type_struct(1, blockcounts, offsets, oldtypes, MPI_WTYPE);
  MPI_Type_commit(MPI_WTYPE);


  //========= creating a HYDRO MPI type =======
  struct HYDRO_MPI _info_hydmpi;

  MPI_Address(&_info_hydmpi.data[0],offsets);
  base=offsets[0];
  MPI_Address(&_info_hydmpi.key,offsets+1);
  MPI_Address(&_info_hydmpi.level,offsets+2);

  offsets[0]=offsets[0]-base;
  offsets[1]=offsets[1]-base;
  offsets[2]=offsets[2]-base;

  oldtypes[0]=*MPI_WTYPE;
  oldtypes[1]=MPI_DOUBLE;
  oldtypes[2]=MPI_INT;

  blockcounts[0]=8;
  blockcounts[1]=1;
  blockcounts[2]=1;

/* Now define structured type and commit it */
  MPI_Type_struct(3, blockcounts, offsets, oldtypes, MPI_HYDRO);
  MPI_Type_commit(MPI_HYDRO);
#endif // WHYDRO2

#ifdef WRAD
  //========= creating a RTYPE MPI type =======

  struct Rtype _info_r;


  MPI_Address(&_info_r.e[0],offsets);
  base=offsets[0];
  MPI_Address(&_info_r.fx[0],offsets+1);
  MPI_Address(&_info_r.fy[0],offsets+2);
  MPI_Address(&_info_r.fz[0],offsets+3);
  MPI_Address(&_info_r.src,offsets+4);

  offsets[0]=offsets[0]-base;
  offsets[1]=offsets[1]-base;
  offsets[2]=offsets[2]-base;
  offsets[3]=offsets[3]-base;
  offsets[4]=offsets[4]-base;

  oldtypes[0]=MPI_REEL;
  oldtypes[1]=MPI_REEL;
  oldtypes[2]=MPI_REEL;
  oldtypes[3]=MPI_REEL;
  oldtypes[4]=MPI_REEL;


  blockcounts[0]=NGRP;
  blockcounts[1]=NGRP;
  blockcounts[2]=NGRP;
  blockcounts[3]=NGRP;
#ifdef WCHEM
  blockcounts[4]=5;
#else
  blockcounts[4]=1;
#endif

#ifdef STARS
  blockcounts[4]++;
#endif

#ifdef HELIUM
  blockcounts[4]+=2;
#endif

  /* Now define structured type and commit it */
  MPI_Type_struct(5, blockcounts, offsets, oldtypes, MPI_RTYPE);
  MPI_Type_commit(MPI_RTYPE);

  //========= creating a RAD MPI type =======

  struct RAD_MPI _info_radmpi;

  MPI_Address(&_info_radmpi.data[0],offsets);
  base=offsets[0];
  MPI_Address(&_info_radmpi.key,offsets+1);
  MPI_Address(&_info_radmpi.level,offsets+2);

  offsets[0]=offsets[0]-base;
  offsets[1]=offsets[1]-base;
  offsets[2]=offsets[2]-base;

  oldtypes[0]=*MPI_RTYPE;
  oldtypes[1]=MPI_DOUBLE;
  oldtypes[2]=MPI_INT;

  blockcounts[0]=8;
  blockcounts[1]=1;
  blockcounts[2]=1;

  /* Now define structured type and commit it */
  MPI_Type_struct(3, blockcounts, offsets, oldtypes, MPI_RAD);
  MPI_Type_commit(MPI_RAD);
#endif // WRAD

  cpu->MPI_PACKET=MPI_PACKET;
#ifdef PIC
  cpu->MPI_PART=MPI_PART;
#endif

#ifdef WHYDRO2
  cpu->MPI_HYDRO=MPI_HYDRO;
#endif

#ifdef WRAD
  cpu->MPI_RAD=MPI_RAD;
#endif

  cpu->comm=MPI_COMM_WORLD;


#if defined(MPIIO) || defined(HDF5)
  cpu->mpiio_ncells = (int*)calloc(cpu->nproc,sizeof(int));
  cpu->mpiio_nparts = (int*)calloc(cpu->nproc,sizeof(int));
#ifdef STARS
  cpu->mpiio_nstars = (int*)calloc(cpu->nproc,sizeof(int));
#endif // STARS
#endif // MPIIO


}
#endif // WMPI
