#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hilbert.h"
#include "prototypes.h"
#include "vector.h"
#include "io.h"
#include "cic.h"
#include "oct.h"
#include "particle.h"
#include "amr.h"
#include "tools.h"
#include "segment.h"
#include "communication.h"
#include <time.h>

// TO COMPILE

//mpicc -lm -DWMPI -DNBUFF=4096 -DNEWASSIGN -DTESTPLUM -DNDUMP=1 -DNSTEP=10 -DLCOARSE=6 -DLMAX=6  -DCIC2 -DNITER=4096 -DALLCELL -DSTRIDE=128 -DDT=1e-4 hilbert.c gilgamesh.c

//gcc -g -lm -DNEWASSIGN -DPART2 -DNDUMP=1 -DNSTEP=10 -DLCOARSE=6 -DLMAX=6  -DCIC2 -DNITER=4096 -DALLCELL -DSTRIDE=128 -DDT=1e-4 hilbert.o gilgamesh.c

//mpicc -lm -DWMPI -DNEWASSIGN -DPART2 -DNDUMP=1 -DNSTEP=10 -DLCOARSE=6 -DLMAX=6  -DCIC2 -DNITER=4096 -DALLCELL -DSTRIDE=128 -DDT=1e-4 hilbert.c gilgamesh.c







//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//==================================================================

//------------------------------------------------------------------------
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------




//------------------------------------------------------------------------

 //------------------------------------------------------------------------
 //------------------------------------------------------------------------
 //------------------------------------------------------------------------


 //------------------------------------------------------------------------
 // the MAIN CODE
 //------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  struct OCT *grid;
  struct OCT **firstoct;
  struct OCT **lastoct;

  int level,levelcoarse,levelmax;
  int ngridmax,ngrid;
  int npartmax;
  int cur,curnext; // flat indexes with boundaries
  int i,il,ichild,icell,inext,ii,ip,j;
  int xp,yp,zp;
  int NBND=1,NBND2=2*NBND;
  float dx;
  int vnei[6],vcell[6]; // arrays to get neighbors
  int vnei2[6],vcell2[6]; // arrays to get neighbors
  int vnei3[6],vcell3[6]; // arrays to get neighbors
  int neip[7]; // contains the index of the six neighbors of the current parent cell +current parent
  int ci,cj,ck;
  int cinext,cjnext,cknext;
  float threshold;
  struct OCT oct;
  struct OCT* nextoct;
  struct OCT* curoct;
  struct OCT* desoct;
  struct CELL * parcell;
  struct CELL * newcell;
  struct CELL * newcell2;
  int tag;
  float dxcur;
  float *dens;
  int firstoct_currl;
  int nxoct;
  int lmap;
  int npart;
  struct PART *part;
  struct PART *nexploc, *curploc;

  struct OCT *endoct; // the very last oct of all the levels;
  struct OCT *newoct;
  int nref=0,ndes=0;

  float xc,yc,zc;
  int stride;
  float **vcomp;
  int ncomp;
  float acc;
  float dt;
  int ntot=0,nlev,noct;
  float ntotd=0.,nlevd=0.;

  float disp,mdisp;
  
  int dir;

  char filename[128]; 
  FILE *fd;
  struct PART *nexp;
  struct PART *nexp2;
  struct PART *curp;

  struct PART *lastpart;

  int curc;

  int nbnd;

  float x,y,z;
  float vx,vy,vz;
  float mass;
  float idx;
  
  unsigned key;

  struct CPUINFO cpu;

  struct PACKET **sendbuffer; 
  struct PACKET **recvbuffer; 

  struct PART_MPI **psendbuffer; 
  struct PART_MPI **precvbuffer; 

  //=========== some initial calls =============
#ifdef WMPI
  MPI_Status stat;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(cpu.nproc));
  MPI_Comm_rank(MPI_COMM_WORLD,&(cpu.rank));

  //========= creating a PACKET MPI type =======
  MPI_Datatype MPI_PACKET,oldtypes[2]; 
  int          blockcounts[2];
  
  /* MPI_Aint type used to be consistent with syntax of */
  /* MPI_Type_extent routine */
  MPI_Aint    offsets[2], extent;
  
  
  /* Setup description of the 8 MPI_FLOAT fields data */
  offsets[0] = 0;
  oldtypes[0] = MPI_FLOAT;
  blockcounts[0] = 8;
  
  /* Setup description of the 2 MPI_INT fields key, level */
  /* Need to first figure offset by getting size of MPI_FLOAT */
  MPI_Type_extent(MPI_FLOAT, &extent);
  offsets[1] = 8 * extent;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 2;

  /* Now define structured type and commit it */
  MPI_Type_struct(2, blockcounts, offsets, oldtypes, &MPI_PACKET);
  MPI_Type_commit(&MPI_PACKET);

  //========= creating a PART MPI type =======
  MPI_Datatype MPI_PART;

  /* Setup description of the 7 MPI_FLOAT fields x,y,z,vx,vy,vz */
  offsets[0] = 0;
  oldtypes[0] = MPI_FLOAT;
  blockcounts[0] = 7;
  
  /* Setup description of the 4 MPI_INT fields idx key level icell*/
  /* Need to first figure offset by getting size of MPI_FLOAT */
  MPI_Type_extent(MPI_FLOAT, &extent);
  offsets[1] = 7 * extent;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 4;

  /* Now define structured type and commit it */
  MPI_Type_struct(2, blockcounts, offsets, oldtypes, &MPI_PART);
  MPI_Type_commit(&MPI_PART);
  
  //============================================


  cpu.MPI_PACKET=&MPI_PACKET;
  cpu.MPI_PART=&MPI_PART;
  cpu.comm=MPI_COMM_WORLD;
#else
  cpu.rank=0;
  cpu.nproc=1;
#endif

  //=========== assigning values =============
  levelcoarse=LCOARSE;
  levelmax=LMAX;

  ngridmax=1000000;
  npartmax=64*64*64*2;
#ifdef PART2
  npart=2;
#else
  npart=64*64*64;
#endif

  threshold=50;
  lmap=LMAX;
  stride=fmax(8,STRIDE);//pow(2,levelcoarse);
  ncomp=8;
  acc=1e-2;
  dt=DT;

  //breakmpi();
  //========== allocations ===================

  if(cpu.rank==0) printf("Allocating %f GB cell=%f GB part=%f GB",(sizeof(struct OCT)*ngridmax+sizeof(struct PART)*npart)/(1024*1024*1024.),sizeof(struct OCT)*ngridmax/(1024*1024*1024.),sizeof(struct PART)*npart/(1024*1024*1024.));

  grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT));
  firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
  lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *));
  part=(struct PART*)calloc(npartmax,sizeof(struct PART));
  cpu.htable=(struct OCT**) calloc(pow(2,3*levelmax-3)/64,sizeof(struct OCT *)); //hashtable assunming a hash function >>6


  vcomp=(float **)calloc(ncomp,sizeof(float*));
  for(i=0;i<ncomp;i++)
    {
      vcomp[i]=(float *)calloc(stride,sizeof(float));
    }
    if(cpu.rank==0) printf("Allocations ok\n");


  //========== setting up the parallel topology ===

  // We segment the oct distributions at levelcoarse 
    cpu.bndoct=NULL;
    cpu.nbuff=NBUFF;
    cpu.allkmin=(int*)calloc(cpu.nproc,sizeof(int));
    cpu.allkmax=(int*)calloc(cpu.nproc,sizeof(int));

    load_balance(levelcoarse,&cpu);

#ifdef WMPI
    MPI_Allgather(&cpu.kmin,1,MPI_INT,cpu.allkmin,1,MPI_INT,cpu.comm);
    MPI_Allgather(&cpu.kmax,1,MPI_INT,cpu.allkmax,1,MPI_INT,cpu.comm);
    MPI_Barrier(cpu.comm);
#else
    cpu.allkmin[0]=cpu.kmin;
    cpu.allkmax[0]=cpu.kmax;
#endif    
    
    

    
  //========== building the initial meshes ===

  if(cpu.rank==0) printf("building initial mesh\n");

  //breakmpi();
  // ZERO WE CREATE A ROOT CELL
  
  struct CELL root;
  root.child=grid;
  

  // FIRST WE POPULATE THE ROOT OCT
  grid->x=0.;
  grid->y=0.;
  grid->z=0.;

  grid->parent=NULL;
  grid->level=1;
  for(i=0;i<6;i++) grid->nei[i]=&root;
  grid->prev=NULL;
  grid->next=NULL;

  // setting the densities in the cells and the index
  for(icell=0;icell<8;icell++){ 
    grid->cell[icell].density=0.;
    grid->cell[icell].pot=0.;
    grid->cell[icell].temp=0.;
    grid->cell[icell].idx=icell;
  }

  grid->cpu=-1;

  // start the creation of the initial amr grid from level 1
  firstoct[0]=grid;
  lastoct[0]=grid;
  int noct2;
  int segok;

  newoct=grid+1;
  for(level=1;level<levelcoarse;level++){ // sweeping the levels from l=1 to l=levelcoarse
    dxcur=1./pow(2,level);
    nextoct=firstoct[level-1];
    noct2=0;
    if(nextoct==NULL) continue;
    do // sweeping level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++){ // sweeping the cells

	  segok=segment_cell(curoct,icell,&cpu,levelcoarse);// the current cell will be splitted according to a segmentation condition
	  if(segok==1){ 
	    //if(level==levelcoarse-1) printf(" segok=%d\n",segok);

	    noct2++;
	    // the newoct is connected to its mother cell
	    curoct->cell[icell].child=newoct;
	    
	    // a newoct is created
	    newoct->parent=&(curoct->cell[icell]);
	    newoct->level=curoct->level+1;
	    newoct->x=curoct->x+( icell   %2)*dxcur;
	    newoct->y=curoct->y+((icell/2)%2)*dxcur;
	    newoct->z=curoct->z+( icell   /4)*dxcur;

	    // filling the cells
	    for(ii=0;ii<8;ii++){
	      newoct->cell[ii].marked=0;
	      newoct->cell[ii].child=NULL;
	      newoct->cell[ii].density=0.;
	      newoct->cell[ii].idx=ii;
	      newoct->cell[ii].phead=NULL;
	    }
	    
	    //the neighbours
	    getcellnei(icell, vnei, vcell);
	    for(ii=0;ii<6;ii++){
	      if((vnei[ii]!=6)){ 
		newoct->nei[ii]=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
	      }else{
		newoct->nei[ii]=&(curoct->cell[vcell[ii]]);
	      }
	    }
	    
	    // preparing the next creations on level+1
	    newoct->next=NULL;
	    
	    if(firstoct[level]==NULL){
	      firstoct[level]=newoct;
	      newoct->prev=NULL;
	    }
	    else{
	      newoct->prev=lastoct[level];
	      lastoct[level]->next=newoct;
	    }
	    lastoct[level]=newoct;


	    // next oct ready
	    newoct++; 
	  }
 	}
      }while(nextoct!=NULL);
    if(cpu.rank==0) printf("level=%d noct=%d\n",level,noct2);
  }


 // ==================================== assigning CPU number to levelcoarse OCTS // filling the hash table // Setting up the MPI COMMS

  setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax);

  // allocating the communication buffers

  sendbuffer=(struct PACKET **)(calloc(cpu.nnei,sizeof(struct PACKET*)));
  recvbuffer=(struct PACKET **)(calloc(cpu.nnei,sizeof(struct PACKET*)));
  for(i=0;i<cpu.nnei;i++) {
    sendbuffer[i]=(struct PACKET *) (calloc(cpu.nbuff,sizeof(struct PACKET)));
    recvbuffer[i]=(struct PACKET *) (calloc(cpu.nbuff,sizeof(struct PACKET)));
  }

  psendbuffer=(struct PART_MPI **)(calloc(cpu.nnei,sizeof(struct PART_MPI*)));
  precvbuffer=(struct PART_MPI **)(calloc(cpu.nnei,sizeof(struct PART_MPI*)));
  for(i=0;i<cpu.nnei;i++) {
    psendbuffer[i]=(struct PART_MPI *) (calloc(cpu.nbuff,sizeof(struct PART_MPI)));
    precvbuffer[i]=(struct PART_MPI *) (calloc(cpu.nbuff,sizeof(struct PART_MPI)));
  }


  //===================================================================================================
  
  // ==== some initial dump

  sprintf(filename,"data/levstart.%05d.p%05d",0,cpu.rank);
  dumpcube(lmap,firstoct,0,filename);
  sprintf(filename,"data/cpustart.%05d.p%05d",0,cpu.rank);
  dumpcube(lmap,firstoct,3,filename);

  // =====  computing the memory location of the last oct 

  endoct=lastoct[0];
  for(i=0;i<levelcoarse;i++) {
    if(lastoct[i]>endoct) endoct=lastoct[i];
  }


#if 1  // ==================================== assigning particles to cells

  if(cpu.rank==0) printf("==> starting part\n");
  firstoct_currl=0;
  for(il=1;il<levelcoarse;il++) firstoct_currl+=pow(pow(2,il-1),3); // the index of the first oct of current level
 
  // initialisation of particles
  

#ifdef PART2

  int ir,nr=2;
  ip=0;
  for(ir=0;ir<nr;ir++) {
    // first we read the position etc... (eventually from the file)
    if(ir==0){
      x=0.5;
      y=0.5;
      z=0.5;
      
      vx=0.;
      vy=0.;
      vz=0.;
      
      mass=0.999;
    }
    else if(ir==1){
      x=0.5+0.2;
      y=0.5;
      z=0.5;
      
      vx=0.;
      vy=sqrt(0.999/0.2);
      vz=0.;
      
      mass=0.001;
    }
    
    // periodic boundary conditions
    
    x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
    y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
    z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 
    
    // it it belongs to the current cpu, we proceed and assign the particle to the particle array
    if(segment_part(x,y,z,&cpu,levelcoarse)){
      part[ip].x=x;
      part[ip].y=y;
      part[ip].z=z;
      
      part[ip].vx=vx;
      part[ip].vy=vy;
      part[ip].vz=vz;
      
      part[ip].mass=mass;
      lastpart=part+ip;
      ip++;
    }
  }
  
  npart=ip; // we compute the localnumber of particle

#endif

#ifdef TESTPLUM
  int dummy;
  float dummyf;
  int npartf;

  //breakmpi();
  fd=fopen("utils/data.inp","r");
  if(fd==NULL) {
    printf("Error while reading particle file ABORT\n");
    abort();
  }
  fscanf(fd,"%d",&dummy);
  fscanf(fd,"%d",&npartf);
  fscanf(fd,"%f",&dummyf);

  ip=0.;
  for(i=0;i<npartf;i++)
    {
      //fscanf(fd,"%d %f %f %f %f %f %f %f",&part[i].idx,&part[i].mass,&(part[i].x),&(part[i].y),&(part[i].z),&(part[i].vx),&(part[i].vy),&(part[i].vz));
      fscanf(fd,"%d %f %f %f %f %f %f %f",&dummy,&mass,&x,&y,&z,&vx,&vy,&vz);
      
      x+=0.5;
      y+=0.5;
      z+=0.5;
      // periodic boundary conditions
    
      x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
      y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
      z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 

      // it it belongs to the current cpu, we proceed and assign the particle to the particle array
      if(segment_part(x,y,z,&cpu,levelcoarse)){
	part[ip].x=x;
	part[ip].y=y;
	part[ip].z=z;
	
	part[ip].vx=vx;
	part[ip].vy=vy;
	part[ip].vz=vz;
	
	part[ip].mass=mass;
	lastpart=part+ip;
	ip++;
      }
      
    }
  fclose(fd);
  npart=ip; // we compute the localnumber of particle

#endif  

  // assigning particles to cells in coarse octs (assuming octs are aligned)

  if(cpu.rank==0) printf("start populating coarse grid with particles\n");
  struct PART* lastp[8]; // will contain the last particle of the 8 cells in each oct

  // FIRST WE CONSIDER THE LEVEL 1
  for(ii=0;ii<8;ii++) lastp[ii]=NULL; // we initialise the last part of each sub cell
  dxcur=0.5;
  for(i=0;i<npart;i++)
    {
      curc=(int)((part[i].x-grid[0].x)/dxcur)+(int)((part[i].y-grid[0].y)/dxcur)*2+(int)((part[i].z-grid[0].z)/dxcur)*4;
      
      if(grid[0].cell[curc].phead==NULL){
	grid[0].cell[curc].phead=&part[i];
	lastp[curc]=&part[i];
      }
      else{
	lastp[curc]->next=&part[i];
	part[i].prev=lastp[curc];
	lastp[curc]=&part[i];
      }
    }
  if(cpu.rank==0) printf("Part assigned root level ok\n");
  
  // WE POPULATE THE NEXT LEVELS BY SUBDIVISIONS
  int np=0;
  for(level=1;level<=levelcoarse-1;level++) // we stop at level coarse -1 because it will be assigned from levelcoarse-1
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  dxcur=1./pow(2,level+1); // size of a CELL at level +1
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if(curoct->cell[icell].child!=NULL){ // a child has been detected so we split the particle in 8 cells
		for(ii=0;ii<8;ii++) lastp[ii]=NULL; // we initialise the last part of each sub cell
		newoct=curoct->cell[icell].child; // temp oct for practical reasons
		nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
		if(nexp!=NULL){ 
		  do{  
		    curp=nexp; 
		    nexp=curp->next; 
		    
		    //curc is the index of the new cell at level+1
		    curc=(int)((curp->x-newoct->x)/dxcur)+(int)((curp->y-newoct->y)/dxcur)*2+(int)((curp->z-newoct->z)/dxcur)*4;
		    if(lastp[curc]==NULL){
		      // first particle in the current subcell
		      newoct->cell[curc].phead=curp;
		      curp->prev=NULL;
		      curp->next=NULL;
		      lastp[curc]=curp;
		    }
		    else{
		      // the current particle is linked to the last one in the current part
		      curp->prev=lastp[curc];
		      lastp[curc]->next=curp;
		      curp->next=NULL;
		      lastp[curc]=curp;
		    }
		  }while(nexp!=NULL); 
		  
		  // we empty the mother cell from particles
		  curoct->cell[icell].phead=NULL;
		  
		}
	      }
	    }
	}while(nextoct!=NULL);
    }


#endif
#if 1
  // ==================================== Check the number of particles and octs

  multicheck(firstoct,npart,levelmax,cpu.rank);

  sprintf(filename,"data/parstart.%05d.p%05d",0,cpu.rank);
  dumppart(firstoct,filename,npart,levelcoarse,levelmax);

#endif	

#if 1
  // ==================================== performing the CIC assignement

  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#endif

#ifdef WMPI
    // ==================================== performing the CIC BOUNDARY CORRECTION 

  mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);

  // ======================================= Density boundary mpi update 

  mpi_exchange(&cpu,sendbuffer,recvbuffer,1);

#endif



  sprintf(filename,"data/denstart.%05d.p%05d",0,cpu.rank);
  printf("%s\n",filename);
  dumpcube(lmap,firstoct,1,filename);





  //================================================================================
  //================================================================================
  //================================================================================
  //
  //          AT THIS STAGE THE INITIAL SETUP HAS BEEN COMPLETED
  //
  //================================================================================
  //================================================================================
  //================================================================================


  int nsteps;
  int pass;
  int smark;
  int ismooth,nsmooth=2;
  int marker;

  for(nsteps=0;nsteps<=NSTEP;nsteps++){
    
    if(cpu.rank==0) printf("============== STEP %d ================\n",nsteps);
    //printf("endoct=%p\n",endoct);
#if 1
    // ==================================== marking the cells
    //if(nsteps==1)  breakmpi();
    
    sprintf(filename,"data/dmstart.%05d.p%05d",nsteps+1,cpu.rank);
    printf("%s\n",filename);
    dumpcube(lmap,firstoct,1,filename);

    mark_cells(levelcoarse,levelmax,firstoct,nsmooth,threshold,&cpu,sendbuffer,recvbuffer);

    sprintf(filename,"data/markstart.%05d.p%05d",nsteps+1,cpu.rank);
    printf("%s\n",filename);
    dumpcube(lmap,firstoct,4,filename);
    
    
  
    
    // ==================================== refining (and destroying) the octs

    
    refine_cells(levelcoarse,levelmax,firstoct,lastoct,endoct,&cpu);
    //if(nsteps==1) breakmpi();
    // cleaning the marks
    for(level=1;level<=levelmax;level++) // looping over levels
      {
	float maxd=0.,mind=1e30,avg=0.;
	int ncell=0;
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
/* 	      if(curoct->cell[icell].marked>50){ */
/* 		printf("ouhl %f\n",curoct->cell[icell].marked); */
/* 		//abort(); */
/* 	      } */
	      curoct->cell[icell].marked=0.;
	      ncell++;
	      avg+=curoct->cell[icell].density;
	      if(curoct->cell[icell].density>maxd) maxd=curoct->cell[icell].density;
	      if(curoct->cell[icell].density<mind) mind=curoct->cell[icell].density;
	    }
	}while(nextoct!=NULL);

      printf("level=%d avg=%f mind=%f maxd=%f\n",level,avg/ncell,mind,maxd);
      }
    
  // recomputing the last oct;
  
  printf("==> memory state\n");
  printf("endoct=%p\n",endoct);
  
  endoct=lastoct[0];
  for(i=0;i<levelmax;i++) {
    if(lastoct[i]>endoct) endoct=lastoct[i];
    //printf("i=%d %p %p\n",i+1,firstoct[i],lastoct[i]);
  }
  //printf("endoct=%p\n",endoct);

#endif

  sprintf(filename,"data/levstart.%05d.p%05d",nsteps+1,cpu.rank);
  dumpcube(lmap,firstoct,0,filename);



#ifdef WMPI
  // ==================================== after refinement we should remap the boundary cells
  setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax);
#endif

  sprintf(filename,"data/cpustart.%05d.p%05d",nsteps+1,cpu.rank);
  dumpcube(lmap,firstoct,3,filename);


#if 1
  // ==================================== performing the CIC assignement
  //breakmpi();

  call_cic(levelmax,levelcoarse,firstoct,&cpu);

#ifdef WMPI
    // ==================================== performing the CIC BOUNDARY CORRECTION 

  mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);

  // ======================================= Density boundary mpi update 

  mpi_exchange(&cpu,sendbuffer,recvbuffer,1);
#endif

#endif

  sprintf(filename,"data/denstart.%05d.p%05d",nsteps+1,cpu.rank);
  printf("%s\n",filename);
  dumpcube(lmap,firstoct,1,filename);


#if 1
  // ==================================== Check the number of particles and octs
  multicheck(firstoct,npart,levelmax,cpu.rank);
#endif	 

#if 1
// ==================================== POISSON Testing the jacobi iteration

#ifdef TIME_JAC
  FILE *ft;
  sprintf(filename,"data/timejac.%05d.p%05d",nsteps+1,cpu.rank);
  ft=fopen(filename,"w");
  
  double tg1,tg2;
  double tl1,tl2;
  double ta1,ta2;
  double ts1,ts2;
  double tm1,tm2;
  double t1,t2;

#endif

  

  printf("==> Poisson Start \n");
  int icomp,iter,niter=NITER;
  float norm_d;
  int nread;
  for(level=levelcoarse;level<=levelmax;level++)
    {
      // COARSE LEVEL TREATMENT ============================================
      if(level==levelcoarse){
	norm_d=0.;
	//fixbound_reg(grid,levelcoarse,NBND);
	for(iter=0;iter<niter;iter++){

	  if((iter%64==0)&&(cpu.rank==0)) printf("iter=%d ",iter);

#ifdef WMPI
	  t1=MPI_Wtime();
#endif 
	  nextoct=firstoct[level-1];
	  if(nextoct!=NULL){
	    dx=pow(0.5,level);
	    do{ 
	      curoct=nextoct;
	    
#ifdef WMPI
	  tg1=MPI_Wtime();
#endif 
	      // First we gather the potential in all neighbors
	      for(icomp=0;icomp<=6;icomp++){
		memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
		nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
	      }
	      
	      // Second we gather the local density
	      memset(vcomp[7],0,stride*sizeof(float)); // reset the vcomp;
	      nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);
#ifdef WMPI
	  tg2=MPI_Wtime();
#endif 

#ifdef WMPI
	  ta1=MPI_Wtime();
#endif 

	      // 2.5 we contrast the density by removing the average density value
	      remove_avg(vcomp[7],stride,1.);

	      // we compute the square of the norm of the density (first iteration only)
	      if(iter==0) norm_d+=square(vcomp[7],nread);
#ifdef WMPI
	  ta2=MPI_Wtime();
#endif 
	  
#ifdef WMPI
	  tl1=MPI_Wtime();
#endif 

	      // Third we perform the calculation (eventually on GPU)
	      laplacian(vcomp,stride,dx);

#ifdef WMPI
	  tl2=MPI_Wtime();
#endif 


#ifdef WMPI
	  ts1=MPI_Wtime();
#endif 
	      // Fourth we scatter back the potential estimation to the temp position
	      nextoct=scattercomp(curoct, vcomp[6], 6, 2, stride,&cpu);
#ifdef WMPI
	  ts2=MPI_Wtime();
#endif 

	    }while(nextoct!=NULL);
	  }
	  
	  // 4.5 we copy the result in the temp position to the potential
	  nextoct=firstoct[level-1];
	  if(nextoct!=NULL){
	    do{ 
	      curoct=nextoct;
	    
	      memset(vcomp[0],0,stride*sizeof(float)); // reset the vcomp;

	      nextoct=gathercomp(curoct, vcomp[0], 6, 2, stride,&cpu,&nread); // getting the data in the temp field
	      nextoct=scattercomp(curoct, vcomp[0], 6, 1, stride,&cpu);

	    }while(nextoct!=NULL);
	  }

#ifdef WMPI
	  tm1=MPI_Wtime();
	  mpi_exchange(&cpu,sendbuffer,recvbuffer,2);
	  if(iter==0) MPI_Allreduce(MPI_IN_PLACE,&norm_d,1,MPI_FLOAT,MPI_SUM,cpu.comm);
	  tm2=MPI_Wtime();

#endif

#ifdef WMPI
	  t2=MPI_Wtime();
#endif 

	    // Fifth we compute the residuals

	    nextoct=firstoct[level-1];
	    float res=0.;
	    if(nextoct!=NULL){
	      dx=pow(0.5,level);
	      do{ 
		curoct=nextoct;
	      
		// First we gather the potential in all neighbors + local
		for(icomp=0;icomp<=6;icomp++){
		  memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
		  nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
		}

		// Second we gather the local density
		memset(vcomp[7],0,stride*sizeof(float)); // reset the vcomp;
		nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);

		// 2.5 we contrast the density by removing the average value
		remove_avg(vcomp[7],stride,1.);
	      
		// Third we perform the square of the residual
		res+=square_res(vcomp,nread,dx);

	      }while(nextoct!=NULL);
	    }
	    
#ifdef WMPI
	    // reducing the residuals
	    MPI_Barrier(cpu.comm);
	    float restot;
	    //if(iter%64==0) printf("res = %f on rank=%d\n",res,cpu.rank);
	    MPI_Allreduce(&res,&restot,1,MPI_FLOAT,MPI_SUM,cpu.comm);
	    res=restot;
#endif
	    if((iter%64==0)&&(cpu.rank==0)) printf("dens=%e res=%e relative residual=%e\n ",sqrt(norm_d),sqrt(res),sqrt(res/norm_d));
	    if(sqrt(res/norm_d)<acc) break;
#ifdef WMPI
	    fprintf(ft,"%d %e %e %e %e %e\n",level,tg2-tg1,ta2-ta1,tl2-tl1,ts2-ts1,tm2-tm1);
#endif

	}
	
      }
      // FINE LEVEL TREATMENT ============================================
      else{
	norm_d=0.;
	
	//initial guess from parent cell
	nextoct=firstoct[level-1];
	//printf("first=%p next=%p\n",nextoct,nextoct->next);
	if(nextoct==NULL){
	  continue; // we skip to next level if the firstoct is empty
	}
	else{
	  do{ 
	    curoct=nextoct;
	    nextoct=curoct->next;
	    if(curoct->cpu!=cpu.rank) continue;
	    for(icell=0;icell<8;icell++){
	      curoct->cell[icell].pot=curoct->parent->pot;
	    }
	  }while(nextoct!=NULL);
	}
#ifdef WMPI
	mpi_exchange(&cpu,sendbuffer,recvbuffer,2);
#endif
	printf("guess ok\n");
#if 1
	// fine level relaxation
	for(iter=0;iter<niter;iter++){
	  nextoct=firstoct[level-1];
	  if((iter%16==0)&&(cpu.rank==0)) printf("level=%d iter=%d ",level,iter);
	  int ncell=0;
	  if(nextoct!=NULL){
	    dx=pow(0.5,level);
	    do{ 
	      ncell++;
	      curoct=nextoct;
	      // First we gather the potential in all neighbors
	      for(icomp=0;icomp<=6;icomp++){
		memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
		nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
	      }

	      // Second we gather the local density
	      memset(vcomp[7],0,stride*sizeof(float)); // reset the vcomp;
	      nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);
	      
	      // 2.5 we contrast the density by removing the average value
	      remove_avg(vcomp[7],stride,1.);
	      
	      // we compute the square of the norm of the density (first iteration only)
	      if(iter==0) norm_d+=square(vcomp[7],nread);
	      
	      // Third we perform the calculation (eventually on GPU)
	      laplacian(vcomp,stride,dx);
	      
	      // Fourth we scatter back the potential estimation
	      nextoct=scattercomp(curoct, vcomp[6], 6, 1, stride,&cpu);

	      
	    }while(nextoct!=NULL);

#ifdef WMPI
	    mpi_exchange(&cpu,sendbuffer,recvbuffer,2);
	    if(iter==0) MPI_Allreduce(MPI_IN_PLACE,&norm_d,1,MPI_FLOAT,MPI_SUM,cpu.comm);
#endif
	    // Fifth we compute the residuals
	    nextoct=firstoct[level-1];
	    float res=0.;
	    if(nextoct!=NULL){
	      dx=pow(0.5,level);
	      do{ 
		curoct=nextoct;
	      
		// First we gather the potential in all neighbors + local
		for(icomp=0;icomp<=6;icomp++){
		  memset(vcomp[icomp],0,stride*sizeof(float)); // reset the vcomp;
		  nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
		}

		// Second we gather the local density
		memset(vcomp[7],0,stride*sizeof(float)); // reset the vcomp;
		nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);

		// 2.5 we contrast the density by removing the average value
		remove_avg(vcomp[7],nread,1.);
	      
		// Third we perform the square of the residual
		res+=square_res(vcomp,nread,dx);
	      
	      }while(nextoct!=NULL);
	    }

#ifdef WMPI
	    // reducing the residuals
	    MPI_Barrier(cpu.comm);
	    float restot;
	    //if(iter%64==0) printf("res = %f on rank=%d\n",res,cpu.rank);
	    MPI_Allreduce(&res,&restot,1,MPI_FLOAT,MPI_SUM,cpu.comm);
	    res=restot;
#endif
	    if((iter%16==0)&&(cpu.rank==0)) printf("dens=%e res=%e relative residual=%e\n ",sqrt(norm_d),sqrt(res),sqrt(res/norm_d));
	    if(sqrt(res/norm_d)<acc) break;
	  }
	}
#endif
      }
    }
#endif

#ifdef TIME_JAC
  fclose(ft);
#endif

  sprintf(filename,"data/potstart.%05d.p%05d",nsteps+1,cpu.rank);
  printf("%s\n",filename);
  dumpcube(lmap,firstoct,2,filename);
	    

  // ==================================== Force calculation and velocity update   // corrector step


  if(nsteps!=0){
    forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dt*0.5,&cpu,sendbuffer,recvbuffer);
  }

  // ==================================== DUMP AFTER SYNCHRONIZATION
#if 1
  if(nsteps%NDUMP==0){
    // ===== Casting rays to fill a map

  /*   sprintf(filename,"data/level.%05d",nsteps); */
  /*   dumpmap(lmap,firstoct,0,filename,0.,1.); */
  /*   sprintf(filename,"data/dens.%05d",nsteps); */
  /*   dumpmap(lmap,firstoct,1,filename,0.,1.); */
  /*   sprintf(filename,"data/pot3d.%05d",nsteps); */
  /*   dumpcube(lmap,firstoct,2,filename); */
  /*   sprintf(filename,"data/lev3d.%05d",nsteps); */
  /*   dumpcube(lmap,firstoct,0,filename); */

  
  /* //==== Gathering particles for dump */

  /*   sprintf(filename,"data/part.%05d",nsteps); */
  /*   dumppart(firstoct,filename,npart,levelcoarse,levelmax); */

  }
#endif
  
  // ==================================== Force calculation and velocity update   // predictor step
  
  forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dt*0.5,&cpu,sendbuffer,recvbuffer);
  

  printf("Moving particles\n");



  // ==================================== Moving Particles + Oct management
  

  // Computing displacement (predictor)

  movepart(levelcoarse,levelmax,firstoct,dt);


  // Moving particles through cells (3 passes)

  partcellreorg(levelcoarse,levelmax,firstoct);

#ifdef WMPI

  // Communication of particles
  int deltan;
  deltan=mpi_exchange_part(&cpu,psendbuffer,precvbuffer,&lastpart);

  // Recounting particles
  npart=npart+deltan;
#endif

  //==== Gathering particles for dump

  sprintf(filename,"data/partstart.%05d.p%05d",nsteps+1,cpu.rank);
  dumppart(firstoct,filename,npart,levelcoarse,levelmax);

#if 1
  // ==================================== Check the number of particles and octs
  multicheck(firstoct,npart,levelmax,cpu.rank);
#endif	 

  //==================================== timestep completed, looping

#ifdef WMPI
  if(nsteps==1){
    printf("ABORTING !!\n");
    MPI_Barrier(cpu.comm);
    MPI_Abort(cpu.comm,42);
  }
#else
    abort();
#endif

  }
  return 0;
}
      
