
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


#ifdef WGPU
#include "interface.h"
#include "vector_gpu.h"
#endif


 //------------------------------------------------------------------------
 // the MAIN CODE
 //------------------------------------------------------------------------

#ifdef TESTCOSMO
float f_aexp(float aexp, float omegam, float omegav)
{
  return 1./sqrtf(omegam/aexp+omegav*aexp*aexp);
}
#endif

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
  float tsim;
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
  float dtnew;
  int nbnd;

  float x,y,z;
  float vx,vy,vz;
  float mass,mtot;
  float idx;
  float faexp, faexp2;
  unsigned key;

  struct CPUINFO cpu;

  struct PACKET **sendbuffer; 
  struct PACKET **recvbuffer; 

  struct PART_MPI **psendbuffer; 
  struct PART_MPI **precvbuffer; 

  struct RUNPARAMS param;

  size_t rstat;

  unsigned int *cartdict;


  //========== TEST ZONE (IF REQUIRED)==========

/*   printf("size =%d\n",sizeof(struct CELL)); */
/*   printf("size =%d\n",sizeof(struct OCT)); */
/*   abort(); */


  //=========== some initial calls =============
  GetParameters(argv[1],&param); // reading the parameters file
    

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
  levelcoarse=param.lcoarse;
  levelmax=param.lmax;

  ngridmax=3000000;
  npartmax=128*128*128;
#ifdef PART2
  npart=2;
#else
  npart=128*128*128;
#endif

  threshold=param.amrthresh;
  lmap=param.levelmap;
  stride=fmax(8,param.stride);//pow(2,levelcoarse);
  ncomp=10;
  acc=param.poissonacc;
  dt=param.dt;
  cpu.maxhash=param.maxhash;
  //breakmpi();
  //========== allocations ===================
  
  //  if(cpu.rank==0) printf("Allocating %f GB cell=%f GB part=%f GB book=%f",(sizeof(struct OCT)*ngridmax+sizeof(struct PART)*npart+cpu.maxhash*sizeof(struct OCT*)+stride*ncomp*sizeof(float))/(1024*1024*1024.),sizeof(struct OCT)*ngridmax/(1024*1024*1024.),sizeof(struct PART)*npart/(1024*1024*1024.),(cpu.maxhash*sizeof(struct OCT*)+stride*ncomp*sizeof(float))/(1024.*1024.*1024.));

  int memsize=0.;
  grid=(struct OCT*)calloc(ngridmax,sizeof(struct OCT)); memsize+=ngridmax*sizeof(struct OCT);// the oct grid
  firstoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *)); memsize+=levelmax*sizeof(struct OCT *);// the firstoct of each level
  lastoct=(struct OCT **)calloc(levelmax,sizeof(struct OCT *)); memsize+=levelmax*sizeof(struct OCT *);// the last oct of each level
  part=(struct PART*)calloc(npartmax,sizeof(struct PART)); memsize+=npartmax*sizeof(struct PART);// the particle array
  cpu.htable=(struct OCT**) calloc(cpu.maxhash,sizeof(struct OCT *)); memsize+=cpu.maxhash*sizeof(struct OCT*);// the htable keys->oct address
  cpu.noct=(int *)calloc(levelmax,sizeof(int)); memsize+=levelmax*sizeof(int);// the number of octs per level

  lastpart=part-1; // the last particle points before the first at the very beginning

  vcomp=(float **)calloc(ncomp,sizeof(float*));
  for(i=0;i<ncomp;i++)
    {
      vcomp[i]=(float *)calloc(stride,sizeof(float));
    }

  memsize+=ncomp*stride*sizeof(float);


  //===================================================================================================

  // allocating the vectorized tree
  
  float *vecpot; //contains the potential in "stride" octs
  float *vecpotnew; //contains the potential in "stride" octs
  float *vecden; //contains the density   in "stride" octs
  int *vecnei;//contains the cell neighbors of the octs
  int *vecl; // contains the level of the octs
  int *vecicoarse; // contains the level of the octs
  
#ifdef NEWJACK
  vecpot=(float*)calloc(stride*8,sizeof(float));
  vecpotnew=(float*)calloc(stride*8,sizeof(float));
  vecden=(float*)calloc(stride*8,sizeof(float));
  vecnei=(int *)calloc(stride*6,sizeof(int));
  vecl=(int *)calloc(stride,sizeof(int));
  vecicoarse=(int *)calloc(stride,sizeof(int));
  memsize+= stride*32*4;

#endif 

#ifdef WGPU
  float *vecden_d;
  int *vecl_d;
  float *vec2_d;
  float *vecsum_d;

  float *vecpot_d; 
  float *vecpotnew_d; 
  int *vecnei_d;
  int *vecicoarse_d; 

  cudaMalloc((void **)&vecl_d,sizeof(int)*stride);
  cudaMalloc((void **)&vecden_d,sizeof(float)*stride*8);
  cudaMalloc((void **)&vec2_d,sizeof(float)*stride*8);
  cudaMalloc((void **)&vecsum_d,sizeof(float)*stride*8);

  cudaMalloc((void **)&vecpot_d,sizeof(float)*stride*8);
  cudaMalloc((void **)&vecpotnew_d,sizeof(float)*stride*8);

  cudaMalloc((void **)&vecnei_d,sizeof(int)*stride*6);
  cudaMalloc((void **)&vecicoarse_d,sizeof(int)*stride);
#endif
    
  if(cpu.rank==0) printf("Allocations %f GB done\n",memsize/(1024.*1024*1024));

  //========== setting up the parallel topology ===

  // We segment the oct distributions at levelcoarse 
    cpu.bndoct=NULL;
    cpu.mpinei=NULL;
    cpu.dict=NULL;

    cpu.nbuff=param.nbuff;
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
  grid->vecpos=-1;

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

	    // vector data
	    newoct->vecpos=-1;
	    
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

  int newloadb=1;
  setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax,newloadb); // out of WMPI to compute the hash table
  newloadb=0;

#ifdef WMPI
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
#endif


  //===================================================================================================
  
  // ==== some initial dump

  /* sprintf(filename,"data/levstart.%05d.p%05d",0,cpu.rank); */
  /* dumpcube(lmap,firstoct,0,filename,0.); */
  /* sprintf(filename,"data/cpustart.%05d.p%05d",0,cpu.rank); */
  //dumpcube(lmap,firstoct,3,filename);

  // =====  computing the memory location of the last oct 

  endoct=lastoct[0];
  for(i=0;i<levelcoarse;i++) {
    if(lastoct[i]>endoct) endoct=lastoct[i];
  }


#if 1  // ==================================== assigning particles to cells
  //breakmpi();
  if(cpu.rank==0) printf("==> starting part\n");
  firstoct_currl=0;
  for(il=1;il<levelcoarse;il++) firstoct_currl+=pow(pow(2,il-1),3); // the index of the first oct of current level
 
  // initialisation of particles
  

#ifdef PART2

  int ir,nr=2;
  ip=0;
  float dxcell=1./pow(2.,levelcoarse);
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
      part[ip].idx=ir;
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

  float omegam,omegav,Hubble;
#ifdef TESTCOSMO
  int dummy;
  float dummyf;
  int npartf;

  int nploc;
  float munit;
  float ainit;
  float lbox;

  fd=fopen("utils/IC.PM.0","rb");

  rstat=fread(&dummy,sizeof(dummy),1,fd);
  rstat=fread(&nploc,sizeof(int),1,fd);
  rstat=fread(&munit,sizeof(float),1,fd);
  rstat=fread(&ainit,sizeof(float),1,fd);
  rstat=fread(&lbox,sizeof(float),1,fd);
  rstat=fread(&omegam,sizeof(float),1,fd);
  rstat=fread(&omegav,sizeof(float),1,fd);
  rstat=fread(&Hubble,sizeof(float),1,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);

  if(cpu.rank==0) printf("%f %d %f %f\n",ainit,nploc,omegav,Hubble);


  float *pos;
  float *vel;
  
  pos=(float *)malloc(sizeof(float)*3*nploc);
  vel=(float *)malloc(sizeof(float)*3*nploc);

  rstat=fread(&dummy,sizeof(dummy),1,fd); 
  rstat=fread(pos,sizeof(float),nploc,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);

  rstat=fread(&dummy,sizeof(dummy),1,fd);
  rstat=fread(pos+nploc,sizeof(float),nploc,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);

  rstat=fread(&dummy,sizeof(dummy),1,fd);
  rstat=fread(pos+2*nploc,sizeof(float),nploc,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);
  
  rstat=fread(&dummy,sizeof(dummy),1,fd);
  rstat=fread(vel,sizeof(float),nploc,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);

  rstat=fread(&dummy,sizeof(dummy),1,fd);
  rstat=fread(vel+nploc,sizeof(float),nploc,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);

  rstat=fread(&dummy,sizeof(dummy),1,fd);
  rstat=fread(vel+2*nploc,sizeof(float),nploc,fd);
  rstat=fread(&dummy,sizeof(dummy),1,fd);

  fclose(fd);

  mass=munit;
  tsim=ainit;

  ip=0.;
  for(i=0;i<nploc;i++)
    {
      x=pos[i];
      y=pos[i+nploc];
      z=pos[i+2*nploc];

      vx=vel[i];
      vy=vel[i+nploc];
      vz=vel[i+2*nploc];
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
	part[ip].idx=ip;
	lastpart=part+ip;
	ip++;
      }
      
    }
  
  npart=ip;
  free(pos);
  free(vel);
  if (cpu.rank==0) printf("cosmo readpart done with munit=%e\n",munit);
#endif


  // we set all the "remaining" particles mass to -1
  for(ii=npart;ii<npartmax;ii++) part[ii].mass=-1.0;


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

  mtot=multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);

#endif	

  // ==================================== performing the CIC assignement
#ifndef WGPU
  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#else
  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#endif

#ifdef WMPI
    // ==================================== performing the CIC BOUNDARY CORRECTION 

  mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);

  // ======================================= Density boundary mpi update 

  mpi_exchange(&cpu,sendbuffer,recvbuffer,1,1);

#endif

  mtot=multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);

  sprintf(filename,"data/partstart.%05d.p%05d",0,cpu.rank);
  dumppart(firstoct,filename,npart,levelcoarse,levelmax,tsim);
  
  sprintf(filename,"data/denstart.%05d.p%05d",0,cpu.rank);
  dumpcube(lmap,firstoct,1,filename,tsim);
  //  abort();

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

  float tmax;
#ifdef TESTCOSMO
  tmax=1.;
#else
  tmax=1000.;
#endif

  FILE *fegy;

  //breakmpi();
  for(nsteps=0;(nsteps<=param.nsteps)*(tsim<=tmax);nsteps++){
    
#ifdef TESTCOSMO
    if(cpu.rank==0) printf("\n============== STEP %d aexp=%e z=%f ================\n",nsteps,tsim,1./tsim-1.);
#else
    if(cpu.rank==0) printf("\n============== STEP %d tsim=%e ================\n",nsteps,tsim);
#endif

#if 1
    // ==================================== marking the cells
    mark_cells(levelcoarse,levelmax,firstoct,nsmooth,threshold,&cpu,sendbuffer,recvbuffer);

    // ==================================== refining (and destroying) the octs

    
    refine_cells(levelcoarse,levelmax,firstoct,lastoct,endoct,&cpu);

    
  // recomputing the last oct;
  
  printf("==> memory state\n");
  printf("endoct=%p\n",endoct);
  
  endoct=lastoct[0];
  for(i=0;i<levelmax;i++) {
    if(lastoct[i]>endoct) endoct=lastoct[i];
  }

#endif

/*   sprintf(filename,"data/levstart.%05d.p%05d",nsteps+1,cpu.rank); */
/*   dumpcube(lmap,firstoct,0,filename); */



#ifdef WMPI
  // ==================================== after refinement we should remap the boundary cells
  setup_mpi(&cpu,firstoct,levelmax,levelcoarse,ngridmax,newloadb);
#endif



  // ==================================== performing the CIC assignement

#ifndef WGPU
  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#else
  call_cic(levelmax,levelcoarse,firstoct,&cpu);
#endif

#ifdef WMPI
    // ==================================== performing the CIC BOUNDARY CORRECTION 

  mpi_cic_correct(&cpu,sendbuffer,recvbuffer,0);

  // ======================================= Density boundary mpi update 

  mpi_exchange(&cpu,sendbuffer,recvbuffer,1,1);
#endif

    // cleaning the marks
    for(level=1;level<=levelmax;level++) // looping over levels
      {
	/* float maxd=0.,mind=1e30,avg=0.; */
	/* int ncell=0; */
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		curoct->cell[icell].marked=0.;
		/* ncell++; */
		/* avg+=curoct->cell[icell].density; */
		/* if(curoct->cell[icell].density>maxd) maxd=curoct->cell[icell].density; */
		/* if(curoct->cell[icell].density<mind) mind=curoct->cell[icell].density; */
	      }
	  }while(nextoct!=NULL);

	//printf("level=%d avg=%e mind=%e maxd=%e\n",level,avg/ncell,mind,maxd);
      }



  /*   sprintf(filename,"data/denstart.%05d.p%05d",nsteps+1,cpu.rank); */
  /* printf("%s\n",filename); */
  /* dumpcube(lmap,firstoct,1,filename); */


  // ==================================== Check the number of particles and octs
  
    mtot=multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);
#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&mtot,1,MPI_FLOAT,MPI_SUM,cpu.comm);
#endif


// ==================================== POISSON Testing the jacobi iteration

#ifdef TIME_JAC
  FILE *ft;
  //  sprintf(filename,"data/timejac.%05d.p%05d",nsteps+1,cpu.rank);
  //ft=fopen(filename,"w");
  
  double tg1,tg2;
  double tl1,tl2;
  double ta1,ta2;
  double ts1,ts2;
  double tm1,tm2;
  double t1,t2;
  double tc1;
  double tg,tl,ts,ta;
#endif



  if(cpu.rank==0){
    printf("=======================================\n");
    printf("==> Poisson Start with ncell coarse =%d \n",cpu.noct[levelcoarse-1]*8);
    printf("=======================================\n");
  }
  int icomp,iter,niter=param.niter;
  float norm_d;
  int nread;
  char compnei[8];
  char compvar[8];
  int ic;

  for(level=levelcoarse;level<=levelmax;level++)
    {
      norm_d=0.;
	
      /* 	// FINE LEVEL ONLY initial guess from parent cell */
      if(level>levelcoarse){
	nextoct=firstoct[level-1];
	if(nextoct==NULL){
	  printf("Proc %d skipping\n",cpu.rank);
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
	mpi_exchange(&cpu,sendbuffer,recvbuffer,2,1);
#endif
      }


#ifndef NEWJACK
      // looping over iterations
      for(iter=0;iter<niter;iter++){
	  
	int vpass=0;
	  
	nextoct=firstoct[level-1];
	float res=0.; // the square of the residual

	double tt[10];

	//start
	if(nextoct!=NULL){
	  dx=pow(0.5,level);
	  do{ 
	    curoct=nextoct;
	      
	  
	    // ==== We gather vector data if this the first iteration or at each iteration if the stride is too small
	    //First we gather the potential in all neighbors
#ifdef WMPI
	    tt[0]=MPI_Wtime();
#endif
	    /* for(icomp=0;icomp<6;icomp++){ */
	    /*   nextoct=gathercompempty(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread); */
	    /* } */
#ifdef WMPI	      
	    tt[1]=MPI_Wtime();
#endif
	    for(icomp=0;icomp<6;icomp++){
	      nextoct=gathercomp(curoct, vcomp[icomp], icomp, 1, stride,&cpu,&nread);
	    }
#ifdef WMPI	    
	    tt[2]=MPI_Wtime();
#endif
	    if((stride<8*cpu.noct[level-1])||(iter==0)){
	      // Second we gather the local density and the local potential
	      nextoct=gathercomp(curoct, vcomp[6], 6, 1, stride,&cpu,&nread);
	      nextoct=gathercomp(curoct, vcomp[7], 6, 0, stride,&cpu,&nread);
	    }
#ifdef WMPI
	    tt[3]=MPI_Wtime();
#endif
	  
	    // 2.5 ==== we contrast the density by removing the average density value 
	    if((stride<8*cpu.noct[level-1])||(iter==0)) remove_avg(vcomp[7],nread,1.);

	    // we compute the square of the norm of the density (first iteration only)
#ifdef TESTCOSMO
	    if(iter==0) norm_d+=square(vcomp[7],nread)*pow(1.5*omegam/tsim,2);
#else
	    if(iter==0) norm_d+=square(vcomp[7],nread)*pow(4*M_PI,2);
#endif

#ifdef WMPI
	    tt[4]=MPI_Wtime();
#endif
	    // Third we perform the calculation (eventually on GPU) also returns the residual
	    float dummy;
	    if(stride<8*cpu.noct[level-1]){
#ifdef TESTCOSMO
	      dummy=laplaciancosmo(vcomp,nread,dx,8,omegam,tsim);
#else
	      dummy=laplacian(vcomp,nread,dx,8);
#endif
	    }
	    else{
#ifdef TESTCOSMO
	      dummy=laplaciancosmo(vcomp,nread,dx,6,omegam,tsim);
#else
	      dummy=laplacian(vcomp,nread,dx,6);
#endif
	    }
	    res+=dummy;
#ifdef WMPI
	    tt[5]=MPI_Wtime();
#endif
	    if(stride<8*cpu.noct[level-1]){
	      // Fourth we scatter back the potential estimation to the temp position 
	      nextoct=scattercomp(curoct, vcomp[8], 6, 2, stride,&cpu);
	    }

#ifdef WMPI
	    tt[6]=MPI_Wtime();
	    //fprintf(ft,"%d %d %e %e %e %e %e %e %e\n",level,stride,tt[6]-tt[0],tt[1]-tt[0],tt[2]-tt[1],tt[3]-tt[2],tt[4]-tt[3],tt[5]-tt[4],tt[6]-tt[5]);
#endif
	  }while(nextoct!=NULL);
	}

	// ===============  we copy the result in the temp position to the potential 
	  
	if(stride<8*cpu.noct[level-1]){
	  nextoct=firstoct[level-1];
	  if(nextoct!=NULL){
	    do{ 
	      curoct=nextoct;
	      nextoct=gathercomp(curoct, vcomp[0], 6, 2, stride,&cpu,&nread); // getting the data in the temp field
	      nextoct=scattercomp(curoct, vcomp[0], 6, 1, stride,&cpu);
		
	    }while(nextoct!=NULL);
	  }
	}
	else{
	  nextoct=firstoct[level-1];
	  if(nextoct!=NULL){
	    do{ 
	      curoct=nextoct;
	      nextoct=scattercomp(curoct, vcomp[6], 6, 1, stride,&cpu);
	    }while(nextoct!=NULL);
	  }
	}
	  
	    
#ifdef WMPI
	// ====We communicate the buffer
	mpi_exchange(&cpu,sendbuffer,recvbuffer,2,(iter==0));
	if(iter==0) MPI_Allreduce(MPI_IN_PLACE,&norm_d,1,MPI_FLOAT,MPI_SUM,cpu.comm);

	// reducing the residuals
	float restot;
	//if(iter%64==0) printf("res = %f on rank=%d\n",res,cpu.rank);
	MPI_Allreduce(&res,&restot,1,MPI_FLOAT,MPI_SUM,cpu.comm);
	res=restot;
#endif
	//if((iter%64==0)&&(cpu.rank==0)) printf("level=%2d iter=%4d dens=%e res=%e relative residual=%e\n ",level,iter,sqrt(norm_d),sqrt(res),sqrt(res/norm_d));
	// if the level is absent on all processors we skip */
	if((res==0.)&&(norm_d==0.)){
	    if(cpu.rank==0) printf("Level %d skipped\n",level);
	    break;
	}

	if((iter%128==0)&&(cpu.rank==0)) printf("level=%2d iter=%4d relative residual=%e res=%e den=%e \n ",level,iter,sqrt(res/norm_d),sqrt(res),sqrt(norm_d));

	// convergence achieved
	if(sqrt(res/norm_d)<acc) {
	  break;
	}

      }
#else

      // cleaning the vector positions
      clean_vec(levelmax,firstoct);

      // some checks

      if((stride<pow(2,levelcoarse))){
	  printf("still debugging\n");
	  abort();
	}

      memset(vecpot,0,sizeof(float)*stride*8);
      memset(vecpotnew,0,sizeof(float)*stride*8);
      memset(vecden,0,sizeof(float)*stride*8);
      memset(vecnei,0,sizeof(int)*stride*6);
      memset(vecl,0,sizeof(int)*stride);
      memset(vecicoarse,0,sizeof(int)*stride);

      // setting neighbors to -1
      for(i=0;i<stride*6;i++){vecnei[i]=-1;}

      dx=pow(0.5,level);

      curoct=firstoct[level-1];

      // gathering the neighbors

#ifdef NEWJACK2
      nextoct=gathervecnei2(curoct,vecnei,vecpot,1,vecl,stride,&cpu,&nread);
#else
      nextoct=gathervecnei(curoct,vecnei,vecpot,1,vecl,stride,&cpu,&nread);
#endif


      // gathering the values
      

#ifdef NEWJACK2
      if(level>levelcoarse){
	//gathering the data from coarse levels
	nextoct=gathervec2(firstoct[level-2],vecden,0,vecl,vecicoarse,stride,&cpu,&nread); // density
	nextoct=gathervec2(firstoct[level-2],vecpot,1,vecl,vecicoarse,stride,&cpu,&nread); // density
      }

      nextoct=gathervec2(curoct,vecden,0,vecl,vecicoarse,stride,&cpu,&nread); // density
      nextoct=gathervec2(curoct,vecpot,1,vecl,vecicoarse,stride,&cpu,&nread); // density

#else
      nextoct=gathervec(curoct,vecden,0,vecl,stride,&cpu,&nread); // density
      nextoct=gathervec(curoct,vecpot,1,vecl,stride,&cpu,&nread); // potential
#endif

      //printf("nread=%d\n",nread);
      // we contrast the density by removing the average value
#ifndef WGPU      
      remove_valvec(vecden,nread,stride,1.,level,vecl);
#else
      CPU2GPU(vecden_d,vecden,sizeof(float)*stride*8);
      CPU2GPU_INT(vecl_d  ,vecl,sizeof(int)*stride);
      remove_valvec_GPU(vecden_d,nread,stride,1.,level,vecl_d);
      GPU2CPU(vecden,vecden_d,sizeof(float)*stride*8);
#endif

      // we square and sum the density
#ifndef WGPU
      // CPU CASE
#ifdef TESTCOSMO
      norm_d+=square_vec(vecden,nread,stride,level,vecl)*pow(1.5*omegam/tsim,2);
#else
      norm_d+=square_vec(vecden,nread,stride,level,vecl)*pow(4.0*M_PI,2);
#endif

#else
      // GPU CASE
#ifdef TESTCOSMO
      norm_d+=square_vec_GPU(vecden_d,nread,stride,level,vecl_d,vec2_d,vecsum_d)*pow(1.5*omegam/tsim,2);
#else
      norm_d+=square_vec_GPU(vecden_d,nread,stride,level,vecl_d,vec2_d,vecsum_d)*pow(4.0*M_PI,2);
#endif
      
#endif


#ifdef WGPU
      // sending data to GPU
      CPU2GPU(vecpotnew_d,vecpotnew,sizeof(float)*stride*8);
      CPU2GPU(vecpot_d,vecpot,sizeof(float)*stride*8);
      CPU2GPU_INT(vecnei_d,vecnei,sizeof(int)*stride*6);
      CPU2GPU_INT(vecicoarse_d,vecicoarse,sizeof(int)*stride);
#endif

      //======================== looping over iterations

      float res;

      for(iter=0;iter<niter;iter++){
	
	// computing the laplacian
#ifdef NEWJACK2
#ifdef WGPU
	cudaMemset(vec2_d,0,sizeof(float)*stride*8);
	cudaMemset(vecsum_d,0,sizeof(float)*stride*8);
	res=laplacian_vec2_GPU(vecden_d,vecpot_d,vecpotnew_d,vecnei_d,vecl_d, vecicoarse_d,level,nread,stride,dx,omegam,tsim,vec2_d,vecsum_d);
#else
	res=laplacian_vec2(vecden,vecpot,vecpotnew,vecnei,vecl,vecicoarse,level,nread,stride,dx,omegam,tsim);
#endif
#else
	res=laplacian_vec(vecden,vecpot,vecpotnew,vecnei,vecl,level,nread,stride,dx,omegam,tsim);
#endif



	// exchange potential evaluations
#ifdef WGPU
	GPU2GPU(vecpot_d,vecpotnew_d,sizeof(float)*stride*8);
#else
	memcpy(vecpot,vecpotnew,sizeof(float)*stride*8);
#endif

	// skip level if it does not exist
	if((res==0.)&&(norm_d==0.)){
	    if(cpu.rank==0) printf("Level %d skipped\n",level);
	    break;
	}

	// some verbosity
	if((iter%128==0)&&(cpu.rank==0)) printf("level=%2d iter=%4d relative residual=%e res=%e den=%e \n ",level,iter,sqrt(res/norm_d),sqrt(res),sqrt(norm_d));
	
	// breaking condition
	if(sqrt(res/norm_d)<acc) {
	  break;
	}
      }

#ifdef WGPU
      GPU2CPU(vecpot,vecpot_d,sizeof(float)*8*stride);
#endif
      //scatter back the result
      nextoct=scattervec(curoct,vecpot,1,stride,&cpu,nread); // density      



#endif
    }      
  
#ifdef TIME_JAC
  //fclose(ft);
#endif


  // ==================================== Force calculation and velocity update   // corrector step


  if(nsteps!=0){
    printf("Corrector\n");
#ifdef TESTCOSMO
    faexp=f_aexp(tsim,omegam,omegav);
#else
    faexp=1.0;
#endif

    forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dt*0.5*faexp,&cpu,sendbuffer,recvbuffer);
  }

#ifdef EGYCSV
  // ==================================== Energy Conservation Test
  
  float egy;
  egy=egypart(levelcoarse,levelmax,firstoct,&cpu);
  
#ifdef WMPI
  MPI_Allreduce(MPI_IN_PLACE,&egy,1,MPI_FLOAT,MPI_SUM,cpu.comm);
#endif

  if(cpu.rank==0){
    printf("== Energy check etot=%e\n",egy);
    if(nsteps==0){
      fegy=fopen("egy.txt","w");
    }
    fprintf(fegy,"%e %e\n",tsim,egy);
  }

#endif  

  // ==================================== DUMP AFTER SYNCHRONIZATION
#if 1
  if(nsteps%(param.ndumps)==0){
    // ===== Casting rays to fill a map

    sprintf(filename,"data/pot3d.%05d.p%05d",nsteps,cpu.rank);
    dumpcube(lmap,firstoct,2,filename,tsim);
    sprintf(filename,"data/den3d.%05d.p%05d",nsteps,cpu.rank);
    dumpcube(lmap,firstoct,1,filename,tsim);
    sprintf(filename,"data/lev3d.%05d.p%05d",nsteps,cpu.rank);
    dumpcube(lmap,firstoct,0,filename,tsim);

#ifdef WMPI
    sprintf(filename,"data/cpu3d.%05d.p%05d",nsteps,cpu.rank);
    dumpcube(lmap,firstoct,3,filename,tsim);
#endif
  
  //==== Gathering particles for dump

    sprintf(filename,"data/part.%05d.p%05d",nsteps,cpu.rank);
    dumppart(firstoct,filename,npart,levelcoarse,levelmax,tsim);

  }
#endif
  
  // ==================================== Force calculation and velocity update   // predictor step
  
#ifdef TESTCOSMO
  faexp=f_aexp(tsim,omegam,omegav);
#else
  faexp=1.0;
#endif
  
  printf("Predictor\n");
  forcevel(levelcoarse,levelmax,firstoct,vcomp,stride,dt*0.5*faexp,&cpu,sendbuffer,recvbuffer);

  // ==================================== Moving Particles + Oct management
  

  printf("Moving particles\n");
  // Computing displacement (predictor)

#ifdef TESTCOSMO
  faexp2=f_aexp(tsim+dt*0.5,omegam,omegav)/((tsim+dt*0.5)*(tsim+dt*0.5));
#else
  faexp2=1.0;
#endif

  dtnew=movepart(levelcoarse,levelmax,firstoct,dt*faexp2,&cpu);
  dt=dtnew/faexp2;
  printf("dt=%e\n",dt);
  // Moving particles through cells (3 passes)

#ifdef WGPU
  partcellreorg(levelcoarse,levelmax,firstoct);
#else
  partcellreorg(levelcoarse,levelmax,firstoct);
#endif

#ifdef WMPI

  // Communication of particles
  int deltan;
  deltan=mpi_exchange_part(&cpu,psendbuffer,precvbuffer,&lastpart);

  // Recounting particles
  npart=npart+deltan;
#endif

  //==== Gathering particles for dump

  //  sprintf(filename,"data/partstart.%05d.p%05d",nsteps+1,cpu.rank);
  //dumppart(firstoct,filename,npart,levelcoarse,levelmax);

#if 1
  // ==================================== Check the number of particles and octs
  multicheck(firstoct,npart,levelcoarse,levelmax,cpu.rank,cpu.noct);
#endif	 

  //==================================== timestep completed, looping
  tsim+=dt;

  }


#ifdef NEWJACK
      free(vecpot); //contains the potential in "stride" octs
      free(vecpotnew); //contains the potential in "stride" octs
      free(vecden); //contains the density   in "stride" octs
      free(vecnei);//contains the cell neighbors of the octs
      free(vecl); // contains the level of the octs
      free(vecicoarse); // contains the level of the octs

#ifdef WGPU
      cudaFree(vecl_d);
      cudaFree(vecden_d);
#endif

#endif


#ifdef WMPI
  MPI_Barrier(cpu.comm);
  
#endif
  if(cpu.rank==0){
    printf("Done .....\n");
  }

#ifdef EGYCSV
  fclose(fegy);
#endif

#ifdef WMPI
  MPI_Barrier(cpu.comm);
  //breakmpi();
  MPI_Finalize();
#endif
  return 0;
}
      
