#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef WMPI
#include <mpi.h>
//#include "communication.h"
#endif
#include "Interface.h"
#include "timestep.h"
#include "Io.h"
#include "GPU.h"
#include "Allocation.h"

extern void topo_cartesian(int rank, int *dims, int *coords, int *neigh);
//**********************************************************
//**********************************************************

float *egy,*egy_new;
float *flx,*flx_new;

float *dedd;
float *src0;
int *src0pos;


float *temperature,*density,*xion;
float *buff;

float t=0.;
float dt;
int nstep,ntsteps;
float aexp;
float zexp;


int nalist;
float *alist;
float *tlist;
//**********************************************************
//**********************************************************

float *cuegy,*cuegy_new;
float *cuflx,*cuflx_new;

float *cudedd;
float *cusrc0;
int *cusrc0pos;

float *cutemperature,*cuxion;
float *cudensity;

float *cubuff;
//**********************************************************
//**********************************************************
int ncells;
float dx;
float tmax;
float courantnumber;
int ndumps;
int ndisp;
int nmax;
char rootname[256];
int verbose;
int interactive;
int nrestart;

float c;
float c_r=2.99792458e8;
float kb=1.3806e-23;
float Tcmb=2.726;
float effective_speed_of_light;
float unit_length;
float unit_time;
float unit_number;

float egy_min;
float defdens;
float deftemp;
int boundary;
int cooling;
float fudgecool;

int ncvgcool;

int nsource;
float fesc;
float aboost;
float kboost;
float a0;
float clump;

float astart=0.;
float Hubble0;
float omegam;
float omegav;

int fieldlist;
char fieldname[256];
char filealist[256];




//**********************************************************
//**********************************************************

void logo(void)
{
  puts("");
  puts("******************************************");
  puts("               CUDaton V.0.2               *");
  puts("******************************************");
  puts("");
  
}

//**********************************************************
//**********************************************************

int GetParameters(char *fparam)
{
  FILE *buf;
  char stream[256];
  buf=fopen(fparam,"r");
  if(buf==NULL)
    {
      printf("ERROR : cannot open %s, please check\n",fparam);
      return 1;
    }
  else
    {
      fscanf(buf,"%s %d",stream,&verbose);
      fscanf(buf,"%s %s",stream,rootname);
      fscanf(buf,"%s %d",stream,&nrestart);

      fscanf(buf,"%s %d",stream,&ncells);
      if(ncells!=NCELLS)
	{
	  printf("file ncells different from the Hard coded value ncells=%d NCELLS=%d\n",ncells,NCELLS);
	  abort();
	}

      fscanf(buf,"%s %f",stream,&dx);

      fscanf(buf,"%s %f",stream,&tmax);
      fscanf(buf,"%s %f",stream,&courantnumber);
      fscanf(buf,"%s %d",stream,&ndumps);
      fscanf(buf,"%s %d",stream,&ndisp);
      fscanf(buf,"%s %d",stream,&nmax);
      //printf("ndumps=%d\n",ndumps);


      fscanf(buf,"%s %f",stream,&unit_length);
      //printf("unit_length=%e\n",unit_length);

      fscanf(buf,"%s %f",stream,&unit_time);
      //printf("unit_time=%e\n",unit_time);

      fscanf(buf,"%s %f",stream,&unit_number);
      //printf("unit_time=%e\n",unit_time);

      fscanf(buf,"%s %f",stream,&effective_speed_of_light);
      //printf("sol=%e\n",effective_speed_of_light);


      fscanf(buf,"%s %f",stream,&egy_min);
      fscanf(buf,"%s %f",stream,&defdens);
      fscanf(buf,"%s %f",stream,&deftemp);
      //printf("egymin=%e\n",egy_min);

      fscanf(buf,"%s %d",stream,&boundary);
      //      printf("bound=%d\n",boundary);


      fscanf(buf,"%s %d",stream,&cooling);
      fscanf(buf,"%s %f",stream,&fudgecool);
      
      fscanf(buf,"%s %d",stream,&ncvgcool);

      fscanf(buf,"%s %d",stream,&nsource);
      fscanf(buf,"%s %f",stream,&fesc);
      fscanf(buf,"%s %f",stream,&aboost);
      fscanf(buf,"%s %f",stream,&kboost);

      fscanf(buf,"%s %f",stream,&Hubble0);
      fscanf(buf,"%s %f",stream,&omegam);
      fscanf(buf,"%s %f",stream,&omegav);

      fscanf(buf,"%s %d",stream,&fieldlist);
      //printf("%d",fieldlist);
      if(fieldlist)
	{
	  fscanf(buf,"%s %s",stream,fieldname);
	  fscanf(buf,"%s %s",stream,filealist);
	}

      fclose(buf);
      
#ifdef COSMO
      if((omegam+omegav)!=1.)
	{
	  printf(" Error omegam+omegav= %f !=1 . Only flat models please !\n",omegam+omegav);
	}
#endif

      //printf("I dx=%e ul=%e\n",dx,unit_length);
      dx=dx*unit_length;
      //printf("II dx=%e\n",dx);
      dt=dt*unit_time;
      tmax=tmax*unit_time;

    }
  return 0;
}

//**********************************************************

int main(int argc, char *argv[])
{
  int code;



#ifdef WMPI
  int mpi_rank,mpi_size,ic_rank;
  int ierr;
  int ndevice;
  int mpi_coords[3];
  int mpi_neighbors[6];
  int dims[3];

// Number of proc per dimension
  dims[0]=NGPUX;
  dims[1]=NGPUY;
  dims[2]=NGPUZ;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  ndevice=countdevices(mpi_rank);
  if(mpi_rank==0) printf("Starting with %d GPUs\n",mpi_size);

#ifdef TITANE
  //----  Some Stuff specific to Titane to set the devices properly START

  char hostname[10];
  int hostnum;
  int *vhostnum;

  gethostname(hostname, 10);
  sscanf(hostname,"curie%04d",&hostnum); // getting the host numbers
  vhostnum=(int*)calloc(mpi_size,sizeof(int));
  vhostnum[mpi_rank]=hostnum;

  //printf("gather start on %d\n",mpi_rank);
  MPI_Allgather(&hostnum,1,MPI_INT,vhostnum,1,MPI_INT,MPI_COMM_WORLD);
  //printf("on proc %d : %d %d %d %d %d %d %d %d\n",mpi_rank,vhostnum[0],vhostnum[1],vhostnum[2],vhostnum[3],vhostnum[4],vhostnum[5],vhostnum[6],vhostnum[7]);
  //MPI_Abort(MPI_COMM_WORLD,2);

  //printf("gather stop on %d\n",mpi_rank);

  initlocaldeviceTITANE(mpi_rank,ndevice,vhostnum,mpi_size);

  //printf("init ok on %d\n",mpi_rank);

  //----  Some Stuff specific to Titane to set the devices properly END
#else
  initlocaldevice(mpi_rank,2);
#endif

  topo_cartesian(mpi_rank,dims,mpi_coords,mpi_neighbors);
  ic_rank=mpi_coords[0]+mpi_coords[1]*dims[0]+mpi_coords[2]*dims[0]*dims[1];

  if(mpi_rank==0)
    {
#else 
  initlocaldevice(1,1);
#endif
      
      if(argc<2)
	{
	  puts("usage: ./caton file.param");
	  return 1;
	}

      logo();

#ifdef WMPI
    }
#endif

#ifdef TESTWRITE
  FILE *ftest;
  char fnt[256];
  int dum=1;
  
  // WRITE TEST
  sprintf(fnt,"/ccc/cont005/home/gch0003/aubertd/noerase/tata/test_%05d_0_%05d",hostnum,mpi_rank);
  ftest=fopen(fnt,"wb");
  fwrite(&dum,sizeof(int),1,ftest);
  fclose(ftest);
#endif
  
  code=GetParameters(argv[1]);
  if(code!=0) return 1;
  if(verbose) puts("Parameters ok");


#ifdef TESTWRITE
  // WRITE TEST
  sprintf(fnt,"/ccc/cont005/home/gch0003/aubertd/noerase/tata/test_%05d_1_%05d",hostnum,mpi_rank);
  ftest=fopen(fnt,"wb");
  fwrite(&dum,sizeof(int),1,ftest);
  fclose(ftest);
#endif

#ifdef WMPI
  if(mpi_rank==0) printf("proc %d : astart=%f\n",mpi_rank,astart);
#endif
  
  Allocation();

#ifdef TESTWRITE
  // WRITE TEST
  sprintf(fnt,"/ccc/cont005/home/gch0003/aubertd/noerase/tata/test_%05d_2_%05d",hostnum,mpi_rank);
  ftest=fopen(fnt,"wb");
  fwrite(&dum,sizeof(int),1,ftest);
  fclose(ftest);
#endif

  cuAllocation();
  //  if(verbose) puts("Allocation ok");
  
#ifdef TESTWRITE

  // WRITE TEST
  sprintf(fnt,"/ccc/cont005/home/gch0003/aubertd/noerase/tata/test_%05d_3_%05d",hostnum,mpi_rank);
  ftest=fopen(fnt,"wb");
  fwrite(&dum,sizeof(int),1,ftest);
  fclose(ftest);
  
#endif

#ifdef WMPI
#ifndef TEST_STROMGREN
  cuGetIC(nrestart,ic_rank);
#else
  int middle=0;
  middle=(mpi_coords[0]==dims[0]/2)*(mpi_coords[1]==dims[1]/2)*(mpi_coords[2]==dims[2]/2);
  cuGetIC(nrestart,middle);
#endif
#else
  cuGetIC(nrestart,0);
#endif

  nstep=0;
  ntsteps=nrestart;
  if(mpi_rank==0) printf("astart=%e\n",astart);

  //printf("Initial Conditions ok on mpi_rank=%d ic_rank=%d\n",mpi_rank,ic_rank);

#ifdef TESTWRITE
  // WRITE TEST
  sprintf(fnt,"/ccc/cont005/home/gch0003/aubertd/noerase/tata/test_%05d_4_%05d",hostnum,mpi_rank);
  ftest=fopen(fnt,"wb");
  fwrite(&dum,sizeof(int),1,ftest);
  fclose(ftest);
#endif

  if(fieldlist)
    {
#ifdef WMPI
      if(mpi_rank==0) puts("Start field setup");
      getalist(mpi_rank);
      if(mpi_rank==0) printf("astart in list=%e\n",alist[0]);
      if(mpi_rank==0) puts("field setup ok");
#else
      puts("Start field setup");
      getalist(0);
      if(mpi_rank==0) printf("astart in list=%e\n",alist[0]);
      puts("field setup ok");
#endif
    }
  else
    {
      //puts("no fieldlist");
    }

#ifdef TESTWRITE
  sprintf(fnt,"/ccc/cont005/home/gch0003/aubertd/noerase/tata/test_%05d_5_%05d",hostnum,mpi_rank);
  ftest=fopen(fnt,"wb");
  fwrite(&dum,sizeof(int),1,ftest);
  fclose(ftest);
#endif

#ifdef WMPI
  MPI_Barrier(MPI_COMM_WORLD);
  Mainloop(mpi_rank,mpi_coords,mpi_neighbors,ic_rank);
  MPI_Finalize();
#else
  int mpi_coords[3];
  int mpi_neighbors[6];
  Mainloop(0,mpi_coords,mpi_neighbors,0);
#endif

  return 0;
}


