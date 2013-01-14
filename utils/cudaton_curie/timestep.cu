#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cutil.h>
#include "params.h"
#include "common.h"
#include "bnd.h"
#include "GPU.h"
#include "Boundary.h"
#include "cosmo.h"
#include "Allocation.h"
#include "Io.h"
#include "Explicit.h"
#include "Atomic.h"
#ifdef WMPI
#include "communication.h"
#include "Interface.h"
#endif



//**********************************************************
//**********************************************************

extern "C" int Mainloop(int rank, int *pos, int *neigh, int ic_rank);

//**********************************************************
//**********************************************************


#define CUERR() //printf("\n %s on %d \n",cudaGetErrorString(cudaGetLastError()),ic_rank)

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define NCELLS3 (NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)

#define N_INT 2048
#define A_INT_MAX 0.166667


//**********************************************************
//**********************************************************


int Mainloop(int rank, int *pos, int *neigh, int ic_rank)
{

#ifdef TIMINGS
    FILE *timefile;
    char ftime[256];
    sprintf(ftime,"time.out.%05d",rank);
    if(rank==0) timefile=fopen(ftime,"w");
  
#endif

  if(rank==0) printf("Mainloop entered by proc %d\n",rank);

  float tnext;


  dim3 blockion(NCELLX);           // USED BY IONISATION
  dim3 gridion(NCELLY,NCELLZ);

  dim3 bcool(BLOCKCOOL);           // USED BY COOLING
  dim3 gcool(GRIDCOOLX,GRIDCOOLY);
  
  dim3 blocksimple(NCELLX);        // USED BY ADVECTION THREADS
  dim3 gridsimple(NCELLY,NCELLZ);


#ifdef SDISCRETE
  int nthreadsource=min(nsource,128);
  dim3 gridsource((int)(round((float)(nsource)/float(nthreadsource))));
  dim3 blocksource(nthreadsource);
#endif

#ifndef WMPI

  dim3 blockboundx(NCELLY);
  dim3 gridboundx(NCELLZ);

  dim3 blockboundy(NCELLX);
  dim3 gridboundy(NCELLZ);

  dim3 blockboundz(NCELLX);
  dim3 gridboundz(NCELLY);

for (int igrp=0;igrp<NGRP;igrp++)
	{
  if(boundary==0) // transmissive boundary conditions
    {
      cusetboundarytrans_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundarytrans_yp<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundarytrans_zp<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundarytrans_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundarytrans_ym<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundarytrans_zm<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }
  else if(boundary==1) // reflexive boundary conditions
    {
      cusetboundaryref_zp<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_zm<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_yp<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_ym<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }
  else if(boundary==2) // Periodic boundary conditions
    {
      cusetboundaryper_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_yp<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_zp<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_ym<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_zm<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }
  else if(boundary==3) // Mixed boundary conditions
    {
      cusetboundarytrans_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_yp  <<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_zp  <<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundarytrans_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_ym  <<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_zm  <<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }
  }
#else
  dim3 blockboundx(NCELLY);
  dim3 gridboundx(NCELLZ);

  dim3 blockboundy(NCELLX);
  dim3 gridboundy(NCELLZ);

  dim3 blockboundz(NCELLX);
  dim3 gridboundz(NCELLY);

  if(neigh[5]!=rank)  
    {  
      for (int igrp=0;igrp<NGRP;igrp++){
	cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
	mpisynch();
	exchange_zp(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[2]%2);
      }
      for (int igrp=0;igrp<NGRP;igrp++){
	cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
	mpisynch();
	exchange_zm(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[2]%2);
      }
    }
  else
    {
      for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_zp<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_zm<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }
  
  if(neigh[3]!=rank)
    {
      for (int igrp=0;igrp<NGRP;igrp++){
	cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
	mpisynch();
	exchange_yp(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[1]%2);
      }
      for (int igrp=0;igrp<NGRP;igrp++){
	cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
	mpisynch();
	exchange_ym(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[1]%2);
      }
    }
  else
    {
      for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_yp<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_ym<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }

  if(neigh[1]!=rank)
    {
      for (int igrp=0;igrp<NGRP;igrp++){
	cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
	mpisynch();
	exchange_xp(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[0]%2);
      }
      for (int igrp=0;igrp<NGRP;igrp++){
	cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
	mpisynch();
	exchange_xm(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[0]%2);
      }

    }
  else
    {
      for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }

  if(boundary==0)
    {
      if(pos[0]==0)   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      if(pos[1]==0)   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_ym<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      if(pos[2]==0)   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_zm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
	  
      if(pos[0]==(NGPUX-1))   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      if(pos[1]==(NGPUY-1))   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_yp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      if(pos[2]==(NGPUZ-1))   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_zp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3) ;
	  
    }
    
#endif


#ifndef COSMO
  dt=courantnumber*dx/3./c;
  if(rank==0) printf("dx=%e cfl=%e dt=%e\n",dx,courantnumber,dt);
  tnext=t;//+ndumps*dt;
#else
  aexp=astart;
#ifndef FLAT_COSMO
  t=a2tgen(aexp,omegam,omegav,Hubble0);// Hubble0 in sec-1
#else
  t=a2t(aexp,omegav,Hubble0);// Hubble0 in sec-1
#endif
  
  tnext=t;
  float tstart=t;
  if(rank==0) printf("aexp= %f tstart=%f tmax=%f\n",aexp,t/unit_time,tmax/unit_time);

#ifndef FLAT_COSMO
  if(rank==0) printf("Building Expansion factor table");

  float da=(A_INT_MAX-aexp)/N_INT;
  float a_int[N_INT],t_int[N_INT];
  for(int i_int=0;i_int<N_INT;i_int++)
    {
      a_int[i]=aexp+i_int*da;
      t_int[i]=a2tgen(a_int[i],omegam,omegav,Hubble0); // Hubble0 in sec-1
    }

  int n_int=0;

#endif


#endif
  
  // some variables for field update
  int changefield=0;
  int forcedump;
  int ifield=0; // 2 because tfield stores the NEXT field
  float tfield;
  if(fieldlist){
    while(aexp>alist[ifield])
      {
	ifield++;
	if(rank==0) printf("a=%e afield=%e ifield=%d\n",t,alist[ifield],ifield);
      }
    tfield=tlist[ifield];
    if(rank==0) printf("ICs (astart=%f) between field #%d (a=%f) and field #%d (a=%f)\n",aexp,ifield-1,alist[ifield-1],ifield,alist[ifield]);
    if(rank==0) printf("starting with NEXT field #%d @ afield =%f with astart=%f\n",ifield,alist[ifield],aexp);// -1 because tfield stores the NEXT field
  }

  // log file
  /* FILE *logfile; */
  /* char flog[256]; */
  /* sprintf(flog,"log.out.%05d",rank); */
  /* logfile=fopen(flog,"w"); */



  //float ft=1./powf(2.,20);
  float ft=1.;
#ifdef COSMO
  float factfesc=1.;
#endif

float *factgrp;
factgrp=(float*)malloc(NGRP*sizeof(float));
FACTGRP;

  unsigned int timer;
  float q0=0.,q1=0.,q3;
#ifdef TIMINGS
  float q4,q7,q8,q9,q10,q11;
  double time_old,time_new;
#endif  
  if(rank==0)
    {
      cutCreateTimer(&timer);
      cutStartTimer(timer);

    }

  
  // MAIN LOOP STARTS HERE ======================================================>>>>
  // ============================================================================>>>>
  // ============================================================================>>>>
  // ============================================================================>>>>
  // ============================================================================>>>>
  // ============================================================================>>>>

  cudaThreadSynchronize();
#ifdef WMPI	  
  mpisynch();
#endif
  
  //  cuDumpResults(0,t,aexp,0);

  while(t<=tmax)
    {  
      

      cudaThreadSynchronize();
#ifdef WMPI	  
      get_elapsed(&time_old);
      mpisynch();
#endif
      if(rank==0)
	{
	  q3=q1-q0;
	  q0=cutGetTimerValue(timer);
	}
      

#ifndef COSMO
      dt=courantnumber*dx/3./c*ft;
      if(((nstep%ndisp)==0)&&(rank==0))
	{
	  printf(" ------------------ \n");
	  printf(" Step= %d Time= %f dt=%f tnext=%f cgpu (msec)=%f\n",nstep,t/unit_time,dt/unit_time,tnext/unit_time,q3);
	  printf(" ------------------ \n");
	}
#else
      dt=courantnumber*dx/3./c*ft;

      // computing boost factor
      float boost;
      if(aboost==0.){
	boost=1.;
      }
      else{
	boost=max(1.,aboost*exp(kboost/(t/unit_time)));
      }

      if(((nstep%ndisp)==0)&&(rank==0))
	{
	  printf(" ------------------------------\n");
	  printf(" Step= %d Time= %f Elapsed= %f dt= %f aexp=%f z=%f fesc=%f clump= %f boost=%e Next tfield=%f cgpu=%f\n",nstep,t/unit_time,(t-tstart)/unit_time,dt/unit_time,aexp,1./aexp-1.,factfesc*fesc,clump,boost,tfield/unit_time,q3);
	  printf(" ----------------------------- \n");
	  //fprintf(logfile,"%d %f %f %f %f %f %f %f\n",nstep,t/unit_time,(t-tstart)/unit_time,dt/unit_time,aexp,1./aexp-1.,tfield/unit_time,q3);
	}
#endif
      

      if(fieldlist)
	{
	  // we must not go further than the next field
	  if(dt>=tfield-t)
	    {
#ifdef WMPI
	      if(rank==0) printf("last timestep with field #%d : next field= %f t=%f t+dt=%f\n",ifield,tfield/unit_time,t/unit_time,(t+dt)/unit_time);

	      if(((tfield-t)/unit_time)==0.)
		{
		  if(rank==0) printf("WARNING FIELD DT=O -> switch immediatly to next field\n"); 
		  cuGetField(ifield,ic_rank);
		  changefield=0;
		  ifield++;
		  tfield=tlist[ifield];
		  ft=1./powf(2.,20);
		}
	      else
		{
		  changefield=1;
		  dt=tfield-t;
		  if(rank==0) printf("dt set to %f\n",dt/unit_time);
		}
#else
	      if(rank==0) printf("last timestep with field #%d : next field= %f t=%f t+dt=%f\n",ifield,tfield/unit_time,t/unit_time,(t+dt)/unit_time);

	      if(((tfield-t)/unit_time)==0.)
		{
		  if(rank==0) printf("WARNING FIELD DT=O -> switch immediatly to next field\n"); 
		  cuGetField(ifield,ic_rank);
		  changefield=0;
		  ifield++;
		  tfield=tlist[ifield];
		  ft=1./powf(2.,20);
		}
	      else
		{
		  changefield=1;
		  dt=tfield-t;
		  if(rank==0) printf("dt set to %f\n",dt/unit_time);
		}
#endif
	    }
	}

      //================================== UNSPLIT 3D SCHEME=============================
      
#ifdef LOGTOUCH
      {
	FILE *logfile;
	char flog[256];
	sprintf(flog,"log.out.1.n%d.%05d",nstep,rank);
	logfile=fopen(flog,"w");
	fclose(logfile);
      }
#endif
	for (int igrp=0;igrp<NGRP;igrp++)
		{

		  //if(rank==0) printf("igrp=%d fact=%e)\n",igrp,factgrp[igrp]);
		#ifdef COSMO
		  cuComputeELF<<<gridsimple,blocksimple>>>(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cuegy_new+igrp*NCELLS3, c, dx, dt, nstep,aexp,egy_min*factgrp[igrp]);
#else
		  cuComputeELF<<<gridsimple,blocksimple>>>(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cuegy_new+igrp*NCELLS3, c, dx, dt, nstep,1.,egy_min*factgrp[igrp]);
		#endif
		  
		  cudaThreadSynchronize();
		  CUERR();
		  if(verbose) puts("Hyperbolic Egy ok");
		  
#ifdef COSMO
		  cuComputeF_TOTAL_LF<<<gridsimple,blocksimple>>>(cuflx+igrp*NCELLS3*3,cuflx_new+igrp*NCELLS3*3,c,dx,dt,nstep,cuegy+igrp*NCELLS3, aexp);
#else
		  cuComputeF_TOTAL_LF<<<gridsimple,blocksimple>>>(cuflx+igrp*NCELLS3*3,cuflx_new+igrp*NCELLS3*3,c,dx,dt,nstep,cuegy+igrp*NCELLS3,1.);
#endif
		  cudaThreadSynchronize();
		  CUERR();
	
		#ifdef SDISCRETE
		#ifdef COSMO
		      if(kf!=0.) factfesc=exp(kf-powf(aexp/a0,af));
		      cuAddSource<<<gridsource,blocksource>>>(cuegy_new+igrp*NCELLS3,cuflx_new+igrp*NCELLS3*3,cusrc0,cusrc0pos,dt*fesc*factfesc*factgrp[igrp],dx,nsource,aexp,c);
		#else
		      cuAddSource<<<gridsource,blocksource>>>(cuegy_new+igrp*NCELLS3,cuflx_new+igrp*NCELLS3*3,cusrc0,cusrc0pos,dt*fesc*factgrp[igrp],dx,nsource,1.,c);
		#endif
		
		      CUERR();
		      if(verbose) puts("Add Source ok");
		#endif
		
		      if(verbose) puts("Hyperbolic Flux ok");
		
		      cudaThreadSynchronize();

		}

#ifdef TIMINGS     
#ifdef WMPI	  
      mpisynch();
#endif
      if(rank==0)
	{
	  q11=cutGetTimerValue(timer);
	}
#endif

	
#ifdef TESTCOOL  
#ifdef COSMO
      cuComputeIon<<<gridion,blockion>>>(cuegy_new, cuflx_new, cuxion, cudensity, cutemperature, dt/cooling, c, egy_min,unit_number,aexp);
#else
      cuComputeIon<<<gridion,blockion>>>(cuegy_new, cuflx_new, cuxion, cudensity, cutemperature, dt/cooling, c, egy_min,unit_number,1.);
#endif
#endif
      CUERR();
      if(verbose) puts("Chemistry     ok");
      cudaThreadSynchronize();
#ifdef WMPI
      mpisynch();
#endif

#ifdef TIMINGS
      if(rank==0)
	{
	  q4=cutGetTimerValue(timer);
	}
#endif

#ifdef LOGTOUCH
      {
	FILE *logfile;
	char flog[256];
	sprintf(flog,"log.out.2.n%d.%05d",nstep,rank);
	logfile=fopen(flog,"w");
	fclose(logfile);
      }
#endif 

	  // Here cuegy is used to store the temperature
#ifdef COSMO
      float hubblet=Hubble0*sqrtf(omegam/aexp+omegav*(aexp*aexp))/aexp;
      cuComputeTemp<<<gcool,bcool>>>( cuxion, cudensity, cutemperature, cuegy_new, fudgecool, c, dt/cooling, unit_number, ncvgcool, aexp, hubblet, cuflx_new, clump,egy_min,fesc,boost,cusrc0);
#else
      cuComputeTemp<<<gcool,bcool>>>( cuxion, cudensity, cutemperature, cuegy_new, fudgecool, c, dt/cooling, unit_number, ncvgcool, 1.,   0., cuflx_new, clump,egy_min,fesc,boost,cusrc0);
#endif
      CUERR();
      if(verbose) puts("Cooling  ok");
      cudaThreadSynchronize();
#ifdef WMPI	  
      mpisynch();
#endif

#ifdef TIMINGS
      cudaThreadSynchronize();
#ifdef WMPI
      mpisynch();
#endif
      if(rank==0)
	{
	  q8=cutGetTimerValue(timer);
	}
#endif

#ifdef LOGTOUCH
      {
	FILE *logfile;
	char flog[256];
	sprintf(flog,"log.out.3.n%d.%05d",nstep,rank);
	logfile=fopen(flog,"w");
	fclose(logfile);
      }
#endif
      
      cudaMemcpy(cuegy,cuegy_new,NCELLS3*sizeof(float)*NGRP,cudaMemcpyDeviceToDevice);
      cudaMemcpy(cuflx,cuflx_new,NCELLS3*sizeof(float)*3*NGRP,cudaMemcpyDeviceToDevice);

      /* cudaMemcpy(egy,cuegy,NCELLS3*sizeof(float)*NGRP,cudaMemcpyDeviceToHost); */
      /* printf("Rank=%d cuegy in =%e\n",rank,egy[2395605]); */
#ifdef LOGTOUCH
      {
	FILE *logfile;
	char flog[256];
	sprintf(flog,"log.out.4.n%d.%05d",nstep,rank);
	logfile=fopen(flog,"w");
	fclose(logfile);
      }
#endif 

#ifdef TIMINGS
      cudaThreadSynchronize();
#ifdef WMPI
      mpisynch();
#endif
      if(rank==0)
	{
	  q10=cutGetTimerValue(timer);
	}
#endif


      if(verbose) puts("Dealing with boundaries");


#ifndef WMPI
for (int igrp=0;igrp<NGRP;igrp++)
	{
  if(boundary==0) // transmissive boundary conditions
    {
      cusetboundarytrans_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NGRP*3);
      cusetboundarytrans_yp<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NGRP*3);
      cusetboundarytrans_zp<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NGRP*3);
      cusetboundarytrans_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NGRP*3);
      cusetboundarytrans_ym<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NGRP*3);
      cusetboundarytrans_zm<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NGRP*3);
    }
  else if(boundary==1) // reflexive boundary conditions
    {
      cusetboundaryref_zp<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_zm<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_yp<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_ym<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }
  else if(boundary==2) // Periodic boundary conditions
    {
      cusetboundaryper_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_yp<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_zp<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_ym<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryper_zm<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }
  else if(boundary==3) // Mixed boundary conditions
    {
      cusetboundarytrans_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_yp  <<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_zp  <<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundarytrans_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_ym  <<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
      cusetboundaryref_zm  <<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
    }
  
}
#else

 if(neigh[5]!=rank)  
   {  
     for (int igrp=0;igrp<NGRP;igrp++) {
	cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
	mpisynch();

       exchange_zp(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[2]%2);
     }
     for (int igrp=0;igrp<NGRP;igrp++){
       cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
       mpisynch();

       exchange_zm(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[2]%2);
     }
   }
 else
   {
     for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_zp<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
     for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_zm<<<gridboundz,blockboundz>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
   }
  
 if(neigh[3]!=rank)
   {
     for (int igrp=0;igrp<NGRP;igrp++){
       cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
       mpisynch();

       exchange_yp(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[1]%2);
     }
     for (int igrp=0;igrp<NGRP;igrp++){
       cudaMemset(cubuff,0,sizeof(float)*NBUFF);
	memset(buff,0,NBUFF*sizeof(float));
      mpisynch();

       exchange_ym(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[1]%2);
     }
   }
 else
   {
     for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_yp<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
     for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_ym<<<gridboundy,blockboundy>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
   }

 if(neigh[1]!=rank)
   {
     for (int igrp=0;igrp<NGRP;igrp++){
       cudaMemset(cubuff,0,sizeof(float)*NBUFF);
       memset(buff,0,NBUFF*sizeof(float)); 
       mpisynch();

       exchange_xp(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[0]%2);
     }
     for (int igrp=0;igrp<NGRP;igrp++){
       cudaMemset(cubuff,0,sizeof(float)*NBUFF);
       memset(buff,0,NBUFF*sizeof(float));
       mpisynch();
       
       exchange_xm(cuegy+igrp*NCELLS3, cuflx+igrp*NCELLS3*3, cubuff, buff, neigh, pos[0]%2);
     }
   }
 else
   {
     for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
     for (int igrp=0;igrp<NGRP;igrp++) cusetboundaryper_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
   }

 if(boundary==0)
   {
     if(pos[0]==0)   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_xm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
     if(pos[1]==0)   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_ym<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
     if(pos[2]==0)   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_zm<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
	  
     if(pos[0]==(NGPUX-1))   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_xp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
     if(pos[1]==(NGPUY-1))   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_yp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3);
     if(pos[2]==(NGPUZ-1))   for (int igrp=0;igrp<NGRP;igrp++) cusetboundarytrans_zp<<<gridboundx,blockboundx>>>(cuegy+igrp*NCELLS3, cuxion, cudensity, cutemperature, cuflx+igrp*NCELLS3*3) ;
	  
   }
#endif

  cudaThreadSynchronize(); 
#ifdef WMPI
  mpisynch();
#endif

#ifdef TIMINGS
  if(rank==0)
    {
      q7=cutGetTimerValue(timer);
    }
#endif
  
  //printf("proc %d ready to dump\n",ic_rank);

#ifdef SYNCHDUMPFIELD
  if((((nstep%ndumps)==0)||(forcedump))||(changefield))
#else
  if(((nstep%ndumps)==0)||(forcedump))
#endif
	{
	  ntsteps=ntsteps+1;
	  forcedump=0;
#ifdef COSMO
#ifdef FLAT_COSMO
	  float aexpdump=t2a(t+dt,omegav,Hubble0);
#else
	  if(t+dt>t_int_max)
	    {
	      aexpdump=(a_int[int_step+2]-a_int[int_step+1])/(t_int[int_step+2]-t_int[int_step+1])*(t+dt-t_int[int_step+1]);
	    }
	  else
	    {
	      aexpdump=(a_int[int_step+1]-a_int[int_step])/(t_int[int_step+1]-t_int[int_step])*(t+dt-t_int[int_step]);
	    }
#endif
	  cuDumpResults(ntsteps,t+dt,aexpdump,ic_rank);
#else
	  cuDumpResults(ntsteps,t+dt,0.,ic_rank);
#endif
	  tnext=tnext+ndumps*dt/ft;
	  if(rank==0) printf("tnext=%f\n",tnext/unit_time);
#ifdef WMPI
       mpisynch();
#endif
	}

      //--------------------------------------------------------------------
      // Dealing with fieldlists
      //--------------------------------------------------------------------

      ft=fminf(ft*2.,1.);
      
      if(fieldlist)
	{
	  if(changefield)
	    {
	    int ercode;
#ifdef WMPI
	      ercode=cuGetField(ifield,ic_rank);
#else
	      ercode=cuGetField(ifield,0);
#endif
	      if(ercode==38)
		{
		  //fclose(logfile);
		  if(rank==0) fclose(timefile);
		  abort();
		}
	      forcedump=0;
	      changefield=0;
	      ifield++;
	      tfield=tlist[ifield];
	      ft=1./powf(2.,20);
#ifdef WMPI	  
         mpisynch();
#endif

	      //ft=1.;
	    }
	}


      // UPDATING VARIABLES

      t=t+dt;
      if(t>tmax)
	{
	  puts("t > tmax -----> run will be terminated");
	}
#ifdef COSMO

#ifdef FLAT_COSMO
      aexp=t2a(t,omegav,Hubble0); // A CHANGER PAR INTERPOLATION
#else
      if(t>t_int_max)
	{
	  int_step++;
	}
      aexp=(a_int[int_step+1]-a_int[int_step])/(t_int[int_step+1]-t_int[int_step])*(t-t_int[int_step]);
#endif


      c=c_r/aexp;
#endif       
      
      cudaThreadSynchronize();
#ifdef WMPI
      mpisynch();
#endif
      if(rank==0)
	{
	  q1=cutGetTimerValue(timer);
	}


      nstep++;
      if(nstep==nmax) {
	if(rank==0) puts("Max number of steps achieved: STOP");
	break;
      }

      cudaThreadSynchronize();
#ifdef WMPI
      get_elapsed(&time_new);
      time_new=time_new-time_old;
      mpireducemax(&time_new);
      mpisynch();
#endif

#ifdef TIMINGS
      if(rank==0){
	q9=cutGetTimerValue(timer);
	printf("transport=%f chem=%f cool=%f update=%f bound=%f IO=%f,grand total=%f time_new=%lf\n",q11-q0,q4-q11,q8-q4,q10-q8,q7-q10,q9-q7,q9-q0,time_new);
	fprintf(timefile,"%d %f %f %f %f %f %f %f\n",nstep-1,q11-q0,q4-q11,q8-q4,q10-q8,q7-q10,q9-q7,q9-q0,time_new);
      }


#endif

      cudaThreadSynchronize();
#ifdef WMPI	  
      mpisynch();
#endif

    }

  //fclose(logfile);
#ifdef TIMINGS
  if(rank==0) fclose(timefile);
#endif
  return 0;
}

