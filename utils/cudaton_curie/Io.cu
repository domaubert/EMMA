#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"
#include "common.h"
#include "bnd.h"
#include "cosmo.h"
#include "GPU.h"
#include "Atomic.h"

extern "C" void getalist(int);
extern "C" void cuGetIC(int,int);

//*************************************************************************
//*************************************************************************

#define NCELLS3 (NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)

//*************************************************************************
//*************************************************************************
void getalist(int rank)
{
  FILE *fp;
  int i,idx;
  float dummy,dummy2;
  fp=fopen(filealist,"r");
  fscanf(fp,"%d",&nalist);
  //printf("# of outputs in the list=%d\n",nalist);
  alist=(float*)(calloc(nalist,sizeof(float)));
  tlist=(float*)(calloc(nalist,sizeof(float)));
  for(i=0;i<nalist;i++)
    {
      fscanf(fp,"%d %f %f",&idx,&dummy,&dummy2);
      alist[i]=dummy;
      if(dummy>0){
	tlist[i]=a2t(alist[i],omegav,Hubble0);
      }
      else{
	tlist[i]=0.;
      }
      if (rank==0) printf("(%f %f) ",alist[i],tlist[i]/unit_time);
    }
  //  printf("\n");
  //  printf("\n");
}


//*************************************************************************
//*************************************************************************
void cuDumpResults(int i, float time, float aexp, int rank)
{
  char fname[256];
  char fmt[256];
  FILE *fp;
  int nc=ncells+NBOUND2;
  float tt=time/unit_time;
  
#ifndef WMPI
  strcpy(fmt,rootname);
  strcat(fmt,".%05d");
  sprintf(fname,fmt,i);
#else
  strcpy(fmt,rootname);
  strcat(fmt,".%05d.p%05d");
  sprintf(fname,fmt,i,rank);
#endif
  
  if(rank==0) printf("Writing output in %s on proc %d\n",fname,rank);


  cudaMemcpy(egy,cuegy,NGRP*NCELLS3*sizeof(float),cudaMemcpyDeviceToHost);

#ifdef DUMPFLUX
  cudaMemcpy(flx,cuflx,NGRP*NCELLS3*sizeof(float)*3,cudaMemcpyDeviceToHost);
#endif

#ifndef SDISCRETE
  cudaMemcpy(src0,cusrc0,NCELLS3*sizeof(float),cudaMemcpyDeviceToHost);
#else
  cudaMemcpy(src0,cusrc0,nsource*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(src0pos,cusrc0pos,3*nsource*sizeof(int),cudaMemcpyDeviceToHost);
#endif

  cudaMemcpy(xion,cuxion,NCELLS3*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(density,cudensity,NCELLS3*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(temperature,cutemperature,NCELLS3*sizeof(float),cudaMemcpyDeviceToHost);
 
  fp=fopen(fname,"wb");

#ifndef DUMPGRID
	fwrite(&nc,sizeof(int),1,fp);
	fwrite(&nsource,sizeof(int),1,fp);
	fwrite(&tt,sizeof(float),1,fp);
	//fwrite(&ngrp,sizeof(int),1,fp);
	fwrite(egy,sizeof(float),NGRP*NCELLS3,fp);
	#ifdef DUMPFLUX
		fwrite(flx,sizeof(float),3*NGRP*NCELLS3,fp);
	#endif
	fwrite(xion,sizeof(float),NCELLS3,fp);
	fwrite(temperature,sizeof(float),NCELLS3,fp);
	#ifndef SDISCRETE
	//fwrite(src0,sizeof(float),NCELLS3,fp);
	#else
		fwrite(src0,sizeof(float),nsource,fp);
		fwrite(src0pos,sizeof(int),3*nsource,fp);
	#endif
		//fwrite(density,sizeof(float),NCELLS3,fp);
	fwrite(&aexp,sizeof(float),1,fp); // fin
#else
	fwrite(&nc,sizeof(int),1,fp);
	fwrite(&tt,sizeof(float),1,fp);
	fwrite(xion,sizeof(float),NCELLS3,fp);
	fwrite(temperature,sizeof(float),NCELLS3,fp);
#endif

  fclose(fp);
  if(rank==0) printf("Done on proc #%d\n",rank);

}

//*************************************************************************
//*************************************************************************

//*************************************************************************
//*************************************************************************

void cuGetIC(int iic,int rank)
{

  int i;
  int factgrp[NGRP];
  FACTGRP;


#ifndef TEST_STROMGREN
  char fname[256];
  char fmt[256];
  FILE *fp;
  int nc;
  int ns;

#ifndef WMPI
  if(iic!=0)
    {
      strcpy(fmt,rootname);
      strcat(fmt,".%05d");
      sprintf(fname,fmt,iic);
    }
  else
    {
      strcpy(fmt,rootname);
      strcat(fmt,".ic");
      sprintf(fname,fmt,iic);
    }
#else
  if(iic!=0)
    {
      strcpy(fmt,rootname);
      strcat(fmt,".%05d.p%05d");
      sprintf(fname,fmt,iic,rank);
    }
  else
    {
      strcpy(fmt,rootname);
      strcat(fmt,".ic.p%05d");
      sprintf(fname,fmt,rank);
    }
#endif
  
  if(rank==0) printf("Reading ICs in %s on proc #%d\n",fname,rank);
  
  fp=fopen(fname,"rb");
  if(fp==NULL) {
    printf("ERROR : IC file does not exist !\n");
    abort();
  }     

  fread(&nc,sizeof(int),1,fp);

  if(nc!=ncells+NBOUND2)
    {
      puts("ERROR while reading ICs : cell number is inconsistent !");
      printf("nc=%d ncells+NBOUND2=%d\n",nc,ncells+NBOUND2);
      abort();
    }
  
  fread(&ns,sizeof(int),1,fp);

#ifdef SDISCRETE
  if(ns!=nsource)
    {
      printf("ERROR while reading ICs : source number is inconsistent ns=%d nsource=%d!\n",ns,nsource);
      abort();
    }
#endif
  fread(&t,sizeof(float),1,fp);
  if(rank==0) printf("nc=%d ns=%d t=%e\n",nc,ns,t);
  if((iic!=0)&&(NGRP>1)){ 
    fread(egy,sizeof(float),NCELLS3*NGRP,fp);
  }
  else{
    fread(egy,sizeof(float),NCELLS3,fp);
  }

  if((iic!=0)&&(NGRP>1)){ 
    fread(flx,sizeof(float),3*NCELLS3*NGRP,fp);
  }
  else{
    fread(flx,sizeof(float),3*NCELLS3,fp);
  }

  fread(xion,sizeof(float),NCELLS3,fp);
  fread(temperature,sizeof(float),NCELLS3,fp);
  
#ifdef SDISCRETE
  fread(src0,sizeof(float),nsource,fp);
  fread(src0pos,sizeof(int),3*nsource,fp);
#else
  fread(src0,sizeof(float),NCELLS3,fp);
#endif
  fread(density,sizeof(float),NCELLS3,fp);
  fread(&astart,sizeof(float),1,fp);
  fclose(fp);
  if(rank==0) printf("astart=%e\n",astart);


#ifdef SAFE
  for(i=0;i<NCELLS3;i++){
    egy[i]+=1e-33;
    /* flx[i]=0.; */
    /* xion[i]=1e-4; */
    /* temperature[i]=1e2; */
    /* src0[i]*=1.;//7e10; */
    /* density[i]+=1e-2;//1e24+1e-6; */
  }
#endif

#else

  // STROMGREN TEST CASE
  if(rank==0)
    {
      printf("Self-consistent generation for Stromgren Sphere\n");
    }

  int 	ii,jj,kk,igrp;

  union
  {
    float f;
    unsigned char b[4];
  } dat1, dat2;

#define swap(X)  dat1.f = X;  dat2.b[0] = dat1.b[3];  dat2.b[1] = dat1.b[2];  dat2.b[2] = dat1.b[1];  dat2.b[3] = dat1.b[0];  X=dat2.f;



  float	Z,dummy,density_temp;

  FILE *rho;
  
  rho=fopen("../sources/density.bin","rb");
  fseek(rho,4,SEEK_CUR);
  fread(&Z,4,1,rho);
  swap(Z);
	fseek(rho,8,SEEK_CUR);
  fread(&dummy,sizeof(float),1,rho);
  swap(dummy);
	fseek(rho,4,SEEK_CUR);

  for(kk=0;kk<NCELLZ;kk++)
    {
	fseek(rho,4,SEEK_CUR);
      for(jj=0;jj<NCELLY;jj++)
	{	
	  for(ii=0;ii<NCELLX;ii++)
	    {
	      int idx=(ii+NBOUND)+(jj+NBOUND)*(NCELLX+2*NBOUND)+(kk+NBOUND)*(NCELLX+2*NBOUND)*(NCELLY+2*NBOUND);
		for (igrp=0;igrp<NGRP;igrp++)
			{
		      egy[idx+igrp*NCELLS3]=0.;
		      flx[idx+0*(NCELLX+2*NBOUND)*(NCELLY+2*NBOUND)*(NCELLZ+2*NBOUND)+igrp*NCELLS3*3]=0.;
		      flx[idx+1*(NCELLX+2*NBOUND)*(NCELLY+2*NBOUND)*(NCELLZ+2*NBOUND)+igrp*NCELLS3*3]=0.;
		      flx[idx+2*(NCELLX+2*NBOUND)*(NCELLY+2*NBOUND)*(NCELLZ+2*NBOUND)+igrp*NCELLS3*3]=0.;
			}
#ifndef COOLING
	      xion[idx]=1.2e-3;
	      temperature[idx]= 1e4; //100K 
#else 
	      xion[idx]=1e-6;
	      temperature[idx]= 1e2;
#endif

//		density[idx]=1000.;

/*		int size=8;
		if( ((ii>NCELLY*3./4-size) && (ii<NCELLY*3./4+size)) && ((jj>NCELLZ*3./4-size) && (jj<NCELLZ*3./4+size)) )
		     { density[idx]=2000.;}
		else { density[idx]=1000.;}	*/

		fread(&density_temp,sizeof(float),1,rho);
		swap(density_temp);
		density[idx]=density_temp*1e6;		
	    }
	}
	fseek(rho,4,SEEK_CUR);
    }
  fclose(rho);

  astart=1.;


#ifndef WMPI
/*
  src0[0]=5e48/dx/dx/dx/8.;
  src0pos[0]=0;
  src0pos[1]=0;
  src0pos[2]=0;

*/
int isrc,jsrc;
int src0pos[16*3];
float src0[16];

FILE * src;
src=fopen("../sources/sources.dat","r");

for (isrc=0;isrc<16;isrc++)
	{
	for (jsrc=0;jsrc<3;jsrc++)
		{
		fscanf(src,"%i",&src0pos[isrc+jsrc*16]);
		}
	fscanf(src,"%f",&src0[isrc]);
	src0[isrc] *= 1e52/dx/dx/dx;
	}
fclose(src);


#else
  if(rank==0)
    {
      src0[0]=0.;
      src0pos[0]=NCELLX/2;
      src0pos[1]=NCELLY/2;
      src0pos[2]=NCELLZ/2;
    }
  else
    {
      src0[0]=5e48/dx/dx/dx;
      src0pos[0]=0;
      src0pos[1]=0;
      src0pos[2]=0;

    }
#endif


  
  
#endif


#ifndef COSMO
  c=effective_speed_of_light*c_r;
#else
  c=effective_speed_of_light*c_r/astart;
  Hubble0=Hubble0/(9.7776e9*3.155815e7); // H0 in sec-1
#endif
  
#ifdef SDISCRETE
  if((nsource!=0)&&(rank==0))
    {
      printf("%d sources found\n",nsource);
#ifndef WMPI
      for(i=0;i<nsource;i++) printf(" %d %d %d %e\n",src0pos[i],src0pos[i+nsource],src0pos[i+2*nsource],src0[i]);
#endif
      printf("tstart= %e\n",t);
    }
#endif

  
  if(iic==0){
    for(i=0;i<NCELLS3;i++){
      if(density[i]<0){
	density[i]=defdens;
	temperature[i]=deftemp;
      }
   
#ifdef FORCET
      temperature[i]=deftemp;
#endif
 
      for(int j=NGRP-1;j>=0;j--) egy[i+j*NCELLS3]=fmaxf(egy_min*factgrp[j],egy[i+j*NCELLS3]);
    }  
  }

  cudaMemcpy(cuegy,egy,NCELLS3*sizeof(float)*NGRP,cudaMemcpyHostToDevice);
  cudaMemcpy(cuflx,flx,NCELLS3*sizeof(float)*3*NGRP,cudaMemcpyHostToDevice);

  
  int odx=15684351;
  int ii,jj,kk;

  kk=odx/((NCELLX+NBOUND2)*(NCELLY+NBOUND2));
  jj=(odx-kk*(NCELLX+NBOUND2)*(NCELLY+NBOUND2))/(NCELLX+NBOUND2);
  ii=odx-kk*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)-jj*(NCELLX+NBOUND2);

  
  //  printf("Rank=%d flx =%e %e %e egy=%e i=%d j=%d k=%d\n",rank,flx[odx-1],flx[odx],flx[odx+1],egy[odx],ii-NBOUND,jj-NBOUND,kk-NBOUND);
  
  cudaMemcpy(cuxion,xion,NCELLS3*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(cudensity,density,NCELLS3*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(cutemperature,temperature,NCELLS3*sizeof(float),cudaMemcpyHostToDevice);
  
#ifdef SDISCRETE
  cudaMemcpy(cusrc0,src0,nsource*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(cusrc0pos,src0pos,3*nsource*sizeof(int),cudaMemcpyHostToDevice);
#else
  cudaMemcpy(cusrc0,src0,NCELLS3*sizeof(float),cudaMemcpyHostToDevice);
#endif      
  
  if(rank==0)  printf("in read astart=%e\n",astart);
  if(rank==0) printf("Device Memory allocated on proc #%d\n",rank);
}


//*************************************************************************

int cuGetField(int iic, int rank)
{
  char fname[256];
  char fmt[256];
  FILE *fp;
  int nc;
  int ns;
  float tloc;

  strcpy(fmt,fieldname);
#ifdef WMPI
  strcat(fmt,"_%03d.p%05d");
  sprintf(fname,fmt,iic,rank);
#else
  strcat(fmt,"_%03d");
  sprintf(fname,fmt,iic);
#endif
      
  if(rank==0) printf("Reading Field in %s\n",fname);
  
  fp=fopen(fname,"rb");
  if(fp==NULL) {
    printf("ERROR : IC file does not exist !\n");
    return 38;
  }     
  
  fread(&nc,sizeof(int),1,fp);
  if(nc!=ncells+NBOUND2)
    {
      puts("ERROR while reading Field : cell number is inconsistent !");
      abort();
    return 38;
    }
  
  fread(&ns,sizeof(int),1,fp);
  if(ns!=nsource)
    {
      puts("ERROR while reading Field : source number is inconsistent !");
      abort();
    return 38;
    }

#ifndef LIGHTFIELD
  // regular snapshot format for input fields 
  fread(&tloc,sizeof(float),1,fp);
  fread(egy,sizeof(float),NCELLS3,fp);
  fread(flx,sizeof(float),3*NCELLS3,fp);
  fread(xion,sizeof(float),NCELLS3,fp);
  fread(temperature,sizeof(float),NCELLS3,fp);
  
#ifdef SDISCRETE
  fread(src0,sizeof(float),nsource,fp);
#ifdef RAND_SRC
  srand(rank);
  for(int i=0;i<nsource;i++)
    {
      src0[i]*=(float)(rand())/(float)(RAND_MAX);
    }
#endif
  fread(src0pos,sizeof(int),3*nsource,fp);
#else
  fread(src0,sizeof(float),NCELLS3,fp);
#endif
  fread(density,sizeof(float),NCELLS3,fp);
  fclose(fp);
#else
  
  // light format for input fields
  fread(&tloc,sizeof(float),1,fp);
#ifdef SDISCRETE
  fread(src0,sizeof(float),nsource,fp);
  fread(src0pos,sizeof(int),3*nsource,fp);
#else
  fread(src0,sizeof(float),NCELLS3,fp);
#endif
  fread(density,sizeof(float),NCELLS3,fp);
  fclose(fp);

#endif      

  for(int i=0;i<NCELLS3;i++){
    if(density[i]<0){
      density[i]=defdens;
    }
  }  

  // sending data to GPU
  cudaMemcpy(cudensity,density,NCELLS3*sizeof(float),cudaMemcpyHostToDevice);

#ifdef SDISCRETE
  cudaMemcpy(cusrc0,src0,nsource*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(cusrc0pos,src0pos,3*nsource*sizeof(int),cudaMemcpyHostToDevice);
#else
  cudaMemcpy(cusrc0,src0,NCELLS3*sizeof(float),cudaMemcpyHostToDevice);
#endif      
      
  if(rank==0) puts("Fields updated");
  return 0;
}

//==========================================================
//==========================================================
