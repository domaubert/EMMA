#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../src/prototypes.h"


// utilisation
// ./utils/part2cic data/part.00035 7 32 0 0 0 0
// ./utils/part2cic fname level ncpu x0 y0 z0 isstar?

float ReadQuartz(float *map, char *inputfile, int nproc, int nmap, float x0, float y0, float z0, float dxmap, int isstar)
{
  int j,dummy;
  float pos[3], mass;
  int ifile;
  char currfile[256];
  char tempi[256];
  int nparttotal=0,npart;
  float time;
  float x;
  float y;
  float z;
  float dum;
  int ip=0;
  FILE *fp;


  /* for(ifile=0;ifile<nfile;ifile++){ */
  /*   strcpy(tempi,inputfile); */
  /*   strcat(tempi,".p%05d"); */
  /*   sprintf(currfile,tempi,ifile); */
  /*   //printf("%s\n",currfile); */
  /*   fp=fopen(currfile,"r"); */
  /*   fread(&npart,sizeof(int),1,fp); */
  /*   fclose(fp); */
  /*   nparttotal+=npart; */
  /* } */

  //kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
  for(ifile=0;ifile<nproc;ifile++){
    strcpy(tempi,inputfile);
    strcat(tempi,".p%05d");
    sprintf(currfile,tempi,ifile);
    fp=fopen(currfile,"r");
    printf("Reading %s\n",currfile);
    fread(&npart,sizeof(int),1,fp);
    fread(&time,sizeof(float),1,fp);

    for(j=0;j<npart;j++){
      fread(&x,sizeof(float),1,fp);
      fread(&y,sizeof(float),1,fp);
      fread(&z,sizeof(float),1,fp);
      fread(&dum,sizeof(float),7,fp);
      if(isstart)       fread(&dum,sizeof(float),1,fp);


      // CIC ASSIGNEMENT
      float dx,dy,dz;
      float tx,ty,tz;
      int xc,yc,zc;

      xc=(int)((x-x0)/dxmap);
      yc=(int)((y-y0)/dxmap);
      zc=(int)((z-z0)/dxmap);

      dx=x/dxmap-xc;
      dy=y/dxmap-yc;
      dz=z/dxmap-zc;

      tx=1.-dx;
      ty=1.-dy;
      tz=1.-dz;
      
      /* printf("%e %e %e\n",x,y,z); */
      /* printf("%d %d %d\n",xc,yc,zc); */
      /* printf("%e %e %e\n",x0,y0,z0); */
      int xp1,yp1,zp1;
      
      xp1=(xc+1)%nmap;
      yp1=(yc+1)%nmap;
      zp1=(zc+1)%nmap;

      map[xc+yc*nmap+zc*nmap*nmap]+=tx*ty*tz;

      map[(xp1)+yc*nmap+zc*nmap*nmap]+=dx*ty*tz;
      map[xc  +(yp1)*nmap+zc*nmap*nmap]+=tx*dy*tz;
      map[xc  +(yc  )*nmap+(zp1)*nmap*nmap]+=tx*ty*dz;

      map[xp1+(yp1)*nmap+zc*nmap*nmap]+=dx*dy*tz;
      map[xc  +(yp1)*nmap+(zp1)*nmap*nmap]+=tx*dy*dz;
      map[xp1+(yc  )*nmap+(zp1)*nmap*nmap]+=dx*ty*dz;

      map[xp1+(yp1)*nmap+(zp1)*nmap*nmap]+=dx*dy*dz;

      ip++;
    }

    fclose(fp);

  }

  return time;
}

int main(int argc, char *argv[]){
  
  int lmap;
  int nmap;
  float dxmap;

  float *map;
  char inputfile[256];
  char format[256];
  int ncpu;
  int isstar;

  //getting the resolution
  sscanf(argv[2],"%d",&lmap);
  nmap=pow(2,lmap);
  printf("nmap=%d\n",nmap);
  map=calloc(sizeof(float),nmap*nmap*nmap);

  dxmap=1./nmap;
  
  // getting the number of CPUs
  sscanf(argv[3],"%d",&ncpu);

  // building the output file name
  strcpy(inputfile,argv[1]);


  float x0=0.;
  float y0=0.;
  float z0=0.;

  // getting the number of CPUs
  sscanf(argv[4],"%f",&x0);
  // getting the number of CPUs
  sscanf(argv[5],"%f",&y0);
  // getting the number of CPUs
  sscanf(argv[6],"%f",&z0);

  
  float time=ReadQuartz(map,inputfile,ncpu,nmap,x0,y0,z0,dxmap,isstar);

  char fname2[256];
  FILE *fp;
  int I;
  int status;
  
  strcpy(fname2,inputfile);
  strcat(fname2,".cic");

  printf("dumping %s with nmap=%d hohoh\n",fname2,nmap*nmap*nmap);
  fp=fopen(fname2,"wb");
  fwrite(&nmap,1,sizeof(int),fp);
  fwrite(&nmap,1,sizeof(int),fp);
  fwrite(&nmap,1,sizeof(int),fp);
  fwrite(&time,1,sizeof(float),fp);
  for(I=0;I<nmap;I++) fwrite(map+I*nmap*nmap,nmap*nmap,sizeof(float),fp);
  status=ferror(fp);
  fclose(fp);
  

  return 0;
}


