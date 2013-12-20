

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "silo.h"
#include "prototypes.h"
//=================================================================================================



/*

 mpicc  -I./silo/include -DPIC -DWGRAV -DTESTCOSMO -DFASTGRAV -DGRAFIC -DONFLYRED -DGPUAXL -o ./part2silo ./part2silo.c ./silo/lib/libsilo.a

 */
int main(int argc, char *argv[])
{
  //(int lmap,struct OCT **firstoct,int field,char filename[],REAL tsim)
  int status;
  int I;
  FILE *fp;
  REAL tsim=0.;
  float tsimf;
  char fname[256];
  char format[256];
  char fname2[256];
  int ncpu;
  int icpu;

  float *x;
  float *y;
  float *z;
  float *vx;
  float *vy;
  float *vz;
  int npart;
  int npartout;
  int samp;

  float xp,yp,zp;
  int id;
  int ip;

  if(argc!=4){
    printf("USAGE: /a.out input nproc sampling\n");
    return 0;
  }

  // getting the number of CPUs
  sscanf(argv[2],"%d",&ncpu);

  // getting the sampling
  sscanf(argv[3],"%d",&samp);


  // scanning the cpus
  for(icpu=0;icpu<ncpu;icpu++){
    
    // building the input file name
    strcpy(format,argv[1]);
    strcat(format,".p%05d");
    sprintf(fname,format,icpu);

    // building the output file name
    strcpy(format,argv[1]);

    sprintf(fname2,format,icpu);


    fp=fopen(fname,"rb");
  
    fread(&npart,sizeof(int),1,fp);
    fread(&tsimf,sizeof(float),1,fp);
    
    if(icpu==0){
      npartout=npart/samp;
      x=(float *)malloc(npartout*sizeof(float));
      y=(float *)malloc(npartout*sizeof(float));
      z=(float *)malloc(npartout*sizeof(float));
      vx=(float *)malloc(npartout*sizeof(float));
      vy=(float *)malloc(npartout*sizeof(float));
      vz=(float *)malloc(npartout*sizeof(float));
    }
    

    if(fp==NULL){
      printf("--ERROR -- file not found\n");
      return 1;
    }
    else{
      printf("Reading %d particles from %s\n",npart,fname);
    }
    
    int ipout=0;
    for(ip=0;ip<npart;ip++){

      fread(x+ip,sizeof(float),1,fp);
      fread(y+ip,sizeof(float),1,fp);
      fread(z+ip,sizeof(float),1,fp);
      fread(vx+ip,sizeof(float),1,fp);
      fread(vy+ip,sizeof(float),1,fp);
      fread(vz+ip,sizeof(float),1,fp);
      fread(&id,sizeof(int),1,fp);
      
      x[ip]=tsimf*(x[ip]-0.5);
      y[ip]=tsimf*(y[ip]-0.5);
      z[ip]=tsimf*(z[ip]-0.5);
    }

    printf("tsimf=%f\n",tsimf);
    

    fclose(fp);
    
    //============= dump
    
    DBfile *dbfile=NULL;
    dbfile=DBCreate(strcat(fname2,".silo"),DB_CLOBBER, DB_LOCAL,"silo file created by Quartz",DB_PDB);
      
    if(dbfile==NULL){
      printf("Error while writing file");
    }

      
    float *coords[] = {(float*)x, (float*)y, (float*)z};
    int ndims=3;
    /* Write a point mesh. */
    DBPutPointmesh(dbfile, "pointmesh", ndims, coords, npart, DB_FLOAT, NULL);
    DBPutPointvar1(dbfile, "vx", "pointmesh", vx, npart, DB_FLOAT, NULL);
    DBPutPointvar1(dbfile, "vy", "pointmesh", vy, npart, DB_FLOAT, NULL);
    DBPutPointvar1(dbfile, "vz", "pointmesh", vz, npart, DB_FLOAT, NULL);


    DBClose(dbfile);
    
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
  }
}
