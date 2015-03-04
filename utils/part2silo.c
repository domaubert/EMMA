

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
  
  size_t dummy;

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
  float *mass;
  float *age;
  float ep,ek;
  
  int isstar=0;
  if(argc!=4){
    printf("USAGE: /a.out input nproc isstar\n");
    return 0;
  }

  // getting the number of CPUs
  sscanf(argv[2],"%d",&ncpu);

  // getting the sampling
  sscanf(argv[3],"%d",&isstar);


  // scanning the cpus
  npartout=0;
  for(icpu=0;icpu<ncpu;icpu++){
    
    // building the input file name
    strcpy(format,argv[1]);
    strcat(format,".p%05d");
    sprintf(fname,format,icpu);

    // building the output file name
    strcpy(format,argv[1]);

    sprintf(fname2,format,icpu);


    fp=fopen(fname,"rb");
  
    dummy=fread(&npart,sizeof(int),1,fp);
    dummy=fread(&tsimf,sizeof(float),1,fp);
    npartout+=npart;
    fclose(fp);
  }

  printf("Looking for a grand total of %d\n",npartout);
  printf("tsimf=%f\n",tsimf);

  // scanning the cpus
  int ip=0;
  int ipp;
  for(icpu=0;icpu<ncpu;icpu++){
    
    // building the input file name
    strcpy(format,argv[1]);
    strcat(format,".p%05d");
    sprintf(fname,format,icpu);

    // building the output file name
    strcpy(format,argv[1]);

    sprintf(fname2,format,icpu);


    fp=fopen(fname,"rb");
  
    dummy=fread(&npart,sizeof(int),1,fp);
    dummy=fread(&tsimf,sizeof(float),1,fp);
    
    if(icpu==0){
      x=(float *)malloc(npartout*sizeof(float));
      y=(float *)malloc(npartout*sizeof(float));
      z=(float *)malloc(npartout*sizeof(float));
      vx=(float *)malloc(npartout*sizeof(float));
      vy=(float *)malloc(npartout*sizeof(float));
      vz=(float *)malloc(npartout*sizeof(float));
      mass=(float *)malloc(npartout*sizeof(float));
      if(isstar) age=(float *)malloc(npartout*sizeof(float));
    }
    

    if(fp==NULL){
      printf("--ERROR -- file not found\n");
      return 1;
    }
    else{
      printf("Reading %d particles from %s\n",npart,fname);
    }
    
    for(ipp=0;ipp<npart;ipp++){

       dummy=fread(x+ip,sizeof(float),1,fp);
       dummy=fread(y+ip,sizeof(float),1,fp);
       dummy=fread(z+ip,sizeof(float),1,fp);
       dummy=fread(vx+ip,sizeof(float),1,fp);
       dummy=fread(vy+ip,sizeof(float),1,fp);
       dummy=fread(vz+ip,sizeof(float),1,fp);
       dummy=fread(&id,sizeof(int),1,fp);
       dummy=fread(mass+ip,sizeof(float),1,fp);
       dummy=fread(&ep,sizeof(int),1,fp);
       dummy=fread(&ek,sizeof(int),1,fp);
       if(isstar)dummy=fread(age+ip,sizeof(float),1,fp);
       ip++;
    }

    

    fclose(fp);
  }
    
    //============= dump
    
    DBfile *dbfile=NULL;
    dbfile=DBCreate(strcat(argv[1],".silo"),DB_CLOBBER, DB_LOCAL,"silo file created by EMMA",DB_PDB);
    npart=npartout;
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
    DBPutPointvar1(dbfile, "mass", "pointmesh", mass, npart, DB_FLOAT, NULL);
    if(isstar) DBPutPointvar1(dbfile, "age", "pointmesh", age, npart, DB_FLOAT, NULL);


    DBClose(dbfile);
    
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
    free(mass);
    if(isstar) free(age);
}
