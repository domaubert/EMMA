
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "silo.h"
//=================================================================================================
#define max(a,b) (a>=b?a:b)

int main(int argc, char *argv[])
{
  char fileout[256];
  char filein[256];
  size_t dummy;
  // building the output file name
  
  strcpy(filein,argv[1]);
  strcpy(fileout,filein);
  strcat(fileout,".silo");
  printf("cube2silo streams to  %s\n",fileout);

  // opening the input cube

  FILE *fp;
  fp=fopen(filein,"rb");
  int nmapx,nmapy,nmapz;
  dummy=fread(&nmapx,1,sizeof(int),fp);
  dummy=fread(&nmapy,1,sizeof(int),fp);
  dummy=fread(&nmapz,1,sizeof(int),fp);
  printf("allocating %dx%dx%d cells\n",nmapx,nmapy,nmapz);
  float *map;
  map=(float *)calloc(nmapx*nmapy*nmapz,sizeof(float));
  int I;
  for(I=0;I<nmapz;I++) dummy=fread(map+I*nmapx*nmapy,nmapx*nmapy,sizeof(float),fp);
  fclose(fp);
  printf("cube data ok\n");

  // dumping into Silo file
  
  printf("dumping silo in %s\n",fileout);
  DBfile *dbfile=NULL;
  dbfile=DBCreate(fileout,DB_CLOBBER, DB_LOCAL,"silo file created by EMMA",DB_PDB);
  if(dbfile==NULL){
    printf("Error while writing file");
  }

  float *x;
  float *y;
  float *z;
  // note we assume cubical data
  int dims[]={nmapx+1,nmapx+1,nmapx+1};
  int ndims=3;
  
  x=(float*)malloc(sizeof(float)*(nmapx+1));
  y=(float*)malloc(sizeof(float)*(nmapx+1));
  z=(float*)malloc(sizeof(float)*(nmapx+1));
  
  float *coords[]={x,y,z};
  
  int i;
  for(i=0;i<=nmapx;i++){
    x[i]=((i*1.0)/nmapx-0.);
    y[i]=((i*1.0)/nmapx-0.);
    z[i]=((i*1.0)/nmapx-0.);
  }
  
  int dimsvar[]={nmapx,nmapx,nmapx};
  int ndimsvar=3;
  
  
  DBPutQuadmesh(dbfile,"quadmesh",NULL,coords,dims,ndims,DB_FLOAT,DB_COLLINEAR,NULL);
  DBPutQuadvar1(dbfile,"monvar","quadmesh",map,dimsvar,ndimsvar,NULL,0,DB_FLOAT,DB_ZONECENT,NULL);

  DBClose(dbfile);
  free(map);
  printf("Done\n");
  return 0;
}
