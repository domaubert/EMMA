#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../src/prototypes.h"
#include "../src/hilbert.h"


int main(){

  FILE *fx = NULL;
  FILE *fy = NULL;
  FILE *fz = NULL;

  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
  int dummy;
  char filename[256];

  int level=7;

  
  printf("test\n");
  int ifile;

  for(ifile=0;ifile<7;ifile++){

    char field[256];


    switch(ifile){
    case 0:
      strcpy(field,"ic_deltab");
      break;
    case 1:
      strcpy(field,"ic_velcx");
      break;
    case 2:
      strcpy(field,"ic_velcy");
      break;
    case 3:
      strcpy(field,"ic_velcz");
      break;
    case 4:
      strcpy(field,"ic_velbx");
      break;
    case 5:
      strcpy(field,"ic_velby");
      break;
    case 6:
      strcpy(field,"ic_velbz");
      break;
    }


    
    sprintf(filename,"./level_%03d/%s",level,field);
    fx=fopen(filename,"rb");
    if(fx == NULL) {
    
      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/%s",level,field);
      printf("trying %s\n", filename);
    
      fx=fopen(filename,"rb");
      if(fx == NULL) {
	printf("Cannot open %s\n", filename);
	abort();
      }
    }
  
  
    printf("Doing %s\n",filename);
    // reading the headers
    size_t outf;

    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&np1,1,4,fx);
    outf=fread(&np2,1,4,fx);
    outf=fread(&np3,1,4,fx);
    outf=fread(&dx,1,4,fx);
    //printf("DX=%e\n",dx);
    outf=fread(&x1o,1,4,fx);
    outf=fread(&x2o,1,4,fx);
    outf=fread(&x3o,1,4,fx);
    outf=fread(&astart,1,4,fx);
    outf=fread(&om,1,4,fx);
    outf=fread(&ov,1,4,fx);
    outf=fread(&h0,1,4,fx);
    outf=fread(&dummy,1,sizeof(dummy),fx);


    float *velx;
    int icpu;
    int nproc=8;
    int *cpu_plan;
    int *cpu_open;
    int *cpu_plan_min_x;
    int *cpu_plan_max_x;
    int *cpu_plan_min_y;
    int *cpu_plan_max_y;
    int *cpu_n;

    velx=(float*)malloc(np1*np2*sizeof(float));
    cpu_plan=(int*)malloc(nproc*sizeof(int)); 
    cpu_n=(int*)calloc(nproc,sizeof(int));
    cpu_open=(int*)calloc(nproc,sizeof(int));
    cpu_plan_min_x=(int*)malloc(nproc*sizeof(int));
    cpu_plan_max_x=(int*)malloc(nproc*sizeof(int));
    cpu_plan_min_y=(int*)malloc(nproc*sizeof(int));
    cpu_plan_max_y=(int*)malloc(nproc*sizeof(int));


    unsigned long long keyrange=pow(2,3*(level-1))/nproc;
    unsigned long long key;


    FILE **fpout;
    fpout=(FILE **)calloc(nproc,sizeof(FILE*));
    



    // scanning the z-plane
    int pz;

    for(pz=0;pz<np3;pz++){
      
      outf=fread(&dummy,1,sizeof(dummy),fx);
      outf=fread(velx,np1*np2,sizeof(float),fx);
      outf=fread(&dummy,1,sizeof(dummy),fx);

      int py,px;
      int cpuloc;

      memset(cpu_plan,0,nproc*sizeof(int));
      memset(cpu_plan_min_x,0,nproc*sizeof(int));
      memset(cpu_plan_max_x,0,nproc*sizeof(int));
      memset(cpu_plan_min_y,0,nproc*sizeof(int));
      memset(cpu_plan_max_y,0,nproc*sizeof(int));

      for(py=0;py<np2;py+=1){
      for(px=0;px<np1;px+=1){
	  bitmask_t c[8*sizeof(bitmask_t)]; // integer coordinates of the oct
	  c[0]=px/2;
	  c[1]=py/2;
	  c[2]=pz/2;
	  key=(unsigned long long)(hilbert_c2i(3,level,c));
	    
	  cpuloc=key/keyrange;
	 
	  //	  if((px>(np1/2))&&(py<(np2/2))&&(pz>(np3/2))) printf("cpu=%d key=%ld \n",cpuloc,key);
   
	  if(cpu_plan[cpuloc]==0){
	    cpu_plan_min_x[cpuloc]=px;
	    cpu_plan_min_y[cpuloc]=py;
	    cpu_plan[cpuloc]=1;
	  }
	    
	  if(px>cpu_plan_max_x[cpuloc])cpu_plan_max_x[cpuloc]=px;
	  if(py>cpu_plan_max_y[cpuloc])cpu_plan_max_y[cpuloc]=py;
	    
	}
      }
      
      for(icpu=0;icpu<nproc;icpu++){
	if(cpu_plan[icpu]!=0){
	  int nxloc=(cpu_plan_max_x[icpu]-cpu_plan_min_x[icpu]+1);
	  int nyloc=(cpu_plan_max_y[icpu]-cpu_plan_min_y[icpu]+1);
	  int nzloc=np1*np2*np3/nproc/nxloc/nyloc;

	  cpu_n[icpu]+=nxloc*nyloc;

	  if(cpu_open[icpu]==0){
	    // file found -> open the file and write header
	    char fname[256];
	    
	    sprintf(fname,"level_%03d/%s.p%05d",level,field,icpu);
	    
	    fpout[icpu]=fopen(fname,"wb");


	    outf=fwrite(&dummy,1,sizeof(dummy),fpout[icpu]);
	    outf=fwrite(&nxloc,1,4,fpout[icpu]);
	    outf=fwrite(&nyloc,1,4,fpout[icpu]);
	    outf=fwrite(&nzloc,1,4,fpout[icpu]);
	    outf=fwrite(&dx,1,4,fpout[icpu]);

	    printf("nx=%d ny=%d nz=%d\n",nxloc,nyloc,nzloc);

	    float x1oloc,x2oloc,x3oloc;

	    x1oloc=cpu_plan_min_x[icpu]*dx;
	    x2oloc=cpu_plan_min_y[icpu]*dx;
	    x3oloc=pz*dx;

	    outf=fwrite(&x1oloc,1,4,fpout[icpu]);
	    outf=fwrite(&x2oloc,1,4,fpout[icpu]);
	    outf=fwrite(&x3oloc,1,4,fpout[icpu]);
	    outf=fwrite(&astart,1,4,fpout[icpu]);
	    outf=fwrite(&om,1,4,fpout[icpu]);
	    outf=fwrite(&ov,1,4,fpout[icpu]);
	    outf=fwrite(&h0,1,4,fpout[icpu]);
	    outf=fwrite(&dummy,1,sizeof(dummy),fpout[icpu]);

	    cpu_open[icpu]=1;
	  }

	  int ix,iy;
	  outf=fwrite(&dummy,1,sizeof(dummy),fpout[icpu]);
	  for(iy=cpu_plan_min_y[icpu];iy<=cpu_plan_max_y[icpu];iy++){
	    for(ix=cpu_plan_min_x[icpu];ix<=cpu_plan_max_x[icpu];ix++){
	      fwrite(&velx[ix+iy*np1],1,4,fpout[icpu]);
	    }
	  }
	  outf=fwrite(&dummy,1,sizeof(dummy),fpout[icpu]);

	  
	  
	}
      }

    }

    for(icpu=0;icpu<nproc;icpu++){
      printf("cpu=%d ncell=%d\n",icpu,cpu_n[icpu]);
      fclose(fpout[icpu]);
    }

    fclose(fx);

  }

  
  return 0;
}
