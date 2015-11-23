#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

int main(int argc, char *argv[]){


// DOCUMENTATION
  if(argc!=10){
    printf("USAGE : ./a.out snapdirectory fielddirectory xc yc zc dx dy dz lmap\n");
    printf("will produce a cube of resolution lmap centered around [xc,yc,zc]\n");
    printf("with side length 2x [dx,dy,dz]\n");
    printf("EXAMPLE :\n");
    printf("./a.out data/00013 grid_xion 0.5 0.5 0.5 0.2 0.2 0.2 10\n");
    printf("produces a cube \n");
    return 0;

  }

  // PARAMETERS
  char repname[256]="data/00013/";
  char ffname[256]="grid_xion/";
  float xc=0.5,yc=0.5,zc=0.5;
  float Dx=0.2,Dy=0.2,Dz=0.2;
  int lmap=10;


  strcpy(repname,argv[1]);
  strcat(repname,"/");
  strcpy(ffname,argv[2]);
  strcat(ffname,"/");
  sscanf(argv[3],"%f",&xc);
  sscanf(argv[4],"%f",&yc);
  sscanf(argv[5],"%f",&zc);
  sscanf(argv[6],"%f",&Dx);
  sscanf(argv[7],"%f",&Dy);
  sscanf(argv[8],"%f",&Dz);
  sscanf(argv[9],"%d",&lmap);

  // TEMP VARIABLES
  char xfname[256];
  char yfname[256];
  char zfname[256];
  char lfname[256];
  char fname[256];
  float dmap=pow(0.5,lmap);

  // ALLOC GRID
  float xmin,xmax;
  float ymin,ymax;
  float zmin,zmax;

  xmin=xc-Dx;
  ymin=yc-Dy;
  zmin=zc-Dz;

  xmax=xc+Dx;
  ymax=yc+Dy;
  zmax=zc+Dz;

  // integerize

  xmin=(int)(xmin/dmap)*dmap;
  xmax=(int)(xmax/dmap)*dmap;
  int nx=(xmax-xmin)/dmap;
  
  ymin=(int)(ymin/dmap)*dmap;
  ymax=(int)(ymax/dmap)*dmap;
  int ny=(ymax-ymin)/dmap;

  zmin=(int)(zmin/dmap)*dmap;
  zmax=(int)(zmax/dmap)*dmap;
  int nz=(zmax-zmin)/dmap;



  printf("Casting in rays in %s\n",repname);
  printf("xmin=%e xmax=%e nx=%d\n",xmin,xmax,nx);
  printf("ymin=%e ymax=%e ny=%d\n",ymin,ymax,ny);
  printf("zmin=%e zmax=%e nz=%d\n",zmin,zmax,nz);
  
  // allocation of the grid

  float *grid;
  grid=(float *)calloc(nx*ny*nz,sizeof(float));

  // Preparing File names

  strcpy(xfname,repname);
  strcpy(yfname,repname);
  strcpy(zfname,repname);
  strcpy(lfname,repname);
  strcpy(fname,repname);

  strcat(xfname,"/grid_x/");
  strcat(yfname,"/grid_y/");
  strcat(zfname,"/grid_z/");
  strcat(lfname,"/grid_l/");
  strcat(fname,ffname);


  FILE *fp;
  FILE *xfp;
  FILE *yfp;
  FILE *zfp;
  FILE *lfp;

  // First we scan the directory to detect the number of processors
  int nproc;
  struct dirent **namelist;
  nproc=scandir(lfname,&namelist,NULL,alphasort);
  // extract the snap number
  int snapnum;
  sscanf(namelist[2]->d_name,"l.%d.p00000",&snapnum);

  // additional scan to switch onto the field files
  nproc=scandir(fname,&namelist,NULL,alphasort);


  printf("Found %d files in %s for snap number %d\n",nproc-2,fname,snapnum); // -2 to remove '.' and '..'


  // scanning files
  int i;
  for(i=2;i<nproc;i++){
    FILE *flocf;
    char fnameloc[256];

    // build local fname
    strcpy(fnameloc,fname);
    strcat(fnameloc,namelist[i]->d_name);

    // extract domain
    //printf("%s\n",fnameloc);
    flocf=fopen(fnameloc,"rb");
    int ncell;
    float bounds[6];

    fread(&ncell,sizeof(int),1,flocf);
    fread(bounds,sizeof(float),6,flocf);
    
    // should we consider this subfile ?
    int nox=(bounds[0]>xmax)||(bounds[1]<xmin);
    int noy=(bounds[2]>ymax)||(bounds[3]<ymin);
    int noz=(bounds[4]>zmax)||(bounds[5]<zmin);
    
    if((nox||noy)||noz){
      //printf("SKIP ");
    }
    else{
      // we are good
      
      // let us open the FILES

      char format[256];
      float bounds[6];

      FILE *flocx;
      strcpy(format,xfname);
      strcat(format,"x.%05d.p%05d");
      sprintf(fnameloc,format,snapnum,i-2);
      flocx=fopen(fnameloc,"rb");
      fread(&ncell,sizeof(int),1,flocx);
      fread(bounds,sizeof(float),6,flocx);


      FILE *flocy;
      strcpy(format,yfname);
      strcat(format,"y.%05d.p%05d");
      sprintf(fnameloc,format,snapnum,i-2);
      flocy=fopen(fnameloc,"rb");
      fread(&ncell,sizeof(int),1,flocy);
      fread(bounds,sizeof(float),6,flocy);

      FILE *flocz;
      strcpy(format,zfname);
      strcat(format,"z.%05d.p%05d");
      sprintf(fnameloc,format,snapnum,i-2);
      flocz=fopen(fnameloc,"rb");
      fread(&ncell,sizeof(int),1,flocz);
      fread(bounds,sizeof(float),6,flocz);


      FILE *flocl;
      strcpy(format,lfname);
      strcat(format,"l.%05d.p%05d");
      sprintf(fnameloc,format,snapnum,i-2);
      flocl=fopen(fnameloc,"rb");
      fread(&ncell,sizeof(int),1,flocl);
      fread(bounds,sizeof(float),6,flocl);
      

      float *x;
      float *y;
      float *z;
      float *l;
      float *f;

      
      x=(float *)calloc(ncell,sizeof(float));
      y=(float *)calloc(ncell,sizeof(float));
      z=(float *)calloc(ncell,sizeof(float));
      l=(float *)calloc(ncell,sizeof(float));
      f=(float *)calloc(ncell,sizeof(float));


      fread(x,sizeof(float),ncell,flocx);
      fread(y,sizeof(float),ncell,flocy);
      fread(z,sizeof(float),ncell,flocz);
      fread(l,sizeof(float),ncell,flocl);
      fread(f,sizeof(float),ncell,flocf);


      fclose(flocx);
      fclose(flocy);
      fclose(flocz);
      fclose(flocl);

      // DATA IS AVAILABLE, LET'S CAST
      
      int pcell;
      int icell,jcell,kcell;
      for(pcell=0;pcell<ncell;pcell++){
	
	       // computing local indexes
	       icell=(x[pcell]-xmin)/dmap;
	       jcell=(y[pcell]-ymin)/dmap;
	       kcell=(z[pcell]-zmin)/dmap;
	
	         // is the current cell within the domain ?
	       int condx=(icell<0)||(icell>=nx);
	       int condy=(jcell<0)||(jcell>=ny);
	       int condz=(kcell<0)||(kcell>=nz);

	       if((condx||condy)||condz){
	       }
	       else{
           
           if(l[pcell]>=lmap){
            grid[icell+jcell*nx+kcell*nx*ny]+=f[pcell]*pow(0.5,3*(l[pcell]-lmap));
           }
           else{
            int iimax=pow(2,lmap-l[pcell]);
            int ii,jj,kk;
            for(kk=kcell;kk<kcell+iimax;kk++){
              if(kk>=nz) continue;
              for(jj=jcell;jj<jcell+iimax;jj++){
                if(jj>=ny) continue;
                for(ii=icell;ii<icell+iimax;ii++){
                  if(ii>=nx) continue;
                  grid[ii+jj*nx+kk*nx*ny]+=f[pcell];
                }
              }
            }
          }  


	         
	       }
      }

      // FREE DATA
      free(x);
      free(y);
      free(z);
      free(l);
      free(f);


    }

    //    printf("ncells=%d bounds= %e %e | %e %e | %e %e\n",ncell,bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
    
    

    fclose(flocf);

  }


  // DUMPING RESULTS
  FILE *dumpfile;
  dumpfile=fopen("dump.dat","wb");
  fwrite(&nx,sizeof(int),1,dumpfile);
  fwrite(&ny,sizeof(int),1,dumpfile);
  fwrite(&nz,sizeof(int),1,dumpfile);
  fwrite(&xmin,sizeof(float),1,dumpfile);
  fwrite(&xmax,sizeof(float),1,dumpfile);
  fwrite(&ymin,sizeof(float),1,dumpfile);
  fwrite(&ymax,sizeof(float),1,dumpfile);
  fwrite(&zmin,sizeof(float),1,dumpfile);
  fwrite(&zmax,sizeof(float),1,dumpfile);
  fwrite(grid,sizeof(float),nx*ny*nz,dumpfile);
  fclose(dumpfile);

  return 0;

}


  



