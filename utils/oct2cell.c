

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "silo.h"
#include "prototypes.h"
//=================================================================================================

#define NREAL 5

/*

 mpicc  -I./silo/include -DPIC -DWGRAV -DTESTCOSMO -DFASTGRAV -DGRAFIC -DONFLYRED -DGPUAXL -o ./oct2grid ./oct2grid.c ./silo/lib/libsilo.a

avconv -i mpeg_link%04d.jpeg -r 24 -b 65536k video.mp4

 */
int main(int argc, char *argv[])
{
  //(int lmap,struct OCT **firstoct,int field,char filename[],REAL tsim)
  int status;
  int I;
  int lmap;
  float *map;
  int imap,jmap,kmap;
  int nmap;
  REAL dxmap,dxcur;
  int icur,ii,jj,kk,ic,icell;
  int level;
  struct LOCT * nextoct;
  struct LOCT oct;
  FILE *fp;
  REAL xc,yc,zc;
  int field;
  REAL tsim=0.;
  float tsimf;
  char fname[256];
  char format[256];
  char fname2[256];
  int ncpu;
  long int nc;
  int icpu;
  struct OCT *zerooct;
#ifdef WRAD
  struct UNITS unit;
#endif

  size_t dummy;

  if(argc<7){
    printf("USAGE: /a.out input level field output nproc\n");
    printf("Ex : data/grid.00002 7 101 data/grid.den 8\n");
    printf("field values :\n");
    printf(" case 0 oct.level;\n case 1 density;\n case 2 pot;\n case 3 cpu;\n case 4 marked;\n case 5 temp;\n case 101 field.d;\n case 102 field.u;\n case 103 field.v;\n case 104 field.w;\n case 105 field.p;\n");
    return 0;
  }

  // getting the field
  sscanf(argv[3],"%d",&field);

  //getting the resolution
  sscanf(argv[2],"%d",&lmap);
  nmap=pow(2,lmap);
  dxmap=1./nmap;

  // getting the number of CPUs
  sscanf(argv[5],"%d",&ncpu);
  printf("ncpu=%d\n",ncpu);

  // silo file
  int dumpsilo=0;
  sscanf(argv[6],"%d",&dumpsilo);
  REAL zmin,zmax;
  int nmapz;
  int k0;

  int mono;
  sscanf(argv[7],"%d",&mono);
  printf("mono=%d\n",mono);

  if(argc>8) {
    sscanf(argv[8],"%f",&zmin);
    sscanf(argv[9],"%f",&zmax);
  }
  else{
    zmin=0;
    zmax=1.;
  }

  nmapz=(zmax-zmin)/dxmap;
  k0=zmin/dxmap;

  if((dumpsilo==1)&&(nmapz!=nmap)) {
    printf("Silo can only handle cubic data !\n");
    abort();
  }

  // scanning the cpus

  int icpustart,icpustop;
  if(mono<0){
    icpustart=0;
    icpustop=ncpu;
  }
  else{
    icpustart=mono;
    icpustop=mono+1;
    printf("Cast single cpu #%d\n",mono) ;
  }


  nc=0;
  for(icpu=icpustart;icpu<icpustop;icpu++){


    // building the input file name
    strcpy(format,argv[1]);
    strcat(format,".p%05d");
    sprintf(fname,format,icpu);

    // building the output file name
    strcpy(format,argv[4]);

    sprintf(fname2,format,icpu);

    // counting cells

    fp=fopen(fname,"rb");

    if(fp==NULL){
      printf("--ERROR -- file not found\n");
      return 1;
    }

    // reading the time
    dummy=fread(&tsim,sizeof(REAL),1,fp);
#ifdef WRAD
    dummy=fread(&unit.unit_l,sizeof(REAL),1,fp);
    printf("unit=%e\n",unit.unit_l);
    dummy=fread(&unit.unit_v,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_t,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_n,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_mass,sizeof(REAL),1,fp);
#endif

    // reading the zerooct
    dummy=fread(&zerooct,sizeof(struct OCT *),1,fp);


    dummy=fread(&oct,sizeof(struct LOCT),1,fp);
    while(!feof(fp)){
      for(icell=0;icell<8;icell++) // looping over cells in oct
	{
	  if(oct.cell[icell].child==0)
	    {
	      if((zc<zmin)||(zc>zmax)){
		continue;
	      }

	      int flag;
	      if(mono>=0){
		flag=(oct.cpu>=0);
	      }
	      else{
		flag=(oct.cpu==icpu);
	      }

	      if(flag) {
		//printf("helo\n,");
		nc++;
		//printf("nc=%d\n",nc);
	      }
	    }
	}
      dummy=fread(&oct,sizeof(struct OCT),1,fp); //reading next oct
    }
    fclose(fp);
  }
  printf("Found %ld leaf cells\n",nc);

  map=(float *)calloc(nc*NREAL,sizeof(float)); // NREAL=5 3 for pos, 1 for level, 1 for data

  ic=0;
  for(icpu=icpustart;icpu<icpustop;icpu++){


    // building the input file name
    strcpy(format,argv[1]);
    strcat(format,".p%05d");
    sprintf(fname,format,icpu);

    // building the output file name
    strcpy(format,argv[4]);

    sprintf(fname2,format,icpu);

    // getting data
    printf("Looking for data in %s\n",fname);
    fp=fopen(fname,"rb");

    if(fp==NULL){
      printf("--ERROR -- file not found\n");
      return 1;
    }


    // reading the time
    dummy=fread(&tsim,sizeof(REAL),1,fp);
#ifdef WRAD
    dummy=fread(&unit.unit_l,sizeof(REAL),1,fp);
    printf("unit=%e\n",unit.unit_l);
    dummy=fread(&unit.unit_v,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_t,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_n,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_mass,sizeof(REAL),1,fp);
#endif

    printf("tsim=%e\n",tsim);
    tsimf=tsim;
    // reading the zerooct
    dummy=fread(&zerooct,sizeof(struct OCT *),1,fp);
    //printf("0oct=%p\n",zerooct);

    //printf("size of the OCT=%ld\n",sizeof(struct OCT));


    dummy=fread(&oct,sizeof(struct OCT),1,fp);
    while(!feof(fp)){
      dxcur=1./pow(2.,oct.level);
      for(icell=0;icell<8;icell++) // looping over cells in oct
	{
	  if(oct.cell[icell].child==0)
	    {
	      xc=oct.x+( icell   %2)*dxcur+0.5*dxcur;
	      yc=oct.y+((icell/2)%2)*dxcur+0.5*dxcur;
	      zc=oct.z+( icell/4   )*dxcur+0.5*dxcur;
	      if((zc<zmin)||(zc>zmax)){
		continue;
	      }

	      int flag;
	      if(mono>=0){
		flag=(oct.cpu>=0);
	      }
	      else{
		flag=(oct.cpu==icpu);
	      }
	      if(flag) {
		map[ic*NREAL+0]=xc;
		map[ic*NREAL+1]=yc;
		map[ic*NREAL+2]=zc;
		map[ic*NREAL+3]=(float)oct.level;
		switch(field){
		case 0:
		  map[ic*NREAL+4]=oct.level;
		  break;
#ifdef WGRAV
		case 1:
		  map[ic*NREAL+4]=oct.cell[icell].den;
		  break;
#ifdef PIC
		case -1:
		  map[ic*NREAL+4]=oct.cell[icell].density;
		  break;
#endif
		case 2:
		  map[ic*NREAL+4]=oct.cell[icell].pot;
		  break;
		case 6:
		  map[ic*NREAL+4]+=oct.cell[icell].res;
		  break;
		case 201:
		  map[ic*NREAL+4]=oct.cell[icell].f[0];
		  break;
		case 202:
		  map[ic*NREAL+4]=oct.cell[icell].f[1];
		  break;
		case 203:
		  map[ic*NREAL+4]=oct.cell[icell].f[2];
		  break;
#endif
		case 3:
		  map[ic*NREAL+4]=oct.cpu;
		  break;
		case 4:
		  map[ic*NREAL+4]=oct.cell[icell].marked;
		  break;
#ifdef WHYDRO2
		case 101:
		  map[ic*NREAL+4]=oct.cell[icell].d;
		  //printf("%f\n",oct.cell[icell].d);
		  break;
		case 102:
		  map[ic*NREAL+4]=oct.cell[icell].u;
		  break;
		case 103:
		  map[ic*NREAL+4]=oct.cell[icell].v;
		  break;
		case 104:
		  map[ic*NREAL+4]=oct.cell[icell].w;
		  break;
		case 105:
		  map[ic*NREAL+4]=oct.cell[icell].p;
		  break;
#endif

#ifdef WRAD
		case 701:
		  map[ic*NREAL+4]=log10(oct.cell[icell].e[0]);
		  break;
		case 702:
		  map[ic*NREAL+4]=oct.cell[icell].fx[0];
		  break;
		case 703:
		  map[ic*NREAL+4]=oct.cell[icell].fy[0];
		  break;
		case 704:
		  map[ic*NREAL+4]=oct.cell[icell].fz[0];
		  break;
		case 705:
		  map[ic*NREAL+4]=oct.cell[icell].src;
		  break;
#ifdef WCHEM
		case 706:
		  map[ic*NREAL+4]=oct.cell[icell].xion;
		  break;
		case 707:
		  map[ic*NREAL+4]=oct.cell[icell].temp;
		  break;
		case 708:
		  map[ic*NREAL+4]=oct.cell[icell].d/pow(unit.unit_l,3.);
		  break;
#endif
#endif

		}
		ic++;
	      }
	    }
	}
      dummy=fread(&oct,sizeof(struct OCT),1,fp); //reading next oct
    }
    fclose(fp);
  }
  printf("done with %d cells\n",ic);

  //============= dump

    if(!dumpsilo){
      printf("dumping %s with nleaf=%d\n",fname2,nc);
      fp=fopen(fname2,"wb");
      fwrite(&nc,1,sizeof(int),fp);
      fwrite(&tsimf,1,sizeof(float),fp);
      for(I=0;I<nc;I++) fwrite(map+I*NREAL,sizeof(float),NREAL,fp);
      status=ferror(fp);
      fclose(fp);
    }
    else{
       DBfile *dbfile=NULL;
      dbfile=DBCreate(strcat(fname2,".silo"),DB_CLOBBER, DB_LOCAL,"silo file created by Quartz",DB_PDB);
      if(dbfile==NULL){
	printf("Error while writing file");
      }


      float *x;
      float *y;
      float *z;
      int dims[]={nmap+1,nmap+1,nmap+1};
      int ndims=3;

      x=(float*)malloc(sizeof(float)*(nmap+1));
      y=(float*)malloc(sizeof(float)*(nmap+1));
      z=(float*)malloc(sizeof(float)*(nmap+1));

      float *coords[]={x,y,z};

      int i;
      for(i=0;i<=nmap;i++){
	x[i]=((i*1.0)/nmap-0.);
	y[i]=((i*1.0)/nmap-0.);
	z[i]=((i*1.0)/nmap-0.);
      }

      int dimsvar[]={nmap,nmap,nmap};
      int ndimsvar=3;

      DBPutQuadmesh(dbfile,"quadmesh",NULL,coords,dims,ndims,DB_FLOAT,DB_COLLINEAR,NULL);
      DBPutQuadvar1(dbfile,"monvar","quadmesh",map,dimsvar,ndimsvar,NULL,0,DB_FLOAT,DB_ZONECENT,NULL);

      DBClose(dbfile);
    }
  printf("status=%d\n",status);
  free(map);

}
