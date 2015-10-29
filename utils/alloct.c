

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "silo.h"
#include "../src/prototypes.h"
//=================================================================================================

/*

 mpicc  -I./silo/include -DPIC -DWGRAV -DTESTCOSMO -DFASTGRAV -DGRAFIC -DONFLYRED -DGPUAXL -o ./oct2grid ./oct2grid.c ./silo/lib/libsilo.a

avconv -i mpeg_link%04d.jpeg -r 24 -b 65536k video.mp4

 */

struct cdata{
  float x;
  float y;
  float z;
  float level;
  float data;
};


float assign_cube(int field, int icell,struct LOCT *oct, struct UNITS *unit);


int main(int argc, char *argv[])
{
  int status;
  int I;
  int lmap;
  //REAL *map;
  float data;
  REAL dxmap,dxcur;
  int icur,ii,jj,kk,ic,icell;
  int level;
  struct LOCT * nextoct;
  struct LOCT oct;
  struct LOCT *boct;
  int ioct;
  FILE *fp;
  FILE *fp2;
  REAL xc,yc,zc;
  int field;
  REAL tsim=0.;
  float tsimf;
  char fname[256];
  char format[256];
  char fname2[256];
  int ncpu;
  int icpu=0;
  struct OCT *zerooct;
  struct UNITS unit;
  size_t dummy;
  struct cdata CD;
  int flagdata;

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
  dxmap=pow(0.5,lmap);

  // getting the number of CPUs
  sscanf(argv[5],"%d",&ncpu);


  // silo file
  int dumpsilo=0;
  sscanf(argv[6],"%d",&dumpsilo);
  float zmin,zmax,xmin,xmax,ymin,ymax;

  int mono;
  sscanf(argv[7],"%d",&mono);

  if(argc>8) {
    sscanf(argv[8],"%f",&xmin);
    sscanf(argv[9],"%f",&xmax);

    sscanf(argv[10],"%f",&ymin);
    sscanf(argv[11],"%f",&ymax);

    sscanf(argv[12],"%f",&zmin);
    sscanf(argv[13],"%f",&zmax);
  }
  else{
    xmin=0;
    xmax=1.;
    ymin=0;
    ymax=1.;
    zmin=0;
    zmax=1.;
  }


  // scanning the cpus

  int icpustart,icpustop;
  if(mono<0){
    icpustart=0;
    icpustop=ncpu;
    printf("Asked to parse %d files\n",ncpu);
  }
  else{
    icpustart=mono;
    icpustop=mono+1;
    printf("Cast single cpu #%d\n",mono) ;
  }
  // building the output file name
  strcpy(format,argv[4]);

  sprintf(fname2,format,icpu);
  int ndummy=0;
  fp2=fopen(fname2,"wb");
  fwrite(&ndummy,1,sizeof(int),fp2);
  fwrite(&tsim,1,sizeof(float),fp2);



  int idat=0;
  for(icpu=icpustart;icpu<icpustop;icpu++){


    // building the input file name
    strcpy(format,argv[1]);
    strcat(format,".p%05d");
    sprintf(fname,format,icpu);


    printf("Looking for %s\n",fname);
    fp=fopen(fname,"rb");

    if(fp==NULL){
      printf("--ERROR -- file not found\n");
      return 1;
    }


    // reading the number of octs

    int noct;
    fseek(fp, -sizeof(int), SEEK_END);
    dummy=fread(&noct,sizeof(int),1,fp);
    printf("noct in file=%d\n",noct);
    rewind(fp);

    // alloc the data
    boct=(struct LOCT*)malloc(sizeof(struct LOCT)*noct);



    // reading the time
    dummy=fread(&tsim,sizeof(REAL),1,fp);

#ifdef WRAD
    dummy=fread(&unit.unit_l,sizeof(REAL),1,fp);
    printf("unit=%e\n",unit.unit_l);
    dummy=fread(&unit.unit_v,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_t,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_n,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_d,sizeof(REAL),1,fp);
    dummy=fread(&unit.unit_N,sizeof(REAL),1,fp);
#endif

    printf("tsim=%e\n",tsim);
    tsimf=tsim;
    // reading the zerooct
    dummy=fread(&zerooct,sizeof(struct OCT *),1,fp);
    printf("0oct=%p\n",zerooct);

    printf("size of the OCT=%ld\n",sizeof(struct LOCT));

    ic=0;

    //    dummy=fread(&oct,sizeof(struct LOCT),1,fp);
    //while(!feof(fp)){
    dummy=fread(boct,sizeof(struct LOCT),noct,fp);

    for(ioct=0;ioct<noct;ioct++){
      memcpy(&oct,boct+ioct,sizeof(struct LOCT));
      if(oct.level<=lmap){
	dxcur=pow(0.5,oct.level);
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    if(((oct.cell[icell].child==0)||(oct.level==lmap)))
	      {
		xc=oct.x+( icell   %2)*dxcur;//+0.5*dxcur;
		yc=oct.y+((icell/2)%2)*dxcur;//+0.5*dxcur;
		zc=oct.z+( icell/4   )*dxcur;//+0.5*dxcur;
		if((zc<zmin)||(zc>zmax)){
		  continue;
		}
		if((yc<ymin)||(yc>ymax)){
		  continue;
		}
		if((xc<xmin)||(xc>xmax)){
		  continue;
		}

		int flag=0;
		if(mono>=0){
		  flag=(oct.cpu>=0);
		}
		else{
		  flag=(oct.cpu==icpu);
		}

		if(flag) {
		  data=assign_cube(field,icell,&oct,&unit);
		  idat++;
		  CD.x=(float)xc;CD.y=(float)yc;CD.z=(float)zc;
		  CD.level=(float)oct.level;
		  CD.data=(float)data;
		  fwrite(&CD,1,sizeof(struct cdata),fp2);
		}

	      }
	  }
      }
      //dummy=fread(&oct,sizeof(struct LOCT),1,fp); //reading next oct
    }
    printf("icpu=%d icell=%d\n",icpu,idat);
    fclose(fp);
    free(boct);
  }


  //============= dump
  printf("done with %d cells %e\n",idat,tsim);
  rewind(fp2);
  fwrite(&idat,1,sizeof(int),fp2);
  fwrite(&tsimf,1,sizeof(float),fp2);
  status=ferror(fp2);
  fclose(fp2);

  printf("status=%d\n",status);
}

// ================================================================
#if 1


float assign_cube(int field, int icell,struct LOCT *oct, struct UNITS *unit){

  float res;

  switch(field){
  case 0:
    res=oct->level;
    break;
#ifdef WGRAV
  case 1:
    res=oct->cell[icell].den;
    break;
#ifdef PIC
  case -1:
    res=oct->cell[icell].density;
    break;
#endif
  case 2:
    res=oct->cell[icell].pot;
    break;
  case 6:
    res+=oct->cell[icell].res;
    break;
  case 201:
    res=oct->cell[icell].f[0];
    break;
  case 202:
    res=oct->cell[icell].f[1];
    break;
  case 203:
    res=oct->cell[icell].f[2];
    break;
#endif
  case 3:
    res=oct->cpu;
    break;
  case 4:
    res=oct->cell[icell].marked;
    break;
#ifdef WHYDRO2
  case 101:
    res=oct->cell[icell].d;
    break;
  case 102:
    res=oct->cell[icell].u;
    break;
  case 103:
    res=oct->cell[icell].v;
    break;
  case 104:
    res=oct->cell[icell].w;
    break;
  case 105:
    res=oct->cell[icell].p;
    break;
  case 106:
    res=SQRT(
        oct->cell[icell].u*oct->cell[icell].u +
        oct->cell[icell].v*oct->cell[icell].v +
        oct->cell[icell].w*oct->cell[icell].w );
    break;
  case 107:
{
	
// radial velocity

    REAL vx= oct->cell[icell].u;
    REAL vy= oct->cell[icell].v;
    REAL vz= oct->cell[icell].w;

    REAL dx = POW(0.5,oct->level);

    REAL x=(oct->x+( icell    & 1)*dx+dx*0.5)-0.5;
    REAL y=(oct->y+((icell>>1)& 1)*dx+dx*0.5)-0.5;
    REAL z=(oct->z+( icell>>2    )*dx+dx*0.5)-0.5;

    res= (x*vx+y*vy+z*vz)/SQRT(x*x+y*y+z*z);
}
    break;

#endif

#ifdef WRAD
  case 701:
    res=oct->cell[icell].e[0];
    break;
  case 702:
    res=oct->cell[icell].fx[0];
    break;
  case 703:
    res=oct->cell[icell].fy[0];
    break;
  case 704:
    res=oct->cell[icell].fz[0];
    break;
  case 711:
    res=oct->cell[icell].e[1];
    break;
  case 712:
    res=oct->cell[icell].fx[1];
    break;
  case 713:
    res=oct->cell[icell].fy[1];
    break;
  case 714:
    res=oct->cell[icell].fz[1];
    break;
  case 709:
    res=oct->cell[icell].snfb;
    break;


  case 721:
    res=oct->cell[icell].e[2];
    break;
  case 722:
    res=oct->cell[icell].fx[2];
    break;
  case 723:
    res=oct->cell[icell].fy[2];
    break;
  case 724:
    res=oct->cell[icell].fz[2];
    break;

  case 705:
    res=log10(oct->cell[icell].src[0]+1e-15);
    break;
#ifdef WCHEM
  case 706:
    res=oct->cell[icell].xion;
    break;
#ifdef WRADHYD
  case 1006:
    res=oct->cell[icell].dX;
    break;
#ifdef HELIUM
  case 1008:
    res=oct->cell[icell].dXHE;
    break;
  case 1009:
    res=oct->cell[icell].dXXHE;
    break;
  case 708:
    res=oct->cell[icell].xHE;
    break;
  case 709:
    res=oct->cell[icell].xxHE;
    break;
#endif
#endif
  case 707:
    res=oct->cell[icell].temp;
    break;
#endif
#endif

  }

  return res;

}

// ===============================================================
// ===============================================================




#endif
