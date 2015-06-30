

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "silo.h"
#include "../src/prototypes.h"
//=================================================================================================
#define max(a,b) (a>=b?a:b)

/*

 mpicc  -I./silo/include -DPIC -DWGRAV -DTESTCOSMO -DFASTGRAV -DGRAFIC -DONFLYRED -DGPUAXL -o ./oct2grid ./oct2grid.c ./silo/lib/libsilo.a

avconv -i mpeg_link%04d.jpeg -r 24 -b 65536k video.mp4

 */


void assign_cube(int field, int icell, float *map, int imap, int jmap, int kmap, int ii, int jj, int kk, int i0, int j0, int k0, int nmapx, int nmapy, int nmapz, struct LOCT *oct, struct UNITS *unit);
void assign_zmap(int field, int icell, float *map, int imap, int jmap, int kmap, int ii, int jj, int kk, int i0, int j0, int k0, int nmapx, int nmapy, int nmapz, struct LOCT *oct, struct UNITS *unit);


int main(int argc, char *argv[])
{
  //(int lmap,struct LOCT **firstoct,int field,char filename[],REAL tsim)
  int status;
  int I;
  int lmap;
  //REAL *map;
  float *map;
  float data;
  int imap,jmap,kmap;
  int nmap;
  REAL dxmap,dxcur;
  int icur,ii,jj,kk,ic,icell;
  int level;
  struct LOCT * nextoct;
  struct LOCT oct;
  struct LOCT *boct;
  int ioct;
  FILE *fp;
  REAL xc,yc,zc;
  int field;
  REAL tsim=0.;
  float tsimf;
  char fname[256];
  char format[256];
  char fname2[256];
  int ncpu;
  int icpu;
  struct OCT *zerooct;
  struct UNITS unit;
  size_t dummy;
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
  nmap=pow(2,lmap);
  printf("nmap=%d\n",nmap);
  dxmap=1./nmap;

  // getting the number of CPUs
  sscanf(argv[5],"%d",&ncpu);


  // silo file
  // Note dumpsilo is kept for consistency
  // but is deprecated --> should be removed

  int dumpsilo=0;
  sscanf(argv[6],"%d",&dumpsilo);
  float zmin,zmax,xmin,xmax,ymin,ymax;
  int nmapz,nmapx,nmapy;
  int i0,j0,k0;

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

  int proj=0;
  if(argc>14) {
    sscanf(argv[14],"%d",&proj);
  }

  nmapx=(xmax-xmin)/dxmap;
  i0=xmin/dxmap;
  nmapy=(ymax-ymin)/dxmap;
  j0=ymin/dxmap;
  nmapz=(zmax-zmin)/dxmap;
  k0=zmin/dxmap;

  // should we dump cell data ?
  int dumpcell=0;
  if(argc>15){
    sscanf(argv[14],"%d",&dumpcell);
  }




  /// allocating the cube or the map
  if(proj==0){
    map=(float *)calloc(nmapx*nmapy*nmapz,sizeof(float));
    printf("nmapx=%d nmapy=%d nmapz=%d\n",nmapx,nmapy,nmapz);
  }
  else if(proj==3){
    map=(float *)calloc(nmapx*nmapy,sizeof(float));
  }
  else{
    printf("proj other than 0 or 3 not implemented yet");
    abort();
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

  for(icpu=icpustart;icpu<icpustop;icpu++){


    // building the input file name
    strcpy(format,argv[1]);
    strcat(format,".p%05d");
    sprintf(fname,format,icpu);

    // building the output file name
    strcpy(format,argv[4]);

    /* if(!dumpsilo){ */
    /*   strcat(format,".p%05d"); */
    /* } */
    sprintf(fname2,format,icpu);

    printf("Looking for %s\n",fname);
    fp=fopen(fname,"rb");

    if(fp==NULL){
      printf("--ERROR -- file not found\n");
      return 1;
    }
    else{
      printf("Casting Rays on %dx%dx%d cube from %s\n",nmapx,nmapy,nmapz,fname);
    }

    printf("size= %ld\n",nmapx*nmapy*nmapz*sizeof(float)+sizeof(int)*2);

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

    printf("size of the OCT=%ld\n",sizeof(struct OCT));

    ic=0;
    int idat=0;

    //    dummy=fread(&oct,sizeof(struct LOCT),1,fp);
    //while(!feof(fp)){
    dummy=fread(boct,sizeof(struct LOCT),noct,fp);

    for(ioct=0;ioct<noct;ioct++){
      memcpy(&oct,boct+ioct,sizeof(struct LOCT));
      if(oct.level<=lmap){
	ic++;
	dxcur=1./pow(2.,oct.level);

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

		imap=xc/dxmap;
		jmap=yc/dxmap;
		kmap=zc/dxmap;
		flagdata=0;
		for(kk=0;kk<pow(2,lmap-oct.level);kk++)
		  {
		    if((kmap+kk)>=(nmapz+k0)){
		      continue;
		    }
		    for(jj=0;jj<pow(2,lmap-oct.level);jj++)
		      {
			if((jmap+jj)>=(nmapy+j0)){
			  continue;
			}
			for(ii=0;ii<pow(2,lmap-oct.level);ii++)
			  {
			    if((imap+ii)>=(nmapx+i0)){
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
			      if(proj==0){
				assign_cube(field,icell,map,imap,jmap,kmap,ii,jj,kk,i0,j0,k0,nmapx,nmapy,nmapz,&oct,&unit);
			      }
			      else if(proj==3){
				assign_zmap(field,icell,map,imap,jmap,kmap,ii,jj,kk,i0,j0,k0,nmapx,nmapy,nmapz,&oct,&unit);
			      }
			    }



			  }
		      }
		  }

	      }
	  }
      }
      //dummy=fread(&oct,sizeof(struct LOCT),1,fp); //reading next oct
    }
    printf("icpu=%d ic=%d\n",icpu,ic);
    fclose(fp);
    free(boct);
  }
  printf("done with %d octs\n",ic);

  //============= dump


  // simple array style
  if(!dumpsilo){
    if(proj==3){
      nmapz=1;
    }
    printf("dumping %s with nmap=%d hohoh\n",fname2,nmap*nmapx*nmapy);
    fp=fopen(fname2,"wb");
    fwrite(&nmapx,1,sizeof(int),fp);
    fwrite(&nmapy,1,sizeof(int),fp);
    fwrite(&nmapz,1,sizeof(int),fp);
    fwrite(&tsimf,1,sizeof(float),fp);
    for(I=0;I<nmapz;I++) fwrite(map+I*nmapx*nmapy,nmapx*nmapy,sizeof(float),fp);
    status=ferror(fp);
    fclose(fp);
  }


#if 0
  else{
    // silo style

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
#endif

  printf("status=%d\n",status);
  free(map);
}

// ================================================================
#if 1


void assign_cube(int field, int icell, float *map, int imap, int jmap, int kmap, int ii, int jj, int kk, int i0, int j0, int k0, int nmapx, int nmapy, int nmapz, struct LOCT *oct, struct UNITS *unit){

  switch(field){
  case 0:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->level;
    break;
#ifdef WGRAV
  case 1:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].den;
    break;
#ifdef PIC
  case -1:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].density;
    break;
#endif
  case 2:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].pot;
    break;
  case 6:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]+=oct->cell[icell].res;
    break;
  case 201:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].f[0];
    break;
  case 202:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].f[1];
    break;
  case 203:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].f[2];
    break;
#endif
  case 3:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cpu;
    break;
  case 4:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].marked;
    break;
#ifdef WHYDRO2
  case 101:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].d;
    break;
  case 102:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].u;
    break;
  case 103:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].v;
    break;
  case 104:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].w;
    break;
  case 105:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].p;
    break;
  case 106:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=SQRT(
        oct->cell[icell].u*oct->cell[icell].u +
        oct->cell[icell].v*oct->cell[icell].v +
        oct->cell[icell].w*oct->cell[icell].w );
    break;
#endif

#ifdef WRAD
  case 701:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].e[0];
    break;
  case 702:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fx[0];
    break;
  case 703:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fy[0];
    break;
  case 704:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fz[0];
    break;
  case 711:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].e[1];
    break;
  case 712:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fx[1];
    break;
  case 713:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fy[1];
    break;
  case 714:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fz[1];
    break;
  case 709:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].snfb;
    break;


  case 721:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].e[2];
    break;
  case 722:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fx[2];
    break;
  case 723:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fy[2];
    break;
  case 724:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].fz[2];
    break;

  case 705:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=log10(oct->cell[icell].src+1e-15);
    break;
#ifdef WCHEM
  case 706:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].xion;
    break;
#ifdef WRADHYD
  case 1006:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].dX;
    break;
#ifdef HELIUM
  case 1009:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].dXXHE;
    break;
  case 1008:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].dXHE;
    break;
  case 708:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].xHE;
    break;
  case 709:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].xxHE;
    break;
#endif
#endif
  case 707:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx+(kmap+kk-k0)*nmapx*nmapy]=oct->cell[icell].temp;
    break;
#endif
#endif

  }


}

// ===============================================================
// ===============================================================

void assign_zmap(int field, int icell, float *map, int imap, int jmap, int kmap, int ii, int jj, int kk, int i0, int j0, int k0, int nmapx, int nmapy, int nmapz, struct LOCT *oct, struct UNITS *unit){

  switch(field){
  case 0:
    //map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->level;
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]=max(oct->level, map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]);
    break;
#ifdef WGRAV
  case 1:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].den;
    break;
#ifdef PIC
  case -1:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].density;
    break;
#endif
  case 2:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].pot;
    break;
  case 6:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].res;
    break;
  case 201:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].f[0];
    break;
  case 202:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].f[1];
    break;
  case 203:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].f[2];
    break;
#endif
  case 3:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cpu;
    break;
  case 4:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].marked;
    break;
#ifdef WHYDRO2
  case 101:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].d;
    break;
  case 102:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].u;
    break;
  case 103:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].v;
    break;
  case 104:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].w;
    break;
  case 105:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].p;
    break;

  case 106:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*SQRT(
        oct->cell[icell].u*oct->cell[icell].u +
        oct->cell[icell].v*oct->cell[icell].v +
        oct->cell[icell].w*oct->cell[icell].w );
    break;
#endif

#ifdef WRAD
  case 701:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].e[0];
    break;
  case 702:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fx[0];
    break;
  case 703:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fy[0];
    break;
  case 704:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fz[0];
    break;
  case 711:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].e[1];
    break;
  case 712:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fx[1];
    break;
  case 713:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fy[1];
    break;
  case 714:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fz[1];
    break;

  case 721:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].e[2];
    break;
  case 722:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fx[2];
    break;
  case 723:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fy[2];
    break;
  case 724:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].fz[2];
    break;

  case 705:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].src;
    break;
#ifdef WCHEM
  case 706:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].xion;
    break;
  case 707:
    map[(imap+ii-i0)+(jmap+jj-j0)*nmapx]+=(1./nmapz)*oct->cell[icell].temp;
    break;
#endif
#endif

  }
}



#endif
