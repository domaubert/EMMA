

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "silo.h"
#include "prototypes.h"
//=================================================================================================



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
  REAL *map;
  int imap,jmap,kmap;
  int nmap;
  REAL dxmap,dxcur;
  int icur,ii,jj,kk,ic,icell;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
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


  // silo file
  int dumpsilo=0;
  sscanf(argv[6],"%d",&dumpsilo);
  REAL zmin,zmax;
  int nmapz;
  int k0;

  int mono;
  sscanf(argv[7],"%d",&mono);

  if(argc>8) {
    sscanf(argv[8],"%lf",&zmin);
    sscanf(argv[9],"%lf",&zmax);
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
 
  map=(REAL *)calloc(nmap*nmap*nmapz,sizeof(REAL));

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
      printf("Casting Rays on %dx%dx%d cube from %s\n",nmap,nmap,nmapz,fname);
    }
    
    printf("size= %ld\n",nmap*nmap*nmapz*sizeof(REAL)+sizeof(int)*2);
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
    printf("0oct=%p\n",zerooct);

    printf("size of the OCT=%ld\n",sizeof(struct OCT));

    ic=0;
    
    dummy=fread(&oct,sizeof(struct OCT),1,fp);
    while(!feof(fp)){
      if(oct.level<=lmap){
	ic++;
	dxcur=1./pow(2.,oct.level);
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if(((oct.cell[icell].child==NULL)||(oct.level==lmap)))
		{
		  xc=oct.x+( icell   %2)*dxcur;//+0.5*dxcur;
		  yc=oct.y+((icell/2)%2)*dxcur;//+0.5*dxcur;
		  zc=oct.z+( icell/4   )*dxcur;//+0.5*dxcur;
		  if((zc<zmin)||(zc>zmax)){
		    //printf("zc=%e SKIP\n",zc);
		    continue;
		  }

		  imap=xc/dxmap;
		  jmap=yc/dxmap;
		  kmap=zc/dxmap;
	      
		  for(kk=0;kk<pow(2,lmap-oct.level);kk++)
		    {
		      if((kmap+kk)>=(nmapz+k0)){
			//printf("zc=%e kmap=%d kk=%d nmpaz=%d k0=%d SKIP2\n",zc,kmap,kk,nmapz,k0);
			continue;
		      }
		      for(jj=0;jj<pow(2,lmap-oct.level);jj++)
			{
		      for(ii=0;ii<pow(2,lmap-oct.level);ii++)
			{

			  int flag;
			  if(mono>=0){
			    flag=(oct.cpu>=0);
			  }
			  else{
			    flag=(oct.cpu==icpu);
			  }
			  if(flag) {
			    switch(field){
			    case 0:
			      if(oct.level==10) printf("hello\n");
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.level;
			      break;
#ifdef WGRAV
			    case 1:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].gdata.d;
			      break;
#ifdef PIC
			    case -1:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].density;
			      break;
#endif
			    case 2:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].gdata.p;
			      break;
			    case 6:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]+=oct.cell[icell].res;
			      break;
			    case 201:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].f[0];
			      break;
			    case 202:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].f[1];
			      break;
			    case 203:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].f[2];
			      break;
#endif
			    case 3:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cpu;
			      break;
			    case 4:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].marked;
			      break;
#ifdef WHYDRO2
			    case 101:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].field.d;
			      //printf("%f\n",oct.cell[icell].field.d);
			      break;
			    case 102:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].field.u;
			      break;
			    case 103:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].field.v;
			      break;
			    case 104:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].field.w;
			      break;
			    case 105:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].field.p;
			      break;
			    case 106:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].field.E;
			      break;
			    case 107:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=(oct.cell[icell].field.E-0.5*oct.cell[icell].field.d*(oct.cell[icell].field.u*oct.cell[icell].field.u+oct.cell[icell].field.v*oct.cell[icell].field.v+oct.cell[icell].field.w*oct.cell[icell].field.w))*(GAMMA-1.); // alternative pressure
			      break;
			    case 110:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=sqrt(pow(oct.cell[icell].field.w,2)+pow(oct.cell[icell].field.v,2)+pow(oct.cell[icell].field.u,2))/oct.cell[icell].field.a;
			      break;
#ifdef WRADHYD
			    case 108:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].field.X;
			      break;
#endif
#endif			  

#ifdef WRAD
			    case 701:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.e[0];
			      break;
			    case 702:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fx[0];
			      break;
			    case 703:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fy[0];
			      break;
			    case 704:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fz[0];
			      break;
			    case 711:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.e[1];
			      break;
			    case 712:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fx[1];
			      break;
			    case 713:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fy[1];
			      break;
			    case 714:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fz[1];
			      break;

			    case 721:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.e[2];
			      break;
			    case 722:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fx[2];
			      break;
			    case 723:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fy[2];
			      break;
			    case 724:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.fz[2];
			      break;

			    case 705:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.src;
			      break;
#ifdef WCHEM
			    case 706:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.xion;
			      break;
			    case 707:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.temp;
			      break;
			    case 708:
			      map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk-k0)*nmap*nmap]=oct.cell[icell].rfield.nh/pow(unit.unit_l,3.);
			      break;
#endif
#endif
			      
			    }
			  }
			}
			}
		    } 
		}
	    }
      }
      dummy=fread(&oct,sizeof(struct OCT),1,fp); //reading next oct
    }    
  }
    fclose(fp);
    printf("done with %d octs\n",ic);
      
    //============= dump
    
    if(!dumpsilo){
      printf("dumping %s with nmap=%d\n",fname2,nmap*nmap*nmap);
      fp=fopen(fname2,"wb");
      fwrite(&nmap,1,sizeof(int),fp);
      fwrite(&nmapz,1,sizeof(int),fp);
      fwrite(&tsimf,1,sizeof(float),fp);
      for(I=0;I<nmapz;I++) fwrite(map+I*nmap*nmap,nmap*nmap,sizeof(REAL),fp);
      status=ferror(fp);
      fclose(fp);
    }
    else{
       DBfile *dbfile=NULL;
      dbfile=DBCreate(strcat(fname2,".silo"),DB_CLOBBER, DB_LOCAL,"silo file created by Quartz",DB_PDB);
      if(dbfile==NULL){
	printf("Error while writing file");
      }


      REAL *x;
      REAL *y;
      REAL *z;
      int dims[]={nmap+1,nmap+1,nmap+1};
      int ndims=3;
      
      x=(REAL*)malloc(sizeof(REAL)*(nmap+1));
      y=(REAL*)malloc(sizeof(REAL)*(nmap+1));
      z=(REAL*)malloc(sizeof(REAL)*(nmap+1));

      REAL *coords[]={x,y,z};
      
      int i;
      for(i=0;i<=nmap;i++){
	x[i]=((i*1.0)/nmap-0.);
	y[i]=((i*1.0)/nmap-0.);
	z[i]=((i*1.0)/nmap-0.);
      }

      int dimsvar[]={nmap,nmap,nmap};
      int ndimsvar=3;
      
      if(sizeof(REAL)==4){
	DBPutQuadmesh(dbfile,"quadmesh",NULL,coords,dims,ndims,DB_FLOAT,DB_COLLINEAR,NULL);
	DBPutQuadvar1(dbfile,"monvar","quadmesh",map,dimsvar,ndimsvar,NULL,0,DB_FLOAT,DB_ZONECENT,NULL);
      }
      else if(sizeof(REAL)==8){
	DBPutQuadmesh(dbfile,"quadmesh",NULL,coords,dims,ndims,DB_DOUBLE,DB_COLLINEAR,NULL);
	DBPutQuadvar1(dbfile,"monvar","quadmesh",map,dimsvar,ndimsvar,NULL,0,DB_DOUBLE,DB_ZONECENT,NULL);
      }

      DBClose(dbfile);
    }
  printf("status=%d\n",status);
  free(map);
  
}
