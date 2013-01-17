

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "silo.h"
#include "prototypes.h"
//=================================================================================================



/*

 mpicc  -I./silo/include -DPIC -DWGRAV -DTESTCOSMO -DFASTGRAV -DGRAFIC -DONFLYRED -DGPUAXL -o ./oct2grid ./oct2grid.c ./silo/lib/libsilo.a

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


  if(argc!=7){
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
  map=(float *)calloc(nmap*nmap*nmap,sizeof(float));
  dxmap=1./nmap;

  // getting the number of CPUs
  sscanf(argv[5],"%d",&ncpu);


  // silo file
  int dumpsilo=0;
  sscanf(argv[6],"%d",&dumpsilo);

  // scanning the cpus
  for(icpu=0;icpu<ncpu;icpu++){
    
    memset(map,0,nmap*nmap*nmap*sizeof(float));
    
    // building the input file name
    strcpy(format,argv[1]);
    strcat(format,".p%05d");
    sprintf(fname,format,icpu);

    // building the output file name
    strcpy(format,argv[4]);

    if(!dumpsilo){
      strcat(format,".p%05d");
    }
    sprintf(fname2,format,icpu);


    fp=fopen(fname,"rb");
  
    if(fp==NULL){
      printf("--ERROR -- file not found\n");
      return 1;
    }
    else{
      printf("Casting Rays on %dx%dx%d cube from %s\n",nmap,nmap,nmap,fname);
    }
    
    printf("size= %ld\n",nmap*nmap*nmap*sizeof(float)+sizeof(int)*2);
    // reading the time
    fread(&tsim,sizeof(REAL),1,fp);
    printf("tsim=%e\n",tsim);
    tsimf=tsim;
    // reading the zerooct
    fread(&zerooct,sizeof(struct OCT *),1,fp);
    printf("0oct=%p\n",zerooct);

    ic=0;
    
    fread(&oct,sizeof(struct OCT),1,fp);
    while(!feof(fp)){
      if(oct.level<=lmap){
      ic++;
      dxcur=1./pow(2.,oct.level);
      for(icell=0;icell<8;icell++) // looping over cells in oct
	{
	  if((oct.cell[icell].child==NULL)||(oct.level==lmap))
	    {
	      xc=oct.x+( icell   %2)*dxcur;//+0.5*dxcur;
	      yc=oct.y+((icell/2)%2)*dxcur;//+0.5*dxcur;
	      zc=oct.z+( icell/4   )*dxcur;//+0.5*dxcur;
	      imap=xc/dxmap;
	      jmap=yc/dxmap;
	      kmap=zc/dxmap;
	      
	      for(kk=0;kk<pow(2,lmap-oct.level);kk++)
		{
		  for(jj=0;jj<pow(2,lmap-oct.level);jj++)
		    {
		      for(ii=0;ii<pow(2,lmap-oct.level);ii++)
			{
			
			  switch(field){
			  case 0:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.level;
			    break;
#ifdef WGRAV
 			  case 1:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].gdata.d;
			    break;
			  case 2:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.cell[icell].gdata.p;
			    break;
			  case 6:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.cell[icell].res;
			    break;
			  case 201:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.cell[icell].f[0];
			    break;
			  case 202:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.cell[icell].f[1];
			    break;
			  case 203:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.cell[icell].f[2];
			    break;
#endif
			  case 3:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cpu;
			    break;
			  case 4:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].marked;
			    break;
			  case 5:
			    //map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].fx;
			    break;
#ifdef WHYDRO2
			  case 101:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].field.d;
			    //printf("%f\n",oct.cell[icell].field.d);
			    break;
			  case 102:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].field.u;
			    break;
			  case 103:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].field.v;
			    break;
			  case 104:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].field.w;
			    break;
			  case 105:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].field.p;
			    break;
			  case 106:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].field.E;
			    break;
			  case 107:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=(oct.cell[icell].field.E-0.5*oct.cell[icell].field.d*(oct.cell[icell].field.u*oct.cell[icell].field.u+oct.cell[icell].field.v*oct.cell[icell].field.v+oct.cell[icell].field.w*oct.cell[icell].field.w))*(GAMMA-1.); // alternative pressure
			    break;
#endif			  

#ifdef WRAD
			  case 701:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].rfield.e[0];
			    break;
			  case 702:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].rfield.fx[0];
			    break;
			  case 703:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].rfield.fy[0];
			    break;
			  case 704:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].rfield.fz[0];
			    break;
			  case 705:
			    map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].rfield.src[0];
			    break;
#endif


			  }
			}
		    }
		} 
	    }
	}
      }
      fread(&oct,sizeof(struct OCT),1,fp); //reading next oct
      
    }
    
    fclose(fp);
    printf("done with %d octs\n",ic);
      
    //============= dump
    
    if(!dumpsilo){
      printf("dumping %s with nmap=%d\n",fname2,nmap*nmap*nmap);
      fp=fopen(fname2,"wb");
      fwrite(&nmap,1,sizeof(int),fp);
      fwrite(&tsimf,1,sizeof(float),fp);
      for(I=0;I<nmap;I++) fwrite(map+I*nmap*nmap,nmap*nmap,sizeof(float),fp);
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
  }
  printf("status=%d\n",status);
  free(map);
  
}
