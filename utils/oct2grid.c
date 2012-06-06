

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "prototypes.h"

//=================================================================================================

int main(int argc, char *argv[])
{
  //(int lmap,struct OCT **firstoct,int field,char filename[],float tsim)
  int lmap;
  float *map;
  int imap,jmap,kmap;
  int nmap;
  float dxmap,dxcur;
  int icur,ii,jj,kk,ic,icell;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  FILE *fp;
  float xc,yc,zc;
  int field;
  float tsim=0.;
  char fname[256];

  fp=fopen(argv[1],"rb");
  
  if(fp==NULL){
    printf("--ERROR -- file not found\n");
  }
  else{
    sscanf(argv[2],"%d",&lmap);
    nmap=pow(2,lmap);
    //printf("Casting Rays on %dx%dx%d cube from %s\n",nmap,nmap,nmap,argv[1]);
  }

  strcpy(fname,argv[1]);
  sscanf(argv[3],"%d",&field);


  map=(float *)calloc(nmap*nmap*nmap,sizeof(float));
  dxmap=1./nmap;


  
  ic=0;
    
  fread(&oct,sizeof(struct OCT),1,fp);
  while(!feof(fp)){
    ic++;
    dxcur=1./pow(2,oct.level);
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
	    		case 1:
	    		  map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].density;
	    		  break;
	    		case 2:
	    		  map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]+=oct.cell[icell].pot;
	    		  break;
	    		case 3:
	    		  map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cpu;
	    		  break;
	    		case 4:
	    		  map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].marked;
	    		  break;
	    		case 5:
	    		  map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].temp;
	    		  break;
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
			  
	    		}
	    	      }
	    	  }
	      } 
	  }
      }
    
    fread(&oct,sizeof(struct OCT),1,fp); //reading next oct

  }

  fclose(fp);
  printf("done with %d cells\n",ic);

  strcat(fname,".den");


  
  //============= dump
  
  //printf("dumping %s\n",filename);
  fp=fopen(fname,"wb");
  fwrite(&nmap,1,sizeof(int),fp);
  fwrite(&tsim,1,sizeof(float),fp);
  fwrite(map,nmap*nmap*nmap,sizeof(float),fp);
  fclose(fp);
  free(map);

}
