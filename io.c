#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"

//------------------------------------------------------------------------
void dumpmap(int lmap,struct OCT **firstoct,int field,char filename[],float zmin, float zmax)
{
  float *map;
  int imap,jmap;
  int nmap=pow(2,lmap);
  float dxmap=1./nmap,dxcur;
  int icur,ii,jj,ic,icell;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  FILE *fp;
  float xc,yc,zc;

  map=(float *)calloc(nmap*nmap,sizeof(float));

  //printf("==>  start map \n");
  for(level=1;level<=lmap;level++) // looping over octs
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];

      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
	  nextoct=oct.next;
	  //	  printf("%f %f %f\n",oct.x,oct.y,oct.z);
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if((oct.cell[icell].child==NULL)||(oct.level==lmap))
		{
		  xc=oct.x+( icell   %2)*dxcur;//+0.5*dxcur;
		  yc=oct.y+((icell/2)%2)*dxcur;//+0.5*dxcur;
		  zc=oct.z+( icell/4   )*dxcur;//+0.5*dxcur;
		  imap=xc/dxmap;
		  jmap=yc/dxmap;
		  
		  if((zc>zmax)||(zc<zmin)) continue;

		  //if(grid[icur].level==lmap) printf("%d %f %f %d %d %f\n",icell,xc,yc,imap,jmap,grid[icur].dens[icell]);
		  for(jj=0;jj<pow(2,lmap-oct.level);jj++)
		    {
		      for(ii=0;ii<pow(2,lmap-oct.level);ii++)
			{

			  switch(field){
			  case 0:
			    map[(imap+ii)+(jmap+jj)*nmap]=fmax(oct.level,map[(imap+ii)+(jmap+jj)*nmap]);
			    break;
			  case 1:
			    map[(imap+ii)+(jmap+jj)*nmap]+=oct.cell[icell].density*pow(2,lmap-oct.level);
			    break;
			  case 2:
			    map[(imap+ii)+(jmap+jj)*nmap]+=oct.cell[icell].pot*pow(2,lmap-oct.level);
			    break;
			  }
			}
		    }
		}
	    }
	}while(nextoct!=NULL);
    }
  
  
  //============= dump

  //printf("dumping %s\n",filename);
  fp=fopen(filename,"wb");
  fwrite(&nmap,1,sizeof(int),fp);
  fwrite(map,nmap*nmap,sizeof(float),fp);
  fclose(fp);
  free(map);

}

void dumpcube(int lmap,struct OCT **firstoct,int field,char filename[],float tsim)
{
  float *map;
  int imap,jmap,kmap;
  int nmap=pow(2,lmap);
  float dxmap=1./nmap,dxcur;
  int icur,ii,jj,kk,ic,icell;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  FILE *fp;
  float xc,yc,zc;

  map=(float *)calloc(nmap*nmap*nmap,sizeof(float));

  //printf("==> start map \n");
  for(level=1;level<=lmap;level++) // looping over octs
    {
      //printf("level=%d\n",level);
      // setting the first oct

      nextoct=firstoct[level-1];

      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
	  nextoct=oct.next;
	  //	  printf("%f %f %f\n",oct.x,oct.y,oct.z);
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
		  
		  //if(grid[icur].level==lmap) printf("%d %f %f %d %d %f\n",icell,xc,yc,imap,jmap,grid[icur].dens[icell]);
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
#ifdef AXLFORCE
			      case 6:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].fx;
				break;
			      case 7:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].fy;
				break;
			      case 8:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].fz;
				break;
#endif
			      }
			    }
			}
		    }
		}
	    }
	}while(nextoct!=NULL);
    }
  
  //============= dump
  
  //printf("dumping %s\n",filename);
  fp=fopen(filename,"wb");
  fwrite(&nmap,1,sizeof(int),fp);
  fwrite(&tsim,1,sizeof(float),fp);
  fwrite(map,nmap*nmap*nmap,sizeof(float),fp);
  fclose(fp);
  free(map);

}
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
void dumppart(struct OCT **firstoct,char filename[],int npart, int levelcoarse, int levelmax, float tsim){

  FILE *fp;
  float val;
  int vali;
  int ipart=0;
  int level;
  struct OCT *nextoct;
  struct OCT oct;
  float dxcur;
  struct PART *nexp;
  struct PART *curp;
  int icell;

  fp=fopen(filename,"wb");
  fwrite(&npart,1,sizeof(int),fp);
  fwrite(&tsim,1,sizeof(float),fp);
  for(level=levelcoarse;level<=levelmax;level++) // looping over levels
    {
      //printf("level=%d\n",level);
      // setting the first oct
      
      nextoct=firstoct[level-1];
      
      do // sweeping through the octs of level
	{
	  if(nextoct==NULL) continue; // in case the level is empty
	  oct=(*nextoct);
	  dxcur=1./pow(2,oct.level);
	  nextoct=oct.next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=oct.cell[icell].phead; //sweeping the particles of the current cell
	      if(nexp!=NULL){
		do{ 
		  curp=nexp; 
		  nexp=curp->next; 
		  val=curp->x;fwrite(&val,1,sizeof(float),fp);
		  val=curp->y;fwrite(&val,1,sizeof(float),fp);
		  val=curp->z;fwrite(&val,1,sizeof(float),fp);
#ifndef PARTN
		  val=curp->vx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vy;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vz;fwrite(&val,1,sizeof(float),fp);
#else
		  val=curp->fx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fz;fwrite(&val,1,sizeof(float),fp);
#endif
		  val=(float)(curp->idx);fwrite(&val,1,sizeof(float),fp);
		  ipart++;
		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }
  fclose(fp);
  printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

//------------------------------------------------------------------------

void GetParameters(char *fparam, struct RUNPARAMS *param)
{
  FILE *buf; 
  char stream[256];
  size_t rstat;

  buf=fopen(fparam,"r");
  if(buf==NULL)
    {
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
    }
  else
    {
      rstat=fscanf(buf,"%s %d",stream,&param->npartmax);
      rstat=fscanf(buf,"%s %d",stream,&param->ngridmax);
      rstat=fscanf(buf,"%s %d",stream,&param->nbuff);
      rstat=fscanf(buf,"%s %d",stream,&param->ndumps);
      rstat=fscanf(buf,"%s %d",stream,&param->nsteps);
      rstat=fscanf(buf,"%s %d",stream,&param->lcoarse);
      rstat=fscanf(buf,"%s %d",stream,&param->lmax);
      rstat=fscanf(buf,"%s %d",stream,&param->levelmap);
      rstat=fscanf(buf,"%s %d",stream,&param->niter);
      rstat=fscanf(buf,"%s %f",stream,&param->poissonacc);
      rstat=fscanf(buf,"%s %d",stream,&param->mgridlmin);
      rstat=fscanf(buf,"%s %d",stream,&param->nvcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->stride);
      rstat=fscanf(buf,"%s %f",stream,&param->dt);
      rstat=fscanf(buf,"%s %f",stream,&param->amrthresh);
      fclose(buf);
    }

  // computing the maxhash
  int val=(pow(2,param->lmax-1)<64?pow(2,param->lmax-1):64); // limit to 2097152 octs in hash table i.e. 16e6 cells
  param->maxhash=pow(val,3);
  printf("maxhash=%d\n",param->maxhash);
}
