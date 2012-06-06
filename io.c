#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "friedmann.h"
#include "segment.h"

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


void dumpgrid(int lmap,struct OCT **firstoct,int field,char filename[],float tsim)
{

  int icur,ii,jj,kk,ic,icell;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  FILE *fp;
  float xc,yc,zc;


  fp=fopen(filename,"wb");
  
  ic=0;
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
	  nextoct=oct.next;
	  
	  fwrite(&oct,sizeof(struct OCT),1,fp);
	  if(oct.level==lmap) ic+=8;
	}while(nextoct!=NULL);
    }

 
  fclose(fp);
 
}


//=================================================================================================

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

#ifdef WHYDRO2
			      case 101:
				map[(imap+ii)+(jmap+jj)*nmap+(kmap+kk)*nmap*nmap]=oct.cell[icell].field.d;
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
      rstat=fscanf(buf,"%s %d",stream,&param->ngridmax);
      rstat=fscanf(buf,"%s %d",stream,&param->npartmax);
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
      rstat=fscanf(buf,"%s %d",stream,&param->nrelax);
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

//==================================================================================
//==================================================================================

struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, float *munit, float *ainit, float *omegam, float *omegav, float *Hubble, int *npart, float omegab){
  
  FILE *fx, *fy, *fz;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0;
  int dummy;
  struct PART *lastpart;
  int ip;

  fx=fopen("utils/grafic_src/ic_velcx","rb");
  fy=fopen("utils/grafic_src/ic_velcy","rb");
  fz=fopen("utils/grafic_src/ic_velcz","rb");
  

  // reading the headers

  fread(&dummy,1,sizeof(dummy),fx);
  fread(&np1,1,4,fx);
  fread(&np2,1,4,fx);
  fread(&np3,1,4,fx);
  fread(&dx,1,4,fx);
  fread(&x1o,1,4,fx);
  fread(&x2o,1,4,fx);
  fread(&x3o,1,4,fx);
  fread(&astart,1,4,fx);
  fread(&om,1,4,fx);
  fread(&ov,1,4,fx);
  fread(&h0,1,4,fx);
  fread(&dummy,1,sizeof(dummy),fx);

  fread(&dummy,1,sizeof(dummy),fy);
  fread(&np1,1,4,fy);
  fread(&np2,1,4,fy);
  fread(&np3,1,4,fy);
  fread(&dx,1,4,fy);
  fread(&x1o,1,4,fy);
  fread(&x2o,1,4,fy);
  fread(&x3o,1,4,fy);
  fread(&astart,1,4,fy);
  fread(&om,1,4,fy);
  fread(&ov,1,4,fy);
  fread(&h0,1,4,fy);
  fread(&dummy,1,sizeof(dummy),fy);

  fread(&dummy,1,sizeof(dummy),fz);
  fread(&np1,1,4,fz);
  fread(&np2,1,4,fz);
  fread(&np3,1,4,fz);
  fread(&dx,1,4,fz);
  fread(&x1o,1,4,fz);
  fread(&x2o,1,4,fz);
  fread(&x3o,1,4,fz);
  fread(&astart,1,4,fz);
  fread(&om,1,4,fz);
  fread(&ov,1,4,fz);
  fread(&h0,1,4,fz);
  fread(&dummy,1,sizeof(dummy),fz);

  printf("============================================\n");
  printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
  printf("om=%f ov=%f h0=%f\n",om,ov,h0);
  printf("dx=%f np1*dx=%f\n",dx,np1*dx);
  printf("astart=%f zstart=%f\n",astart,1./astart-1.);
  printf("============================================\n");


  // computing Zeldovich displacement quantities
  
  double vfact;
  vfact=fomega(astart,om,ov)*h0*dladt(astart,om,ov)/astart;
  printf("vfact=%f\n",vfact);
  // reading the grafic planes

  float *velx;
  float *vely;
  float *velz;

  float x;
  float y;
  float z;

  float vx;
  float vy;
  float vz;

  float x0,y0,z0;

  int i1,i2,i3;

  float mass;

#ifdef WHYDRO
  mass=(1.-omegab/om)/(np1*np2*np3);
#else
  mass=1./(np1*np2*np3);
#endif

  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);


  ip=0;
  for(i3=1;i3<=np3;i3++){

    fread(&dummy,1,sizeof(dummy),fx);
    fread(velx,np1*np2,sizeof(float),fx);
    fread(&dummy,1,sizeof(dummy),fx);

    fread(&dummy,1,sizeof(dummy),fy);
    fread(vely,np1*np2,sizeof(float),fy);
    fread(&dummy,1,sizeof(dummy),fy);

    fread(&dummy,1,sizeof(dummy),fz);
    fread(velz,np1*np2,sizeof(float),fz);
    fread(&dummy,1,sizeof(dummy),fz);

    z0=(i3-0.5)*dx;

    for(i2=1;i2<=np2;i2++){
      y0=(i2-0.5)*dx;
      for(i1=1;i1<=np1;i1++){
	x0=(i1-0.5)*dx;
	// computing the displacements
	x=(x0+velx[(i1-1)+(i2-1)*np1]/vfact)/(np1*dx);
	y=(y0+vely[(i1-1)+(i2-1)*np1]/vfact)/(np2*dx);
	z=(z0+velz[(i1-1)+(i2-1)*np1]/vfact)/(np3*dx);

	// periodic boundary conditions

	x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x); 
	y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y); 
	z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z); 

	// computing the velocities
	vx=velx[(i1-1)+(i2-1)*np1]*astart/(np1*dx*h0)/(sqrt(om)*0.5);
	vy=vely[(i1-1)+(i2-1)*np1]*astart/(np2*dx*h0)/(sqrt(om)*0.5);
	vz=velz[(i1-1)+(i2-1)*np1]*astart/(np3*dx*h0)/(sqrt(om)*0.5);

	// if it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,cpu,cpu->levelcoarse)){
	  
	  part[ip].x=x;
	  part[ip].y=y;
	  part[ip].z=z;
	  
	  part[ip].vx=vx;
	  part[ip].vy=vy;
	  part[ip].vz=vz;
	  
	  part[ip].mass=mass;
	  part[ip].idx=(i1-1)+(i2-1)*np1+(i3-1)*np1*np2;
	  lastpart=part+ip;
	  ip++;
	}
      }
    }
  }

  free(velx);
  free(vely);
  free(velz);
  fclose(fx);
  fclose(fy);
  fclose(fz);

  *munit=mass;
  *ainit=astart;
  *omegam=om;
  *omegav=ov;
  *Hubble=h0;
  *npart=ip;


  printf("Grafic Particle Read ok\n");
  return lastpart;
}

//==================================================================================
//==================================================================================

#ifdef WHYDRO
int read_grafic_hydro(struct CPUINFO *cpu, float omegab, float omegam){
  
  FILE *fx;
  FILE *fy;
  FILE *fz;
  FILE *fdx;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0;
  int dummy;
  struct PART *lastpart;
  int ip;

  fdx=fopen("utils/grafic_src/ic_deltab","rb");
  fx=fopen("utils/grafic_src/ic_velcx","rb");
  fy=fopen("utils/grafic_src/ic_velcy","rb");
  fz=fopen("utils/grafic_src/ic_velcz","rb");
  
  
  // reading the headers

  fread(&dummy,1,sizeof(dummy),fdx);
  fread(&np1,1,4,fdx);
  fread(&np2,1,4,fdx);
  fread(&np3,1,4,fdx);
  fread(&dx,1,4,fdx);
  fread(&x1o,1,4,fdx);
  fread(&x2o,1,4,fdx);
  fread(&x3o,1,4,fdx);
  fread(&astart,1,4,fdx);
  fread(&om,1,4,fdx);
  fread(&ov,1,4,fdx);
  fread(&h0,1,4,fdx);
  fread(&dummy,1,sizeof(dummy),fdx);

  fread(&dummy,1,sizeof(dummy),fx);
  fread(&np1,1,4,fx);
  fread(&np2,1,4,fx);
  fread(&np3,1,4,fx);
  fread(&dx,1,4,fx);
  fread(&x1o,1,4,fx);
  fread(&x2o,1,4,fx);
  fread(&x3o,1,4,fx);
  fread(&astart,1,4,fx);
  fread(&om,1,4,fx);
  fread(&ov,1,4,fx);
  fread(&h0,1,4,fx);
  fread(&dummy,1,sizeof(dummy),fx);

  fread(&dummy,1,sizeof(dummy),fy);
  fread(&np1,1,4,fy);
  fread(&np2,1,4,fy);
  fread(&np3,1,4,fy);
  fread(&dx,1,4,fy);
  fread(&x1o,1,4,fy);
  fread(&x2o,1,4,fy);
  fread(&x3o,1,4,fy);
  fread(&astart,1,4,fy);
  fread(&om,1,4,fy);
  fread(&ov,1,4,fy);
  fread(&h0,1,4,fy);
  fread(&dummy,1,sizeof(dummy),fy);

  fread(&dummy,1,sizeof(dummy),fz);
  fread(&np1,1,4,fz);
  fread(&np2,1,4,fz);
  fread(&np3,1,4,fz);
  fread(&dx,1,4,fz);
  fread(&x1o,1,4,fz);
  fread(&x2o,1,4,fz);
  fread(&x3o,1,4,fz);
  fread(&astart,1,4,fz);
  fread(&om,1,4,fz);
  fread(&ov,1,4,fz);
  fread(&h0,1,4,fz);
  fread(&dummy,1,sizeof(dummy),fz);

  printf("============================================\n");
  printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
  printf("om=%f ov=%f h0=%f\n",om,ov,h0);
  printf("dx=%f np1*dx=%f\n",dx,np1*dx);
  printf("astart=%f zstart=%f\n",astart,1./astart-1.);
  printf("============================================\n");


  // reading the grafic planes

  float *deltab;
  float *velz;
  float *vely;
  float *velx;

  int i1,i2,i3;
  int icx,icy,icz,icell;
  unsigned long key;
  struct OCT *curoct;
  struct OCT *nextoct;
  int hidx;
  int found;
  float z0,y0,x0;
  int ifound=0;


  deltab=(float*)malloc(sizeof(float)*np1*np2);
  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);

  double rhob,pressure;
  double Gnewt=6.674e-11; // Newton G
  double H0=h0*1e3/3.08568025e22; // Hubble constant (s-1)
  double rhoc=3.*H0*H0/(8.*M_PI*Gnewt); // comoving critical density (kg/m3)
  double mp=1.67262158e-27; // proton mass (kg)
  double mu=0.59; // mean molecular weight
  double kboltz=1.3806503e-23; // boltzmann constant SI
  double zstart=1./astart-1.;
  double temp=550.*((1.0+zstart)*(1.0+zstart)/(201.*201.)); // baryon temperature (to check) in K
  //double temp=1e4;

  // supercomoving unit values
  double rhostar;
  double rstar;
  double vstar;
  double tstar;
  double tstar2;
  double pstar;
  double mpc=3.08568025e22; // Mpc in m

  rstar= np1*dx*mpc; // box size in m
  rhostar=rhoc*om;
  tstar=2./H0/sqrt(om); // sec
  tstar2=2./h0/sqrt(om); // Mpc sec / km
  vstar=rstar/tstar;
  pstar=rhostar*vstar*vstar;

  printf("rhoc=%e temperature=%lf rstar=%e pstar=%e tstar=%e vstar=%e rhostar=%e\n",rhoc,temp,rstar,pstar,tstar,vstar,rhostar);

  for(i3=0;i3<np3;i3++){

    fread(&dummy,1,sizeof(dummy),fdx);
    fread(deltab,np1*np2,sizeof(float),fdx);
    fread(&dummy,1,sizeof(dummy),fdx);

    fread(&dummy,1,sizeof(dummy),fx);
    fread(velx,np1*np2,sizeof(float),fx);
    fread(&dummy,1,sizeof(dummy),fx);

    fread(&dummy,1,sizeof(dummy),fy);
    fread(vely,np1*np2,sizeof(float),fy);
    fread(&dummy,1,sizeof(dummy),fy);

    fread(&dummy,1,sizeof(dummy),fz);
    fread(velz,np1*np2,sizeof(float),fz);
    fread(&dummy,1,sizeof(dummy),fz);

    z0=(i3*1.0)/(np3);
    for(i2=0;i2<np2;i2++){
      y0=(i2*1.0)/(np2);
      for(i1=0;i1<np1;i1++){
	x0=(i1*1.0)/(np1);
	
	key=pos2key(x0,y0,z0,cpu->levelcoarse);

	// first we compute the adress from the hashfunction
	hidx=hfun(key,cpu->maxhash);
	nextoct=cpu->htable[hidx];

	// looking for the oct
	found=0;
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,cpu->levelcoarse)==key)&&(curoct->level==cpu->levelcoarse));
	  }while((nextoct!=NULL)&&(!found));
	}

	// filling the cell

	if(found){
	  icx=i1%2;
	  icy=i2%2;
	  icz=i3%2;
	  
	  icell=icx+icy*2+icz*4;
	  /* if(i1>115) printf("x0=%f y0=%f z0=%f key=%d oct=%p deltab=%f icell=%d ix=%f %d\n",x0,y0,z0,key,curoct,deltab[i1+np1*i2],icell,(x0/pow(0.5,6)),(int)(x0*pow(2,6))); */

	  if(curoct->cell[icell].d!=0){
	    abort();
	  }
	
	  rhob=(deltab[i1+i2*np1]+1.0)*omegab*rhoc/(astart*astart*astart); // physical baryon density in kg/m3
	  pressure=(GAMMA-1.0)*1.5*rhob*kboltz*temp/(mu*mp); // physical pressure

	  // filling the cells using supercomoving values
	  curoct->cell[icell].d=(deltab[i1+i2*np1]+1.0)*omegab/om;
	  curoct->cell[icell].u=(velx[i1+i2*np1]      )*astart/(np1*dx*h0)/(sqrt(om)*0.5);
	  curoct->cell[icell].v=(vely[i1+i2*np1]      )*astart/(np1*dx*h0)/(sqrt(om)*0.5);
	  curoct->cell[icell].w=(velz[i1+i2*np1]      )*astart/(np1*dx*h0)/(sqrt(om)*0.5);
	  curoct->cell[icell].p=pressure*pow(astart,5)/pstar;
	  ifound++;
	}



      }
    }
  }

  fclose(fdx);
  fclose(fx);
  fclose(fy);
  fclose(fz);

  free(deltab);
  free(velx);
  free(vely);
  free(velz);

  printf("Grafic hydro read ok\n");
  return ifound;
}

#endif
