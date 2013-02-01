#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "friedmann.h"
#include "segment.h"
#include <string.h>


//====================================================================================================
//====================================================================================================

void dumpgrid(int levelmax,struct OCT **firstoct, char filename[],REAL tsim)
{

  int icur,ii,jj,kk;
  int level;
  struct OCT * nextoct;
  struct OCT oct;
  FILE *fp;


  fp=fopen(filename,"wb");

  //printf("tsim=%f\n",tsim);
  fwrite(&tsim,sizeof(REAL),1,fp); 

  //printf("==> start map \n");
  
  /* // dumping the offsets and pointers informations */

  /* fwrite(&levelmax,sizeof(int),1,fp); */

  /* for(level=1;level<=levelmax;level++){ */
  /*   fwrite(firstoct+level-1,sizeof(struct OCT*),1,fp); */
  /*   printf("%p ",firstoct[level-1]); */
  /* }   */
  /* printf("\n"); */

  /* for(level=1;level<=levelmax;level++){ */
  /*   fwrite(lastoct+level-1,sizeof(struct OCT*),1,fp); */
  /*   printf("%p ",lastoct[level-1]); */
  /* } */
  /* printf("\n"); */

  // dumping the zero oct

  //printf("%p\n",firstoct[0]);
  fwrite(&(firstoct[0]),sizeof(struct OCT*),1,fp);
  printf("size of the OCT= %ld\n",sizeof(struct OCT));

  for(level=1;level<=levelmax;level++) // looping over octs
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
	}while(nextoct!=NULL);
    }

  

 
  fclose(fp);
}


//====================================================================================================
//=================================================================================================

  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
#ifdef PIC
void dumppart(struct OCT **firstoct,char filename[],int npart, int levelcoarse, int levelmax, REAL tsim){

  FILE *fp;
  float val;
  int vali;
  int ipart=0;
  int level;
  struct OCT *nextoct;
  struct OCT oct;
  REAL dxcur;
  struct PART *nexp;
  struct PART *curp;
  int icell;
  float tsimf=tsim;

  fp=fopen(filename,"wb");
  fwrite(&npart,1,sizeof(int),fp);
  fwrite(&tsimf,1,sizeof(float),fp);
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
#ifdef PART_EGY
		  val=curp->ekin+curp->epot;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;fwrite(&val,1,sizeof(float),fp);
#else
		  val=curp->vx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vy;fwrite(&val,1,sizeof(float),fp);
		  val=curp->vz;fwrite(&val,1,sizeof(float),fp);
#endif
#else
		  val=curp->fx;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fy;fwrite(&val,1,sizeof(float),fp);
		  val=curp->fz;fwrite(&val,1,sizeof(float),fp);
#endif
		  val=(REAL)(curp->idx);fwrite(&val,1,sizeof(float),fp);
		  ipart++;
		}while(nexp!=NULL);
	      }
	    }
	}while(nextoct!=NULL);
    }
  fclose(fp);
  printf("wrote %d particles (%d expected) in %s\n",ipart,npart,filename);
}
#endif

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

//------------------------------------------------------------------------

void GetParameters(char *fparam, struct RUNPARAMS *param)
{
  FILE *buf; 
  char stream[256];
  size_t rstat;
  REAL dummyf;

  buf=fopen(fparam,"r");
  if(buf==NULL)
    {
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
    }
  else
    {
      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ngridmax);
      rstat=fscanf(buf,"%s %d",stream,&param->npartmax);
      rstat=fscanf(buf,"%s %d",stream,&param->nbuff);
      
      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ndumps);
      rstat=fscanf(buf,"%s %d",stream,&param->nsteps);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->dt=dummyf;

      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->lcoarse);
      rstat=fscanf(buf,"%s %d",stream,&param->lmax);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->amrthresh=dummyf;

      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->niter);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->poissonacc=dummyf;
      rstat=fscanf(buf,"%s %d",stream,&param->mgridlmin);
      rstat=fscanf(buf,"%s %d",stream,&param->nvcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nrelax);

      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->nrestart);

      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->stride);
      rstat=fscanf(buf,"%s %d",stream,&param->nsubcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nthread);
      rstat=fscanf(buf,"%s %d",stream,&param->nstream);

#ifdef WRAD
      fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->clight=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->srcthresh=dummyf;
      rstat=fscanf(buf,"%s %lf",stream,&dummyf);param->srcint=dummyf;
#endif
      fclose(buf);
    }

  // computing the maxhash
  int val=(pow(2,param->lmax-1)<64?pow(2,param->lmax-1):64); // limit to 2097152 octs in hash table i.e. 16e6 cells
  param->maxhash=pow(val,3);
  printf("maxhash=%d\n",param->maxhash);
}

//==================================================================================
//==================================================================================
#ifdef PIC
#ifdef GRAFIC
struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param)
{
  FILE *fx, *fy, *fz;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
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


  if(np1*np2*np3>param->npartmax){
    printf("Error : Number of particles greater than npartmax Np(grafic)=%d Npartmax=%d!\n",np1*np2*np3,param->npartmax);
    abort();
  }

  //setting omegab

  ob=OMEGAB;


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

#ifdef WHYDRO2
  mass=(1.-ob/om)/(np1*np2*np3);
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
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=ip;


  printf("Grafic Particle Read ok\n");
  return lastpart;
}
#endif
#endif
//==================================================================================
//==================================================================================

#ifdef WHYDRO2
int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, struct RUNPARAMS *param){
  
  FILE *fx;
  FILE *fy;
  FILE *fz;
  FILE *fdx;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,ob,h0;
  int dummy;
  struct PART *lastpart;
  int ip;
  struct Wtype W;

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


  // setting baryon density parameter
  ob=OMEGAB;
  
  printf("============================================\n");
  printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
  printf("om=%f ov=%f ob=%f h0=%f\n",om,ov,ob,h0);
  printf("dx=%f np1*dx=%f\n",dx,np1*dx);
  printf("astart=%f zstart=%f\n",astart,1./astart-1.);
  printf("============================================\n");


  if(np1!=(int)pow(2,cpu->levelcoarse)){
    printf("ERROR !ABORT! Grafic file not compliant with parameter file : ngrafic=%d nquartz=%d\n",np1,(int)pow(2,cpu->levelcoarse));
    abort();
  }


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

  // ---------- Setting the initial temperature ---- //

  //  double temp=550.*((1.0+zstart)*(1.0+zstart)); // baryon temperature (to check) in K
  //double temp=2.7*(1+zstart);
  //double temp=1e4;
  //double temp=170.*(1.+zstart)*(1.+zstart)/10000.;
  double temp=0.0874545+0.0302621*zstart+0.00675076*zstart*zstart; // recfast ob fit


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
  vstar=rstar/tstar; //m/s
  pstar=rhostar*vstar*vstar;

  printf("rhoc=%e temperature=%lf rstar=%e(%e) pstar=%e tstar=%e vstar=%e rhostar=%e\n",rhoc,temp,rstar,np1*dx,pstar,tstar,vstar,rhostar);

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
	
	  rhob=(deltab[i1+i2*np1]+1.0)*ob*rhoc/pow(astart,3); // physical baryon density in kg/m3
	  pressure=(GAMMA-1.0)*1.5*(rhob/(mu*mp))*kboltz*temp; // physical pressure

	  // filling the cells using supercomoving values
	  
	  //abort();

	  W.d=(deltab[i1+i2*np1]+1.0)*ob/om;
	  W.u=(velx[i1+i2*np1]*1e3)*astart/vstar; // vstar is expressed in m/s and grafic vel are in km/s
	  W.v=(vely[i1+i2*np1]*1e3)*astart/vstar;
	  W.w=(velz[i1+i2*np1]*1e3)*astart/vstar;
	  W.p=pressure/pstar*pow(astart,5);
	  W.a=sqrt(GAMMA*W.p/W.d);
	  getE(&W);

	  memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	  ifound++;
	}
	else{
	  printf("euh pas trouve!");
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

  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  param->cosmo->H0=h0;

#ifdef WRAD
  param->unit.unit_l=rstar;
  param->unit.unit_v=vstar;
  param->unit.unit_t=param->unit.unit_l/param->unit.unit_v;
  param->unit.unit_n=1.;
  param->unit.unit_mass=rhostar*pow(param->unit.unit_l,3);
#endif

  
  printf("Grafic hydro read ok\n");
  return ifound;
}

#endif


