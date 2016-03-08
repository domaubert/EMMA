#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prototypes.h"
#include "friedmann.h"
#include "segment.h"
#include "hydro_utils.h"


#ifdef PIC
//==================================================================================
#ifdef GRAFIC
//==================================================================================
struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param,int level)
{
  FILE *fx = NULL;
  FILE *fy = NULL;
  FILE *fz = NULL;

#ifndef EMMAZELDO
  FILE *fpx = NULL;
  FILE *fpy = NULL;
  FILE *fpz = NULL;
#endif

  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
  int dummy;
  struct PART *lastpart=NULL;
  int ip;
  char filename[256];

  int res_level = level+param->DM_res;
  if(cpu->rank==0){

    sprintf(filename,"./level_%03d/ic_velbx",level);
    fx=fopen(filename,"rb");
    if(fx == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_velcx",level);
      printf("trying %s\n", filename);

      fx=fopen(filename,"rb");
      if(fx == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }

    sprintf(filename,"./level_%03d/ic_velby",level);
    fy=fopen(filename,"rb");
    if(fy == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_velcy",level);
      printf("trying %s\n", filename);

      fy=fopen(filename,"rb");
      if(fy == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }

    sprintf(filename,"./level_%03d/ic_velbz",level);
    fz=fopen(filename,"rb");
    if(fz == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_velcz",level);
      printf("trying %s\n", filename);

      fz=fopen(filename,"rb");
      if(fz == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }


#ifndef EMMAZELDO

    printf("ICS: READING DISPLACEMENTS FROM FILE\n");

    sprintf(filename,"./level_%03d/ic_poscz",level);
    fpz=fopen(filename,"rb");
    if(fpz == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_poscz",level);
      printf("trying %s\n", filename);

      fpz=fopen(filename,"rb");
      if(fpz == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }

    sprintf(filename,"./level_%03d/ic_poscz",level);
    fpz=fopen(filename,"rb");
    if(fpz == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_poscz",level);
      printf("trying %s\n", filename);

      fpz=fopen(filename,"rb");
      if(fpz == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }

    sprintf(filename,"./level_%03d/ic_poscy",level);
    fpy=fopen(filename,"rb");
    if(fpy == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_poscy",level);
      printf("trying %s\n", filename);

      fpy=fopen(filename,"rb");
      if(fpy == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }

    sprintf(filename,"./level_%03d/ic_poscx",level);
    fpx=fopen(filename,"rb");
    if(fpx == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_poscx",level);
      printf("trying %s\n", filename);

      fpx=fopen(filename,"rb");
      if(fpx == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }
#else
    printf("ICS: APPLYING ZELDOVICH DISPLACEMENT\n");
#endif


    // reading the headers
    size_t outf;

    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&np1,1,4,fx);
    outf=fread(&np2,1,4,fx);
    outf=fread(&np3,1,4,fx);
    outf=fread(&dx,1,4,fx);
    printf("DX=%e\n",dx);
    outf=fread(&x1o,1,4,fx);
    outf=fread(&x2o,1,4,fx);
    outf=fread(&x3o,1,4,fx);
    outf=fread(&astart,1,4,fx);
    outf=fread(&om,1,4,fx);
    outf=fread(&ov,1,4,fx);
    outf=fread(&h0,1,4,fx);
    outf=fread(&dummy,1,sizeof(dummy),fx);

    outf=fread(&dummy,1,sizeof(dummy),fy);
    outf=fread(&np1,1,4,fy);
    outf=fread(&np2,1,4,fy);
    outf=fread(&np3,1,4,fy);
    outf=fread(&dx,1,4,fy);
    outf=fread(&x1o,1,4,fy);
    outf=fread(&x2o,1,4,fy);
    outf=fread(&x3o,1,4,fy);
    outf=fread(&astart,1,4,fy);
    outf=fread(&om,1,4,fy);
    outf=fread(&ov,1,4,fy);
    outf=fread(&h0,1,4,fy);
    outf=fread(&dummy,1,sizeof(dummy),fy);

    outf=fread(&dummy,1,sizeof(dummy),fz);
    outf=fread(&np1,1,4,fz);
    outf=fread(&np2,1,4,fz);
    outf=fread(&np3,1,4,fz);
    outf=fread(&dx,1,4,fz);
    outf=fread(&x1o,1,4,fz);
    outf=fread(&x2o,1,4,fz);
    outf=fread(&x3o,1,4,fz);
    outf=fread(&astart,1,4,fz);
    outf=fread(&om,1,4,fz);
    outf=fread(&ov,1,4,fz);
    outf=fread(&h0,1,4,fz);
    outf=fread(&dummy,1,sizeof(dummy),fz);

#ifndef EMMAZELDO
    outf=fread(&dummy,1,sizeof(dummy),fpx);
    outf=fread(&np1,1,4,fpx);
    outf=fread(&np2,1,4,fpx);
    outf=fread(&np3,1,4,fpx);
    outf=fread(&dx,1,4,fpx);
    outf=fread(&x1o,1,4,fpx);
    outf=fread(&x2o,1,4,fpx);
    outf=fread(&x3o,1,4,fpx);
    outf=fread(&astart,1,4,fpx);
    outf=fread(&om,1,4,fpx);
    outf=fread(&ov,1,4,fpx);
    outf=fread(&h0,1,4,fpx);
    outf=fread(&dummy,1,sizeof(dummy),fpx);

    outf=fread(&dummy,1,sizeof(dummy),fpy);
    outf=fread(&np1,1,4,fpy);
    outf=fread(&np2,1,4,fpy);
    outf=fread(&np3,1,4,fpy);
    outf=fread(&dx,1,4,fpy);
    outf=fread(&x1o,1,4,fpy);
    outf=fread(&x2o,1,4,fpy);
    outf=fread(&x3o,1,4,fpy);
    outf=fread(&astart,1,4,fpy);
    outf=fread(&om,1,4,fpy);
    outf=fread(&ov,1,4,fpy);
    outf=fread(&h0,1,4,fpy);
    outf=fread(&dummy,1,sizeof(dummy),fpy);

    outf=fread(&dummy,1,sizeof(dummy),fpz);
    outf=fread(&np1,1,4,fpz);
    outf=fread(&np2,1,4,fpz);
    outf=fread(&np3,1,4,fpz);
    outf=fread(&dx,1,4,fpz);
    outf=fread(&x1o,1,4,fpz);
    outf=fread(&x2o,1,4,fpz);
    outf=fread(&x3o,1,4,fpz);
    outf=fread(&astart,1,4,fpz);
    outf=fread(&om,1,4,fpz);
    outf=fread(&ov,1,4,fpz);
    outf=fread(&h0,1,4,fpz);
    outf=fread(&dummy,1,sizeof(dummy),fpz);


#endif


  }

#ifdef WMPI
  // Massive broadcast of grafic header
  MPI_Bcast(&np1,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np2,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np3,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&dx,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x1o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x2o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x3o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&astart,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&om,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&ov,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&h0,1,MPI_FLOAT,0,cpu->comm);
  MPI_Barrier(cpu->comm);
#endif

  if(cpu->rank==0){
    printf("============================================\n");
    printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
    printf("om=%f ov=%f h0=%f\n",om,ov,h0);
    printf("dx=%f np1*dx=%f\n",dx,np1*dx);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }

  if(level==param->lcoarse){
  if((np1*np2*np3)/(cpu->nproc)*1.>param->npartmax){
    printf("Error : Number of particles greater than npartmax Np(grafic, est. )=%d Npartmax=%d!\n",(np1*np2*np3)/(cpu->nproc),param->npartmax);
    abort();
  }
  }
  //setting omegab

  ob=OMEGAB;


  // computing Zeldovich displacement quantities

  double vfact;
  vfact=fomega(astart,om,ov)*h0*dladt(astart,om,ov)/astart;
  if(cpu->rank==0) printf("vfact=%f\n",vfact);
  // reading the grafic planes

  float *velx;
  float *vely;
  float *velz;

#ifndef EMMAZELDO
  float *dispx;
  float *dispy;
  float *dispz;

#endif

  double x;
  double y;
  double z;

  double vx;
  double vy;
  double vz;

  double x0,y0,z0;

  int i1,i2,i3;
  int offidx=0;
  int keep;
  double mass;
  double veloff=0.;

#ifdef WHYDRO2
  mass=(1.-ob/om)/(np1*np2*np3);
#else
  mass=1./(np1*np2*np3);
#endif

#ifdef ZOOM
  if(level>param->lcoarse){
    offidx=POW(2,3*(level-param->lcoarse)); // for ids of particles
  }
#endif

  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);

#ifndef EMMAZELDO
  dispx=(float*)malloc(sizeof(float)*np1*np2);
  dispy=(float*)malloc(sizeof(float)*np1*np2);
  dispz=(float*)malloc(sizeof(float)*np1*np2);
#endif

  //REAL rmin=2.;

  ip=0;
  size_t outf;
  for(i3=1;i3<=np3;i3++){

    if(cpu->rank==0){
      printf("\r %f percent done",(i3*1.0)/np3*100.);
      outf=fread(&dummy,1,sizeof(dummy),fx);
      outf=fread(velx,np1*np2,sizeof(float),fx);
      outf=fread(&dummy,1,sizeof(dummy),fx);

      outf=fread(&dummy,1,sizeof(dummy),fy);
      outf=fread(vely,np1*np2,sizeof(float),fy);
      outf=fread(&dummy,1,sizeof(dummy),fy);

      outf=fread(&dummy,1,sizeof(dummy),fz);
      outf=fread(velz,np1*np2,sizeof(float),fz);
      outf=fread(&dummy,1,sizeof(dummy),fz);

#ifndef EMMAZELDO
      outf=fread(&dummy,1,sizeof(dummy),fpx);
      outf=fread(dispx,np1*np2,sizeof(float),fpx);
      outf=fread(&dummy,1,sizeof(dummy),fpx);

      outf=fread(&dummy,1,sizeof(dummy),fpy);
      outf=fread(dispy,np1*np2,sizeof(float),fpy);
      outf=fread(&dummy,1,sizeof(dummy),fpy);

      outf=fread(&dummy,1,sizeof(dummy),fpz);
      outf=fread(dispz,np1*np2,sizeof(float),fpz);
      outf=fread(&dummy,1,sizeof(dummy),fpz);
#endif


    }


#ifdef WMPI
    MPI_Barrier(cpu->comm);
    MPI_Bcast(velx,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(vely,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velz,np1*np2,MPI_FLOAT,0,cpu->comm);

#ifndef EMMAZELDO
    MPI_Bcast(dispx,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(dispy,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(dispz,np1*np2,MPI_FLOAT,0,cpu->comm);
#endif

#endif


    z0=(i3-0.5)*dx;
    for(i2=1;i2<=np2;i2++){
      y0=(i2-0.5)*dx;
      for(i1=1;i1<=np1;i1++){
	x0=(i1-0.5)*dx;
	// computing the displacements
#ifdef EMMAZELDO
	x=(x0+velx[(i1-1)+(i2-1)*np1]/vfact)/(np1*dx);
	y=(y0+vely[(i1-1)+(i2-1)*np1]/vfact)/(np2*dx);
	z=(z0+velz[(i1-1)+(i2-1)*np1]/vfact)/(np3*dx);
#else
	x=(x0+dispx[(i1-1)+(i2-1)*np1]/(h0/100.))/(np1*dx);
	y=(y0+dispy[(i1-1)+(i2-1)*np1]/(h0/100.))/(np2*dx);
	z=(z0+dispz[(i1-1)+(i2-1)*np1]/(h0/100.))/(np3*dx);
#endif
	// periodic boundary conditions

	x+=((x<0.)-(x>1.))*1.;
	y+=((y<0.)-(y>1.))*1.;
	z+=((z<0.)-(z>1.))*1.;

	// ugly fix for huge config in SINGLE FLOAT precision
	// generally affects a tiny fraction of particle (like 1 over 1e7)
 	if(x>0.99999) x=0.;
	if(y>0.99999) y=0.;
	if(z>0.99999) z=0.;

	// computing the velocities
	vx=(velx[(i1-1)+(i2-1)*np1]+veloff/astart*1e-3)*astart/(np1*dx*h0)/(sqrt(om)*0.5);
	vy=(vely[(i1-1)+(i2-1)*np1])*astart/(np2*dx*h0)/(sqrt(om)*0.5);
	vz=(velz[(i1-1)+(i2-1)*np1])*astart/(np3*dx*h0)/(sqrt(om)*0.5);

	// if it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,cpu,cpu->levelcoarse)){

	  keep=1;
#ifdef ZOOM
	  // is the current particle at the correct level?
	  int lzoom;
	  lzoom=pos2levelzoom(x,y,z,param);
	  REAL rloc;
	  rloc=sqrt(((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)));

	  if(lzoom!=level){
	    keep=0;
	  }
 #endif

	  if(keep) {
 	    part[ip].x=x;
	    part[ip].y=y;
	    part[ip].z=z;

 	    //rmin=(rloc<rmin?rloc:rmin);

	    part[ip].vx=vx;
	    part[ip].vy=vy;
	    part[ip].vz=vz;
	    part[ip].level=level;

	    part[ip].mass=mass;
	    part[ip].idx=(i1-1)+(i2-1)*np1+(i3-1)*np1*np2+offidx;
	    lastpart=part+ip;
	    ip++;
	  }
	}
      }
    }
  }

  free(velx);
  free(vely);
  free(velz);

#ifndef EMMAZELDO
  free(dispx);
  free(dispy);
  free(dispz);
#endif

  if(cpu->rank==0){
    fclose(fx);
    fclose(fy);
    fclose(fz);
  }


  // supercomoving unit values
  double rhostar;
  double rstar;
  double vstar;
  double tstar;
  double tstar2;
  double pstar;
  double mpc=3.08568025e22; // Mpc in m
  double H0=h0*1e3/3.08568025e22; // Hubble constant (s-1)

  rstar= np1*dx*mpc; // box size in m


  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  param->cosmo->H0=h0;
  param->cosmo->unit_l=rstar;
  *npart=ip;



#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif

  if(cpu->rank==RANK_DISP)
  {
    printf("\nGrafic Particle Read ok\n");
  }
  return lastpart;
}
#endif // GRAFIC

//==================================================================================
#ifdef ZELDOVICH
//==================================================================================
struct PART * read_zeldovich_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param, struct OCT **firstoct)
{
  FILE *fd;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,h0,ob;
  int dummy;
  struct PART *lastpart;
  int ip;

  int nploc;
  float munit_z;
  float lbox;
  float ainit_z;
  REAL mass;
  size_t outf;

  fd=fopen("utils/grafic_src/ZEL.PM.0","rb");

  outf=fread(&dummy,sizeof(dummy),1,fd);
  outf=fread(&nploc,sizeof(int),1,fd);
  outf=fread(&munit_z,sizeof(float),1,fd);
  outf=fread(&ainit_z,sizeof(float),1,fd);
  outf=fread(&lbox,sizeof(float),1,fd);
  outf=fread(&om,sizeof(float),1,fd);
  outf=fread(&ov,sizeof(float),1,fd);
  outf=fread(&h0,sizeof(float),1,fd);
  outf=fread(&dummy,sizeof(dummy),1,fd);

  astart=ainit_z;


  if(cpu->rank==0){
    printf("============================================\n");
    printf("ntot%d\n",nploc);
    printf("om=%f ov=%f h0=%f\n",om,ov,h0);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }

  if((nploc)/(cpu->nproc)*1.>param->npartmax){
    printf("Error : Number of particles greater than npartmax Np(grafic, est. )=%d Npartmax=%d!\n",(nploc)/(cpu->nproc),param->npartmax);
    abort();
  }

  //setting omegab

  ob=OMEGAB;


#ifdef WHYDRO2
  mass=(1.-ob/om)/(nploc);
#else
  mass=1./(nploc);
#endif


  // reading the grafic planes
  float *pos;
  float *vel;
  int nread=nploc; // quick fixes
  int npatch=1.; // quick fixes
  int ipatch;
  int i;
  pos=(float *)malloc(sizeof(REAL)*3*nread);
  vel=(float *)malloc(sizeof(REAL)*3*nread);

  double x;
  double y;
  double z;

  double vx;
  double vy;
  double vz;
  size_t rstat;

  int pstart=ftell(fd);

  ip=0.;
  for(ipatch=0;ipatch<npatch;ipatch++) {
    //    rstat=outf=fread(&dummy,sizeof(dummy),1,fd);
    //    fseek(fd,pstart,SEEK_SET);
    fseek(fd,pstart+(0*nploc+ipatch*nread)*sizeof(float)+1*sizeof(dummy),SEEK_SET);
    outf=fread(pos,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(1*nploc+ipatch*nread)*sizeof(float)+3*sizeof(dummy),SEEK_SET);
    outf=fread(pos+nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(2*nploc+ipatch*nread)*sizeof(float)+5*sizeof(dummy),SEEK_SET);
    outf=fread(pos+2*nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(3*nploc+ipatch*nread)*sizeof(float)+7*sizeof(dummy),SEEK_SET);
    outf=fread(vel,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(4*nploc+ipatch*nread)*sizeof(float)+9*sizeof(dummy),SEEK_SET);
    outf=fread(vel+nread,sizeof(float),nread,fd);
    //    outf=fread(&dummy,sizeof(dummy),1,fd);

    //    outf=fread(&dummy,sizeof(dummy),1,fd);
    fseek(fd,pstart+(5*nploc+ipatch*nread)*sizeof(float)+11*sizeof(dummy),SEEK_SET);
    outf=fread(vel+2*nread,sizeof(float),nread,fd);
    //outf=fread(&dummy,sizeof(dummy),1,fd);


    for(i=0;i<nread;i++)
      {
	x=pos[i];
	y=pos[i+nread];
	z=pos[i+2*nread];
	vx=vel[i];
	vy=vel[i+nread];
	vz=vel[i+2*nread];
	// periodic boundary conditions

	x+=(x<0)*((int)(-x)+1)-(x>1.)*((int)x);
	y+=(y<0)*((int)(-y)+1)-(y>1.)*((int)y);
	z+=(z<0)*((int)(-z)+1)-(z>1.)*((int)z);
	// it it belongs to the current cpu, we proceed and assign the particle to the particle array
	if(segment_part(x,y,z,cpu,cpu->levelcoarse)){

	  part[ip].x=x;
	  part[ip].y=y;
	  part[ip].z=z;

	  part[ip].vx=vx;
	  part[ip].vy=vy;
	  part[ip].vz=vz;

	  part[ip].mass=mass;
	  part[ip].idx=i;
	  lastpart=part+ip;
	  ip++;
	}
      }
  }


  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=ip;

  if(cpu->rank==RANK_DISP){
    printf("Zeldovich Particle Read ok\n");
  }

#ifdef WHYDRO2

  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;
  REAL ZC=10; // hard coded collapse of redshift
  REAL ZI=1./astart-1.;

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
 		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		W.d=(1.+(1.+ZC)/(1.+ZI)*cos(2.*M_PI*(xc-0.5)))*ob/om;
		W.p=PMIN;
		W.u=-(1.+ZC)/POW(1.+ZI,1.5)*sin(2.*M_PI*(xc-0.5))/(M_PI); // for omegam=1. only
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);
		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

#endif

  return lastpart;
}
#endif // ZELDOVICH


//==================================================================================
#ifdef EDBERT
//==================================================================================
struct PART * read_edbert_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param, struct OCT **firstoct)
{
  float astart,om,ov,h0,ob;
  int dummy;

  struct PART *lastpart;
  int ip,iploc;
  int i,j,k;

  float lbox;
  float ainit_z;
  REAL mass;
  REAL x,y,z,r;
  REAL delta;
  REAL lsphere;
  printf("Start EDBERT\n");

  lsphere=0.05;
  om=0.99999;
  ov=0.00001;
  ob=OMEGAB;
  delta=0.2;
  h0=70.;
  lbox=1.;
  astart=1e-3;

  int n1d=POW(2,param->lcoarse);
  REAL dx=1./n1d;
  iploc=0;
  int nin=0;
  int nout=0;
  REAL m;

#ifdef PIC

  REAL mout,mint;
  REAL vsphere=0.;
  for(k=0;k<n1d;k++)
    {
      for(j=0;j<n1d;j++)
	{
	  for(i=0;i<n1d;i++)
	    {
	      x=(i+0.5)*dx;
	      y=(j+0.5)*dx;
	      z=(k+0.5)*dx;

	      r=sqrt(POW(x-0.5,2)+POW(y-0.5,2)+POW(z-0.5,2));
	      if(r<lsphere){
		nin++;
		mout=-1;
		vsphere+=POW(dx,3);
	      }
	      else{
		nout++;
		mout=1;
	      }

	      if(segment_part(x,y,z,cpu,cpu->levelcoarse)){
		part[iploc].x=x;
		part[iploc].y=y;
		part[iploc].z=z;

		part[iploc].vx=0;
		part[iploc].vy=0;
		part[iploc].vz=0;

		part[iploc].mass=mout;
		part[iploc].idx=-nin;
		lastpart=part+iploc;
		iploc++;
	      }
	    }
	}
    }

  mint=(om-ob)*(1.+delta)*vsphere/nin;
  mout=((om-ob)-mint*nin)/nout;
  printf("mint=%e mout=%e\n",mint,mout);

  for(i=0;i<iploc;i++){
    part[i].mass=(part[i].mass<0?mint:mout);
  }

#endif

  *munit=mass;
  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  //param->cosmo->Hubble=h0;
  *npart=iploc;

  if(cpu->rank==RANK_DISP){
    printf("Edbert Particle Read ok\n");
  }

#ifdef WHYDRO2

  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		rad=sqrt((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5));
		if(rad<lsphere){
		  W.d=(1.+delta)*ob;
		}
		else{
		  W.d=ob*(1.-(1.+delta)*vsphere)/(1.-vsphere);
		}
		W.p=PMIN;
		W.u=0.; // vstar is expressed in m/s and grafic vel are in km/s
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);

		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

#endif // WHYDRO2

  return lastpart;
}
#endif // EDBERT


#endif // PIC

#ifdef WHYDRO2

//==================================================================================
#ifdef TUBE
//==================================================================================

void read_shocktube(struct CPUINFO *cpu, REAL *ainit, struct RUNPARAMS *param, struct OCT **firstoct)
{
  FILE *fd = NULL;

  int level;
  REAL dxcur;
  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;
  REAL xc,yc,zc;
  struct Wtype W;
  REAL rad;


  struct Wtype WL, WR;
  if(cpu->rank==RANK_DISP) printf("Init Hydro\n");

  /* /\*  /\\* // TEST 1 *\\/ *\/ */

  /* WL.d=1.; */
  /* WL.u=0.; */
  /* WL.v=0.; */
  /* WL.w=0.; */
  /* WL.p=1.0; */
  /* WL.a=sqrt(GAMMA*WL.p/WL.d); */
  /* getE(&WL); */

  /* WR.d=0.125; */
  /* WR.u=0.; */
  /* WR.v=0.; */
  /* WR.w=0.; */
  /* WR.p=0.1; */
  /* WR.a=sqrt(GAMMA*WR.p/WR.d); */
  /* getE(&WR); */


  /* REAL X0=0.3125; */

  // SEDOV

  WL.d=1.;
  WL.u=0.;
  WL.v=0.;
  WL.w=0.;
  WL.p=1.0;
  WL.a=sqrt(GAMMA*WL.p/WL.d);
  getE(&WL);

  WR.d=0.125;
#ifdef SED
  WR.d=1.;
#endif
  WR.u=0.;
  WR.v=0.;
  WR.w=0.;
  WR.p=1e5;
  WR.a=sqrt(GAMMA*WR.p/WR.d);
  getE(&WR);

  REAL X0=1./128;

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    if(curoct->cpu!=cpu->rank) continue;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		REAL rc=xc;
#ifdef SED
		rc=SQRT(POW(xc-0.5,2)+POW(yc-0.5,2)+POW(zc-0.5,2));
#endif

		if((rc<X0)*(xc>0.5)*(yc>0.5)*(zc>0.5)){
		  memcpy(&(curoct->cell[icell].field),&WR,sizeof(struct Wtype));
		  printf("coucou\n");
		}
		else{
		  memcpy(&(curoct->cell[icell].field),&WL,sizeof(struct Wtype));
		}
	      }
	  }while(nextoct!=NULL);
      }
}
#endif



//==================================================================================
#ifdef EVRARD
//==================================================================================
int read_evrard_hydro(struct CPUINFO *cpu,struct OCT **firstoct, struct RUNPARAMS *param){

  int level;
  REAL dxcur;
  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  struct Wtype W;
  REAL rad;
  REAL xc,yc,zc;

  //==== parameters of evrard sphere
  REAL R=0.35;
  REAL M=1.;
  REAL rhostar=M/(4./3.*M_PI*R*R*R);
  REAL estar=M/R; //assuming G=1
  REAL pstar=rhostar*estar;
  REAL tstar=sqrt(M_PI*M_PI/8.)*POW(R,1.5)/POW(M,0.5);
  if(cpu->rank==RANK_DISP) printf("Generating Evrard Test Case ts=%e, rhostar=%e\n",tstar,rhostar);

  for(level=param->lcoarse;level<=param->lcoarse;level++) // (levelcoarse only for the moment)
      {
	dxcur=POW(0.5,level);
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
		yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
		zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;

		rad=sqrt((xc-0.5)*(xc-0.5)+(yc-0.5)*(yc-0.5)+(zc-0.5)*(zc-0.5))/R;
		if(rad<1.){
		  W.d=rhostar/rad;
		  W.p=pstar/rad*0.05;
		}
		else{
		  W.d=1e-3;
		  W.p=1e-5;
		}

		W.u=0.; // vstar is expressed in m/s and grafic vel are in km/s
		W.v=0.;
		W.w=0.;
		W.a=sqrt(GAMMA*W.p/W.d);
		getE(&W);

		memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	      }
	  }while(nextoct!=NULL);
      }

  }
#endif


#ifdef TESTCOSMO
//==================================================================================
#ifdef GRAFIC
//==================================================================================
 int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, struct RUNPARAMS *param, int level){

  FILE *fx = NULL;
  FILE *fy = NULL;
  FILE *fz = NULL;
  FILE *fdx = NULL;
  int np1,np2,np3;
  float dx,x1o,x2o,x3o,astart,om,ov,ob,h0;
  int dummy;
  int ip;
  struct Wtype W;
  size_t outf;
  char filename[256];

  // Note only the rank 0 reads the file.

  if(cpu->rank==0){

    sprintf(filename,"./level_%03d/ic_deltab",level);
    fdx=fopen(filename,"rb");		if(fdx == NULL) {printf("Cannot open %s\n", filename); abort();}

    sprintf(filename,"./level_%03d/ic_velbx",level);
    fx=fopen(filename,"rb");
    if(fx == NULL) {
      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_velcx",level);
      printf("trying %s\n", filename);
      fx=fopen(filename,"rb");
      if(fx == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }

    sprintf(filename,"./level_%03d/ic_velby",level);
    fy=fopen(filename,"rb");
    if(fy == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_velcy",level);
      printf("trying %s\n", filename);

      fy=fopen(filename,"rb");
      if(fy == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }

    sprintf(filename,"./level_%03d/ic_velbz",level);
    fz=fopen(filename,"rb");
    if(fz == NULL) {

      printf("Cannot open %s : ", filename);
      sprintf(filename,"./level_%03d/ic_velcz",level);
      printf("trying %s\n", filename);

      fz=fopen(filename,"rb");
      if(fz == NULL) {
        printf("Cannot open %s\n", filename);
        abort();
      }
    }

    /* fdx=fopen("./ic_deltab","rb"); */
    /* fx=fopen("./ic_velcx","rb"); */
    /* fy=fopen("./ic_velcy","rb"); */
    /* fz=fopen("./ic_velcz","rb"); */

    // reading the headers

    outf=fread(&dummy,1,sizeof(dummy),fdx);
    outf=fread(&np1,1,4,fdx);
    outf=fread(&np2,1,4,fdx);
    outf=fread(&np3,1,4,fdx);
    outf=fread(&dx,1,4,fdx);
    outf=fread(&x1o,1,4,fdx);
    outf=fread(&x2o,1,4,fdx);
    outf=fread(&x3o,1,4,fdx);
    outf=fread(&astart,1,4,fdx);
    outf=fread(&om,1,4,fdx);
    outf=fread(&ov,1,4,fdx);
    outf=fread(&h0,1,4,fdx);
    outf=fread(&dummy,1,sizeof(dummy),fdx);

    outf=fread(&dummy,1,sizeof(dummy),fx);
    outf=fread(&np1,1,4,fx);
    outf=fread(&np2,1,4,fx);
    outf=fread(&np3,1,4,fx);
    outf=fread(&dx,1,4,fx);
    outf=fread(&x1o,1,4,fx);
    outf=fread(&x2o,1,4,fx);
    outf=fread(&x3o,1,4,fx);
    outf=fread(&astart,1,4,fx);
    outf=fread(&om,1,4,fx);
    outf=fread(&ov,1,4,fx);
    outf=fread(&h0,1,4,fx);
    outf=fread(&dummy,1,sizeof(dummy),fx);

    outf=fread(&dummy,1,sizeof(dummy),fy);
    outf=fread(&np1,1,4,fy);
    outf=fread(&np2,1,4,fy);
    outf=fread(&np3,1,4,fy);
    outf=fread(&dx,1,4,fy);
    outf=fread(&x1o,1,4,fy);
    outf=fread(&x2o,1,4,fy);
    outf=fread(&x3o,1,4,fy);
    outf=fread(&astart,1,4,fy);
    outf=fread(&om,1,4,fy);
    outf=fread(&ov,1,4,fy);
    outf=fread(&h0,1,4,fy);
    outf=fread(&dummy,1,sizeof(dummy),fy);

    outf=fread(&dummy,1,sizeof(dummy),fz);
    outf=fread(&np1,1,4,fz);
    outf=fread(&np2,1,4,fz);
    outf=fread(&np3,1,4,fz);
    outf=fread(&dx,1,4,fz);
    outf=fread(&x1o,1,4,fz);
    outf=fread(&x2o,1,4,fz);
    outf=fread(&x3o,1,4,fz);
    outf=fread(&astart,1,4,fz);
    outf=fread(&om,1,4,fz);
    outf=fread(&ov,1,4,fz);
    outf=fread(&h0,1,4,fz);
    outf=fread(&dummy,1,sizeof(dummy),fz);
  }

  // setting baryon density parameter
  ob=OMEGAB;

#ifdef WMPI
  // Massive broadcast of grafic header
  MPI_Bcast(&np1,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np2,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&np3,1,MPI_INT,0,cpu->comm);
  MPI_Bcast(&dx,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x1o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x2o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&x3o,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&astart,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&om,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&ov,1,MPI_FLOAT,0,cpu->comm);
  MPI_Bcast(&h0,1,MPI_FLOAT,0,cpu->comm);
  MPI_Barrier(cpu->comm);

#endif

  if(cpu->rank==0){
    printf("============================================\n");
    printf("nx=%d ny=%d nz=%d\n",np1,np2,np3);
    printf("om=%f ov=%f ob=%f h0=%f\n",om,ov,ob,h0);
    printf("dx=%f np1*dx=%f\n",dx,np1*dx);
    printf("astart=%f zstart=%f\n",astart,1./astart-1.);
    printf("============================================\n");
  }



  if(np1!=(int)POW(2,level)){
    printf("ERROR !ABORT! Grafic hydro  file not compliant with parameter file : ngrafic=%d nquartz=%d\n",np1,(int)POW(2,level));
    abort();
  }


  // reading the grafic planes

  float *deltab;
  float *velz;
  float *vely;
  float *velx;

  int i1,i2,i3;
  int icx,icy,icz,icell;
  unsigned long long key;
  struct OCT *curoct;
  struct OCT *nextoct;
  unsigned long hidx;
  int found;
  float z0,y0,x0;
  int ifound=0;


  deltab=(float*)malloc(sizeof(float)*np1*np2);
  velx=(float*)malloc(sizeof(float)*np1*np2);
  vely=(float*)malloc(sizeof(float)*np1*np2);
  velz=(float*)malloc(sizeof(float)*np1*np2);

  double rhob,pressure;
  double H0=h0*1e3/3.08568025e22; // Hubble constant (s-1)
  double rhoc=3.*H0*H0/(8.*M_PI*NEWTON_G); // comoving critical density (kg/m3)
  double zstart=1./astart-1.;

  // ---------- Setting the initial temperature ---- //
  double temp;
  //  double temp=550.*((1.0+zstart)*(1.0+zstart)); // baryon temperature (to check) in K
  //double temp=2.7*(1+zstart);
  //double temp=1e4;
  //double temp=170.*(1.+zstart)*(1.+zstart)/10000.;

  //double temp=0.0874545+0.0302621*zstart+0.00675076*zstart*zstart; // recfast ob fit
#ifdef COOLING
  if(om==1.) {
    temp=33.64/POW(41.,2)*POW(1.+zstart,2);
    if(cpu->rank==RANK_DISP) printf("WARNING: YOU ARE USING SCDM COSMOLOGY\n");
  }
  else{
    if(cpu->rank==RANK_DISP) printf("No temperature law for cosmologies other than SCDM -> F** it\n");
    temp=33.64/POW(41.,2)*POW(1.+zstart,2);
    //    abort();
  }
#else
  temp=1e4;
#endif

  // supercomoving unit values
  double rhostar;
  double rstar;
  double vstar;
  double tstar;
  double tstar2;
  double pstar;
  double mpc=3.08568025e22; // Mpc in m
  double veloff=0.;
#ifdef BULKFLOW
  veloff=VBC; // relative motion in km/s at z=999
#endif

  rstar= np1*dx*mpc; // box size in m
  rhostar=rhoc*om;
  tstar=2./H0/sqrt(om); // sec
  tstar2=2./h0/sqrt(om); // Mpc sec / km
  vstar=rstar/tstar; //m/s
  pstar=rhostar*vstar*vstar;

  if(cpu->rank==RANK_DISP) printf("rhoc=%e temperature=%lf rstar=%e(%e) pstar=%e tstar=%e vstar=%e rhostar=%e\n",rhoc,temp,rstar,np1*dx,pstar,tstar,vstar,rhostar);


  for(i3=0;i3<np3;i3++){

    if(cpu->rank==0){
      printf("\r %f percent done",(i3*1.0)/np3*100.);

      outf=fread(&dummy,1,sizeof(dummy),fdx);
      outf=fread(deltab,np1*np2,sizeof(float),fdx);
      outf=fread(&dummy,1,sizeof(dummy),fdx);

      outf=fread(&dummy,1,sizeof(dummy),fx);
      outf=fread(velx,np1*np2,sizeof(float),fx);
      outf=fread(&dummy,1,sizeof(dummy),fx);

      outf=fread(&dummy,1,sizeof(dummy),fy);
      outf=fread(vely,np1*np2,sizeof(float),fy);
      outf=fread(&dummy,1,sizeof(dummy),fy);

      outf=fread(&dummy,1,sizeof(dummy),fz);
      outf=fread(velz,np1*np2,sizeof(float),fz);
      outf=fread(&dummy,1,sizeof(dummy),fz);
    }

#ifdef WMPI
    MPI_Barrier(cpu->comm);
    MPI_Bcast(deltab,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velx,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(vely,np1*np2,MPI_FLOAT,0,cpu->comm);
    MPI_Bcast(velz,np1*np2,MPI_FLOAT,0,cpu->comm);
#endif

    z0=(i3*1.0)/(np3);
    for(i2=0;i2<np2;i2++){
      y0=(i2*1.0)/(np2);
      for(i1=0;i1<np1;i1++){
	x0=(i1*1.0)/(np1);

	key=pos2key(x0,y0,z0,level);

	// first we compute the adress from the hashfunction
	hidx=hfun(key,cpu->maxhash);
	nextoct=cpu->htable[hidx];

	// looking for the oct
	found=0;
	if(nextoct!=NULL){
	  do{ // resolving collisions
	    curoct=nextoct;
	    nextoct=curoct->nexthash;
	    found=((oct2key(curoct,level)==key)&&(curoct->level==level));
	  }while((nextoct!=NULL)&&(!found));
	}

	// filling the cell

	if(found){
	  icx=i1%2;
	  icy=i2%2;
	  icz=i3%2;

	  icell=icx+icy*2+icz*4;

	  rhob=(deltab[i1+i2*np1]+1.0)*ob*rhoc/POW(astart,3); // physical baryon density in kg/m3
	  pressure=(GAMMA-1.0)*1.5*(rhob/(PROTON_MASS*MOLECULAR_MU))*KBOLTZ*temp; // physical pressure

	  //printf("pres=%e\n",pressure);
	  // filling the cells using supercomoving values

	  //abort();

	  W.d=(deltab[i1+i2*np1]+1.0)*ob/om;
	  W.u=((velx[i1+i2*np1]+veloff)*1e3)*astart/vstar; // vstar is expressed in m/s and grafic vel are in km/s
	  W.v=(vely[i1+i2*np1]*1e3)*astart/vstar;
	  W.w=(velz[i1+i2*np1]*1e3)*astart/vstar;
	  W.p=pressure/pstar*POW(astart,5);
	  W.a=sqrt(GAMMA*W.p/W.d);
	  getE(&W);

	  if(W.p<PMIN) {printf(" YOU SHOULD RECONSIDER PMIN %e %e %e %e\n",W.p,pressure,pstar,PMIN);abort();}

#ifdef WRADHYD
	  // Testing ADVECTION
	  //W.X=(i1/6)%2+((i2+1)/6)%2;
	  W.dX=0.2e-3*W.d*(1.-YHE);
#ifdef HELIUM
	  W.dXHE=0.2e-3*W.d*(YHE)/yHE;
	  W.dXXHE=0.2e-3*W.d*(YHE)/yHE;
#endif
#endif
	  memcpy(&(curoct->cell[icell].field),&W,sizeof(struct Wtype));

	  ifound++;
	}
	else{

	  // this branch corresponds to cell out of the domain
	  //	  printf("euh pas trouve! hidx=%d %p",hidx,cpu->htable[hidx]);
	  //	  abort();
	}
      }
    }
  }

  if(cpu->rank==0){
    fclose(fdx);
    fclose(fx);
    fclose(fy);
    fclose(fz);
  }


  free(deltab);
  free(velx);
  free(vely);
  free(velz);

  *ainit=astart;
  param->cosmo->om=om;
  param->cosmo->ov=ov;
  param->cosmo->ob=ob;
  param->cosmo->H0=h0;
  param->cosmo->unit_l=rstar;


  int setunit=0;

#ifdef WHYDRO2
  setunit=1;
#endif

#ifdef WRAD
  setunit=1;
#endif
  if(setunit){
    param->unit.unit_l=rstar;
    param->unit.unit_v=vstar;
    param->unit.unit_t=param->unit.unit_l/param->unit.unit_v;
    param->unit.unit_n=1.;//(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS)*rhostar; //
    param->unit.unit_mass=rhostar*pow(param->unit.unit_l,3);
    param->unit.unit_d=rhostar; // kg/m3
    param->unit.unit_N=rhostar/PROTON_MASS; // atom/m3
  }

//  REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS*MOLECULAR_MU/param->unit.unit_mass)/POW(param->unit.unit_l,3);
 /*   REAL navg=(param->cosmo->ob/param->cosmo->om)/(PROTON_MASS*MOLECULAR_MU)*rhostar; */

/* if(cpu->rank==RANK_DISP) printf("navg=%e \n",navg); */

#ifdef WMPI
  MPI_Barrier(cpu->comm);
#endif
  if(cpu->rank==RANK_DISP) printf("Grafic hydro read ok\n");
  return ifound;
}
#endif // GRAFIC
#endif // TESTCOSMO
#endif // WHYDRO2

#ifdef WHYDRO2
int init_sedov(struct RUNPARAMS *param, struct OCT **firstoct){
/**
  * initialize the grid for sedov test
  * the grid is fill up with an uniform and unitary medium
  *
  * density is set to 1
  * velocity to 0
  * pressure to 1e-5
  **/

  param->unit.unit_l=0.125;
  param->unit.unit_t=1.;
  //param->unit.unit_v=param->unit.unit_l/param->unit.unit_t;
  param->unit.unit_v=1.;
  param->unit.unit_d=1.;
  param->unit.unit_N=1.;
  param->unit.unit_mass=1.;

  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){
    REAL dxcur=POW(0.5,level);
    struct OCT *nextoct=firstoct[level-1];
    if(nextoct==NULL) continue;
    do{
      struct OCT * curoct=nextoct;
      nextoct=curoct->next;
      int icell;
      for(icell=0;icell<8;icell++){
        struct CELL *curcell= &curoct->cell[icell];
        curcell->field.d=1.0;
        curcell->field.u=0.0;
        curcell->field.v=0.0;
        curcell->field.w=0.0;
        curcell->field.p=1e-5;
        curcell->field.a=SQRT(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d);
        getE(&(curcell->field));
      }
    }while(nextoct!=NULL);
  }
  return 0;
}
#endif // WHYDRO2

#ifdef WHYDRO2
int init_star_test(struct RUNPARAMS *param, struct OCT **firstoct){
/**
  * initialize the grid for star formation test
  * the grid is fill up with an uniform and unitary medium
  *
  * central cells are filled with an overdensity where stars can form
  **/

  param->unit.unit_l=0.125;
  param->unit.unit_t=1.;
  //param->unit.unit_v=param->unit.unit_l/param->unit.unit_t;
  param->unit.unit_v=1.;
  param->unit.unit_d=1.;
  param->unit.unit_N=1.;
  param->unit.unit_mass=1.;

  int level;
  for(level=param->lcoarse;level<=param->lmax;level++){
    REAL dxcur=POW(0.5,level);
    struct OCT *nextoct=firstoct[level-1];
    if(nextoct==NULL) continue;
    do{
      struct OCT * curoct=nextoct;
      nextoct=curoct->next;
      int icell;
      for(icell=0;icell<8;icell++){
        struct CELL *curcell= &curoct->cell[icell];
        curcell->field.d=1.0;
        curcell->field.u=0.0;
        curcell->field.v=0.0;
        curcell->field.w=0.0;
        curcell->field.p=1e-5;
        curcell->field.a=SQRT(GAMMA*curoct->cell[icell].field.p/curoct->cell[icell].field.d);
        getE(&(curcell->field));


        REAL xc=curoct->x+( icell&1)*dxcur+dxcur*0.5;
        REAL yc=curoct->y+((icell>>1)&1)*dxcur+dxcur*0.5;
        REAL zc=curoct->z+((icell>>2))*dxcur+dxcur*0.5;


        REAL r= SQRT((xc-0.5)*(xc-0.5) + (yc-0.5)*(yc-0.5) + (yc-0.5)*(yc-0.5));
        REAL dx = 1./POW2(param->lcoarse);
        if (r<dx) curcell->field.d=100;
      }
    }while(nextoct!=NULL);
  }
  return 0;
}
#endif // WHYDRO2
