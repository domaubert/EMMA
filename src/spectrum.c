#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
REAL cross_section(REAL egy){
/**
  * Compute cross section
  */
  if(egy >= 13.6){
    REAL P=2.963;
    REAL x=egy/4.298e-1;

    REAL a=POW(x-1.,2);
    REAL b=POW(x,P/2.-5.5);
    REAL c=POW(1+SQRT(x/3.288e1),-P);

    return 5.475e4*a*b*c*1e-22;

  }else {
    return 0;
  }
}

REAL eV2lambda(REAL eV){
/**
  * Convert electron Volt to wavelength in meter
  */
	REAL h = KPLANCK ;
	REAL c = LIGHT_SPEED_IN_M_PER_S;
	REAL Ej = eV*ELECTRONVOLT;
	return h*c/Ej;
}

REAL lambda2E(REAL l){
/**
  * Convert wavelength to Joule
  */
  REAL h = KPLANCK ;
	REAL c = LIGHT_SPEED_IN_M_PER_S;
	return h*c/l;
}

void integ_dwave(struct RUNPARAMS *param, REAL* y, REAL space_bound_min, REAL space_bound_max, REAL *out){
/**
  * Integration over wavelength
  */
  *out = 0;
  int ispace;
  for(ispace=0;ispace<param->spectrum.nwave;ispace++){
    if ((param->spectrum.wavelength[ispace]<=space_bound_min)
    &&  (param->spectrum.wavelength[ispace]>space_bound_max)) {

      REAL dx = param->spectrum.wavelength[ispace+1]-param->spectrum.wavelength[ispace];
      REAL y_mean= (y[ispace+1]+y[ispace])/2.;

      *out+=dx*y_mean;
    }
  }
}

void integ_dtime(struct RUNPARAMS *param, REAL* y, REAL time_bound_min, REAL time_bound_max, REAL *out){
/**
  * Integration over time
  */
  *out=0;
  int i_time;
  for(i_time=0;i_time<param->spectrum.ntime;i_time++){
    if ((param->spectrum.time[i_time]>time_bound_min)
    &&  (param->spectrum.time[i_time]<=time_bound_max)){

      REAL dx = param->spectrum.time[i_time+1]-param->spectrum.time[i_time];
      REAL y_mean = (y[i_time+1]+y[i_time])/2.;

//      printf("dx=%e dy=%e out=%e\n", dx, y_mean, *out);

      *out+=dx*y_mean;
    }
  }
}

struct GRP{
  REAL *s_N;
  REAL *s_E;
  REAL *Epho;
  REAL *nphot;
  REAL *factgrp;
};

void get_space_GRP(struct RUNPARAMS *param, int i_time, struct GRP *grp){
/**
  * Compute group caracteristic for a given time
  * Integration over space
  */
  grp->s_N =  (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  grp->s_E =  (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  grp->Epho = (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  grp->nphot= (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  grp->factgrp= (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));

  REAL *s_lambda = (REAL*)calloc(param->spectrum.nwave,sizeof(REAL));
  REAL *N_lambda = (REAL*)calloc(param->spectrum.nwave,sizeof(REAL));
  REAL *tmp_space= (REAL*)calloc(param->spectrum.nwave,sizeof(REAL));

  REAL nphot_tot=0;

  int i_grp_space;
  for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){

    REAL space_bound_min=eV2lambda(param->atomic.space_bound[i_grp_space]);
    REAL space_bound_max=eV2lambda(param->atomic.space_bound[i_grp_space+1]);

    //Eion
    REAL Eion = 0;
    integ_dwave(param, param->spectrum.flux[i_time], space_bound_min, space_bound_max, &Eion);

    //Nphot_per_sec
    REAL Nphot_per_sec=0;
    int ispace;
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      REAL E = lambda2E(param->spectrum.wavelength[ispace]);
      tmp_space[ispace]=param->spectrum.flux[i_time][ispace]/E;
    }
    integ_dwave(param, tmp_space, space_bound_min, space_bound_max, &Nphot_per_sec);

    //s_lambda
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      REAL E = lambda2E(param->spectrum.wavelength[ispace]);
      s_lambda[ispace] = cross_section(E/ELECTRONVOLT);
    }

    //N_lambda
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      REAL E = lambda2E(param->spectrum.wavelength[ispace]);
      N_lambda[ispace] = param->spectrum.flux[i_time][ispace]/E;
    }

    //s_N
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      tmp_space[ispace]=s_lambda[ispace]*N_lambda[ispace];
    }
    integ_dwave(param, tmp_space, space_bound_min, space_bound_max, &grp->s_N[i_grp_space]);
    grp->s_N[i_grp_space]/=Nphot_per_sec;

    //s_E
    for(ispace=0;ispace<param->spectrum.nwave;ispace++){
      tmp_space[ispace]=s_lambda[ispace]*param->spectrum.flux[i_time][ispace];
    }
    integ_dwave(param, tmp_space, space_bound_min, space_bound_max, &grp->s_E[i_grp_space]);
    grp->s_E[i_grp_space]/=Eion;

    grp->Epho[i_grp_space]= Eion/Nphot_per_sec;
    grp->nphot[i_grp_space] = Nphot_per_sec/(1e6*SOLAR_MASS);

    nphot_tot+=grp->nphot[i_grp_space];
  }

  for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){
    grp->factgrp[i_grp_space]=grp->nphot[i_grp_space]/nphot_tot;
  }

  free(s_lambda);
  free(N_lambda);
  free(tmp_space);
}

void get_time_GRP(struct RUNPARAMS *param, struct GRP *grp){
/**
  * Compute group caracteristic for a given wavelength
  * Integration over time
  */

  REAL *tmp_time= (REAL*)calloc(param->spectrum.ntime,sizeof(REAL));

  int i_grp_time;
  for(i_grp_time=0;i_grp_time<param->atomic.ngrp_time;i_grp_time++){

    REAL time_bound_min=param->atomic.time_bound[i_grp_time];
    REAL time_bound_max=param->atomic.time_bound[i_grp_time+1];
    REAL nphot_tot=0;

    int i_grp_space;
    for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){

      int id= i_grp_time*param->atomic.ngrp_space+i_grp_space;
      //printf("id=%d \n", id);

      //Nphot
      int itmp;
      for(itmp=0;itmp<param->spectrum.ntime;itmp++){
        tmp_time[itmp]= grp[itmp].nphot[i_grp_space];
  //      printf("grp[%d].nphot[%d]=%e \n", itmp,i_grp_space,tmp_time[itmp]);
      }
      REAL Nphot=0;
      integ_dtime(param, tmp_time, time_bound_min, time_bound_max, &Nphot);
      printf("Nphot=%e \n", Nphot);

//      nphot_tot+=Nphot;

      //E_tot
      for(itmp=0;itmp<param->spectrum.ntime;itmp++){
        tmp_time[itmp]=grp[itmp].Epho[i_grp_space]*grp[itmp].nphot[i_grp_space];
      }
      integ_dtime(param, tmp_time, time_bound_min, time_bound_max, &param->atomic.hnu[id]);
      param->atomic.hnu[id]/=Nphot;

//    printf("hnu=%e \n", param->atomic.hnu[id]/ELECTRONVOLT);

      //s_N_tot
      for(itmp=0;itmp<param->spectrum.ntime;itmp++){
        tmp_time[itmp]=grp[itmp].s_N[i_grp_space]*grp[itmp].nphot[i_grp_space];
      }
      integ_dtime(param, tmp_time, time_bound_min, time_bound_max, &param->atomic.alphai[id]);
      param->atomic.alphai[id]/=Nphot;
  //   printf("ai=%e \n", param->atomic.alphai[id]);

      //s_E_tot
      for(itmp=0;itmp<param->spectrum.ntime;itmp++){
        tmp_time[itmp]=grp[itmp].s_E[i_grp_space]*grp[itmp].nphot[i_grp_space];
      }
      integ_dtime(param, tmp_time, time_bound_min, time_bound_max, &param->atomic.alphae[id]);
      param->atomic.alphae[id]/=Nphot;
    //  printf("ae=%e \n", param->atomic.alphae[id]);

    }


//			factgrp = np.sum(nphot[igrp])/np.sum(nphot)



    int i_grp_time;
    for(i_grp_time=0;i_grp_time<param->atomic.ngrp_time;i_grp_time++){
      REAL time_bound_min=param->atomic.time_bound[i_grp_time];
      REAL time_bound_max=param->atomic.time_bound[i_grp_time+1];

      REAL nphot_tot=0;

      for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){

        int id= i_grp_time*param->atomic.ngrp_space+i_grp_space;

        int i_time;
        for(i_time=0;i_time<param->spectrum.ntime;i_time++){
          nphot_tot+=grp[i_time].nphot[i_grp_space];
        }
        for(i_time=0;i_time<param->spectrum.ntime;i_time++){
          if ((param->spectrum.time[i_time]>time_bound_min)
          &&  (param->spectrum.time[i_time]<time_bound_max)){
            param->atomic.factgrp[id] += grp[i_time].nphot[i_grp_space]/nphot_tot;
          }
        }
      }
    }
  }
}

void setAtomic(struct RUNPARAMS *param){
/**
  * reduce spectrum to mean energy and cross section, depending of NGRP
  */

  //  param->atomic.ngrp_space=1;   REAL space_bound[] = {13.6,INFINITY};
  param->atomic.ngrp_space=3;   REAL space_bound[] = {13.6, 24.587387, 54.417760,INFINITY};
  //param->atomic.ngrp_time=1;    REAL time_bound[] = {0.,INFINITY};
  param->atomic.ngrp_time=3;    REAL time_bound[] = {0.,10,100,INFINITY};

  int i_time;
  for(i_time=0;i_time<param->atomic.ngrp_time;i_time++){
    time_bound[i_time]*=1e6;
  }

  param->atomic.n = param->atomic.ngrp_space * param->atomic.ngrp_time;
  param->atomic.space_bound= (REAL*)calloc(param->atomic.ngrp_space,sizeof(REAL));
  param->atomic.time_bound= (REAL*)calloc(param->atomic.ngrp_time,sizeof(REAL));

  param->atomic.space_bound= space_bound;
  param->atomic.time_bound= time_bound;

  struct GRP *grp = (struct GRP*)calloc(param->spectrum.ntime,sizeof(struct GRP));

  for(i_time=0;i_time<param->spectrum.ntime;i_time++){
    get_space_GRP(param, i_time,&grp[i_time]);
    //printf("\n");
  }
/*
  for(i_time=0;i_time<param->spectrum.ntime;i_time++){
    int i_grp_space;
    for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){

      printf("hnu=%e \t",grp[i_time].Epho[i_grp_space] /ELECTRONVOLT);
      printf("sN=%e \t", grp[i_time].s_N[i_grp_space]);
      printf("sE=%e \t", grp[i_time].s_E[i_grp_space]);
      printf("fact=%e \n", grp[i_time].factgrp[i_grp_space]);

    }
  }
*/
  //alloc atomic
  param->atomic.hnu=    (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphae= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.alphai= (REAL*)calloc(param->atomic.n,sizeof(REAL));
  param->atomic.factgrp=(REAL*)calloc(param->atomic.n,sizeof(REAL));

  get_time_GRP(param, grp);

  int i_grp_time;
  for(i_grp_time=0;i_grp_time<param->atomic.ngrp_time;i_grp_time++){
      int i_grp_space;
      for(i_grp_space=0;i_grp_space<param->atomic.ngrp_space;i_grp_space++){
        int id= i_grp_time*param->atomic.ngrp_space+i_grp_space;
        printf("hnu=%e \t", param->atomic.hnu[id]/ELECTRONVOLT);
        printf("sE=%e \t", param->atomic.alphae[id]);
        printf("sN=%e \t", param->atomic.alphai[id]);
        printf("fact=%e \n", param->atomic.factgrp[id]);
      }
      printf("\n");
    }
}


void readStarburst99(struct RUNPARAMS *param){
/**
  * Read the Starburst99 spectrum database
  * http://www.stsci.edu/science/starburst99/data/fig1e.dat
  */

  sprintf(param->spectrum.path,"./src/atomic_data/Starburst99_spectrum.dat");

  FILE *buf=fopen(param->spectrum.path,"r");
  if(buf==NULL){
    printf("ERROR : cannot open the parameter file (%s given), please check\n",param->spectrum.path);
    abort();
  }

  //skip human readable header
  int nmax = 2048;
  char stream[nmax];
  fgets(stream, nmax, buf);
  fgets(stream, nmax, buf);
  fgets(stream, nmax, buf);

  //read header
  size_t rstat;
  rstat=fscanf(buf,"%s %d",stream,&param->spectrum.ntime);
  rstat=fscanf(buf,"%s %d",stream,&param->spectrum.nwave);
  printf("%d \n", param->spectrum.ntime);printf("%d \n", param->spectrum.nwave);

  //alloc
  param->spectrum.time=(REAL*)calloc(param->spectrum.ntime,sizeof(REAL));
  param->spectrum.wavelength=(REAL*)calloc(param->spectrum.nwave,sizeof(REAL));

  param->spectrum.flux=(REAL**)calloc(param->spectrum.ntime,sizeof(REAL*));
  int iwave;
  for(iwave=0;iwave<param->spectrum.nwave;iwave++)
    param->spectrum.flux[iwave]=(REAL*)calloc(param->spectrum.nwave,sizeof(REAL));

  //read time
  int itime;
  for(itime=0;itime<param->spectrum.ntime;itime++)
    rstat=fscanf(buf,"%lf",&param->spectrum.time[itime]);

  //read wavelenght and flux
  for(iwave=0;iwave<param->spectrum.nwave; iwave++){
    rstat=fscanf(buf,"%lf",&param->spectrum.wavelength[iwave]);
    for(itime=0;itime<param->spectrum.ntime;itime++){
      rstat=fscanf(buf,"%lf",&param->spectrum.flux[itime][iwave]);
    }
  }

  //set time (Myr in yr)
  for(itime=0;itime<param->spectrum.ntime;itime++){
    param->spectrum.time[itime]*=1e6;
  }

  //set Wavelenght (Angstrom in meter)
  for(iwave=0;iwave<param->spectrum.nwave; iwave++){
    param->spectrum.wavelength[iwave]*=1e-10;
  }

  //set flux (log(Erg/A) in J/m/s)
  for(iwave=0;iwave<param->spectrum.nwave; iwave++){
    for(itime=0;itime<param->spectrum.ntime;itime++){
      //printf("%e\n",param->spectrum.flux[itime][iwave] );
      param->spectrum.flux[itime][iwave] = POW(10.,param->spectrum.flux[itime][iwave]) * 1e-7/1e-10;
      //printf("%e\n",param->spectrum.flux[itime][iwave] );
    }
  }
}
