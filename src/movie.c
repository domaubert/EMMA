// ----------------------------------------------------------
// ----------------------------------------------------------
/** \file movie.c
  * \brief contain functions for movie files output
  * \author Nicolas Deparis
  *
  * Movie mode need the MOVIE preprocessor flag\n
  * It dump a 2D array a each time step containing
  * the projection of some physical field
  *
  * format:\n
  * int x (size along x axis)\n
  * int y (size along y axis)\n
  * float aexp (scale factor)\n
  * x*y float (potential)\n
  * x*y float (density)\n
  * x*y float (ionisation fraction)\n
  * x*y float (temperature)\n
  */
// ----------------------------------------------------------
// ----------------------------------------------------------


#ifdef MOVIE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <sys/stat.h>

#include "prototypes.h"
#include "movie.h"
//=================================================================================================


void dumpMovie(struct RUNPARAMS *param, struct CPUINFO *cpu, float aexp){

  //static int MOVIE_SNAP_NUMBER;
  int MOVIE_SNAP_NUMBER = cpu->nsteps;
	if(cpu->rank==RANK_DISP)  printf("Dumping movie file #%d",MOVIE_SNAP_NUMBER);

// Param------------------------
	const char ffolder[128] = "data/movie/" ;

  const int ndim = 2; //number of dimension


	const int lmap   = param->movie->lmap;
	const REAL xmin  = param->movie->xmin;
	const REAL xmax  = param->movie->xmax;
	const REAL ymin  = param->movie->ymin;
	const REAL ymax  = param->movie->ymax;
	const REAL zmin  = param->movie->zmin;
	const REAL zmax  = param->movie->zmax;
//------------------------------

  const int  nmap  = pow(2,lmap);
  const REAL dxmap = 1./nmap;

  const int nmapx	 =	(xmax-xmin)/dxmap;
  const int nmapy	 =	(ymax-ymin)/dxmap;
  const int nmapz	 =	(zmax-zmin)/dxmap;

  const int i0		 = xmin/dxmap;
  const int j0 		 = ymin/dxmap;
  const int k0		 = zmin/dxmap;

  int ntot;
  if (ndim==2)
    ntot = 	nmapx*nmapy;
  else if (ndim==3)
    ntot = 	nmapx*nmapy*nmapz;

	float * map = param->movie->map;

	float * m1 = map + 0*ntot;
	float * m2 = map + 1*ntot;
	float * m3 = map + 2*ntot;
	float * m4 = map + 3*ntot;

	int ilev;
	for(ilev=param->lcoarse; ilev<=lmap; ilev++){

		const REAL dxcur = pow(2.,    -ilev);
		const int  locN  = pow(2 ,lmap-ilev);

    int iOct;
    for(iOct=0; iOct<cpu->locNoct[ilev-1]; iOct++){
      struct OCT *oct=cpu->octList[ilev-1][iOct];

			int icell;
			for(icell=0;icell<8;icell++){
				struct CELL *cell = &oct->cell[icell];
				if( cell->child==NULL || ilev==lmap ){

					const REAL xc=oct->x+( icell   %2)*dxcur;
					const REAL yc=oct->y+((icell/2)%2)*dxcur;
					const REAL zc=oct->z+( icell/4   )*dxcur;

					if((zc<zmin)||(zc>zmax))	  continue;
					if((yc<ymin)||(yc>ymax))	  continue;
					if((xc<xmin)||(xc>xmax))	  continue;

					const int imap=xc/dxmap;
					const int jmap=yc/dxmap;
					const int kmap=zc/dxmap;

					int ii,jj,kk;
					for(kk=0;kk<locN;kk++){							if((kmap+kk)>=(nmapz+k0)) continue;
						for(jj=0;jj<locN;jj++){						if((jmap+jj)>=(nmapy+j0))	continue;
							for(ii=0;ii<locN;ii++){					if((imap+ii)>=(nmapx+i0)) continue;

								const int x  = imap+ii-i0;
								const int y  = jmap+jj-j0;
                const int z  = kmap+kk-k0;

								int id;
                if (ndim==2)
                  id = x+y*nmapx;
                else if (ndim==3)
                  id = x+y*nmapx+z*nmapx*nmapy;



            float field1=0;
            float field2=0;
            float field3=0;
            float field4=0;


            int set =0;
            switch(set){
              case 0:
                field1 = (float)oct->level;
                field2 = (float)cell->field.d;
                field3 = (float)cell->field.p;
                field4 = (float)oct->cpu;
                break;
              case 1:
                field1 = (float)oct->level;
                field2 = (float)cell->field.d;
                #ifdef WRAD
                field3 = (float)cell->rfield.temp;
                #endif // WRAD

                field4 = (float)oct->cpu;
                break;
            }

            int mode;
            if(strcmp(param->movie->mode_str, "max") ==0){mode=0;}
            if(strcmp(param->movie->mode_str, "avg") ==0){mode=1;}

            switch(mode){

              case 0:
								m1[id] = fmax( m1[id], field1);
								m2[id] = fmax( m2[id], field2);
								m3[id] = fmax( m3[id], field3);
								m4[id] = fmax( m4[id], field4);
                break;

              case 1:
								m1[id] = field1/nmapz;
								m2[id] = field2/nmapz;
								m3[id] = field3/nmapz;
								m4[id] = field4/nmapz;
								break;

              default:
                printf("Couldn't determine movie mode, check movie parameters in param.h\n");
                abort();
            }


							}
						}
					}
				}
			}
    }
	}

  //=======================================
  //============= dump ====================
  //=======================================

	float* mapred = param->movie->map_reduce;
	MPI_Reduce(map, mapred, 4*ntot, MPI_FLOAT, MPI_SUM, 0, cpu->comm);
	if(cpu->rank==RANK_DISP){

    mkdir(ffolder, 0755);
		char fname[128];
		sprintf(fname,"%smovie_%08d",ffolder,MOVIE_SNAP_NUMBER);

		FILE *fp = NULL;
		fp =fopen(fname,"wb");
		if(fp == NULL) printf("Cannot open %s\n", fname);

		fwrite(&nmapx,1,sizeof(int  ),fp);
		fwrite(&nmapy,1,sizeof(int  ),fp);
		fwrite(&aexp, 1,sizeof(float),fp);
		fwrite(mapred,4*ntot,sizeof(float),fp);

		fclose(fp);
	}

  if(cpu->rank==0) printf(" done\n");
}

#endif//MOVIE
