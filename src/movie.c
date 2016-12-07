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
  * x*y float (field)\n
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
#include "io.h" // assign_grid_field

//=================================================================================================


void init_movie(struct RUNPARAMS *param){
	const int n=POW2(param->movie->lmap);
	const int n_field=param->out_grid->n_field_movie;
	param->movie->map=(float*)calloc(n_field*n*n,sizeof(float));
	param->movie->map_reduce = (float*)calloc(n_field*n*n,sizeof(float));
}

void dumpMovie(struct RUNPARAMS *param, struct CPUINFO *cpu, float aexp){

#ifdef WMPI
      MPI_Barrier(cpu->comm);
#endif
  double t_movie_start=MPI_Wtime();



  const int n=POW2(param->movie->lmap);
  const int n_field=param->out_grid->n_field_movie;

  // RAZ
  int i;
  for (i=0; i<n_field*n*n; i++){
    param->movie->map[i]=0;
    param->movie->map_reduce[i]=0;
  }


  const int debug = 0;
  //static int MOVIE_SNAP_NUMBER;
  const int MOVIE_SNAP_NUMBER = cpu->nsteps;
	if(cpu->rank==RANK_DISP)  printf("Dumping movie file #%d\n",MOVIE_SNAP_NUMBER);

  if (debug) printf("n_field=%d\n",n_field);

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

  int mode;
  if(strcmp(param->movie->mode_str, "max") ==0){mode=0;}
  if(strcmp(param->movie->mode_str, "avg") ==0){mode=1;}

  int ntot;
  if (ndim==2)
    ntot = 	nmapx*nmapy;
  else if (ndim==3)
    ntot = 	nmapx*nmapy*nmapz;

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

//============================================================================================================

        // getting the n_field values
        float field[n_field];
        int i,ii2=0;
        for (i=0;i<param->out_grid->n_field_tot; i++){
          if (param->out_grid->field_state_movie[i]){            field[ii2]=assign_grid_field(i,cell);
            ii2++;          }
        }

        switch(mode){
          case 0:
            for(i=0;i<n_field;i++){
              param->movie->map[id+i*ntot] = fmax( param->movie->map[id+i*ntot], field[i]);
            }
            break;

          case 1:
            for(i=0;i<n_field;i++){
              param->movie->map[id+i*ntot] += field[i]/nmapz;
            }
            break;

          default:
            printf("ERROR : Couldn't determine movie mode, check movie parameters in param.run, mode could be 'max' or 'avg' \n");
            abort();
        }

//============================================================================================================
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

#ifdef WMPI
	switch(mode){
    case 0:
      MPI_Reduce(param->movie->map, mapred, param->out_grid->n_field_movie*ntot, MPI_FLOAT, MPI_MAX, RANK_DISP, cpu->comm);
      break;
    case 1:
      MPI_Reduce(param->movie->map, mapred, param->out_grid->n_field_movie*ntot, MPI_FLOAT, MPI_SUM, RANK_DISP, cpu->comm);
      break;
	}
#endif // WMPI


	if(cpu->rank==RANK_DISP){

    int i,ii=0;
    for (i=0;i<param->out_grid->n_field_tot; i++){
      if (param->out_grid->field_state_movie[i]){

        mkdir(ffolder, 0755);

        char fffolder[256];
        sprintf(fffolder,"%s%s/",ffolder,param->out_grid->field_name[i]);

        if(debug) printf("field_name = %s\n", param->out_grid->field_name[i]);

        mkdir(fffolder, 0755);
        char fname[256];
        sprintf(fname,"%smovie_%08d",fffolder,MOVIE_SNAP_NUMBER);

        FILE *fp = NULL;
        fp =fopen(fname,"wb");
        if(fp==NULL) printf("Cannot open %s\n", fname);
        fwrite(&nmapx,1,sizeof(int  ),fp);
        fwrite(&nmapy,1,sizeof(int  ),fp);
        fwrite(&aexp, 1,sizeof(float),fp);
        fwrite(mapred+ii*ntot,ntot,sizeof(float),fp);
        fclose(fp);

        ii++;      }
    }
	}

  if(debug) abort();


#ifdef WMPI
    MPI_Barrier(cpu->comm);
#endif

  	if(cpu->rank==RANK_DISP) printf("==== CPU Movie TOTAL TIME =%e\n",MPI_Wtime()-t_movie_start);

}
#endif//MOVIE
