#ifdef MOVIE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include <sys/types.h>
#include <sys/stat.h>


#include "prototypes.h"
#include "movie.h"
//=================================================================================================

int MOVIE_SNAP_NUMBER = 0;

void dumpMovie(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, int level, float aexp){

	if(cpu->rank==RANK_DISP)  printf("dump movie file ");

// Param------------------------
	const char ffolder[128] = "data/movie/" ;

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

	const int ntot = 	nmapx*nmapy;

	float * map = param->movie->map;

	float * m1 = map + 0*ntot;
	float * m2 = map + 1*ntot;
	float * m3 = map + 2*ntot;
	float * m4 = map + 3*ntot;

	int ilev;
	for(ilev=param->lcoarse; ilev<=lmap; ilev++){

		const REAL dxcur = pow(2.,    -ilev);
		const int  locN  = pow(2 ,lmap-ilev);

		struct OCT  *oct;
		struct OCT  *nextoct=firstoct[ilev-1];
		if(nextoct!=NULL)
		do {
			oct=nextoct;
			nextoct=oct->next;
			if(oct->cpu != cpu->rank) 	continue;

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
								const int id = x+y*nmapx;

#ifdef WGRAV
								m1[id] += (float)cell->gdata.p;
#endif
								m2[id] += (float)cell->field.d;
#ifdef WRADHYD
								m3[id] += (float)cell->field.dX/(float)cell->field.d;
#endif
#ifdef WRAD
								m4[id] += (float)cell->rfield.temp;
#endif
							}
						}
					}
				}
			}
		}while(nextoct!=NULL);
	}

	int ii;
	for(ii=0;ii<4*ntot;ii++) map[ii] /= nmapz * cpu->nproc;

  //=======================================
  //============= dump ====================
  //=======================================

	float* mapred = param->movie->map_reduce;
	MPI_Reduce(map, mapred, 4*ntot, MPI_FLOAT, MPI_SUM, 0, cpu->comm);
	if(cpu->rank==0){


    mkdir(ffolder, 0755);
		char fname[128];
		sprintf(fname,"%smovie_%08d",ffolder,MOVIE_SNAP_NUMBER++);

		FILE *fp = NULL;
		fp =fopen(fname,"wb");
		if(fp == NULL) printf("Cannot open %s\n", fname);

		fwrite(&nmapx,1,sizeof(int  ),fp);
		fwrite(&nmapy,1,sizeof(int  ),fp);
		fwrite(&aexp,1,sizeof(float),fp);
		fwrite(mapred,4*ntot,sizeof(float),fp);

		fclose(fp);
	}

  if(cpu->rank==0) printf("done\n");
}

#endif
