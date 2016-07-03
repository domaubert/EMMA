#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"
#include "oct.h"
#include "cic.h"

#ifdef WMPI
#include <mpi.h>
#endif
#include "segment.h"
#include "communication.h"
#include "particle.h"

#include "hydro_utils.h"
#include "poisson_utils.h"

#ifdef ZOOM
#include "zoom.h"
#endif


#ifdef WHYDRO2
REAL comp_grad_hydro(struct OCT *curoct, int icell){
  REAL gradd[3]={0.,0.,0.};
  REAL gradv[3]={0.,0.,0.};
  REAL gradu[3]={0.,0.,0.};
  REAL gradw[3]={0.,0.,0.};
  REAL gradp[3]={0.,0.,0.};


  int vcell[6],vnei[6];
  struct Wtype W;
  struct Wtype Wi[8];
  int ii;

  REAL ratiod,ratiou,ratiov,ratiow,ratiop,ratio;
  REAL dxcur=POW(0.5,curoct->level);

  getcellnei(icell, vnei, vcell);
  for(ii=0;ii<6;ii++){ // looking for the gradient in 3 directions
    if(vnei[ii]==6){
      memcpy(&W,&(curoct->cell[vcell[ii]].field),sizeof(struct Wtype));

    }
    else{
      // Note that the neibourgh cell may not exist therefore we have to check
      if(curoct->nei[vnei[ii]]->child!=NULL){
	memcpy(&W,&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]].field),sizeof(struct Wtype));

#
      }
      else{
	// the neighbour does not exist we need to interpolate the value at the correct position
	coarse2fine_hydrolin(curoct->nei[vnei[ii]],Wi);
	memcpy(&W,Wi+vcell[ii],sizeof(struct Wtype));

      }

#ifdef TRANSZM
	if(ii==4){
	  if(curoct->z==0.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[icell].field),sizeof(struct Wtype));
	  }
	}
#endif

#ifdef TRANSZP
	if(ii==5){
	  if((curoct->z+2.*dxcur)==1.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[icell].field),sizeof(struct Wtype));
	  }
	}
#endif

#ifdef TRANSXM
	if(ii==0){
	  if(curoct->x==0.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[ii].field),sizeof(struct Wtype));
	  }
	}
#endif

#ifdef TRANSXP
	if(ii==1){
	  if(curoct->x+2.*dxcur==1.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[ii].field),sizeof(struct Wtype));
	  }
	}
#endif

#ifdef TRANSYM
	if(ii==2){
	  if(curoct->y==0.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[icell].field),sizeof(struct Wtype));
	  }
	}
#endif

#ifdef TRANSYP
	if(ii==3){
	  if(curoct->y+2.*dxcur==1.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[icell].field),sizeof(struct Wtype));
	  }
	}
#endif

    }

    int ax=ii/2;
    int fact=((ii%2)==0?-1:1);
    gradd[ax]+=(W.d*fact);
    /* gradu[ax]+=(W.u*fact); */
    /* gradv[ax]+=(W.v*fact);  */
    /* gradw[ax]+=(W.w*fact);    */
    /* gradp[ax]+=(W.p*fact); */

  }

  ratiod=SQRT(POW(gradd[0],2)+POW(gradd[1],2)+POW(gradd[2],2))*0.5/fabs(curoct->cell[icell].field.d+1e-10);
  ratiou=SQRT(POW(gradu[0],2)+POW(gradu[1],2)+POW(gradu[2],2))*0.5/fabs(curoct->cell[icell].field.u+1e-10);
  ratiov=SQRT(POW(gradv[0],2)+POW(gradv[1],2)+POW(gradv[2],2))*0.5/fabs(curoct->cell[icell].field.v+1e-10);
  ratiow=SQRT(POW(gradw[0],2)+POW(gradw[1],2)+POW(gradw[2],2))*0.5/(fabs(curoct->cell[icell].field.w)+1e-3);
  ratiop=SQRT(POW(gradp[0],2)+POW(gradp[1],2)+POW(gradp[2],2))*0.5/fabs(curoct->cell[icell].field.p+1e-10);

  //  if((ratiow>0.1)&&(fabs(curoct->cell[icell].field.w)<1e-15)) abort();

  ratio=ratiod;
  ratio=FMAX(ratio,ratiou);
  ratio=FMAX(ratio,ratiov);
  ratio=FMAX(ratio,ratiow);
  ratio=FMAX(ratio,ratiop);

  return ratio;

}
#endif // WHYDRO2

// =========================================================================================================
// =========================================================================================================

#ifdef WRAD
#ifdef WCHEM

REAL comp_grad_rad(struct OCT *curoct, int icell){
  REAL gradd[3]={0.,0.,0.};
  REAL avgd[3]={0.,0.,0.};
  REAL gradn[3]={0.,0.,0.};
  REAL avgn[3]={0.,0.,0.};
  REAL grade[3]={0.,0.,0.};
  REAL avge[3]={0.,0.,0.};
  REAL val[6];

  int vcell[6],vnei[6];
  struct Rtype W;
  struct Rtype Wi[8];
  int ii;

  REAL ratiox,ratio,ratioe,ration;
  REAL dxcur=POW(0.5,curoct->level);

  getcellnei(icell, vnei, vcell);
  //printf("__\n");
  for(ii=0;ii<6;ii++){ // looking for the gradient in 3 directions
    if(vnei[ii]==6){
      memcpy(&W,&(curoct->cell[vcell[ii]].rfield),sizeof(struct Rtype));
    }
    else{
      // Note that the neibourgh cell may not exist therefore we have to check
      if(curoct->nei[vnei[ii]]->child!=NULL){
	memcpy(&W,&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]].rfield),sizeof(struct Rtype));


      }
      else{
	// the neighbour does not exist we need to interpolate the value at the correct position
	//coarse2fine_radlin(curoct->nei[vnei[ii]],Wi);
	/* coarse2fine_radlin(curoct->nei[vnei[ii]],Wi); */
	/* for(il=0;il<8;il++) memcpy(&Wi[il],&(curoct->nei[vnei[ii]]->rfield),sizeof(struct Rtype));  */

	memcpy(&W,&(curoct->nei[vnei[ii]]->rfield),sizeof(struct Rtype)); // straight injection

      }

      #ifdef TRANSZM
	if(ii==4){
	  if(curoct->z==0.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[icell].rfield),sizeof(struct Rtype));
	  }
	}
#endif

#ifdef TRANSZP
	if(ii==5){
	  if((curoct->z+2.*dxcur)==1.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[icell].rfield),sizeof(struct Rtype));
	  }
	}
#endif

#ifdef TRANSXM
	//printf("clunk3\n");
	if(ii==0){
	  //printf("clunk2\n");
	  if(curoct->x==0.){
	    // the neighbor is a periodic mirror
	    //printf("clunk\n");
	    memcpy(&W,&(curoct->cell[ii].rfield),sizeof(struct Rtype));
	  }
	}
#endif

#ifdef TRANSXP
	if(ii==1){
	  if(curoct->x+2.*dxcur==1.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[ii].rfield),sizeof(struct Rtype));
	  }
	}
#endif

#ifdef TRANSYM
	if(ii==2){
	  if(curoct->y==0.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[icell].rfield),sizeof(struct Rtype));
	  }
	}
#endif

#ifdef TRANSYP
	if(ii==3){
	  if(curoct->y+2.*dxcur==1.){
	    // the neighbor is a periodic mirror
	    memcpy(&W,&(curoct->cell[icell].rfield),sizeof(struct Rtype));
	  }
	}
#endif



    }

    int ax=ii/2;
    int fact=((ii%2)==0?-1:1);
    gradd[ax]+=(W.nhplus*fact);
    avgd[ax]+=W.nhplus*0.5;
    gradn[ax]+=(W.nh*fact);
    avgn[ax]+=W.nh*0.5;
    grade[ax]+=(W.e[0]*fact);
    avge[ax]+=W.e[0]*0.5;
  }


  ratiox=FMAX(fabs(gradd[0]/(avgd[0]+1e-10)),FMAX(fabs(gradd[1]/(avgd[1]+1e-10)),fabs(gradd[2]/(avgd[2]+1e-10))));
  ration=FMAX(fabs(gradn[0]/(avgn[0]+1e-10)),FMAX(fabs(gradn[1]/(avgn[1]+1e-10)),fabs(gradn[2]/(avgn[2]+1e-10))))*0.;
  ratioe=FMAX(fabs(grade[0]/(avge[0]+1e-10)),FMAX(fabs(grade[1]/(avge[1]+1e-10)),fabs(grade[2]/(avge[2]+1e-10))))*0;
  //ratiox=fabs(gradd[0]/(avgd[0]+1e-10));
  /* if((curoct->x==0.)&&(ratiox>1.9)) if(icell%2==0) if(curoct->level==6) { */
  /* 	printf("ratiox=%e %e %e %e\n",ratiox,gradd[0],gradd[1],gradd[2]); */
  /* 	abort(); */
  /*     } */
  //  if((ratiow>0.1)&&(fabs(curoct->cell[icell].field.w)<1e-15)) abort();

  ratio=FMAX(ratiox,ratioe);
  ratio=FMAX(ration,ratio);
  return ratio;

}

#endif
#endif

//========================================================================================================================
//========================================================================================================================

struct OCT * L_refine_cells(int level, struct RUNPARAMS *param, struct OCT **firstoct, struct OCT ** lastoct, struct OCT * freeoct, struct CPUINFO *cpu, struct OCT *limit, REAL aexp)
{
  int nref,ndes;
  struct OCT *newoct;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *desoct;
  REAL dxcur;
  int icell;
  int sump,sum2;
  int ii;
  int vnei[6],vcell[6];
  int ip,xp,yp,zp;
  int ic;

  struct PART *curploc;
  struct PART *nexploc;
  struct PART *nexp;
  struct PART *curp;

#ifdef WHYDRO2
  struct Wtype Wi[8];
#endif

#ifdef WRAD
  struct Rtype Ri[8];
#endif


#ifdef WGRAV
  struct Gtype Gi[8];
#endif



  //if(nsteps==1) abort();
  nref=0;
  ndes=0;

  newoct=freeoct; // the new oct will be the first freeoct

  if(cpu->rank==RANK_DISP) printf("==> start refining on cpu %d lcoarse=%d lmax=%d freeoct=%ld\n",cpu->rank,param->lcoarse,param->lmax,freeoct-firstoct[0]);
  long int dorg=freeoct-firstoct[0];

  dxcur=1./POW(2,level);
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if(newoct->level!=0) {
		printf("Error newoct level=%d\n",newoct->level);
		abort();
	      }

#if 1
	      // destroying octs with level>levelcoarse ==========================
	      if(((curoct->cell[icell].child!=NULL)&&(curoct->cell[icell].marked==0))&&(curoct->level>=param->lcoarse)){


		if(cpu->rank==curoct->cpu) ndes++;


		desoct=curoct->cell[icell].child; // the oct to destroy
#ifdef PIC
		if(curoct->cell[icell].phead!=NULL){
		  printf("non void split cell !\n");
		  abort();
		}
#endif
		// we remove the child from the refined cell
		curoct->cell[icell].child=NULL;

		// we remove the parent from the oct to be destroyed
		desoct->parent=NULL;


		// we remove the oct from the list

		if(desoct->prev!=NULL){
		  desoct->prev->next=desoct->next;
		}
		else{
		  //the oct is the first one at this level
		  if(firstoct[level]!=desoct){
		    printf("oups 1\n");
		    abort();
		  }
		  firstoct[level]=desoct->next;
		}

		if(desoct->next!=NULL){
		  desoct->next->prev=desoct->prev;
		}
		else{
		  //the oct is the last one at this level
		  if(lastoct[level]!=desoct){
		    printf("oups 2\n");
		    abort();
		  }
		  lastoct[level]=desoct->prev;
		}

		desoct->level=0;

		// we fixMak the freeoct list
		desoct->prev=NULL;
		freeoct->prev=desoct;
		desoct->next=freeoct;
		freeoct=desoct;
		// if destruction, newoct should be updated
		newoct=freeoct;
#ifdef PIC
		///======== dealing with particles
		// we gather the particles of the oct to be destroyed in the parent cell

		curploc=findlastpart(curoct->cell[icell].phead); // we get the last particle from the parent cell (should be nil in principle)
		sump=0;
		sum2=0;
		for(ii=0;ii<8;ii++){ // sweeping the cells of desoct
		  nexp=desoct->cell[ii].phead;
		  sump+=countpart(desoct->cell[ii].phead);
		  if(curploc==NULL){
		    curoct->cell[icell].phead=nexp;
		    curploc=findlastpart(nexp);
		  }
		  else if(nexp!=NULL){
		    curploc->next=nexp;
		    nexp->prev=curploc;
		    curploc=findlastpart(nexp);
		  }
		}
		sum2=countpart(curoct->cell[icell].phead);
		if(sum2!=sump){
		  printf("sum2=%d sump=%d\n",sum2,sump);
		  abort();
		}


#endif
	      }
#endif
	      // creation of a new oct ==================
		// gros test

	      if(((curoct->cell[icell].child==NULL)&&(curoct->cell[icell].marked!=0))){


		if(curoct->cell[icell].marked<10){
#ifdef WMPI
		if((curoct->cpu!=cpu->rank)&&(curoct->level>=(param->lcoarse))){
		  int segok;
		  segok=segment_cell(curoct,icell,cpu,param->lcoarse);// the current cell will be splitted according to a segmentation condition
		  if(!segok) {
		    //printf("skipping on proc %d\n",cpu->rank);
		    continue;
		  }
 		}
#endif

 		getcellnei(icell, vnei, vcell);


		// Past here the rule are respected


		//



		if(curoct->cpu==cpu->rank) nref++;

		// a new oct is created
		newoct->parent=&(curoct->cell[icell]);
		newoct->level=curoct->level+1;
		newoct->x=curoct->x+( icell   %2)*dxcur;
		newoct->y=curoct->y+((icell/2)%2)*dxcur;
		newoct->z=curoct->z+( icell   /4)*dxcur;


		// the new oct is connected to parent
		curoct->cell[icell].child=newoct;

		// it shares the same cpu
		newoct->cpu=curoct->cpu;

		// it is not vectorized yet
		/* newoct->vecpos=-1; */
		/* newoct->border=0; */

		//the neighbours
		for(ii=0;ii<6;ii++){
		  newoct->nei[ii]=NULL;//
		  if(vnei[ii]!=6){
		    if(curoct->nei[vnei[ii]]!=NULL){
		      if((curoct->nei[vnei[ii]]->child==NULL)&&(curoct->cpu==cpu->rank)){
#ifdef TRANSXM
			continue;
#endif
#ifndef TRANSXM
#ifndef TRANSXP
#ifndef TRANSYM
#ifndef TRANSYP
#ifndef TRANSZP
#ifndef TRANSZM
			// here we refine too much
			printf("ERROR ouhla rank=%d curoct.cpu=%d\n",cpu->rank,curoct->cpu);
			abort();
#endif
#endif
#endif
#endif
#endif
#endif
		      }
		      // Note that boundary octs are refined but may have missing neighbors
		      if(curoct->nei[vnei[ii]]->child!=NULL) newoct->nei[ii]=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
		    }
		  }else{
		    newoct->nei[ii]=&(curoct->cell[vcell[ii]]);
		  }
		}

#ifdef WHYDRO2
		if(cpu->rank==curoct->cpu){
		  int il;
		  //coarse2fine_hydrolin(&(curoct->cell[icell]),Wi);
		  //coarse2fine_hydro2(&(curoct->cell[icell]),Wi);
		  for(il=0;il<8;il++){
		    //   Ri[il].src=0.; */
		    memcpy(&Wi[il],&curoct->cell[icell].field,sizeof(struct Wtype));
		  }
		}
#endif


#ifdef WRAD
		if(cpu->rank==curoct->cpu){
		    REAL E;
		    REAL F;
		    int il;
		    //curoct->cell[icell].rfield.nhplus=curoct->cell[icell].field.dX/(PROTON_MASS*MOLECULAR_MU/param->unit.unit_mass);

		    //coarse2fine_rad2(&(curoct->cell[icell]),Ri);
  		     for(il=0;il<8;il++){
		       //   Ri[il].src=0.; */
 		       memcpy(&Ri[il],&curoct->cell[icell].rfield,sizeof(struct Rtype));


#ifdef SUPERNOVAE
		       Ri[il].snfb=0.;
#endif

#ifdef WRADHYD
		       // some fix here for RHD values

		       REAL xion=Ri[il].nhplus/Ri[il].nh;
#ifdef HELIUM
		       REAL xhe=Ri[il].nheplus/Ri[il].nh;
		       REAL xxhe=Ri[il].nhepplus/Ri[il].nh;
#endif

		       Wi[il].dX=Wi[il].d*(1.-YHE)*xion;
		       Ri[il].nh= Wi[il].d*(1.-YHE);
		       Ri[il].nhplus= Wi[il].dX;
		       Ri[il].eint= Wi[il].p/(GAMMA-1.);

#ifdef HELIUM
		       Wi[il].dXHE=Wi[il].d*YHE*xhe/yHE;
		       Wi[il].dXXHE=Wi[il].d*YHE*xxhe/yHE;
		       Ri[il].nheplus=Wi[il].dXHE/MHE_OVER_MH;
		       Ri[il].nhepplus=Wi[il].dXXHE/MHE_OVER_MH;
#endif
#endif
		     }

		}

#endif


#ifdef WGRAV
		if(cpu->rank==curoct->cpu){
		  coarse2fine_gravlin(&(curoct->cell[icell]),Gi);
		}
#endif

		// filling the cells
		for(ii=0;ii<8;ii++){
		  newoct->cell[ii].marked=0;
		  newoct->cell[ii].child=NULL;
		  newoct->cell[ii].idx=ii;

#ifdef PIC
		  newoct->cell[ii].density=curoct->cell[icell].density;
		  newoct->cell[ii].phead=NULL;
		  //newoct->cell[ii].temp=0.;
#endif

#ifdef WGRAV
		  for(ic=0;ic<3;ic++) newoct->cell[ii].f[ic]=0.;
#endif

#ifdef WHYDRO2
		  if(cpu->rank==curoct->cpu){
		    memcpy(&(newoct->cell[ii].field),Wi+ii,sizeof(struct Wtype));
		    memcpy(&(newoct->cell[ii].fieldnew),Wi+ii,sizeof(struct Wtype));
		  }
		  else{
		    memset(&(newoct->cell[ii].field),0,sizeof(struct Wtype));
		    memset(&(newoct->cell[ii].fieldnew),0,sizeof(struct Wtype));
		  }
#endif


#ifdef WRAD
		  if(cpu->rank==curoct->cpu){
        newoct->cell[ii].z_first_xion=curoct->cell[icell].z_first_xion;
        newoct->cell[ii].z_last_xion=curoct->cell[icell].z_last_xion;

		    memcpy(&(newoct->cell[ii].rfield),Ri+ii,sizeof(struct Rtype));
		    memcpy(&(newoct->cell[ii].rfieldnew),Ri+ii,sizeof(struct Rtype));
		  }
		  else{
		    memset(&(newoct->cell[ii].rfield),0,sizeof(struct Rtype));
		    memset(&(newoct->cell[ii].rfieldnew),0,sizeof(struct Rtype));
		  }
#endif


#ifdef WGRAV
		  if(cpu->rank==curoct->cpu){
		    memcpy(&(newoct->cell[ii].gdata),Gi+ii,sizeof(struct Gtype));
		    //memset(&(newoct->cell[ii].gdata),0,sizeof(struct Gtype));
		  }
		  else{
		    memset(&(newoct->cell[ii].gdata),0,sizeof(struct Gtype));
		  }
#endif
		}


#ifdef PIC
		// splitting the particles
		nexp=curoct->cell[icell].phead;
		if(nexp!=NULL){
		  do{ //sweeping the particles of the current cell
		    curp=nexp;
		    nexp=curp->next;


		    xp=(int)(2*(curp->x-newoct->x)/dxcur);xp=(xp>1?1:xp);xp=(xp<0?0:xp);xp=(xp==2?1:xp);
		    yp=(int)(2*(curp->y-newoct->y)/dxcur);yp=(yp>1?1:yp);yp=(yp<0?0:yp);yp=(yp==2?1:yp);
		    zp=(int)(2*(curp->z-newoct->z)/dxcur);zp=(zp>1?1:zp);zp=(zp<0?0:zp);zp=(zp==2?1:zp);
		    ip=xp+yp*2+zp*4;

		    //ip=(int)(2*(curp->x-newoct->x)/dxcur)+(int)(2*(curp->y-newoct->y)/dxcur)*2+(int)(2*(curp->z-newoct->z)/dxcur)*4;


		    // actual spliting
		    if(newoct->cell[ip].phead==NULL)
		      {
			// we create a new particle list in the cell
			newoct->cell[ip].phead=curp;
			//			newoct->cell[ip].density=8./(dxcur*dxcur*dxcur);

			// we remove the particle from the parent cell
			if(curp->prev!=NULL) curp->prev->next=curp->next;
			if(curp->next!=NULL) curp->next->prev=curp->prev;
			if(curp->prev==NULL) curoct->cell[icell].phead=curp->next; // we removed the head
			curp->prev=NULL;
			curp->next=NULL;
		      }
		    else{
		      //newoct->cell[ip].density+=8./(dxcur*dxcur*dxcur);
		      // we remove the particle from the parent cell
		      if(curp->prev!=NULL) curp->prev->next=curp->next;
		      if(curp->next!=NULL) curp->next->prev=curp->prev;
		      if(curp->prev==NULL) curoct->cell[icell].phead=curp->next; // we removed the head
		      curp->next=NULL;

		      nexploc=newoct->cell[ip].phead;
		      do{ //sweeping the particles of the new cell
			curploc=nexploc;
			nexploc=curploc->next;
		      }while(nexploc!=NULL);
		      curploc->next=curp;
		      curp->prev=curploc;
		    }
		  }while(nexp!=NULL);
		}


		if(curoct->cell[icell].phead!=NULL){
		  printf("cell not emptied after split !\n");
		  abort();
		};
#endif

		// preparing the next creations on level+1
		//printf("%p %p\n",freeoct,freeoct->next);
		freeoct=newoct->next; // we prepare the new free oct
		freeoct->prev=NULL; // the new free oct has no previous oct

		newoct->next=NULL;

		if(firstoct[level]==NULL){
		  firstoct[level]=newoct;
		  newoct->prev=NULL;
		}
		else{
		  newoct->prev=lastoct[level];
		  lastoct[level]->next=newoct;
		}
		lastoct[level]=newoct;


		newoct=freeoct; //ready to go to the next free oct
		if(newoct>=(limit-1)){
		  printf("\n ERROR === Allocated grid full, please increase ngridmax");
		  abort();
		}



		}
	      }

	    }
	  //printf("nextoct=%p endoct=%p\n",nextoct,endoct);
	}while(nextoct!=NULL);
    }  //printf("level=%d done\n",level);

#ifdef WMPI
  int nreftot;
  int ndestot;
  MPI_Allreduce(&nref,&nreftot,1,MPI_INT,MPI_SUM,cpu->comm);
  nref=nreftot;
  MPI_Allreduce(&ndes,&ndestot,1,MPI_INT,MPI_SUM,cpu->comm);
  ndes=ndestot;
#endif

  setOctList(firstoct[level-1], cpu, param,level);

  if(cpu->rank==RANK_DISP){
    printf("octs created   = %d ",nref);
    printf("octs destroyed = %d freeoctorg=%ld freeoct=%ld\n",ndes,dorg,freeoct-firstoct[0]);
  }

  return freeoct;

}

// ============================================================
// ============================================================

void L_check_rule(int level, struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu)
{
  struct OCT *newoct;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *desoct;
  int icell;
  int ii;
  int vnei[6],vcell[6];
  int ic;
  REAL dxcur=1./POW(2,level);

  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
    do // sweeping level
      {
	curoct=nextoct;
	nextoct=curoct->next;

	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    if(curoct->cell[icell].marked==0) continue;
	    // First we check if we don't violate the refinement rule (adaptative time step)
	    int vrule=0;
	    getcellnei(icell, vnei, vcell);
	    for(ii=0;ii<6;ii++){
#ifdef TRANSXM
	      if((curoct->x==0.)&&(ii==0)) continue;
#endif
#ifdef TRANSYM
	      if((curoct->y==0.)&&(ii==2)) continue;
#endif
#ifdef TRANSZM
	      if((curoct->z==0.)&&(ii==4)) continue;
#endif

#ifdef TRANSXP
	      if((curoct->x+2*dxcur==1.)&&(ii==1)) continue;
#endif
#ifdef TRANSYP
	      if((curoct->y+2*dxcur==1.)&&(ii==2)) continue;
#endif
#ifdef TRANSZP
	      if((curoct->z+2*dxcur==1.)&&(ii==5)) continue;
#endif

	      if((curoct->cpu==cpu->rank)){ // the violation rule is checked only on the current cpu octs
		if(curoct->nei[ii]->child==NULL){
		  // refinement rule is violated so skip
		  vrule=1;
		}
		else{
		  struct OCT *oct;
		  oct=curoct->nei[ii]->child;
		  int ii2;
		  for(ii2=0;ii2<6;ii2++){
		    if(ii2/2==ii/2) continue;

#ifdef TRANSXM
	      if((oct->x==0.)&&(ii2==0)) continue;
#endif
#ifdef TRANSYM
	      if((oct->y==0.)&&(ii2==2)) continue;
#endif
#ifdef TRANSZM
	      if((oct->z==0.)&&(ii2==4)) continue;
#endif

#ifdef TRANSXP
	      if((oct->x+2*dxcur==1.)&&(ii2==1)) continue;
#endif
#ifdef TRANSYP
	      if((oct->y+2*dxcur==1.)&&(ii2==2)) continue;
#endif
#ifdef TRANSZP
	      if((oct->z+2*dxcur==1.)&&(ii2==5)) continue;
#endif


		    if(oct->nei[ii2]->child==NULL){
		      vrule=1;
		    }
		    else{
		      struct OCT *oct2;
		      oct2=oct->nei[ii2]->child;
		      int ii3;
		      for(ii3=0;ii3<6;ii3++){

#ifdef TRANSXM
	      if((oct2->x==0.)&&(ii3==0)) continue;
#endif
#ifdef TRANSYM
	      if((oct2->y==0.)&&(ii3==2)) continue;
#endif
#ifdef TRANSZM
	      if((oct2->z==0.)&&(ii3==4)) continue;
#endif

#ifdef TRANSXP
	      if((oct2->x+2*dxcur==1.)&&(ii3==1)) continue;
#endif
#ifdef TRANSYP
	      if((oct2->y+2*dxcur==1.)&&(ii3==2)) continue;
#endif
#ifdef TRANSZP
	      if((oct2->z+2*dxcur==1.)&&(ii3==5)) continue;
#endif


			if((ii3/2==ii/2)||(ii3/2==ii2/2)) continue;
			if(oct2->nei[ii3]->child==NULL) vrule=1;
		      }
		    }
		  }
		}
	      }
	    }

	    if(vrule) {
	      curoct->cell[icell].marked=10; // the mark is cancelled
	    }

	  }
      }while(nextoct!=NULL);
  }

}


//========================================================================================================================
//========================================================================================================================



//=========================================================================================================================================
//=========================================================================================================================================

void L_mark_cells(int level,struct RUNPARAMS *param, struct OCT **firstoct, int nsmooth, REAL threshold, struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer){

  int nmark;
  int marker;
  int ismooth;
  REAL dx;
  int pass;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *newoct;
  struct OCT *desoct;
  int icell,ii,il,ip;
  int smark;
  int vnei[6],vcell[6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  struct CELL *newcell;
  struct CELL *newcell2;
  struct CELL *newcell3;
  int ichild;
  REAL mcell;
  REAL mmax=0.;
  int stati[3]={0,0,0};
  REAL rin;

  if(cpu->rank==RANK_DISP) printf("==> start marking");
  //    for(level=levelmax;level>=param->lcoarse;level--) // looping over octs
  marker=0;
  nmark=0;
  for(ismooth=0;ismooth<nsmooth;ismooth++)
    {
	  //printf("level=%d ",level);
	  dx=1./POW(2,level);

#ifdef ZOOM
 	  rin=param->rzoom*POW(param->fzoom,param->lmaxzoom-level-1);
	  //printf("rin=%e\n",rin);
#endif
	  for(pass=0;pass<3;pass++)
	    {
	      marker++;
	      nextoct=firstoct[level-1];
	      if(nextoct==NULL){
		//printf("Mark on level %d skipped by proc %d\n",level,cpu->rank);
	      }
	      else{
		do
		  {
		    curoct=nextoct;
		    nextoct=curoct->next;
		    if(curoct->cpu!=cpu->rank) continue;
		    for(icell=0;icell<8;icell++) // looping over cells in oct
		      {
			if((pass==0)&&(ismooth==0)){
			  if(curoct->cell[icell].marked>0){
          printf("marked>0 check amr.c\n");
          abort();
			  }
			}

			switch(pass){
			  //=========================================================
			case 0: // marking cell already refined or marked marker=1/4
			  smark=0;
			  if((curoct->cell[icell].child!=NULL)&&(curoct->cell[icell].marked==0)){ // CHECK IF ACTING ON ORIGINAL OR COPY
			    newoct=curoct->cell[icell].child;
			    for(ichild=0;ichild<8;ichild++){
			      smark+=((newoct->cell[ichild].marked!=0)||(newoct->cell[ichild].child!=NULL));
			    }
			  }
			  if(smark!=0){
			    curoct->cell[icell].marked=marker;
			    nmark++;
			    stati[0]++;
			  }
			  break;
			  //=========================================================
			case 1: //marking neighbors marker=2/5
			  if((curoct->cell[icell].marked<marker)&&(curoct->cell[icell].marked>0)){
			    getcellnei(icell, vnei, vcell);
			    for(ii=0;ii<6;ii++){ // marking the 6 cardinal neighbors
			      if(vnei[ii]==6){
				newcell=&(curoct->cell[vcell[ii]]);
				if(curoct->cell[vcell[ii]].marked==0){
				  curoct->cell[vcell[ii]].marked=marker;
				  nmark++;stati[1]++;
				}
			      }
			      else{
				// Note that the neibourgh cell may not exist therefore we have to check
				if(curoct->nei[vnei[ii]]->child!=NULL){
#ifdef TRANSXM
				  if((ii==0)&&(curoct->x==0.)){
				    newcell=NULL;
				    continue;
				  }
#endif

#ifdef TRANSXP
				  if((ii==1)&&((curoct->x+2.*dx)==1.)){
				    newcell=NULL;
				    continue;
				  }
#endif

#ifdef TRANSYM
				  if((ii==2)&&(curoct->y==0.)){
				    newcell=NULL;
				    continue;
				  }
#endif

#ifdef TRANSYP
				  if((ii==3)&&(curoct->y+2.*dx==1.)){
				    newcell=NULL;
				    continue;
				  }
#endif

#ifdef TRANSZM
				  if((ii==4)&&(curoct->z==0.)){
				    newcell=NULL;
				    continue;
				  }
#endif

#ifdef TRANSZP
				  if((ii==5)&&(curoct->z+2.*dx==1.)){
				    newcell=NULL;
				    continue;
				  }
#endif

				  newcell=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
				  if(curoct->nei[vnei[ii]]->child->cell[vcell[ii]].marked==0) {
				    curoct->nei[vnei[ii]]->child->cell[vcell[ii]].marked=marker;
				    nmark++;stati[1]++;
				  }
				}
				else{
				  newcell=NULL;
				}
			      }

			      // each of the 6 cardinal neighbors will mark 4 side neighbors
			      if(newcell!=NULL){
				newoct=cell2oct(newcell); // we get the parent oct
				getcellnei(newcell->idx, vnei2, vcell2); //we get the neighbors
				for(il=0;il<6;il++){
				  if((il/2)==(ii/2)) continue;
				  if(vnei2[il]==6){
				    newcell2=&(newoct->cell[vcell2[il]]);
				    if(newoct->cell[vcell2[il]].marked==0){
				      newoct->cell[vcell2[il]].marked=marker;
				      nmark++;stati[1]++;
				    }
				  }
				  else{
				    if(newoct->nei[vnei2[il]]->child!=NULL){
#ifdef TRANSXM
				      if((il==0)&&(newoct->x==0.)){
					newcell2=NULL;
					continue;
				      }
#endif

#ifdef TRANSXP
				      if((il==1)&&((newoct->x+2.*dx)==1.)){
					newcell2=NULL;
					continue;
				      }
#endif

#ifdef TRANSYM
				      if((il==2)&&(newoct->y==0.)){
					newcell2=NULL;
					continue;
				      }
#endif

#ifdef TRANSYP
				      if((il==3)&&(newoct->y+2.*dx==1.)){
					newcell2=NULL;
					continue;
				      }
#endif

#ifdef TRANSZM
				      if((il==4)&&(newoct->z==0.)){
					newcell2=NULL;
					continue;
				      }
#endif

#ifdef TRANSZP
				      if((il==5)&&(newoct->z+2.*dx==1.)){
					newcell2=NULL;
					continue;
				      }
#endif



				      newcell2=&(newoct->nei[vnei2[il]]->child->cell[vcell2[il]]);
				      if(newoct->nei[vnei2[il]]->child->cell[vcell2[il]].marked==0){
					newoct->nei[vnei2[il]]->child->cell[vcell2[il]].marked=marker;
					nmark++;stati[1]++;
				      }
				    }
				    else{
				      newcell2=NULL;
				    }
				  }

				  // ecah of the 4 side neighbors will mark 2 corners
				  if(newcell2!=NULL){
				    desoct=cell2oct(newcell2);
				    getcellnei(newcell2->idx, vnei3, vcell3);
				    for(ip=0;ip<6;ip++){
				      if(((ip/2)==(il/2))||((ip/2)==(ii/2))) continue;
				      if(vnei3[ip]==6){
					if(desoct->cell[vcell3[ip]].marked==0){
					  desoct->cell[vcell3[ip]].marked=marker;
					  nmark++;stati[1]++;
					}
				      }
				      else{
					if(desoct->nei[vnei3[ip]]->child!=NULL){
#ifdef TRANSXM
					  if((ip==0)&&(desoct->x==0.)){
					    continue;
					  }
#endif

#ifdef TRANSXP
					  if((ip==1)&&((desoct->x+2.*dx)==1.)){
					    continue;
					  }
#endif

#ifdef TRANSYM
					  if((ip==2)&&(desoct->y==0.)){
					    continue;
					  }
#endif

#ifdef TRANSYP
					  if((ip==3)&&(desoct->y+2.*dx==1.)){
					    continue;
					  }
#endif

#ifdef TRANSZM
					  if((ip==4)&&(desoct->z==0.)){
					    continue;
					  }
#endif

#ifdef TRANSZP
					  if((ip==5)&&(desoct->z+2.*dx==1.)){
					    continue;
					  }
#endif
					  if(desoct->nei[vnei3[ip]]->child->cell[vcell3[ip]].marked==0){
					    desoct->nei[vnei3[ip]]->child->cell[vcell3[ip]].marked=marker;

					    nmark++;stati[1]++;
					  }
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			  break;
			  //=========================================================
			case 2: // marking cells satisfying user defined criterion marker=3/6
			  if((curoct->level<=param->lmax)&&(ismooth==0)){ // we don't need to test the finest level


#ifdef ZOOM
			    // if within zoomed region the cell is marked in any case
			    int flagzoom=0;
			    if(level<param->lmaxzoom){
			      flagzoom=queryzoom(curoct,icell,dx,rin);
			      if((flagzoom)&&(curoct->cell[icell].marked==0)) {
				curoct->cell[icell].marked=marker;
				nmark++;stati[2]++;
			      }
			    }
#endif

			    REAL den;

#ifdef EVRARD
			    // ===================== EVRARD TEST ================

 			    mcell=comp_grad_hydro(curoct, icell)*(curoct->level>=param->lcoarse);
			    if(mcell>mmax) mmax=mcell;
			    if((mcell>(threshold))&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }

#else
			    // ===================== AMR COSMO ================

			    // First we define the density
#ifdef TESTCOSMO

#ifdef WGRAV
			    den=curoct->cell[icell].gdata.d+1.;
#endif // WGRAV

#ifdef ZELDOVICH
#ifdef WHYDRO2
			    den=curoct->cell[icell].field.d;
#endif // WHYDRO2
#endif // ZELDOVICH

#else // #ifndef TESTCOSMO

#ifdef WGRAV
			    // ->> utilise pour la cosmo // le gaz est utilise
			    den=curoct->cell[icell].gdata.d;
#endif // WGRAV
#endif // TESTCOSMO

			    // Second we apply a criterion

#ifdef PIC
#ifdef EDBERT
			    mcell=den*(curoct->level>=param->lcoarse)*dx*dx*dx;
			    if(mcell>mmax) mmax=mcell;
			    if((mcell>threshold)&&(curoct->cell[icell].marked==0)) {
 			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }

#else // #ifndef EDBERT
#ifdef ZELDOVICH
			    mcell=den*(curoct->level>=param->lcoarse);
			    if(mcell>mmax) mmax=mcell;
			    if((mcell>threshold)&&(curoct->cell[icell].marked==0)) {
 			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }
#else // #ifndef ZELDOVICH

			    // --------------- MAIN AMR COSMO

			    int refarea=1;
#ifdef ZOOM
			    refarea=(curoct->level>=param->lmaxzoom);
#endif // ZOOM

#ifndef AMRPART
			    mcell=den*(curoct->level>=param->lcoarse)*dx*dx*dx*refarea;
#else // #ifdef AMRPART

#ifdef PIC
			    int npart=0;
			    if(curoct->level>=param->lcoarse){
			      countpartDM(&curoct->cell[icell],&npart);
			    }
			    mcell=npart;
			    threshold=param->amrthresh0;

#else // #ifndef PIC
			    printf("AMR on particles SET ON without PIC enabled\n");
			    abort();
#endif // PIC
#endif // AMRPART
			    if(mcell>mmax) mmax=mcell;
			    if((mcell>threshold)&&(curoct->cell[icell].marked==0)) {
  			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }

			    // --------------- MAIN AMR COSMO
#endif // ZELDOVICH
#endif // EDBERT
#else // #ifndef PIC

			    // ===================== AMR NO COSMO ================

#ifdef WGRAV
			    mcell=den*(curoct->level>=param->lcoarse);
			    if((mcell>threshold)&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }
#endif // WGRAV

#ifdef SNTEST
/*
          threshold=1e-5;

			    mcell=comp_grad_hydro(curoct, icell)*(curoct->level>=param->lcoarse);//*(fabs(curoct->y-0.5)<0.05)*(fabs(curoct->z-0.5)<0.05);
			    if(mcell>mmax) mmax=mcell;
			    if((mcell>(threshold))&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }
*/
        //if (cpu->nsteps >= 50)
        {
          REAL epsilon = 1e-1;
          REAL threshold = 1. + epsilon;
          mcell=curoct->cell[icell].field.d;
			    if( (mcell>threshold)&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }
        }



#endif // SNTEST

#ifdef WRAD
#ifdef WRADTEST
			    // == START AMR STRATEGY FOR RAD TESTS
			    mcell=comp_grad_rad(curoct, icell)*(curoct->level>=param->lcoarse);

#ifdef TESTCLUMP
			    REAL den2;
			    den2=curoct->cell[icell].rfield.nh*param->unit.unit_N;
			    den=-1;

			    //mcell=(curoct->cell[icell].rfield.src>0.);
			    if((((den<8e-1)&&(den>1e-1))||(den2>250.))&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }

#else
			    //mcell=0.;
			    den=curoct->cell[icell].rfield.nhplus/curoct->cell[icell].rfield.nh; // xion
#endif // TESTCLUMP


 			    //mcell=(curoct->cell[icell].rfield.src>0.);
			    if(((den<8e-1)&&(den>1e-2))&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }

			    // == END AMR STRATEGY FOR RAD TESTS
#else

#ifdef WCHEM
			    mcell=comp_grad_rad(curoct, icell)*(curoct->level>=param->lcoarse);
			    if((mcell>(threshold))&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }

#endif // WCHEM
#endif // WRADTEST
#endif // WRAD
#endif // PIC

#ifdef TUBE
#ifdef SED
			    mcell=curoct->cell[icell].field.d;
#else
			    mcell=comp_grad_hydro(curoct, icell)*(curoct->level>=param->lcoarse);//*(fabs(curoct->y-0.5)<0.05)*(fabs(curoct->z-0.5)<0.05);
#endif
			    if(mcell>mmax) mmax=mcell;
			    if((mcell>(1.1))&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }
#endif // TUBE
#endif // EVRARD

			  }

			}
		      }
		  }while(nextoct!=NULL);
		//printf("pass=%d nmark=%d\n",pass,nmark);
	      }
#ifdef WMPI
	      MPI_Barrier(cpu->comm);
	      // we correct from the marker diffusion
	      if(marker%3==2) mpi_cic_correct_level(cpu, cpu->sendbuffer, cpu->recvbuffer, 1,level);

#endif

	    }
	  //printf("\n");
    }

 // printf("stat0=%d stat1=%d stat2=%d on rank %d\n",stati[0],stati[1],stati[2],cpu->rank);

#ifdef WMPI
  int stat0,stat1,stat2;
  REAL MMAX;
  MPI_Allreduce(&stati[0],&stat0,1,MPI_INT,MPI_SUM,cpu->comm);
  stati[0]=stat0;
  MPI_Allreduce(&stati[1],&stat1,1,MPI_INT,MPI_SUM,cpu->comm);
  stati[1]=stat1;
  MPI_Allreduce(&stati[2],&stat2,1,MPI_INT,MPI_SUM,cpu->comm);
  stati[2]=stat2;

  MPI_Allreduce(&mmax,&MMAX,1,MPI_REEL,MPI_MAX,cpu->comm);
  mmax=MMAX;
#endif

  if(cpu->rank==RANK_DISP) printf(" STAT MARK 0:%d 1:%d 2:%d mmax=%e thresh=%e\n",stati[0],stati[1],stati[2],mmax,threshold);

}


//=============================================================================

void clean_marks(int levelmax,struct OCT **firstoct){

  int level;
  struct OCT* curoct;
  struct OCT* nextoct;
  int icell;

    for(level=1;level<=levelmax;level++) // looping over levels
      {
	nextoct=firstoct[level-1];
	if(nextoct==NULL) continue;
	do // sweeping level
	  {
	    curoct=nextoct;
	    nextoct=curoct->next;
	    for(icell=0;icell<8;icell++) // looping over cells in oct
	      {
		curoct->cell[icell].marked=0.;
	      }
	  }while(nextoct!=NULL);

      }
}
//=============================================================================
//=============================================================================

void L_clean_marks(int level,struct OCT **firstoct){

  struct OCT* curoct;
  struct OCT* nextoct;
  int icell;

  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
    do // sweeping level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    curoct->cell[icell].marked=0.;
	  }
      }while(nextoct!=NULL);
  }
}
