
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

#include "communication.h"
#include "particle.h"

#include "hydro_utils.h"


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
  REAL dxcur=pow(0.5,curoct->level);

  getcellnei(icell, vnei, vcell);
  for(ii=0;ii<6;ii++){ // looking for the gradient in 3 directions
    if(vnei[ii]==6){
      memcpy(&W,&(curoct->cell[vcell[ii]].field),sizeof(struct Wtype));

    }
    else{
      // Note that the neibourgh cell may not exist therefore we have to check
      if(curoct->nei[vnei[ii]]->child!=NULL){
	memcpy(&W,&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]].field),sizeof(struct Wtype));

#ifdef TRANSZM
	if(ii==4){
	  //if((curoct->nei[vnei[ii]]->child->z-curoct->z)>0.5){
	  if(curoct->z==0.){
	    // the neighbor is a periodic mirror 
	    memcpy(&W,&(curoct->cell[icell].field),sizeof(struct Wtype));
#ifdef REFZM
	    W.w*=-1.0;
	    //dxcur=1./pow(2,curoct->level);
	    //W.p=W.p+GRAV*W.d*dxcur;
#endif
	  }
	}
#endif 

#ifdef TRANSZP
	if(ii==5){
	  //if((curoct->nei[vnei[ii]]->child->z-curoct->z)<0.){
	  if((curoct->z+2.*dxcur)==1.){
	    // the neighbor is a periodic mirror 
	    memcpy(&W,&(curoct->cell[icell].field),sizeof(struct Wtype));
#ifdef REFZP
	    W.w*=-1.0;
	    //dxcur=1./pow(2,curoct->level);
	    //	    W.p=W.p-GRAV*W.d*dxcur;
#endif

	  }
	}
#endif 

#ifdef TRANSXM
	if(ii==0){
	  if((curoct->nei[vnei[ii]]->child->x-curoct->x)>0.5){
	    // the neighbor is a periodic mirror 
	    memcpy(&W,&(curoct->cell[ii].field),sizeof(struct Wtype));
	  }
	}
#endif 

#ifdef TRANSXP
	if(ii==1){
	  if((curoct->nei[vnei[ii]]->child->x-curoct->x)<0.){
	    // the neighbor is a periodic mirror 
	    memcpy(&W,&(curoct->cell[ii].field),sizeof(struct Wtype));
	  }
	}
#endif 

#ifdef TRANSYM
	if(ii==2){
	  if((curoct->nei[vnei[ii]]->child->y-curoct->y)>0.5){
	    // the neighbor is a periodic mirror 
	    memcpy(&W,&(curoct->cell[icell].field),sizeof(struct Wtype));
	  }
	}
    if(vnei[3]>6) abort();

#endif 

#ifdef TRANSYP
	if(ii==3){
	  if((curoct->nei[vnei[ii]]->child->y-curoct->y)<0.){
	    // the neighbor is a periodic mirror 
	    memcpy(&W,&(curoct->cell[icell].field),sizeof(struct Wtype));
	  }
	}
    if(vnei[3]>6) abort();
#endif 
 
      }
      else{
	// the neighbour does not exist we need to interpolate the value at the correct position
	coarse2fine_hydrolin(curoct->nei[vnei[ii]],Wi);
	memcpy(&W,Wi+vcell[ii],sizeof(struct Wtype));
    
      }
    }
    
    int ax=ii/2;
    int fact=((ii%2)==0?-1:1);
    gradd[ax]+=(W.d*fact);   
    /* gradu[ax]+=(W.u*fact); */
    /* gradv[ax]+=(W.v*fact);  */
    /* gradw[ax]+=(W.w*fact);    */
    /* gradp[ax]+=(W.p*fact); */

  }

  ratiod=sqrt(pow(gradd[0],2)+pow(gradd[1],2)+pow(gradd[2],2))*0.5/fabs(curoct->cell[icell].field.d+1e-10);
  ratiou=sqrt(pow(gradu[0],2)+pow(gradu[1],2)+pow(gradu[2],2))*0.5/fabs(curoct->cell[icell].field.u+1e-10);
  ratiov=sqrt(pow(gradv[0],2)+pow(gradv[1],2)+pow(gradv[2],2))*0.5/fabs(curoct->cell[icell].field.v+1e-10);
  ratiow=sqrt(pow(gradw[0],2)+pow(gradw[1],2)+pow(gradw[2],2))*0.5/(fabs(curoct->cell[icell].field.w)+1e-3);
  ratiop=sqrt(pow(gradp[0],2)+pow(gradp[1],2)+pow(gradp[2],2))*0.5/fabs(curoct->cell[icell].field.p+1e-10);

  //  if((ratiow>0.1)&&(fabs(curoct->cell[icell].field.w)<1e-15)) abort();

  ratio=ratiod;
  ratio=fmax(ratio,ratiou);
  ratio=fmax(ratio,ratiov);
  ratio=fmax(ratio,ratiow);
  ratio=fmax(ratio,ratiop);

  return ratio; 
  
}


#endif


//========================================================================================================================
//========================================================================================================================



//========================================================================================================================
//========================================================================================================================

struct OCT * L_refine_cells(int level, struct RUNPARAMS *param, struct OCT **firstoct, struct OCT ** lastoct, struct OCT * freeoct, struct CPUINFO *cpu, struct OCT *limit)
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

#ifdef WGRAV  
  struct Gtype Gi[8];
#endif


		  
  //if(nsteps==1) abort();
  nref=0;
  ndes=0;

  newoct=freeoct; // the new oct will be the first freeoct

  if(cpu->rank==0) printf("==> start refining on cpu %d lcoarse=%d lmax=%d\n",cpu->rank,param->lcoarse,param->lmax);

  dxcur=1./pow(2,level);
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      if(newoct->level!=0) {
		printf("newoct level=%d\n",newoct->level);
		abort();
	      }
	  
#if 1
	      // destroying octs with level>levelcoarse ==========================
	      if(((curoct->cell[icell].child!=NULL)&&(curoct->cell[icell].marked==0))&&(curoct->level>=param->lcoarse)){
		//if(((curoct->cell[icell].child!=NULL)&&(curoct->cell[icell].marked==0))){
		//if((curoct->level==(levelcoarse-1))&&(curoct->cpu==cpu->rank)) continue; // we don't want to destroy the cpu coarse grid

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

		// we cancels some flags
		desoct->vecpos=-1;
		desoct->border=0;
		
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

		// we fix the freeoct list
		freeoct->prev=desoct;
		desoct->next=freeoct;
		freeoct=desoct;
		
#ifdef PIC
		//======== dealing with particles
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
	      if(((curoct->cell[icell].child==NULL)&&(curoct->cell[icell].marked!=0))){


#ifdef WMPI
		if((curoct->cpu!=cpu->rank)&&(curoct->level>=(levelcoarse))){
		  int segok;
		  segok=segment_cell(curoct,icell,cpu,levelcoarse);// the current cell will be splitted according to a segmentation condition
		  if(!segok) continue;
 
		}

#endif

		// First we check if we don't violate the refinement rule (adaptative time step)
		int vrule=0;
		getcellnei(icell, vnei, vcell);
		for(ii=0;ii<6;ii++){
		  if(vnei[ii]!=6){
		    if(curoct->nei[vnei[ii]]!=NULL){	
		      if((curoct->nei[vnei[ii]]->child==NULL)&&(curoct->cpu==cpu->rank)){
			// refinement rule is violated so skip
			vrule=1;
		      }
		    }
		  }
		}
		if(vrule) continue;




		/* #ifdef WMPI */
		/* 		if(curoct->cpu!=cpu->rank){ */
		/* 		  if(curoct->level<levelcoarse){ */
		/* 		    if(curoct->cell[icell].marked==4) continue; */
		/* 		  } */
		/* 		} */
		/* #endif */

		if(curoct->cpu==cpu->rank) nref++;
		
		// a new oct is created
		newoct->parent=&(curoct->cell[icell]);
		newoct->level=curoct->level+1;
		newoct->x=curoct->x+( icell   %2)*dxcur;
		newoct->y=curoct->y+((icell/2)%2)*dxcur;
		newoct->z=curoct->z+( icell   /4)*dxcur;


	
		/* if((newoct->x==0.296875)*(newoct->y==0.25)*(newoct->z==0.)){ */
		/*   printf("SOCT2 FOUND\n"); */
		/*   SOCTX2=newoct; */
		/* } */

		// the new oct is connected to parent
		curoct->cell[icell].child=newoct;
		
		// it shares the same cpu
		newoct->cpu=curoct->cpu;

		// it is not vectorized yet
		newoct->vecpos=-1;
		newoct->border=0;

		//the neighbours
		for(ii=0;ii<6;ii++){
		  if(vnei[ii]!=6){
		    if(curoct->nei[vnei[ii]]!=NULL){	
		      if((curoct->nei[vnei[ii]]->child==NULL)&&(curoct->cpu==cpu->rank)){
			// here we refine too much
			printf("ouhla rank=%d curoct.cpu=%d\n",cpu->rank,curoct->cpu);
			abort();
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
		  coarse2fine_hydro2(&(curoct->cell[icell]),Wi);
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
		  newoct->cell[ii].temp=0.;
#endif

#ifdef WGRAV
		  for(ic=0;ic<3;ic++) newoct->cell[ii].f[ic]=0.;
#endif

#ifdef WHYDRO2
		  if(cpu->rank==curoct->cpu){
		    memcpy(&(newoct->cell[ii].field),Wi+ii,sizeof(struct Wtype)); 
		  }
		  else{
		    memset(&(newoct->cell[ii].field),0,sizeof(struct Wtype));
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

  if(cpu->rank==0){
    printf("octs created   = %d ",nref);
    printf("octs destroyed = %d\n",ndes);
  }

  return freeoct;

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

  if(cpu->rank==0) printf("==> start marking\n");
  //    for(level=levelmax;level>=param->lcoarse;level--) // looping over octs
  marker=0;
  nmark=0;
  for(ismooth=0;ismooth<nsmooth;ismooth++)
    {
	  //printf("level=%d ",level);
	  dx=1./pow(2,level);
	    
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
			  if(curoct->cell[icell].marked>0) abort();
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
				  if((curoct->nei[vnei[ii]]->child->x-curoct->x)>0.5){
				    newcell=NULL;
				    continue;
				  }
#endif

#ifdef TRANSXP
				  if((curoct->nei[vnei[ii]]->child->x-curoct->x)<0.){
				    newcell=NULL;
				    continue;
				  }
#endif

#ifdef TRANSYM
				  if((curoct->nei[vnei[ii]]->child->y-curoct->y)>0.5){
				    newcell=NULL;
				    continue;
				  }
#endif

#ifdef TRANSYP
				  if((curoct->nei[vnei[ii]]->child->y-curoct->y)<0.){
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

			    REAL den;
#ifdef TESTCOSMO
			    den=curoct->cell[icell].gdata.d+1.;
#else
#ifdef WGRAV
			    den=curoct->cell[icell].gdata.d;
#endif
#endif

/* #ifdef PART2 */
/* 			    if((curoct->x<0.5)&&(curoct->x>0.25)){ */
/* 			      if((curoct->y<0.5)&&(curoct->y>0.25)){ */
/* 				if((curoct->z<0.5)&&(curoct->z>0.25)){ */
/* 				  if(curoct->cell[icell].marked==0) { */
/* 				    curoct->cell[icell].marked=marker; */
/* 				    nmark++;stati[2]++; */
/* 				  } */
/* 				} */
/* 			      } */
/* 			    } */
/* #endif */
			    
#ifdef PIC
			    mcell=den*(curoct->level>=param->lcoarse)*dx*dx*dx;
			    if(mcell>mmax) mmax=mcell;
			    if((mcell>threshold)&&(curoct->cell[icell].marked==0)) {
 			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }
#else

 

#ifdef WGRAV
			    mcell=den*(curoct->level>=param->lcoarse);
			    if((mcell>threshold)&&(curoct->cell[icell].marked==0)) {
			      curoct->cell[icell].marked=marker;
			      nmark++;stati[2]++;
			    }
#endif

#endif

/* #ifdef WHYDRO2 */
/* 			    mcell=comp_grad_hydro(curoct, icell)*(curoct->level>=param->lcoarse); */

/* 			    if((mcell>(threshold))&&(curoct->cell[icell].marked==0)) { */
/* 			      curoct->cell[icell].marked=marker; */
/* 			      nmark++;stati[2]++; */
/* 			    } */
/* #endif */
			      

			  }
			  
			}

		      }
		  }while(nextoct!=NULL);
		//printf("pass=%d nmark=%d\n",pass,nmark);
	      }
#ifdef WMPI
	      MPI_Barrier(cpu->comm);
	      // first we correct from the marker diffusion
	      if((marker==2)||(marker==5)) mpi_cic_correct(cpu,sendbuffer,recvbuffer,1);
		
	      // second exchange boundaries
	      //mpi_exchange(cpu,sendbuffer,recvbuffer,3,1);
	      if(level>=(param->lcoarse-1)) mpi_exchange_level(cpu,sendbuffer,recvbuffer,3,1,level);
	      //mpi_exchange_level(cpu,sendbuffer,recvbuffer,3,1,level);
	      //breakmpi();
#endif
	    }
	  //printf("\n");
	}
  printf(" STAT MARK 0:%d 1:%d 2:%d mmax=%e thresh=%e\n",stati[0],stati[1],stati[2],mmax,param->amrthresh);

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
