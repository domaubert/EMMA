
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

void refine_cells(int levelcoarse, int levelmax, struct OCT **firstoct, struct OCT ** lastoct, struct OCT * endoct, struct CPUINFO *cpu)
{
  int nref,ndes;
  struct OCT *newoct;
  struct OCT *nextoct;
  struct OCT *curoct;
  struct OCT *desoct;
  int level;
  float dxcur;
  int icell;
  int sump,sum2;
  int ii;
  int vnei[6],vcell[6];
  int ip,xp,yp,zp;
  struct PART *curploc;
  struct PART *nexploc;
  struct PART *nexp;
  struct PART *curp;

  //if(nsteps==1) abort();
  nref=0;
  ndes=0;

  newoct=endoct+1; // the new oct will be the next to the last one

  if(cpu->rank==0) printf("==> start refining on cpu %d lcoarse=%d lmax=%d\n",cpu->rank,levelcoarse,levelmax);
  //breakmpi();
  //  for(level=levelcoarse;level<levelmax;level++) // looping over levels
  for(level=1;level<levelmax;level++) // looping over levels
    {
      dxcur=1./pow(2,level);
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
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

	      // destroying octs with level>levelcoarse ==========================
	      if(((curoct->cell[icell].child!=NULL)&&(curoct->cell[icell].marked==0))&&(curoct->level>levelcoarse)){
	      	ndes++;
		
		desoct=curoct->cell[icell].child; // the oct to destroy
		
		if(curoct->cell[icell].phead!=NULL){
	      	  printf("non void split cell !\n");
	      	  abort();
	      	}
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
	      }
	      
	      // creation of a new oct ==================
	      if((curoct->cell[icell].child==NULL)&&(curoct->cell[icell].marked!=0)){
		nref++;
		
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

		//the neighbours
		getcellnei(icell, vnei, vcell);
		for(ii=0;ii<6;ii++){
		  if(vnei[ii]!=6){
		
		    if((curoct->nei[vnei[ii]]->child==NULL)&&(curoct->cpu==cpu->rank)){
		      printf("ouhla rank=%d curoct.cpu=%d\n",cpu->rank,curoct->cpu);
		      abort();
		    }
		    // Note that boundary octs are refined but may have missing neighbors
		    if(curoct->nei[vnei[ii]]->child!=NULL) newoct->nei[ii]=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
		  }else{
		    newoct->nei[ii]=&(curoct->cell[vcell[ii]]);
		  }
		}

		// filling the cells
		for(ii=0;ii<8;ii++){
		  newoct->cell[ii].marked=0;
		  newoct->cell[ii].child=NULL;
		  newoct->cell[ii].density=0.;
		  newoct->cell[ii].idx=ii;
		  newoct->cell[ii].phead=NULL;
		}


		// splitting the particles
#if 1
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
		    
		    if((ip<0)||(ip>7)){
		      printf("ip=%d idx=%d\n",ip,curp->idx);
		      abort();
		    }
		    
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
		newoct++;
	      }

	      
	    }
	  //printf("nextoct=%p endoct=%p\n",nextoct,endoct);
	}while(nextoct!=NULL);
      //printf("level=%d done\n",level);
    }

  printf("octs created by rank %d =%d\n",cpu->rank,nref);
  printf("octs destroyed by rank %d=%d\n",cpu->rank,ndes);

}

 //------------------------------------------------------------------------

void mark_cells(int levelcoarse,int levelmax,struct OCT **firstoct, int nsmooth, float threshold, struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer){

  int nmark;
  int level;
  int marker;
  int ismooth;
  float dx;
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

    printf("==> start marking\n");
    //    for(level=levelmax;level>=levelcoarse;level--) // looping over octs
    for(level=levelmax;level>=1;level--) // looping over octs
      {
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
		if(nextoct==NULL) continue;
		
		do
		  {
		    curoct=nextoct;
		    nextoct=curoct->next;
		    
		    for(icell=0;icell<8;icell++) // looping over cells in oct
		      {
			switch(pass){
			  //=========================================================
			case 0: // marking cell already refined or marked marker=1/4
			  smark=0;
			  if(curoct->cell[icell].child!=NULL){ // CHECK IF ACTING ON ORIGINAL OR COPY
			    newoct=curoct->cell[icell].child;
			    for(ichild=0;ichild<8;ichild++){
			      smark+=((newoct->cell[ichild].marked!=0)||(newoct->cell[ichild].child!=NULL));
			    }
			  }
			  if(smark!=0){
			    curoct->cell[icell].marked=marker;
			    nmark++;
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
				  nmark++;
				}
			      }
			      else{
				// Note that the neibourgh cell may not exist therefore we have to check
				if(curoct->nei[vnei[ii]]->child!=NULL){
				  newcell=&(curoct->nei[vnei[ii]]->child->cell[vcell[ii]]);
				  if(curoct->nei[vnei[ii]]->child->cell[vcell[ii]].marked==0) {
				    curoct->nei[vnei[ii]]->child->cell[vcell[ii]].marked=marker;
				    nmark++;
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
				      nmark++;
				    }
				  }
				  else{
				    if(newoct->nei[vnei2[il]]->child!=NULL){
				      newcell2=&(newoct->nei[vnei2[il]]->child->cell[vcell2[il]]);
				      if(newoct->nei[vnei2[il]]->child->cell[vcell2[il]].marked==0){
					newoct->nei[vnei2[il]]->child->cell[vcell2[il]].marked=marker;
					nmark++;
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
					  nmark++;
					}
				      }
				      else{
				  	if(desoct->nei[vnei3[ip]]->child!=NULL){
					  if(desoct->nei[vnei3[ip]]->child->cell[vcell3[ip]].marked==0){
					    desoct->nei[vnei3[ip]]->child->cell[vcell3[ip]].marked=marker;
					    nmark++;
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
			  //if(countpart(curoct->cell[icell].phead)>threshold) curoct->cell[icell].marked=marker;
			  //if(curoct->npart>threshold) curoct->cell[icell].marked=marker;
			  if((curoct->cell[icell].density>threshold)&&(curoct->cell[icell].marked==0)) {
			    curoct->cell[icell].marked=marker;
			    nmark++;
			  }
			  break;
			}
		      }
	      
		  }while(nextoct!=NULL);
		//printf("pass=%d nmark=%d\n",pass,nmark);
#ifdef WMPI

		//breakmpi();
		// first we correct from the marker diffusion
		mpi_cic_correct(cpu,sendbuffer,recvbuffer,1);
		
		// second exchange boundaries
		mpi_exchange(cpu,sendbuffer,recvbuffer,3);
		//breakmpi();
#endif
	      }
	    //printf("\n");
	  }
      }

}
