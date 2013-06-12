#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"

#ifdef PIC

//==================================================================

void part2cell_cic(struct PART *curp, struct OCT *curoct, int icell, char full)
{
  REAL xc,yc,zc;
  //REAL dxcur=1./pow(2,curoct->level);
  REAL dxcur=1./(1<<curoct->level);
  
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  int neip[3];
  REAL tx,ty,tz;
  REAL dx,dy,dz;
  char fullok=0;

  int visit[8]={0,0,0,0,0,0,0,0};
  REAL vcont[8];

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;

  int i1,i2,i3;
  int idx1,idx2; // scan the cic cell being analysed


  xc=curoct->x+( icell   & 1)*dxcur+dxcur*0.5; // coordinates of the cell center 
  yc=curoct->y+((icell>>1)& 1)*dxcur+dxcur*0.5;
  zc=curoct->z+( icell>>2   )*dxcur+dxcur*0.5; 

  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);

  // the CIC weights
  
  tx=(curp->x-xc)/dxcur;tx=(tx>0.?tx:-tx);
  ty=(curp->y-yc)/dxcur;ty=(ty>0.?ty:-ty);
  tz=(curp->z-zc)/dxcur;tz=(tz>0.?tz:-tz);

  dx=1.-tx;
  dy=1.-ty;
  dz=1.-tz;

  if(!full)
    {
      fullok=(tx<=1.)*(ty<=1.)*(tz<=1.);
    } // necessary for boundary cells

  //                   ZYX

  REAL dunit=curp->mass/(dxcur*dxcur*dxcur);
  vcont[0]=dz*dy*dx*dunit; //000
  vcont[1]=dz*dy*tx*dunit; //001
  vcont[2]=dz*ty*dx*dunit; //010
  vcont[4]=tz*dy*dx*dunit; //100

  vcont[3]=dz*ty*tx*dunit; //011
  vcont[6]=tz*ty*dx*dunit; //110
  vcont[5]=tz*dy*tx*dunit; //101

  vcont[7]=tz*ty*tx*dunit; //111


  // contrib to current cell 000 =====

  REAL tot=0;
  REAL tot2=0;
  int ntot=0;
  if(full||fullok){
    curoct->cell[icell].density+=vcont[0];
  }

  if(full){
    // contribs to cardinal neighbors
    getcellnei(icell, vnei, vcell);
    for(i1=0;i1<3;i1++){
      idx1=1<<i1;
      if(vnei[neip[i1]]==6){
	curoct->cell[vcell[neip[i1]]].density+=vcont[idx1];
	newcell=&(curoct->cell[vcell[neip[i1]]]);
      }
      else{
	if(curoct->nei[vnei[neip[i1]]]->child!=NULL){
	  /* if(curoct->nei[vnei[neip[i1]]]->child->cpu!=curoct->cpu){ */
	  /*   printf("hello %d diffuse in %d\n",curoct->cpu,curoct->nei[vnei[neip[i1]]]->child->cpu); */
	  /* } */
	  curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].density+=vcont[idx1];
	  newcell=&(curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]]);
	}
	else{
	  newcell=NULL;
	}
      }
    
      // contrib to 2nd order neighbours
      if(newcell!=NULL){
	for(i2=0;i2<3;i2++){
	  if(i2==i1) continue;
	  idx2=(1<<i1)+(1<<i2);
	  if(visit[idx2]) continue;

	  getcellnei(newcell->idx, vnei2, vcell2);
	  newoct=cell2oct(newcell);
	  if(vnei2[neip[i2]]==6){
	    newoct->cell[vcell2[neip[i2]]].density+=vcont[idx2];
	    newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	  }
	  else{
	    if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	      newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].density+=vcont[idx2];
	      newcell2=&(newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]]);
	      visit[idx2]=1;
	    }
	    else{
	      newcell2=NULL;
	    }
	  }
	}

	// contrib to 3rd order neighbors
	if(newcell2!=NULL){
	  for(i3=0;i3<3;i3++){
	    if((i3==i1)||(i3==i2)) continue;
	    if(visit[7]) continue;
	    getcellnei(newcell2->idx, vnei3, vcell3);
	    newoct2=cell2oct(newcell2);
	    if(vnei3[neip[i3]]==6){
	      newoct2->cell[vcell3[neip[i3]]].density+=vcont[7];
	      visit[7]=1;
	    }
	    else{
	      if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
		newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].density+=vcont[7];
		visit[7]=1;
	      }
	    }
	  }
	}
      }
    }
  }

}

//=====================================================================================
//=====================================================================================


//=======================================================================================================
#ifdef PIC
void cell2part_cic(struct PART *curp, struct OCT *curoct, int icell, REAL dt)
{
  REAL xc,yc,zc;
  //  REAL dxcur=1./pow(2,curoct->level);
  REAL dxcur=1./(1<<curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  int neip[3];
  REAL tx,ty,tz;
  REAL dx,dy,dz;
  REAL contrib;
  struct OCT *curoctlr;
  int ic;
  REAL accel[3]={0.,0.,0.};

  char hres=1; // could be switched to hres=0 if particle is not deep enough

  int visit[8]={0,0,0,0,0,0,0,0};
  REAL vcont[8];
  //REAL vf[8];

  int i1,i2,i3;
  int idx1,idx2; // scan the cic cell being analysed

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;
  

  xc=curoct->x+( icell   & 1)*dxcur+dxcur*0.5; // coordinates of the cell center 
  yc=curoct->y+((icell>>1)& 1)*dxcur+dxcur*0.5;
  zc=curoct->z+( icell>>2   )*dxcur+dxcur*0.5; 
		  
  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);

  // the CIC weights

  tx=(curp->x-xc)/dxcur;tx=(tx>0.?tx:-tx);
  ty=(curp->y-yc)/dxcur;ty=(ty>0.?ty:-ty);
  tz=(curp->z-zc)/dxcur;tz=(tz>0.?tz:-tz);

  /* tx=fabs((curp->x-xc)/dxcur); */
  /* ty=fabs((curp->y-yc)/dxcur); */
  /* tz=fabs((curp->z-zc)/dxcur); */
	  
  dx=1.-tx;
  dy=1.-ty;
  dz=1.-tz;


  vcont[0]=dz*dy*dx; //000
  vcont[1]=dz*dy*tx; //001
  vcont[2]=dz*ty*dx; //010
  vcont[4]=tz*dy*dx; //100

  vcont[3]=dz*ty*tx; //011
  vcont[6]=tz*ty*dx; //110
  vcont[5]=tz*dy*tx; //101

  vcont[7]=tz*ty*tx; //111

  /* if(curp->idx==0){ */
  /*   for(ic=0;ic<8;ic++) printf("%e\n",vcont[ic]); */
  /* } */

  
  // checking the neighbors
  int inei;
  getcellnei(icell, vnei, vcell);
  for(inei=0;inei<6;inei++){
    if(curoct->nei[inei]->child==NULL) hres=0; // particle is not deep enough for this level
  }


  if(hres!=0){ // can be skipped if particle is at border

    // contrib from current cell 000 =====

    contrib=vcont[0];
    if((contrib<=1.)&&(contrib>=0.)){
      for(ic=0;ic<3;ic++) accel[ic]+=curoct->cell[icell].f[ic]*contrib;
    }
    visit[0]=1;

    // contribs to cardinal neighbors
    for(i1=0;i1<3;i1++){
      //idx1=pow(2,i1);
      idx1=(1<<i1);

      contrib=vcont[idx1];
      if(vnei[neip[i1]]==6){
	for(ic=0;ic<3;ic++) accel[ic]+=curoct->cell[vcell[neip[i1]]].f[ic]*contrib;
	newcell=&(curoct->cell[vcell[neip[i1]]]);
	visit[idx1]=1;
	//vf[idx1]=curoct->cell[vcell[neip[i1]]].gdata.p;
      }
      else{
	if(curoct->nei[vnei[neip[i1]]]->child!=NULL){
	  for(ic=0;ic<3;ic++) accel[ic]+=curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].f[ic]*contrib;
	  newcell=&(curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]]);
	  visit[idx1]=1;
	  //vf[idx1]=curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].gdata.p;

	}
	else{
	  // the particle is not deep enough we stop
	  for(ic=0;ic<3;ic++) accel[ic]=0;
	  hres=0;
	  newcell=NULL;
	  break;
	}
      }
    
      // contrib to 2nd order neighbours
      if(newcell!=NULL){
	for(i2=0;i2<3;i2++){
	  //idx2=pow(2,i1)+pow(2,i2);
	  idx2=(1<<i1)+(1<<i2);
	  if(i2==i1) continue;
	  if(visit[idx2]) continue;

	  contrib=vcont[idx2];
	  getcellnei(newcell->idx, vnei2, vcell2);
	  newoct=cell2oct(newcell);
	  if(vnei2[neip[i2]]==6){
	    for(ic=0;ic<3;ic++) accel[ic]+=newoct->cell[vcell2[neip[i2]]].f[ic]*contrib;
	    newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	    //vf[idx2]=newoct->cell[vcell2[neip[i2]]].gdata.p;
	  }
	  else{
	    if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	      for(ic=0;ic<3;ic++) accel[ic]+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].f[ic]*contrib;
	      newcell2=&(newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]]);
	      visit[idx2]=1;
	      //vf[idx2]=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].gdata.p;
	    }
	    else{
	      for(ic=0;ic<3;ic++) accel[ic]=0;
	      hres=0;
	      newcell2=NULL;
	      break;
	    }
	  }
	}

	// contrib to 3rd order neighbors
	if(newcell2!=NULL){
	  for(i3=0;i3<3;i3++){
	    if((i3==i1)||(i3==i2)) continue;
	    if(visit[7]) continue;
	    contrib=vcont[7];
	    getcellnei(newcell2->idx, vnei3, vcell3);
	    newoct2=cell2oct(newcell2);
	    if(vnei3[neip[i3]]==6){
	      for(ic=0;ic<3;ic++) accel[ic]+=newoct2->cell[vcell3[neip[i3]]].f[ic]*contrib;
	      visit[7]=1;
	      //vf[7]=newoct2->cell[vcell3[neip[i3]]].gdata.p;
	    }
	    else{
	      if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
		for(ic=0;ic<3;ic++) accel[ic]+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].f[ic]*contrib;
		visit[7]=1;
		//vf[7]=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].gdata.p;
	      }
	      else{
		for(ic=0;ic<3;ic++) accel[ic]=0;
		hres=0;
		break;
	      }
	    }
	  }
	}
      }
    }
  }


  // we must recompute the force from level-1 if the particle is not deep enough LOWRES
  if(hres!=1){
    
    // Getting the new curoct at low resolution
    curoctlr=cell2oct(curoct->parent);
    icell=curoct->parent->idx;
    dxcur=1./(1<<curoctlr->level);

    // start again

    xc=curoctlr->x+( icell   & 1)*dxcur+dxcur*0.5; // coordinates of the cell center 
    yc=curoctlr->y+((icell>>1)& 1)*dxcur+dxcur*0.5;
    zc=curoctlr->z+( icell>>2   )*dxcur+dxcur*0.5; 

    // here we compute the indexes of the direct neighbors which are involved in the cic
    neip[0]=(curp->x<xc?0:1);
    neip[1]=(curp->y<yc?2:3);
    neip[2]=(curp->z<zc?4:5);

    // the CIC weights
	  
    tx=(curp->x-xc)/dxcur;tx=(tx>0.?tx:-tx);
    ty=(curp->y-yc)/dxcur;ty=(ty>0.?ty:-ty);
    tz=(curp->z-zc)/dxcur;tz=(tz>0.?tz:-tz);

    dx=1.-tx;
    dy=1.-ty;
    dz=1.-tz;

    for(i1=0;i1<8;i1++) visit[i1]=0;

    vcont[0]=dz*dy*dx; //000
    vcont[1]=dz*dy*tx; //001
    vcont[2]=dz*ty*dx; //010
    vcont[4]=tz*dy*dx; //100

    vcont[3]=dz*ty*tx; //011
    vcont[6]=tz*ty*dx; //110
    vcont[5]=tz*dy*tx; //101

    vcont[7]=tz*ty*tx; //111


    // contrib from current cell 000 =====
    /* REAL tot=0; */
    /* int ntot=0; */
    contrib=vcont[0];
    if((contrib<=1.)&&(contrib>=0.)){

      for(ic=0;ic<3;ic++) accel[ic]=curoctlr->cell[icell].f[ic]*contrib;
    }
  
  
    // contribs to cardinal neighbors
    getcellnei(icell, vnei, vcell);
    for(i1=0;i1<3;i1++){
      //idx1=pow(2,i1);
      idx1=(1<<i1);
      contrib=vcont[idx1];
      if(vnei[neip[i1]]==6){
	for(ic=0;ic<3;ic++) accel[ic]+=curoctlr->cell[vcell[neip[i1]]].f[ic]*contrib;
	newcell=&(curoctlr->cell[vcell[neip[i1]]]);
      }
      else{
	if(curoctlr->nei[vnei[neip[i1]]]->child!=NULL){
	  for(ic=0;ic<3;ic++) accel[ic]+=curoctlr->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].f[ic]*contrib;
	  newcell=&(curoctlr->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]]);
	}
      }
    
      // contrib to 2nd order neighbours
      if(newcell!=NULL){
	for(i2=0;i2<3;i2++){
	  //idx2=pow(2,i1)+pow(2,i2);
	  idx2=(1<<i1)+(1<<i2);
	  if(i2==i1) continue;
	  if(visit[idx2]) continue;

	  contrib=vcont[idx2];
	  getcellnei(newcell->idx, vnei2, vcell2);
	  newoct=cell2oct(newcell);
	  if(vnei2[neip[i2]]==6){
	    for(ic=0;ic<3;ic++) accel[ic]+=newoct->cell[vcell2[neip[i2]]].f[ic]*contrib;
	    newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	  }
	  else{
	    if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	      for(ic=0;ic<3;ic++) accel[ic]+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].f[ic]*contrib;
	      newcell2=&(newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]]);
	      visit[idx2]=1;
	    }
	  }
	}

	// contrib to 3rd order neighbors
	if(newcell2!=NULL){
	  for(i3=0;i3<3;i3++){
	    if((i3==i1)||(i3==i2)) continue;
	    if(visit[7]) continue;
	    contrib=vcont[7];
	    getcellnei(newcell2->idx, vnei3, vcell3);
	    newoct2=cell2oct(newcell2);
	    if(vnei3[neip[i3]]==6){
	      for(ic=0;ic<3;ic++) accel[ic]+=newoct2->cell[vcell3[neip[i3]]].f[ic]*contrib;
	      visit[7]=1;
	    }
	    else{
	      if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
		for(ic=0;ic<3;ic++) accel[ic]+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].f[ic]*contrib;
		visit[7]=1;
	      }
	    }
	  }
	}
      }
    }
  }
  


  
  // Once we have the acceleration we can compute the velocity

#ifdef PERFECT
   if(curp->mass<0.5) {
    REAL r=sqrt((curp->x-0.5)*(curp->x-0.5)+(curp->y-0.5)*(curp->y-0.5)+(curp->z-0.5)*(curp->z-0.5));
    accel[0]=(curp->x-0.5)/(r*r*r);
    accel[1]=(curp->y-0.5)/(r*r*r);
    accel[2]=(curp->z-0.5)/(r*r*r);
   }
#endif
   
   curp->vx+=-accel[0]*dt;
   curp->vy+=-accel[1]*dt;
   curp->vz+=-accel[2]*dt;
   
   //printf("idx=%d fx=%e/%e fy=%e/%e fz=%e/%e\n",curp->idx,accel[ic]2,accel[ic],accely2,accely,accelz2,accelz);
   
#ifndef TESTCOSMO
   curp->fx=-accel[0];
   curp->fy=-accel[1];
   curp->fz=-accel[2];
   
   if(curp->idx==1) {
     REAL r=sqrt((curp->x-0.5)*(curp->x-0.5)+(curp->y-0.5)*(curp->y-0.5)+(curp->z-0.5)*(curp->z-0.5));
     REAL f=sqrt((curp->fx)*(curp->fx)+(curp->fy)*(curp->fy)+(curp->fz)*(curp->fz));
     REAL fth=1./(r*r);
     printf("r=%e force=%e fth=%e rap=%e fx=%e fy=%e fz=%e\n",r,f,fth,(f-fth)/fth,curp->fx,curp->fy,curp->fz);
     //for(ic=0;ic<8;ic++) printf("%d %e %e\n",visit[ic],vcont[ic],vf[ic]);
   }
#endif
   //if((curp->idx==1)&&(curp->fx>72)) abort();

}
#endif



//=====================================================================================
//=====================================================================================

REAL cell2part_cic_egy(struct PART *curp, struct OCT *curoct, int icell)
{
  REAL xc,yc,zc;
  REAL dxcur=1./pow(2,curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  int neip[3];
  REAL tx,ty,tz;
  REAL dx,dy,dz;
  REAL contrib;
  struct OCT *curoctlr;
  REAL accel=0.;
  char hres=1; // could be switched to hres=0 if particle is not deep enough

  int visit[8]={0,0,0,0,0,0,0,0};
  REAL vcont[8];

  int i1,i2,i3;
  int idx1,idx2; // scan the cic cell being analysed

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;
  

  xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
  yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
  zc=curoct->z+( icell/4   )*dxcur+dxcur/2; 

		  
  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);
		   


  // the CIC weights
  tx=fabs((curp->x-xc)/dxcur);
  ty=fabs((curp->y-yc)/dxcur);
  tz=fabs((curp->z-zc)/dxcur);
	  
  dx=1.-tx;
  dy=1.-ty;
  dz=1.-tz;


  vcont[0]=dz*dy*dx; //000
  vcont[1]=dz*dy*tx; //001
  vcont[2]=dz*ty*dx; //010
  vcont[4]=tz*dy*dx; //100

  vcont[3]=dz*ty*tx; //011
  vcont[6]=tz*ty*dx; //110
  vcont[5]=tz*dy*tx; //101

  vcont[7]=tz*ty*tx; //111


  // contrib from current cell 000 =====
  REAL tot=0;
  int ntot=0;

  contrib=vcont[0];
  if((contrib<=1.)&&(contrib>=0.)){
    accel+=curoct->cell[icell].gdata.p*contrib;
  }
  

  // contribs to cardinal neighbors
  getcellnei(icell, vnei, vcell);
  for(i1=0;i1<3;i1++){
    idx1=pow(2,i1);
    contrib=vcont[idx1];
    if(vnei[neip[i1]]==6){
      accel+=curoct->cell[vcell[neip[i1]]].gdata.p*contrib;
      newcell=&(curoct->cell[vcell[neip[i1]]]);
      tot+=contrib;
      ntot++;
    }
    else{
      if(curoct->nei[vnei[neip[i1]]]->child!=NULL){
	accel+=curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].gdata.p*contrib;
	newcell=&(curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]]);
	tot+=contrib;
	ntot++;

      }
      else{
	// the particle is not deep enough we stop
	accel=0;
	hres=0;
	newcell=NULL;
	break;
      }
    }
    
    // contrib to 2nd order neighbours
    if(newcell!=NULL){
      for(i2=0;i2<3;i2++){
	idx2=pow(2,i1)+pow(2,i2);
	if(i2==i1) continue;
	if(visit[idx2]) continue;

	contrib=vcont[idx2];
	getcellnei(newcell->idx, vnei2, vcell2);
	newoct=cell2oct(newcell);
	if(vnei2[neip[i2]]==6){
	  accel+=newoct->cell[vcell2[neip[i2]]].gdata.p*contrib;
	  newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	  visit[idx2]=1;
	  tot+=contrib;
	  ntot++;
	}
	else{
	  if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	    accel+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].gdata.p*contrib;
	    newcell2=&(newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	    tot+=contrib;
	    ntot++;
	  }
	  else{
	    accel=0;
	    hres=0;
	    newcell2=NULL;
	    break;
	  }
	}
      }

      // contrib to 3rd order neighbors
      if(newcell2!=NULL){
	for(i3=0;i3<3;i3++){
	  if((i3==i1)||(i3==i2)) continue;
	  if(visit[7]) continue;
	  contrib=vcont[7];
	  getcellnei(newcell2->idx, vnei3, vcell3);
	  newoct2=cell2oct(newcell2);
	  if(vnei3[neip[i3]]==6){
	    accel+=newoct2->cell[vcell3[neip[i3]]].gdata.p*contrib;
	    visit[7]=1;
	    tot+=contrib;
	    ntot++;

	  }
	  else{
	    if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
	      accel+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].gdata.p*contrib;
	      visit[7]=1;
	      tot+=contrib;
	      ntot++;
	    }
	    else{
	      accel=0;
	      hres=0;
	      break;
	    }
	  }
	}
      }
    }
  }

  // we must recompute the force from level-1 if the particle is not deep enough LOWRES
  if(hres!=1){
    
    // Getting the new curoct at low resolution
    curoctlr=cell2oct(curoct->parent);
    icell=curoct->parent->idx;
    dxcur=1./pow(2,curoctlr->level);

    // start again
    
    xc=curoctlr->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
    yc=curoctlr->y+((icell/2)%2)*dxcur+dxcur/2;
    zc=curoctlr->z+( icell/4   )*dxcur+dxcur/2; 
    
    // the CIC weights
    tx=fabs((curp->x-xc)/dxcur);
    ty=fabs((curp->y-yc)/dxcur);
    tz=fabs((curp->z-zc)/dxcur);
	  
    dx=1.-tx;
    dy=1.-ty;
    dz=1.-tz;

    for(i1=0;i1<8;i1++) visit[i1]=0;

    vcont[0]=dz*dy*dx; //000
    vcont[1]=dz*dy*tx; //001
    vcont[2]=dz*ty*dx; //010
    vcont[4]=tz*dy*dx; //100

    vcont[3]=dz*ty*tx; //011
    vcont[6]=tz*ty*dx; //110
    vcont[5]=tz*dy*tx; //101

    vcont[7]=tz*ty*tx; //111


    // contrib from current cell 000 =====
    REAL tot=0;
    int ntot=0;
    contrib=vcont[0];
    if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoctlr->cell[icell].gdata.p*contrib;
    }
  
  
    // contribs to cardinal neighbors
    getcellnei(icell, vnei, vcell);
    for(i1=0;i1<3;i1++){
      idx1=pow(2,i1);
      contrib=vcont[idx1];
      if(vnei[neip[i1]]==6){
	accel+=curoctlr->cell[vcell[neip[i1]]].gdata.p*contrib;
	newcell=&(curoctlr->cell[vcell[neip[i1]]]);
	tot+=contrib;
	ntot++;
      }
      else{
	if(curoctlr->nei[vnei[neip[i1]]]->child!=NULL){
	  accel+=curoctlr->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].gdata.p*contrib;
	  newcell=&(curoctlr->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]]);
	  tot+=contrib;
	  ntot++;

	}
      }
    
      // contrib to 2nd order neighbours
      if(newcell!=NULL){
	for(i2=0;i2<3;i2++){
	  idx2=pow(2,i1)+pow(2,i2);
	  if(i2==i1) continue;
	  if(visit[idx2]) continue;

	  contrib=vcont[idx2];
	  getcellnei(newcell->idx, vnei2, vcell2);
	  newoct=cell2oct(newcell);
	  if(vnei2[neip[i2]]==6){
	    accel+=newoct->cell[vcell2[neip[i2]]].gdata.p*contrib;
	    newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	    tot+=contrib;
	    ntot++;
	  }
	  else{
	    if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	      accel+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].gdata.p*contrib;
	      newcell2=&(newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]]);
	      visit[idx2]=1;
	      tot+=contrib;
	      ntot++;
	    }
	  }
	}

	// contrib to 3rd order neighbors
	if(newcell2!=NULL){
	  for(i3=0;i3<3;i3++){
	    if((i3==i1)||(i3==i2)) continue;
	    if(visit[7]) continue;
	    contrib=vcont[7];
	    getcellnei(newcell2->idx, vnei3, vcell3);
	    newoct2=cell2oct(newcell2);
	    if(vnei3[neip[i3]]==6){
	      accel+=newoct2->cell[vcell3[neip[i3]]].gdata.p*contrib;
	      visit[7]=1;
	      tot+=contrib;
	      ntot++;

	    }
	    else{
	      if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
		accel+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].gdata.p*contrib;
		visit[7]=1;
		tot+=contrib;
		ntot++;
	      }
	    }
	  }
	}
      }
    }
  }
  
  
#ifdef PERFECT
   if(curp->mass<0.5) {
    REAL r=sqrt((curp->x-0.5)*(curp->x-0.5)+(curp->y-0.5)*(curp->y-0.5)+(curp->z-0.5)*(curp->z-0.5));
    accel=-1./r;
   }
#endif


  // we return the energy (stored in accel) 

  return accel;

}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------


  // ==================================== performing the CIC assignement
  
void call_cic(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu){

  int level;
  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  int inei,inei2;
  REAL dxcur;

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;

  int vnei[6],vcell[6];
  int vnei2[6],vcell2[6];
  int il,ip;

  double t[10];

  if(cpu->rank==0) printf("==> start CIC on CPU\n");

  // First we clean the density (or setup the density with hydro values)
  for(level=levelmax;level>=levelcoarse;level--)
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
#ifdef WHYDRO
	  for(icell=0;icell<8;icell++) curoct->cell[icell].density=curoct->cell[icell].d;
#else
	  for(icell=0;icell<8;icell++) curoct->cell[icell].density=0.;
#endif
	}while(nextoct!=NULL);
    }

#ifdef PIC
  //second start CIC
  for(level=levelcoarse;level<=levelmax;level++)
    {
      dxcur=1./pow(2,level);
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  
	  // we skip octs which do not belong to the current CPU (they will be considered through mpi)
	  if(curoct->cpu!=cpu->rank) continue;
	  
	  //==== FIRST WE CONSIDER THE PARTICLES INSIDE THE BOUNDARIES AT LEVEL L
	  for(icell=0;icell<8;icell++) // looping over cells in oct
	    {
	      nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
	      if(nexp!=NULL){ 
		/* if(level==levelcoarse){ */
		/*   printf("ouhla\n"); */
		/*   //		  abort(); */
		/* } */
		do{  
		  curp=nexp; 
		  nexp=curp->next; 
		  part2cell_cic(curp, curoct, icell,1); 
		}while(nexp!=NULL); 
	      }
	      else if(curoct->cell[icell].child!=NULL){
	      // THEN WE LOOK FOR THE PARTICLES IN THE CHILD OCTS
		cic_child(curoct->cell[icell].child,curoct,icell);
	      }
#if 1
	      
	      if(level>levelcoarse){
		//==== SECOND WE CONSIDER THE PARTICLES INSIDE THE NEIGHBORS AT LEVEL L-1

		int visit[27];
		int iv;
		for(inei=0;inei<27;inei++) visit[inei]=0;
		int tot=0;
		// first the cartesian neighbors (6)
		for(inei=0;inei<6;inei++){
		  nexp=curoct->nei[inei]->phead; //sweeping the particles of the neighbour cell at level l-1
		  //		  iv=13+pow(3,inei/2)*(2*(inei%2)-1);
		  iv=13+pow(3,inei>>1)*(2*(inei&1)-1);
		  visit[iv]=1;
		  if(nexp!=NULL){
		    do{ 
		      curp=nexp;
		      nexp=curp->next;
		      part2cell_cic(curp, curoct, icell,0);
		    }while(nexp!=NULL);
		  }
		  // second the fundamental plane (4)
		  // each of the 6 cardinal neighbors will probe 4 neighbors
		
		  newcell=curoct->nei[inei];
		  newoct=cell2oct(newcell); // we get the parent oct;
		  getcellnei(newcell->idx, vnei, vcell); // we get its neighbors
		  for(il=0;il<6;il++){
		    //iv=pow(3,inei/2)*(2*(inei%2)-1)+pow(3,il/2)*(2*(il%2)-1);
		    iv=pow(3,inei>>1)*(2*(inei&1)-1)+pow(3,il>>1)*(2*(il&1)-1);
		    if((il>>1)==(inei>>1)) continue;
		    if(visit[iv]) continue;
		    if(vnei[il]==6){
		      visit[iv]=1;
		      nexp=newoct->cell[vcell[il]].phead;
		      newcell2=&(newoct->cell[vcell[il]]);
		      if(nexp!=NULL){
			do{ 
			  curp=nexp;
			  nexp=curp->next;
			  part2cell_cic(curp, curoct, icell,0);
			}while(nexp!=NULL);
		      }
		    }
		    else{
		      if(newoct->nei[vnei[il]]->child!=NULL){
			nexp=newoct->nei[vnei[il]]->child->cell[vcell[il]].phead;
			newcell2=&(newoct->nei[vnei[il]]->child->cell[vcell[il]]);
			visit[iv]=1;
			if(nexp!=NULL){
			  do{ 
			    curp=nexp;
			    nexp=curp->next;
			    part2cell_cic(curp, curoct, icell,0);
			  }while(nexp!=NULL);
			}
		      }
		      else{
			newcell2=NULL;
		      }
		    }

		    // ecah of the 4 side neighbors will mark 2 corners
		    if(newcell2!=NULL){
		      newoct2=cell2oct(newcell2);
		      getcellnei(newcell2->idx, vnei2, vcell2); // we get its neighbors
		      for(ip=0;ip<6;ip++){
			//iv=pow(3,inei/2)*(2*(inei%2)-1)+pow(3,il/2)*(2*(il%2)-1)+pow(3,ip/2)*(2*(ip%2)-1);
			iv=pow(3,inei>>1)*(2*(inei&1)-1)+pow(3,il>>1)*(2*(il&1)-1)+pow(3,ip>>1)*(2*(ip&1)-1);
			//if(((ip/2)==(il/2))||((ip/2)==(inei/2))) continue;
			if(((ip>>1)==(il>>1))||((ip>>1)==(inei>>1))) continue;
			if(visit[iv]) continue;
			if(vnei2[ip]==6){
			  visit[iv]=1;
			  nexp=newoct2->cell[vcell2[ip]].phead;
			  if(nexp!=NULL){
			    do{ 
			      curp=nexp;
			      nexp=curp->next;
			      part2cell_cic(curp, curoct, icell,0);
			    }while(nexp!=NULL);
			  }
			}
			else{
			  if(newoct2->nei[vnei2[ip]]->child!=NULL){
			    visit[iv]=1;
			    nexp=newoct2->nei[vnei2[ip]]->child->cell[vcell2[ip]].phead;
			    if(nexp!=NULL){
			      do{ 
				curp=nexp;
				nexp=curp->next;
				part2cell_cic(curp, curoct, icell,0);
			      }while(nexp!=NULL);
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }

#endif
	     
	    }
	}while(nextoct!=NULL);
    }
#endif
  //  printf("great total=%f\n",toto);
}

// ====================================================================================
// ====================================================================================

// COULD BE REWRITTEN IN A RECURSIVE FASHION (SEE gatherstencil)

void L_cic(int level,struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu){

  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  struct PART *curp;
  struct PART *nexp;
  int inei,inei2;
  REAL dxcur;

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;

  int vnei[6],vcell[6];
  int vnei2[6],vcell2[6];
  int il,ip;

  double t[10];

  dxcur=1./pow(2,level);
  nextoct=firstoct[level-1];
  do // sweeping level
    {
      curoct=nextoct;
      nextoct=curoct->next;
	  
      // we skip octs which do not belong to the current CPU (they will be considered through mpi)
      if(curoct->cpu!=cpu->rank) continue;
	  
      //==== FIRST WE CONSIDER THE PARTICLES INSIDE THE BOUNDARIES AT LEVEL L
      for(icell=0;icell<8;icell++) // looping over cells in oct
	{
	  nexp=curoct->cell[icell].phead; //sweeping the particles of the current cell */
	  if(nexp!=NULL){ 
	    do{  
	      curp=nexp; 
	      nexp=curp->next; 
	      part2cell_cic(curp, curoct, icell,1); 
	    }while(nexp!=NULL); 
	  }
	  else if(curoct->cell[icell].child!=NULL){
	    // THEN WE LOOK FOR THE PARTICLES IN THE CHILD OCTS
	    cic_child(curoct->cell[icell].child,curoct,icell);
	  }

	  if(level>param->lcoarse){
	    //==== SECOND WE CONSIDER THE PARTICLES INSIDE THE NEIGHBORS AT LEVEL L-1

	    int visit[27];
	    int iv;
	    for(inei=0;inei<27;inei++) visit[inei]=0;
	    int tot=0;
	    // first the cartesian neighbors (6)
	    for(inei=0;inei<6;inei++){
	      nexp=curoct->nei[inei]->phead; //sweeping the particles of the neighbour cell at level l-1
	      iv=13+pow(3,inei>>1)*(2*(inei&1)-1);
	      visit[iv]=1;
	      if(nexp!=NULL){
		do{ 
		  curp=nexp;
		  nexp=curp->next;
		  part2cell_cic(curp, curoct, icell,0);
		}while(nexp!=NULL);
	      }
	      // second the fundamental plane (4)
	      // each of the 6 cardinal neighbors will probe 4 neighbors
		
	      newcell=curoct->nei[inei];
	      newoct=cell2oct(newcell); // we get the parent oct;
	      getcellnei(newcell->idx, vnei, vcell); // we get its neighbors
	      for(il=0;il<6;il++){
		iv=pow(3,inei>>1)*(2*(inei&1)-1)+pow(3,il>>1)*(2*(il&1)-1);
		if((il>>1)==(inei>>1)) continue;
		if(visit[iv]) continue;
		if(vnei[il]==6){
		  visit[iv]=1;
		  nexp=newoct->cell[vcell[il]].phead;
		  newcell2=&(newoct->cell[vcell[il]]);
		  if(nexp!=NULL){
		    do{ 
		      curp=nexp;
		      nexp=curp->next;
		      part2cell_cic(curp, curoct, icell,0);
		    }while(nexp!=NULL);
		  }
		}
		else{
		  if(newoct->nei[vnei[il]]->child!=NULL){
		    nexp=newoct->nei[vnei[il]]->child->cell[vcell[il]].phead;
		    newcell2=&(newoct->nei[vnei[il]]->child->cell[vcell[il]]);
		    visit[iv]=1;
		    if(nexp!=NULL){
		      do{ 
			curp=nexp;
			nexp=curp->next;
			part2cell_cic(curp, curoct, icell,0);
		      }while(nexp!=NULL);
		    }
		  }
		  else{
		    newcell2=NULL;
		  }
		}

		// ecah of the 4 side neighbors will mark 2 corners
		if(newcell2!=NULL){
		  newoct2=cell2oct(newcell2);
		  getcellnei(newcell2->idx, vnei2, vcell2); // we get its neighbors
		  for(ip=0;ip<6;ip++){
		    iv=pow(3,inei>>1)*(2*(inei&1)-1)+pow(3,il>>1)*(2*(il&1)-1)+pow(3,ip>>1)*(2*(ip&1)-1);
		    if(((ip>>1)==(il>>1))||((ip>>1)==(inei>>1))) continue;
		    if(visit[iv]) continue;
		    if(vnei2[ip]==6){
		      visit[iv]=1;
		      nexp=newoct2->cell[vcell2[ip]].phead;
		      if(nexp!=NULL){
			do{ 
			  curp=nexp;
			  nexp=curp->next;
			  part2cell_cic(curp, curoct, icell,0);
			}while(nexp!=NULL);
		      }
		    }
		    else{
		      if(newoct2->nei[vnei2[ip]]->child!=NULL){
			visit[iv]=1;
			nexp=newoct2->nei[vnei2[ip]]->child->cell[vcell2[ip]].phead;
			if(nexp!=NULL){
			  do{ 
			    curp=nexp;
			    nexp=curp->next;
			    part2cell_cic(curp, curoct, icell,0);
			  }while(nexp!=NULL);
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }

	     
	}
    }while(nextoct!=NULL);




}




// ===================================================

void L_clean_dens(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu){

  struct OCT *curoct;
  struct OCT *nextoct;
  int icell;

  // --------------- setting the first oct of the level
  nextoct=firstoct[level-1];
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next; 

      for(icell=0;icell<8;icell++) curoct->cell[icell].density=0.;
      
    }while(nextoct!=NULL);
  }
}


#endif
