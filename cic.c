#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"


//==================================================================

void part2cell_cic(struct PART *curp, struct OCT *curoct, int icell, char full)
{
  float xc,yc,zc;
  //float dxcur=1./pow(2,curoct->level);
  float dxcur=1./(1<<curoct->level);
  
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  int neip[3];
  float tx,ty,tz;
  float dx,dy,dz;
  char fullok=0;

  int visit[8]={0,0,0,0,0,0,0,0};
  float vcont[8];

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;

  int i1,i2,i3;
  int idx1,idx2; // scan the cic cell being analysed

  /* struct OCT **cicoct; */
  /* char ciccell[8]; */
  /* cicoct=(struct OCT **)calloc(8,sizeof(struct OCT *)); */

  xc=curoct->x+( icell   & 1)*dxcur+dxcur*0.5; // coordinates of the cell center 
  yc=curoct->y+((icell>>1)& 1)*dxcur+dxcur*0.5;
  zc=curoct->z+( icell>>2   )*dxcur+dxcur*0.5; 

  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);

  /* // getting the neighbors */
  /* getcellnei(icell, vnei, vcell); */
		  
		  
  /* // here we denote the offset in ZYX */
  /* ciccell[0]=icell;          //cell 000 */
  /* ciccell[1]=vcell[neip[0]]; //cell 001 */
  /* ciccell[2]=vcell[neip[1]]; //cell 010 */
  /* ciccell[4]=vcell[neip[2]]; //cell 100 */
  /* cicoct[0]=curoct; */

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
  float dunit=curp->mass/(dxcur*dxcur*dxcur);
  vcont[0]=dz*dy*dx*dunit; //000
  vcont[1]=dz*dy*tx*dunit; //001
  vcont[2]=dz*ty*dx*dunit; //010
  vcont[4]=tz*dy*dx*dunit; //100

  vcont[3]=dz*ty*tx*dunit; //011
  vcont[6]=tz*ty*dx*dunit; //110
  vcont[5]=tz*dy*tx*dunit; //101

  vcont[7]=tz*ty*tx*dunit; //111


  // contrib to current cell 000 =====

  float tot=0;
  int ntot=0;
  //contrib=vcont[0];
  //  if((contrib<=1.)&&(contrib>=0.)){
    if(full||fullok){
#ifndef NGP
      curoct->cell[icell].density+=vcont[0];
#else
      curoct->cell[icell].density+=1./pow(dxcur,3)*curp->mass;
#endif
      /* tot+=contrib; */
      /* ntot++; */
    }
    //}

#ifndef NGP
  if(full){
    // contribs to cardinal neighbors
    getcellnei(icell, vnei, vcell);
    for(i1=0;i1<3;i1++){
      //idx1=pow(2,i1);
      idx1=1<<i1;
      //contrib=vcont[idx1];
      if(vnei[neip[i1]]==6){
	curoct->cell[vcell[neip[i1]]].density+=vcont[idx1];
	newcell=&(curoct->cell[vcell[neip[i1]]]);
	/* tot+=contrib; */
	/* ntot++; */
      }
      else{
	if(curoct->nei[vnei[neip[i1]]]->child!=NULL){
	  curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].density+=vcont[idx1];
	  newcell=&(curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]]);
	  /* tot+=contrib; */
	  /* ntot++; */

	}
	else{
	  newcell=NULL;
	}
      }
    
      // contrib to 2nd order neighbours
      if(newcell!=NULL){
	for(i2=0;i2<3;i2++){
	  //idx2=pow(2,i1)+pow(2,i2);
	  idx2=(1<<i1)+(1<<i2);
	  if(i2==i1) continue;
	  if(visit[idx2]) continue;

	  //contrib=vcont[idx2];
	  getcellnei(newcell->idx, vnei2, vcell2);
	  newoct=cell2oct(newcell);
	  if(vnei2[neip[i2]]==6){
	    newoct->cell[vcell2[neip[i2]]].density+=vcont[idx2];
	    newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	    /* tot+=contrib; */
	    /* ntot++; */
	  }
	  else{
	    if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	      newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].density+=vcont[idx2];
	      newcell2=&(newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]]);
	      visit[idx2]=1;
	      /* tot+=contrib; */
	      /* ntot++; */
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
	    //contrib=vcont[7];
	    getcellnei(newcell2->idx, vnei3, vcell3);
	    newoct2=cell2oct(newcell2);
	    if(vnei3[neip[i3]]==6){
	      newoct2->cell[vcell3[neip[i3]]].density+=vcont[7];
	      visit[7]=1;
	      /* tot+=contrib; */
	      /* ntot++; */

	    }
	    else{
	      if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
		newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].density+=vcont[7];
		visit[7]=1;
		/* tot+=contrib; */
		/* ntot++; */
	      }
	    }
	  }
	}
      }
    }
  }

  /* if(ntot>8) abort(); */
  /* if((dx+tx)!=1.) abort(); */
  /* if((dy+ty)!=1.) abort(); */
  /* if((dz+tz)!=1.) abort(); */

#endif

}

//=====================================================================================
//=====================================================================================

void cell2part_cic(struct PART *curp, struct OCT *curoct, int icell, char dir, float dt)
{
  float xc,yc,zc;
  float dxcur=1./pow(2,curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  int neip[3];
  float tx,ty,tz;
  float dx,dy,dz;
  float contrib;
  struct OCT *curoctlr;
  float accel=0.;
  char hres=1; // could be switched to hres=0 if particle is not deep enough

  int visit[8]={0,0,0,0,0,0,0,0};
  float vcont[8];

  int i1,i2,i3;
  int idx1,idx2; // scan the cic cell being analysed

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;
  
#ifndef NEWFORCE
  struct OCT **cicoct;
  char ciccell[8];
  cicoct=(struct OCT **)calloc(8,sizeof(struct OCT *));
#endif

  xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
  yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
  zc=curoct->z+( icell/4   )*dxcur+dxcur/2; 

		  
  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);
		   

#ifndef NEWFORCE
  // getting the neighbors
  getcellnei(icell, vnei, vcell);

 // here we denote the offset in ZYX
  ciccell[0]=icell;          //cell 000
  ciccell[1]=vcell[neip[0]]; //cell 001
  ciccell[2]=vcell[neip[1]]; //cell 010
  ciccell[4]=vcell[neip[2]]; //cell 100
  cicoct[0]=curoct;
#endif

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
  float tot=0;
  int ntot=0;

  contrib=vcont[0];
  if((contrib<=1.)&&(contrib>=0.)){
    accel+=curoct->cell[icell].temp*contrib;
  }
  

#ifdef NEWFORCE
  // contribs to cardinal neighbors
  getcellnei(icell, vnei, vcell);
  for(i1=0;i1<3;i1++){
    idx1=pow(2,i1);
    contrib=vcont[idx1];
    if(vnei[neip[i1]]==6){
      accel+=curoct->cell[vcell[neip[i1]]].temp*contrib;
      newcell=&(curoct->cell[vcell[neip[i1]]]);
      tot+=contrib;
      ntot++;
    }
    else{
      if(curoct->nei[vnei[neip[i1]]]->child!=NULL){
	accel+=curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].temp*contrib;
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
	  accel+=newoct->cell[vcell2[neip[i2]]].temp*contrib;
	  newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	  visit[idx2]=1;
	  tot+=contrib;
	  ntot++;
	}
	else{
	  if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	    accel+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].temp*contrib;
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
	    accel+=newoct2->cell[vcell3[neip[i3]]].temp*contrib;
	    visit[7]=1;
	    tot+=contrib;
	    ntot++;

	  }
	  else{
	    if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
	      accel+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].temp*contrib;
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
    float tot=0;
    int ntot=0;
    contrib=vcont[0];
    if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoctlr->cell[icell].temp*contrib;
    }
  
  
    // contribs to cardinal neighbors
    getcellnei(icell, vnei, vcell);
    for(i1=0;i1<3;i1++){
      idx1=pow(2,i1);
      contrib=vcont[idx1];
      if(vnei[neip[i1]]==6){
	accel+=curoctlr->cell[vcell[neip[i1]]].temp*contrib;
	newcell=&(curoctlr->cell[vcell[neip[i1]]]);
	tot+=contrib;
	ntot++;
      }
      else{
	if(curoctlr->nei[vnei[neip[i1]]]->child!=NULL){
	  accel+=curoctlr->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].temp*contrib;
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
	    accel+=newoct->cell[vcell2[neip[i2]]].temp*contrib;
	    newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	    tot+=contrib;
	    ntot++;
	  }
	  else{
	    if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	      accel+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].temp*contrib;
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
	      accel+=newoct2->cell[vcell3[neip[i3]]].temp*contrib;
	      visit[7]=1;
	      tot+=contrib;
	      ntot++;

	    }
	    else{
	      if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
		accel+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].temp*contrib;
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
  
#endif


#ifndef NEWFORCE
  // contrib from 100 cell ===========================================================
      
  if(hres==1){
    contrib=tx*dy*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[0]]==6){
  	accel+=curoct->cell[vcell[neip[0]]].temp*contrib;
  	cicoct[1]=curoct;
      }
      else {
  	if(curoct->nei[vnei[neip[0]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].temp*contrib;
  	  cicoct[1]=curoct->nei[vnei[neip[0]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  // contrib from 010 cell ===========================================================
  if(hres==1){
    contrib=dx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[1]]==6){
  	accel+=curoct->cell[vcell[neip[1]]].temp*contrib;
  	cicoct[2]=curoct;
      }
      else{
  	if(curoct->nei[vnei[neip[1]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].temp*contrib;
  	  cicoct[2]=curoct->nei[vnei[neip[1]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  // contrib from 001 cell ===========================================================

  if(hres==1){
    contrib=dx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[2]]==6){
  	accel+=curoct->cell[vcell[neip[2]]].temp*contrib;
  	cicoct[4]=curoct;
      }
      else{
  	if(curoct->nei[vnei[neip[2]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].temp*contrib;
  	  cicoct[4]=curoct->nei[vnei[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  //contrib from 110 cell 2 paths ======================================================
  
  contrib=tx*ty*dz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){
      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[3]=vcell2[neip[1]];
  	if(vnei2[neip[1]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[1]]].temp*contrib;
  	  cicoct[3]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].temp*contrib;
  	  cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
      

  //contrib from 101 cell 2 paths ======================================================
  contrib=tx*dy*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[5]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[5]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
	
      }
    }
  }

  //contrib from 011 cell 2 paths ======================================================
  contrib=dx*ty*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[2]!=NULL){
  	getcellnei(ciccell[2], vnei2, vcell2);
  	ciccell[6]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[2]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[6]=cicoct[2];
  	}
  	else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}

      }
    }
  }
      
  // contrib from 111 ======================================================
  contrib=tx*ty*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){
		  
      if(cicoct[3]!=NULL){
  	getcellnei(ciccell[3], vnei2, vcell2);
  	ciccell[7]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[3]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[7]=cicoct[3];
  	}
  	else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
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
    
    // getting the neighbors
    getcellnei(icell, vnei, vcell);
    
    // here we compute the indexes of the direct neighbors which are involved in the cic
    neip[0]=(curp->x<xc?0:1);
    neip[1]=(curp->y<yc?2:3);
    neip[2]=(curp->z<zc?4:5);
    
    // here we denote the offset in ZYX
    ciccell[0]=icell;          //cell 000
    ciccell[1]=vcell[neip[0]]; //cell 001
    ciccell[2]=vcell[neip[1]]; //cell 010
    ciccell[4]=vcell[neip[2]]; //cell 100
    cicoct[0]=curoctlr;

    // the CIC weights
    tx=fabs((curp->x-xc)/dxcur);
    ty=fabs((curp->y-yc)/dxcur);
    tz=fabs((curp->z-zc)/dxcur);
    
    dx=1.-tx;
    dy=1.-ty;
    dz=1.-tz;
    
    // contrib from current cell 000 =====
    contrib=dx*dy*dz;
    if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoctlr->cell[icell].temp*contrib;
      //curoct->cell[icell].density+=contrib/pow(dxcur,3);
    }
  

    // contrib from 100 cell ===========================================================
      
    contrib=tx*dy*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[0]]==6){
  	accel+=curoctlr->cell[vcell[neip[0]]].temp*contrib;
  	cicoct[1]=curoctlr;
      }
      else {
  	if(curoctlr->nei[vnei[neip[0]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].temp*contrib;
  	  cicoct[1]=curoctlr->nei[vnei[neip[0]]]->child;
  	}
      }
    }

    // contrib from 010 cell ===========================================================
    contrib=dx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[1]]==6){
  	accel+=curoctlr->cell[vcell[neip[1]]].temp*contrib;
  	cicoct[2]=curoctlr;
      }
      else{
  	if(curoctlr->nei[vnei[neip[1]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].temp*contrib;
  	  cicoct[2]=curoctlr->nei[vnei[neip[1]]]->child;
  	}
      }
    }
    
    // contrib from 001 cell ===========================================================

    contrib=dx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[2]]==6){
  	accel+=curoctlr->cell[vcell[neip[2]]].temp*contrib;
  	cicoct[4]=curoctlr;
      }
      else{
  	if(curoctlr->nei[vnei[neip[2]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].temp*contrib;
  	  cicoct[4]=curoctlr->nei[vnei[neip[2]]]->child;
  	}
      }
    }

    //contrib from 110 cell 2 paths ======================================================
  
    contrib=tx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[3]=vcell2[neip[1]];
  	if(vnei2[neip[1]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[1]]].temp*contrib;
  	  cicoct[3]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].temp*contrib;
  	  cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
  	}
      }
    }
    
    //contrib from 101 cell 2 paths ======================================================
    contrib=tx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[5]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[5]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }

    //contrib from 011 cell 2 paths ======================================================
    contrib=dx*ty*tz;
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[2]!=NULL){
  	getcellnei(ciccell[2], vnei2, vcell2);
  	ciccell[6]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[2]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[6]=cicoct[2];
  	}
  	else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }
      
    // contrib from 111 ======================================================
    contrib=tx*ty*tz;
    if((contrib<=1.)&&(contrib>=0)){
		  
      if(cicoct[3]!=NULL){
  	getcellnei(ciccell[3], vnei2, vcell2);
  	ciccell[7]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[3]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[7]=cicoct[3];
  	}
  	else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }
  }

#endif
  
  // Once we have the acceleration we can compute the velocity
  switch(dir){
  case(0):
    curp->vx+=-accel*dt;
    if(curp->mass<0.5) printf("dir=%d accel=%e vy=%e\n",dir,accel,curp->vx);
    break;
  case(1):
    curp->vy+=-accel*dt;
    if(curp->mass<0.5) printf("dir=%d accel=%e vy=%e\n",dir,accel,curp->vy);
    break;
  case(2):
    curp->vz+=-accel*dt;
    break;
  }
#ifndef NEWFORCE
  free(cicoct);
#endif

}

#ifdef WGPU
void cell2part_cic_GPU(struct PART *curp, struct OCT *curoct, int icell, char dir, float dt)
{
  float xc,yc,zc;
  float xm,ym,zm;
  float dxcur=1./pow(2,curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  int neip[3];
  float tx,ty,tz;
  float dx,dy,dz;
  float contrib;
  struct OCT *curoctlr;
  float accel=0.;
  char hres=1; // could be switched to hres=0 if particle is not deep enough

  int visit[8]={0,0,0,0,0,0,0,0};
  float vcont[8];

  int i1,i2,i3;
  int idx1,idx2; // scan the cic cell being analysed

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;
  
#ifndef NEWFORCE
  struct OCT **cicoct;
  char ciccell[8];
  cicoct=(struct OCT **)calloc(8,sizeof(struct OCT *));
#endif


  xc=curoct->x+( icell   %2)*dxcur; // coordinates of the corner cell 
  yc=curoct->y+((icell/2)%2)*dxcur;
  zc=curoct->z+( icell/4   )*dxcur; 

		  
  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);

  if((curp->x<xc)||(curp->y<yc)){
    printf("x %e %e %e\n",curp->x,xc,curp->x-xc);
    printf("x %e %e %e\n",curp->y,yc,curp->y-yc);
    abort();
  }
		   

#ifndef NEWFORCE
  // getting the neighbors
  getcellnei(icell, vnei, vcell);

 // here we denote the offset in ZYX
  ciccell[0]=icell;          //cell 000
  ciccell[1]=vcell[neip[0]]; //cell 001
  ciccell[2]=vcell[neip[1]]; //cell 010
  ciccell[4]=vcell[neip[2]]; //cell 100
  cicoct[0]=curoct;
#endif

  // the CIC weights
  tx=fabs((curp->x-xc)/dxcur);
  ty=fabs((curp->y-yc)/dxcur);
  tz=fabs((curp->z-zc)/dxcur);

  if(curp->mass<0.5) printf("%e %e %e x=%e xc=%e xo=%e\n",tx,ty,tz,curp->x,xc,curoct->x);
	  
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
  float tot=0;
  int ntot=0;

  contrib=vcont[0];
  if((contrib<=1.)&&(contrib>=0.)){
    accel+=curoct->cell[icell].temp*contrib;
    tot+=contrib;
    ntot++;
    visit[0]=1;
  }
  

#ifdef NEWFORCE
  // contribs to cardinal neighbors
  getcellnei(icell, vnei, vcell);
  for(i1=0;i1<3;i1++){
    idx1=pow(2,i1);
    contrib=vcont[idx1];
    visit[idx1]=1;
    if(vnei[neip[i1]]==6){
      accel+=curoct->cell[vcell[neip[i1]]].temp*contrib;
      newcell=&(curoct->cell[vcell[neip[i1]]]);
      tot+=contrib;
      ntot++;
    }
    else{
      if(curoct->nei[vnei[neip[i1]]]->child!=NULL){
	accel+=curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].temp*contrib;
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
	  accel+=newoct->cell[vcell2[neip[i2]]].temp*contrib;
	  newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	  visit[idx2]=1;
	  tot+=contrib;
	  ntot++;
	}
	else{
	  if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	    accel+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].temp*contrib;
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
	    accel+=newoct2->cell[vcell3[neip[i3]]].temp*contrib;
	    visit[7]=1;
	    tot+=contrib;
	    ntot++;

	  }
	  else{
	    if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
	      accel+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].temp*contrib;
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

  if(ntot!=8){
    int ii;
    printf("tot=%e ntot=%d m=%e v0=%e\n",tot,ntot,curp->mass,vcont[0]);
    for(ii=0;ii<8;ii++) printf("%d ",visit[ii]);
    printf("\n");
    abort();
  }

  // we must recompute the force from level-1 if the particle is not deep enough LOWRES
  if(hres!=1){
    
    // Getting the new curoct at low resolution
    curoctlr=cell2oct(curoct->parent);
    icell=curoct->parent->idx;
    dxcur=1./pow(2,curoctlr->level);

    // start again
    
    xc=curoctlr->x+( icell   %2)*dxcur; // coordinates of the cell center 
    yc=curoctlr->y+((icell/2)%2)*dxcur;
    zc=curoctlr->z+( icell/4   )*dxcur; 
    
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
    float tot=0;
    int ntot=0;
    contrib=vcont[0];
    if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoctlr->cell[icell].temp*contrib;
    }
  
  
    // contribs to cardinal neighbors
    getcellnei(icell, vnei, vcell);
    for(i1=0;i1<3;i1++){
      idx1=pow(2,i1);
      contrib=vcont[idx1];
      if(vnei[neip[i1]]==6){
	accel+=curoctlr->cell[vcell[neip[i1]]].temp*contrib;
	newcell=&(curoctlr->cell[vcell[neip[i1]]]);
	tot+=contrib;
	ntot++;
      }
      else{
	if(curoctlr->nei[vnei[neip[i1]]]->child!=NULL){
	  accel+=curoctlr->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].temp*contrib;
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
	    accel+=newoct->cell[vcell2[neip[i2]]].temp*contrib;
	    newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	    tot+=contrib;
	    ntot++;
	  }
	  else{
	    if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	      accel+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].temp*contrib;
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
	      accel+=newoct2->cell[vcell3[neip[i3]]].temp*contrib;
	      visit[7]=1;
	      tot+=contrib;
	      ntot++;

	    }
	    else{
	      if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
		accel+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].temp*contrib;
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
  
#endif


#ifndef NEWFORCE
  // contrib from 100 cell ===========================================================
      
  if(hres==1){
    contrib=tx*dy*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[0]]==6){
  	accel+=curoct->cell[vcell[neip[0]]].temp*contrib;
  	cicoct[1]=curoct;
      }
      else {
  	if(curoct->nei[vnei[neip[0]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].temp*contrib;
  	  cicoct[1]=curoct->nei[vnei[neip[0]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  // contrib from 010 cell ===========================================================
  if(hres==1){
    contrib=dx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[1]]==6){
  	accel+=curoct->cell[vcell[neip[1]]].temp*contrib;
  	cicoct[2]=curoct;
      }
      else{
  	if(curoct->nei[vnei[neip[1]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].temp*contrib;
  	  cicoct[2]=curoct->nei[vnei[neip[1]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  // contrib from 001 cell ===========================================================

  if(hres==1){
    contrib=dx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[2]]==6){
  	accel+=curoct->cell[vcell[neip[2]]].temp*contrib;
  	cicoct[4]=curoct;
      }
      else{
  	if(curoct->nei[vnei[neip[2]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].temp*contrib;
  	  cicoct[4]=curoct->nei[vnei[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  //contrib from 110 cell 2 paths ======================================================
  
  contrib=tx*ty*dz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){
      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[3]=vcell2[neip[1]];
  	if(vnei2[neip[1]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[1]]].temp*contrib;
  	  cicoct[3]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].temp*contrib;
  	  cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
      

  //contrib from 101 cell 2 paths ======================================================
  contrib=tx*dy*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[5]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[5]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
	
      }
    }
  }

  //contrib from 011 cell 2 paths ======================================================
  contrib=dx*ty*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[2]!=NULL){
  	getcellnei(ciccell[2], vnei2, vcell2);
  	ciccell[6]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[2]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[6]=cicoct[2];
  	}
  	else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}

      }
    }
  }
      
  // contrib from 111 ======================================================
  contrib=tx*ty*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){
		  
      if(cicoct[3]!=NULL){
  	getcellnei(ciccell[3], vnei2, vcell2);
  	ciccell[7]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[3]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[7]=cicoct[3];
  	}
  	else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
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
    
    // getting the neighbors
    getcellnei(icell, vnei, vcell);
    
    // here we compute the indexes of the direct neighbors which are involved in the cic
    neip[0]=(curp->x<xc?0:1);
    neip[1]=(curp->y<yc?2:3);
    neip[2]=(curp->z<zc?4:5);
    
    // here we denote the offset in ZYX
    ciccell[0]=icell;          //cell 000
    ciccell[1]=vcell[neip[0]]; //cell 001
    ciccell[2]=vcell[neip[1]]; //cell 010
    ciccell[4]=vcell[neip[2]]; //cell 100
    cicoct[0]=curoctlr;

    // the CIC weights
    tx=fabs((curp->x-xc)/dxcur);
    ty=fabs((curp->y-yc)/dxcur);
    tz=fabs((curp->z-zc)/dxcur);
    
    dx=1.-tx;
    dy=1.-ty;
    dz=1.-tz;
    
    // contrib from current cell 000 =====
    contrib=dx*dy*dz;
    if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoctlr->cell[icell].temp*contrib;
      //curoct->cell[icell].density+=contrib/pow(dxcur,3);
    }
  

    // contrib from 100 cell ===========================================================
      
    contrib=tx*dy*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[0]]==6){
  	accel+=curoctlr->cell[vcell[neip[0]]].temp*contrib;
  	cicoct[1]=curoctlr;
      }
      else {
  	if(curoctlr->nei[vnei[neip[0]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].temp*contrib;
  	  cicoct[1]=curoctlr->nei[vnei[neip[0]]]->child;
  	}
      }
    }

    // contrib from 010 cell ===========================================================
    contrib=dx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[1]]==6){
  	accel+=curoctlr->cell[vcell[neip[1]]].temp*contrib;
  	cicoct[2]=curoctlr;
      }
      else{
  	if(curoctlr->nei[vnei[neip[1]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].temp*contrib;
  	  cicoct[2]=curoctlr->nei[vnei[neip[1]]]->child;
  	}
      }
    }
    
    // contrib from 001 cell ===========================================================

    contrib=dx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[2]]==6){
  	accel+=curoctlr->cell[vcell[neip[2]]].temp*contrib;
  	cicoct[4]=curoctlr;
      }
      else{
  	if(curoctlr->nei[vnei[neip[2]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].temp*contrib;
  	  cicoct[4]=curoctlr->nei[vnei[neip[2]]]->child;
  	}
      }
    }

    //contrib from 110 cell 2 paths ======================================================
  
    contrib=tx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[3]=vcell2[neip[1]];
  	if(vnei2[neip[1]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[1]]].temp*contrib;
  	  cicoct[3]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].temp*contrib;
  	  cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
  	}
      }
    }
    
    //contrib from 101 cell 2 paths ======================================================
    contrib=tx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[5]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[5]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }

    //contrib from 011 cell 2 paths ======================================================
    contrib=dx*ty*tz;
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[2]!=NULL){
  	getcellnei(ciccell[2], vnei2, vcell2);
  	ciccell[6]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[2]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[6]=cicoct[2];
  	}
  	else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }
      
    // contrib from 111 ======================================================
    contrib=tx*ty*tz;
    if((contrib<=1.)&&(contrib>=0)){
		  
      if(cicoct[3]!=NULL){
  	getcellnei(ciccell[3], vnei2, vcell2);
  	ciccell[7]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[3]->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[7]=cicoct[3];
  	}
  	else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].temp*contrib;
  	  cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }
  }

#endif
  
  // Once we have the acceleration we can compute the velocity
  switch(dir){
  case(0):
    curp->vx+=-accel*dt;
    if(curp->mass<0.5) printf("dir=%d accel=%e vy=%e\n",dir,accel,curp->vx);
    break;
  case(1):
    curp->vy+=-accel*dt;
    if(curp->mass<0.5) printf("dir=%d accel=%e vy=%e\n",dir,accel,curp->vy);
    break;
  case(2):
    curp->vz+=-accel*dt;
    break;
  }
#ifndef NEWFORCE
  free(cicoct);
#endif

}
#endif

//=====================================================================================
//=====================================================================================

float cell2part_cic_egy(struct PART *curp, struct OCT *curoct, int icell)
{
  float xc,yc,zc;
  float dxcur=1./pow(2,curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int vnei3[6],vcell3[6];
  int neip[3];
  float tx,ty,tz;
  float dx,dy,dz;
  float contrib;
  struct OCT *curoctlr;
  float accel=0.;
  char hres=1; // could be switched to hres=0 if particle is not deep enough

  int visit[8]={0,0,0,0,0,0,0,0};
  float vcont[8];

  int i1,i2,i3;
  int idx1,idx2; // scan the cic cell being analysed

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;
  
#ifndef NEWFORCE
  struct OCT **cicoct;
  char ciccell[8];
  cicoct=(struct OCT **)calloc(8,sizeof(struct OCT *));
#endif

  xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
  yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
  zc=curoct->z+( icell/4   )*dxcur+dxcur/2; 

		  
  // here we compute the indexes of the direct neighbors which are involved in the cic
  neip[0]=(curp->x<xc?0:1);
  neip[1]=(curp->y<yc?2:3);
  neip[2]=(curp->z<zc?4:5);
		   

#ifndef NEWFORCE
  // getting the neighbors
  getcellnei(icell, vnei, vcell);

 // here we denote the offset in ZYX
  ciccell[0]=icell;          //cell 000
  ciccell[1]=vcell[neip[0]]; //cell 001
  ciccell[2]=vcell[neip[1]]; //cell 010
  ciccell[4]=vcell[neip[2]]; //cell 100
  cicoct[0]=curoct;
#endif

  // the CIC weights
  tx=fabs((curp->x-xc)/dxcur);
  ty=fabs((curp->y-yc)/dxcur);
  tz=fabs((curp->z-zc)/dxcur);
	  
  dx=1.-tx;
  dy=1.-ty;
  dz=1.-tz;


  vcont[0]=dz*dy*dx*curp->mass; //000
  vcont[1]=dz*dy*tx*curp->mass; //001
  vcont[2]=dz*ty*dx*curp->mass; //010
  vcont[4]=tz*dy*dx*curp->mass; //100

  vcont[3]=dz*ty*tx*curp->mass; //011
  vcont[6]=tz*ty*dx*curp->mass; //110
  vcont[5]=tz*dy*tx*curp->mass; //101

  vcont[7]=tz*ty*tx*curp->mass; //111


  // contrib from current cell 000 =====
  float tot=0;
  int ntot=0;

  contrib=vcont[0];
  if((contrib<=1.)&&(contrib>=0.)){
    accel+=curoct->cell[icell].pot*contrib;
  }
  

#ifdef NEWFORCE
  // contribs to cardinal neighbors
  getcellnei(icell, vnei, vcell);
  for(i1=0;i1<3;i1++){
    idx1=pow(2,i1);
    contrib=vcont[idx1];
    if(vnei[neip[i1]]==6){
      accel+=curoct->cell[vcell[neip[i1]]].pot*contrib;
      newcell=&(curoct->cell[vcell[neip[i1]]]);
      tot+=contrib;
      ntot++;
    }
    else{
      if(curoct->nei[vnei[neip[i1]]]->child!=NULL){
	accel+=curoct->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].pot*contrib;
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
	  accel+=newoct->cell[vcell2[neip[i2]]].pot*contrib;
	  newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	  visit[idx2]=1;
	  tot+=contrib;
	  ntot++;
	}
	else{
	  if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	    accel+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].pot*contrib;
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
	    accel+=newoct2->cell[vcell3[neip[i3]]].pot*contrib;
	    visit[7]=1;
	    tot+=contrib;
	    ntot++;

	  }
	  else{
	    if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
	      accel+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].pot*contrib;
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
    float tot=0;
    int ntot=0;
    contrib=vcont[0];
    if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoctlr->cell[icell].pot*contrib;
    }
  
  
    // contribs to cardinal neighbors
    getcellnei(icell, vnei, vcell);
    for(i1=0;i1<3;i1++){
      idx1=pow(2,i1);
      contrib=vcont[idx1];
      if(vnei[neip[i1]]==6){
	accel+=curoctlr->cell[vcell[neip[i1]]].pot*contrib;
	newcell=&(curoctlr->cell[vcell[neip[i1]]]);
	tot+=contrib;
	ntot++;
      }
      else{
	if(curoctlr->nei[vnei[neip[i1]]]->child!=NULL){
	  accel+=curoctlr->nei[vnei[neip[i1]]]->child->cell[vcell[neip[i1]]].pot*contrib;
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
	    accel+=newoct->cell[vcell2[neip[i2]]].pot*contrib;
	    newcell2=&(newoct->cell[vcell2[neip[i2]]]);
	    visit[idx2]=1;
	    tot+=contrib;
	    ntot++;
	  }
	  else{
	    if(newoct->nei[vnei2[neip[i2]]]->child!=NULL){
	      accel+=newoct->nei[vnei2[neip[i2]]]->child->cell[vcell2[neip[i2]]].pot*contrib;
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
	      accel+=newoct2->cell[vcell3[neip[i3]]].pot*contrib;
	      visit[7]=1;
	      tot+=contrib;
	      ntot++;

	    }
	    else{
	      if(newoct2->nei[vnei3[neip[i3]]]->child!=NULL){
		accel+=newoct2->nei[vnei3[neip[i3]]]->child->cell[vcell3[neip[i3]]].pot*contrib;
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
  
#endif


#ifndef NEWFORCE
  // contrib from 100 cell ===========================================================
      
  if(hres==1){
    contrib=tx*dy*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[0]]==6){
  	accel+=curoct->cell[vcell[neip[0]]].pot*contrib;
  	cicoct[1]=curoct;
      }
      else {
  	if(curoct->nei[vnei[neip[0]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].pot*contrib;
  	  cicoct[1]=curoct->nei[vnei[neip[0]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  // contrib from 010 cell ===========================================================
  if(hres==1){
    contrib=dx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[1]]==6){
  	accel+=curoct->cell[vcell[neip[1]]].pot*contrib;
  	cicoct[2]=curoct;
      }
      else{
  	if(curoct->nei[vnei[neip[1]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].pot*contrib;
  	  cicoct[2]=curoct->nei[vnei[neip[1]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  // contrib from 001 cell ===========================================================

  if(hres==1){
    contrib=dx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[2]]==6){
  	accel+=curoct->cell[vcell[neip[2]]].pot*contrib;
  	cicoct[4]=curoct;
      }
      else{
  	if(curoct->nei[vnei[neip[2]]]->child!=NULL){
  	  accel+=curoct->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].pot*contrib;
  	  cicoct[4]=curoct->nei[vnei[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
  //contrib from 110 cell 2 paths ======================================================
  
  contrib=tx*ty*dz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){
      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[3]=vcell2[neip[1]];
  	if(vnei2[neip[1]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[1]]].pot*contrib;
  	  cicoct[3]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].pot*contrib;
  	  cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
      }
    }
  }
      

  //contrib from 101 cell 2 paths ======================================================
  contrib=tx*dy*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[5]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[5]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}
	
      }
    }
  }

  //contrib from 011 cell 2 paths ======================================================
  contrib=dx*ty*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[2]!=NULL){
  	getcellnei(ciccell[2], vnei2, vcell2);
  	ciccell[6]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[2]->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[6]=cicoct[2];
  	}
  	else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
  	}

      }
    }
  }
      
  // contrib from 111 ======================================================
  contrib=tx*ty*tz;
  if(hres==1){
    if((contrib<=1.)&&(contrib>=0)){
		  
      if(cicoct[3]!=NULL){
  	getcellnei(ciccell[3], vnei2, vcell2);
  	ciccell[7]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[3]->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[7]=cicoct[3];
  	}
  	else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
  	}
  	else{ // the particle is not deep enough
  	  accel=0.;
  	  hres=0;
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
    
    // getting the neighbors
    getcellnei(icell, vnei, vcell);
    
    // here we compute the indexes of the direct neighbors which are involved in the cic
    neip[0]=(curp->x<xc?0:1);
    neip[1]=(curp->y<yc?2:3);
    neip[2]=(curp->z<zc?4:5);
    
    // here we denote the offset in ZYX
    ciccell[0]=icell;          //cell 000
    ciccell[1]=vcell[neip[0]]; //cell 001
    ciccell[2]=vcell[neip[1]]; //cell 010
    ciccell[4]=vcell[neip[2]]; //cell 100
    cicoct[0]=curoctlr;

    // the CIC weights
    tx=fabs((curp->x-xc)/dxcur);
    ty=fabs((curp->y-yc)/dxcur);
    tz=fabs((curp->z-zc)/dxcur);
    
    dx=1.-tx;
    dy=1.-ty;
    dz=1.-tz;
    
    // contrib from current cell 000 =====
    contrib=dx*dy*dz;
    if((contrib<=1.)&&(contrib>=0.)){
      accel+=curoctlr->cell[icell].pot*contrib;
      //curoct->cell[icell].density+=contrib/pow(dxcur,3);
    }
  

    // contrib from 100 cell ===========================================================
      
    contrib=tx*dy*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[0]]==6){
  	accel+=curoctlr->cell[vcell[neip[0]]].pot*contrib;
  	cicoct[1]=curoctlr;
      }
      else {
  	if(curoctlr->nei[vnei[neip[0]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].pot*contrib;
  	  cicoct[1]=curoctlr->nei[vnei[neip[0]]]->child;
  	}
      }
    }

    // contrib from 010 cell ===========================================================
    contrib=dx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[1]]==6){
  	accel+=curoctlr->cell[vcell[neip[1]]].pot*contrib;
  	cicoct[2]=curoctlr;
      }
      else{
  	if(curoctlr->nei[vnei[neip[1]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].pot*contrib;
  	  cicoct[2]=curoctlr->nei[vnei[neip[1]]]->child;
  	}
      }
    }
    
    // contrib from 001 cell ===========================================================

    contrib=dx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){
      if(vnei[neip[2]]==6){
  	accel+=curoctlr->cell[vcell[neip[2]]].pot*contrib;
  	cicoct[4]=curoctlr;
      }
      else{
  	if(curoctlr->nei[vnei[neip[2]]]->child!=NULL){
  	  accel+=curoctlr->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].pot*contrib;
  	  cicoct[4]=curoctlr->nei[vnei[neip[2]]]->child;
  	}
      }
    }

    //contrib from 110 cell 2 paths ======================================================
  
    contrib=tx*ty*dz;
    if((contrib<=1.)&&(contrib>=0)){
      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[3]=vcell2[neip[1]];
  	if(vnei2[neip[1]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[1]]].pot*contrib;
  	  cicoct[3]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].pot*contrib;
  	  cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
  	}
      }
    }
    
    //contrib from 101 cell 2 paths ======================================================
    contrib=tx*dy*tz;
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[1]!=NULL){
  	getcellnei(ciccell[1], vnei2, vcell2);
  	ciccell[5]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[1]->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[5]=cicoct[1];
  	}
  	else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }

    //contrib from 011 cell 2 paths ======================================================
    contrib=dx*ty*tz;
    if((contrib<=1.)&&(contrib>=0)){

      if(cicoct[2]!=NULL){
  	getcellnei(ciccell[2], vnei2, vcell2);
  	ciccell[6]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[2]->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[6]=cicoct[2];
  	}
  	else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }
      
    // contrib from 111 ======================================================
    contrib=tx*ty*tz;
    if((contrib<=1.)&&(contrib>=0)){
		  
      if(cicoct[3]!=NULL){
  	getcellnei(ciccell[3], vnei2, vcell2);
  	ciccell[7]=vcell2[neip[2]];
  	if(vnei2[neip[2]]==6){
  	  accel+=cicoct[3]->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[7]=cicoct[3];
  	}
  	else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
  	  accel+=cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].pot*contrib;
  	  cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
  	}
      }
    }
  }

#endif
  
#ifndef NEWFORCE
  free(cicoct);
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
  float dxcur;

  struct OCT *newoct;
  struct OCT *newoct2;
  
  struct CELL *newcell;
  struct CELL *newcell2;

  int vnei[6],vcell[6];
  int vnei2[6],vcell2[6];
  int il,ip;

  if(cpu->rank==0) printf("==> start CIC\n");

  // First we clean the density
  for(level=levelmax;level>=levelcoarse;level--)
    {
      nextoct=firstoct[level-1];
      if(nextoct==NULL) continue;
      do // sweeping level
	{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  for(icell=0;icell<8;icell++) curoct->cell[icell].density=0.;
	}while(nextoct!=NULL);
    }

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
		do{  
		  curp=nexp; 
		  nexp=curp->next; 
		  part2cell_cic(curp, curoct, icell,1); 
		}while(nexp!=NULL); 
	      }

#ifndef NGP

	      //==== SECOND WE CONSIDER THE PARTICLES INSIDE THE NEIGHBORS AT LEVEL L-1

	      int visit[27];
	      int iv;
	      for(inei=0;inei<27;inei++) visit[inei]=0;
	      int tot=0;
	      // first the cartesian neighbors (6)
	      for(inei=0;inei<6;inei++){
		nexp=curoct->nei[inei]->phead; //sweeping the particles of the neighbour cell at level l-1
		iv=13+pow(3,inei/2)*(2*(inei%2)-1);
		visit[iv]=1;
		if(nexp!=NULL){
		  do{ 
		    curp=nexp;
		    nexp=curp->next;
		    part2cell_cic(curp, curoct, icell,0);
		  }while(nexp!=NULL);
		}
#ifdef NEWCIC
		// second the fundamental plane (4)
		// each of the 6 cardinal neighbors will probe 4 neighbors
		
		newcell=curoct->nei[inei];
		newoct=cell2oct(newcell); // we get the parent oct;
		getcellnei(newcell->idx, vnei, vcell); // we get its neighbors
		for(il=0;il<6;il++){
		  iv=pow(3,inei/2)*(2*(inei%2)-1)+pow(3,il/2)*(2*(il%2)-1);
		  if((il/2)==(inei/2)) continue;
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
		       iv=pow(3,inei/2)*(2*(inei%2)-1)+pow(3,il/2)*(2*(il%2)-1)+pow(3,ip/2)*(2*(ip%2)-1);
		       if(((ip/2)==(il/2))||((ip/2)==(inei/2))) continue;
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
#endif
	      }

#ifndef NEWCIC
	      // second the fundamental plane (4)

	      char dir[4]={2,3,1,0};
	      for(inei=0;inei<4;inei++){
		if(curoct->nei[inei]->child!=NULL){
		  nexp=curoct->nei[inei]->child->nei[dir[inei]]->phead; //sweeping the particles of the neighbour cell at level l-1
		  if(nexp!=NULL){
		    do{ 
		      curp=nexp;
		      nexp=curp->next;
		      part2cell_cic(curp, curoct, icell,0);
		    }while(nexp!=NULL);
		  }
		}
	      }

	      // third the top-bottom cross (8)
	      for(inei=5;inei<6;inei++){
		if(curoct->nei[inei]->child!=NULL){
		  for(inei2=0;inei2<4;inei2++){
		    nexp=curoct->nei[inei]->child->nei[inei2]->phead; //sweeping the particles of the neighbour cell at level l-1
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

	      // fourth the plane of each top/bottom cross returns the corners (8)
	      
	      for(inei=5;inei<6;inei++){
		if(curoct->nei[inei]->child!=NULL){
		  for(inei2=0;inei2<4;inei2++){
		    if(curoct->nei[inei]->child->nei[inei2]->child!=NULL){
		      nexp=curoct->nei[inei]->child->nei[inei2]->child->nei[dir[inei2]]->phead; 
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

#endif
#endif
	      // THIRD WE LOOK FOR THE PARTICLES IN THE CHILD OCTS
	      if(curoct->cell[icell].child!=NULL){
		cic_child(curoct->cell[icell].child,curoct,icell);
	      }
	    }
	}while(nextoct!=NULL);
    }
  //  printf("great total=%f\n",toto);
}
