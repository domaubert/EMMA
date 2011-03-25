#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"


//==================================================================

void part2cell_cic(struct PART *curp, struct OCT *curoct, int icell, char full)
{
  float xc,yc,zc;
  float dxcur=1./pow(2,curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int neip[3];
  float tx,ty,tz;
  float dx,dy,dz;
  float contrib;
  struct OCT **cicoct;
  char ciccell[8];
  char fullok=0;

  cicoct=(struct OCT **)calloc(8,sizeof(struct OCT *));

  xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
  yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
  zc=curoct->z+( icell/4   )*dxcur+dxcur/2; 

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
  cicoct[0]=curoct;

  // the CIC weights
  tx=fabs((curp->x-xc)/dxcur);
  ty=fabs((curp->y-yc)/dxcur);
  tz=fabs((curp->z-zc)/dxcur);
  if(!full)
    {
      fullok=(tx<=1.)*(ty<=1.)*(tz<=1.);
    } // necessary for boundary cells
	  
  dx=1.-tx;
  dy=1.-ty;
  dz=1.-tz;

  //  printf("id=%d %f %f %f %f %f %f \n",curp->idx,dx,dy,dz,tx,ty,tz);

  // contrib to current cell 000 =====
  contrib=dx*dy*dz;
  //printf("000 id=%d contrib=%f\n",curp->idx,contrib);

  if((contrib<=1.)&&(contrib>=0.)){
    if(full||fullok){
      curoct->cell[icell].density+=contrib/pow(dxcur,3)*curp->mass;
      //if(!full) printf("%f %f %f contrib=%f\n",tx,ty,tz,contrib);
    }
  }
  

  if(full) // full=0 for level-1 corrections
    {
      // contrib to 100 cell ===========================================================
  
      contrib=tx*dy*dz;
      //printf("100 id=%d contrib=%f\n",curp->idx,contrib);
      if((contrib<=1.)&&(contrib>=0)){
	if(vnei[neip[0]]==6){
	  curoct->cell[vcell[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	  cicoct[1]=curoct;
	}
	else if(curoct->nei[vnei[neip[0]]]->child!=NULL){
	  curoct->nei[vnei[neip[0]]]->child->cell[vcell[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	  cicoct[1]=curoct->nei[vnei[neip[0]]]->child;
	}
      }
      // contrib to 010 cell ===========================================================
      contrib=dx*ty*dz;
      //printf("010 id=%d contrib=%f\n",curp->idx,contrib);
      if((contrib<=1.)&&(contrib>=0)){
	if(vnei[neip[1]]==6){
	  curoct->cell[vcell[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	  cicoct[2]=curoct;
	}
	else{
	  if(curoct->nei[vnei[neip[1]]]->child!=NULL){
	    curoct->nei[vnei[neip[1]]]->child->cell[vcell[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[2]=curoct->nei[vnei[neip[1]]]->child;
	  }
	}
      }
  
      // contrib to 001 cell ===========================================================
      contrib=dx*dy*tz;
      //printf("001 id=%d contrib=%f\n",curp->idx,contrib);
      if((contrib<=1.)&&(contrib>=0)){
	if(vnei[neip[2]]==6){
	  curoct->cell[vcell[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	  cicoct[4]=curoct;
	}
	else{
	  if(curoct->nei[vnei[neip[2]]]->child!=NULL){
	    curoct->nei[vnei[neip[2]]]->child->cell[vcell[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[4]=curoct->nei[vnei[neip[2]]]->child;
	  }
	}
      }
      //contrib to 110 cell 2 paths ======================================================
  
      contrib=tx*ty*dz;
      //printf("110 id=%d contrib=%f\n",curp->idx,contrib);
      if((contrib<=1.)&&(contrib>=0)){
	if(cicoct[1]!=NULL){
	  getcellnei(ciccell[1], vnei2, vcell2);
	  ciccell[3]=vcell2[neip[1]];
	  if(vnei2[neip[1]]==6){
	    cicoct[1]->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[3]=cicoct[1];
	  }
	  else if(cicoct[1]->nei[vnei2[neip[1]]]->child!=NULL){
	    cicoct[1]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[3]=cicoct[1]->nei[vnei2[neip[1]]]->child;
	  }
	}
  
	if(cicoct[3]==NULL){
	  if(cicoct[2]!=NULL){
	    getcellnei(ciccell[2], vnei2, vcell2);
	    ciccell[3]=vcell2[neip[0]];
	    if(vnei2[neip[0]]==6){
	      cicoct[2]->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[3]=cicoct[2];
	    }
	    else if(cicoct[2]->nei[vnei2[neip[0]]]->child!=NULL){
	      cicoct[2]->nei[vnei2[neip[0]]]->child->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[3]=cicoct[2]->nei[vnei2[neip[0]]]->child;
	    }
	  }
	}
      }
  
      //contrib to 101 cell 2 paths ======================================================
      contrib=tx*dy*tz;
      //printf("101 id=%d contrib=%f\n",curp->idx,contrib);
      if((contrib<=1.)&&(contrib>=0)){

	if(cicoct[1]!=NULL){
	  getcellnei(ciccell[1], vnei2, vcell2);
	  ciccell[5]=vcell2[neip[2]];
	  if(vnei2[neip[2]]==6){
	    cicoct[1]->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[5]=cicoct[1];
	  }
	  else if(cicoct[1]->nei[vnei2[neip[2]]]->child!=NULL){
	    cicoct[1]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[5]=cicoct[1]->nei[vnei2[neip[2]]]->child;
	  }
	}
  
	if(cicoct[5]==NULL){
	  if(cicoct[4]!=NULL){
	    getcellnei(ciccell[4], vnei2, vcell2);
	    ciccell[5]=vcell2[neip[0]];
	    if(vnei2[neip[0]]==6){
	      cicoct[4]->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[5]=cicoct[4];
	    }
	    else if(cicoct[4]->nei[vnei2[neip[0]]]->child!=NULL){
	      cicoct[4]->nei[vnei2[neip[0]]]->child->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[5]=cicoct[4]->nei[vnei2[neip[0]]]->child;
	    }
	  }
	}
      }
		  
      //contrib to 011 cell 2 paths ======================================================
      contrib=dx*ty*tz;
      //printf("011 id=%d contrib=%f\n",curp->idx,contrib);
      if((contrib<=1.)&&(contrib>=0)){

	if(cicoct[2]!=NULL){
	  getcellnei(ciccell[2], vnei2, vcell2);
	  ciccell[6]=vcell2[neip[2]];
	  if(vnei2[neip[2]]==6){
	    cicoct[2]->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[6]=cicoct[2];
	  }
	  else if(cicoct[2]->nei[vnei2[neip[2]]]->child!=NULL){
	    cicoct[2]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[6]=cicoct[2]->nei[vnei2[neip[2]]]->child;
	  }
	}
		  
	if(cicoct[6]==NULL){
	  if(cicoct[4]!=NULL){
	    getcellnei(ciccell[4], vnei2, vcell2);
	    ciccell[6]=vcell2[neip[1]];
	    if(vnei2[neip[1]]==6){
	      cicoct[4]->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[6]=cicoct[4];
	    }
	    else if(cicoct[4]->nei[vnei2[neip[1]]]->child!=NULL){
	      cicoct[4]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[6]=cicoct[4]->nei[vnei2[neip[1]]]->child;
	    }
	  }
	}
      }		  
      // contrib to 111 ======================================================
      contrib=tx*ty*tz;
      //printf("111 id=%d contrib=%f\n",curp->idx,contrib);
      if((contrib<=1.)&&(contrib>=0)){
		  
	if(cicoct[3]!=NULL){
	  getcellnei(ciccell[3], vnei2, vcell2);
	  ciccell[7]=vcell2[neip[2]];
	  if(vnei2[neip[2]]==6){
	    cicoct[3]->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[7]=cicoct[3];
	  }
	  else if(cicoct[3]->nei[vnei2[neip[2]]]->child!=NULL){
	    cicoct[3]->nei[vnei2[neip[2]]]->child->cell[vcell2[neip[2]]].density+=contrib/pow(dxcur,3)*curp->mass;
	    cicoct[7]=cicoct[3]->nei[vnei2[neip[2]]]->child;
	  }
	}
		  
	if(cicoct[7]==NULL){
	  if(cicoct[5]!=NULL){
	    getcellnei(ciccell[5], vnei2, vcell2);
	    ciccell[7]=vcell2[neip[1]];
	    if(vnei2[neip[1]]==6){
	      cicoct[5]->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[7]=cicoct[5];
	    }
	    else if(cicoct[5]->nei[vnei2[neip[1]]]->child!=NULL){
	      cicoct[5]->nei[vnei2[neip[1]]]->child->cell[vcell2[neip[1]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[7]=cicoct[5]->nei[vnei2[neip[1]]]->child;
	    }
	  }
	}
		  
	if(cicoct[7]==NULL){
	  if(cicoct[6]!=NULL){
	    getcellnei(ciccell[6], vnei2, vcell2);
	    ciccell[7]=vcell2[neip[0]];
	    if(vnei2[neip[0]]==6){
	      cicoct[6]->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[7]=cicoct[6];
	    }
	    else if(cicoct[6]->nei[vnei2[neip[0]]]->child!=NULL){
	      cicoct[6]->nei[vnei2[neip[0]]]->child->cell[vcell2[neip[0]]].density+=contrib/pow(dxcur,3)*curp->mass;
	      cicoct[7]=cicoct[6]->nei[vnei2[neip[0]]]->child;
	    }
	  }
	}
      }
    }

#if 0
  //check total
  float tot=0.;
  for(ii=0;ii<8;ii++)
    {
      float xcc,ycc,zcc;
      tot+=cicoct[ii]->cell[ciccell[ii]].density;
      xcc=cicoct[ii]->x+( ciccell[ii]   %2)*dxcur+dxcur/2; // coordinates of the cell center 
      ycc=cicoct[ii]->y+((ciccell[ii]/2)%2)*dxcur+dxcur/2;
      zcc=cicoct[ii]->z+( ciccell[ii]/4   )*dxcur+dxcur/2; 
      
      printf("%f %f %f %d %f\n",xcc,ycc,zcc,ciccell[ii],cicoct[ii]->cell[ciccell[ii]].density*pow(dxcur,3));
    }
  printf("tot=%f\n",tot*pow(dxcur,3));
#endif

  free(cicoct);
}



void cell2part_cic(struct PART *curp, struct OCT *curoct, int icell, char dir, float dt)
{
  float xc,yc,zc;
  float dxcur=1./pow(2,curoct->level);
  int vnei [6],vcell [6];
  int vnei2[6],vcell2[6];
  int neip[3];
  float tx,ty,tz;
  float dx,dy,dz;
  float contrib;
  struct OCT **cicoct;
  struct OCT *curoctlr;
  char ciccell[8];
  char fullok=0;
  float accel=0.;
  char hres=1; // could be switched to hres=0 if particle is not deep enough

  cicoct=(struct OCT **)calloc(8,sizeof(struct OCT *));

  xc=curoct->x+( icell   %2)*dxcur+dxcur/2; // coordinates of the cell center 
  yc=curoct->y+((icell/2)%2)*dxcur+dxcur/2;
  zc=curoct->z+( icell/4   )*dxcur+dxcur/2; 

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
  cicoct[0]=curoct;

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
      accel+=curoct->cell[icell].temp*contrib;
      //curoct->cell[icell].density+=contrib/pow(dxcur,3);
  }
  

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

  
  // Once we have the acceleration we can compute the velocity
  switch(dir){
  case(0):
    curp->vx+=-accel*dt;
    break;
  case(1):
    curp->vy+=-accel*dt;
    break;
  case(2):
    curp->vz+=-accel*dt;
    break;
  }

  free(cicoct);
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
  for(level=levelmax;level>=levelcoarse;level--)
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

#if 1
	      //==== SECOND WE CONSIDER THE PARTICLES INSIDE THE NEIGHBORS AT LEVEL L-1
	      // first the cartesian neighbors (6)
	      for(inei=0;inei<6;inei++){
		nexp=curoct->nei[inei]->phead; //sweeping the particles of the neighbour cell at level l-1
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
		  if((il/2)==(inei/2)) continue;
		  if(vnei[il]==6){
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
		       if(((ip/2)==(il/2))||((ip/2)==(inei/2))) continue;
		       if(vnei2[ip]==6){
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
#endif
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
	      // THIRD WE LOOK FOR THE PARTICLES IN THE CHILD OCTS
	      if(curoct->cell[icell].child!=NULL){
		cic_child(curoct->cell[icell].child,curoct,icell);
	      }
	    }

	}while(nextoct!=NULL);
    }
  //  printf("great total=%f\n",toto);
}
