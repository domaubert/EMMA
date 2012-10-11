#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>
#include <mpi.h>

#ifdef GPUAXL
#include "poisson_utils_gpu.h"
#endif


#ifdef WGRAV
//================================================================
// ===============================================================
// ===============================================================

void interpminmodgrav(struct Gtype *U0, struct Gtype *Up, struct Gtype *Dx, struct Gtype *Dy, struct Gtype *Dz,REAL dx,REAL dy,REAL dz){
  
  Up->d =U0->d  +dx*Dx->d  +dy*Dy->d  +dz*Dz->d;
  Up->p =U0->p  +dx*Dx->p  +dy*Dy->p  +dz*Dz->p;
}

//===============================================
void minmod2grav(struct Gtype *Um, struct Gtype *Up, struct Gtype *Ur){
  REAL r;
  REAL xi;
  REAL w=0.;
  REAL beta=1.0;
  // SLOPE LIMITER

  if(Up->d==0.){
    xi=0.;}
  else{
    r=Um->d/Up->d;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xi=1.;
  Ur->d=(0.5*(1.+w)*Um->d+0.5*(1.-w)*Up->d)*xi;

  if(Up->p==0.){
    xi=0.;}
  else{
    r=Um->p/Up->p;
    if(r<=0.){
      xi=0.;
    }
    else if(r<=1.){
      xi=r;
    }
    else{
      xi=fmin(1.0,2.0*beta/(1.-w+(1.+w)*r));
    }
  }
  xi=1.;
  Ur->p=(0.5*(1.+w)*Um->p+0.5*(1.-w)*Up->p)*xi;
}

// ================== performs the difference between two Us

void diffUgrav(struct Gtype *U2, struct Gtype *U1, struct Gtype *UR){
  
  UR->d =U2->d - U1->d;
  UR->p=U2->p- U1->p;
}


// ===============================================================
// ==============================================

void coarse2fine_grav(struct CELL *cell, struct Gtype *Wi){ 

  
  struct OCT * oct;
	  
  struct Gtype U0;
  struct Gtype Up;
  struct Gtype Um;
  struct Gtype Dp,Dm;
  struct Gtype D[3];
  struct Gtype *W;
  int inei2;
  int vcell[6],vnei[6];
  int dir;
  REAL dxcur;
  
  oct=cell2oct(cell);
  getcellnei(cell->idx, vnei, vcell); // we get the neighbors
  dxcur=pow(0.5,oct->level);
  
  W=&(cell->gdata);
  U0=*W;
  // Limited Slopes
  for(dir=0;dir<3;dir++){
    
    inei2=2*dir;
    if(vnei[inei2]==6){
      W=&(oct->cell[vcell[inei2]].gdata);
    }
    else{
      W=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].gdata);
      
#ifdef TRANSXM
      if((oct->x==0.)&&(inei2==0)){
	W=&(cell->gdata);
      }
#endif
      
#ifdef TRANSYM
      if((oct->y==0.)&&(inei2==2)){
	W=&(cell->gdata);
      }
#endif
      
#ifdef TRANSZM
      if((oct->z==0.)&&(inei2==4)){
	W=&(cell->gdata);
      }
#endif
      
      
    }
    Um=*W;
    
    inei2=2*dir+1;
    if(vnei[inei2]==6){
      W=&(oct->cell[vcell[inei2]].gdata);
    }
    else{
      W=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].gdata);
      
#ifdef TRANSXP
      if(((oct->x+2.*dxcur)==1.)&&(inei2==1)){
	W=&(cell->gdata);
      }
#endif
      
#ifdef TRANSYP
      if(((oct->y+2.*dxcur)==1.)&&(inei2==3)){
	W=&(cell->gdata);
      }
#endif
      
#ifdef TRANSZP
      //if((oct->nei[vnei[inei2]]->child->z-oct->z)<0.){
      if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
	W=&(cell->gdata);
      }
#endif
      
    }
    Up=*W;
    
    diffUgrav(&Up,&U0,&Dp); 
    diffUgrav(&U0,&Um,&Dm); 
    
    
    minmod2grav(&Dm,&Dp,D+dir);
}
  
  // Interpolation
  int ix,iy,iz;
  int icell;
  
for(iz=0;iz<2;iz++){
    for(iy=0;iy<2;iy++){
      for(ix=0;ix<2;ix++){
	icell=ix+iy*2+iz*4;
	interpminmodgrav(&U0,&Up,D,D+1,D+2,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation
	memcpy(Wi+icell,&Up,sizeof(struct Gtype));
      }
    }
  }

}
//================================================================
void oldrecursive_neighbor_gather_oct_grav(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct HGRID *stencil,char *visit){

  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;
  static int face[8]={0,1,2,3,4,5,6,7};
  REAL dxcur;

  struct Gtype Wi[8];
  struct OCT *oct;
  struct OCT *neioct;
  struct CELL *neicell;

  if(order==1){
    ineiloc=inei;
  }
  else if(order==2){
    ineiloc=inei2;
  }
  else if(order==3){
    ineiloc=inei3;

  }

  if(cell->child!=NULL){
    // the oct at the right level exists
    neicell=cell->child->nei[ineiloc];
  }
  else{
    getcellnei(cell->idx, vnei, vcell); // we get the neighbors
    oct=cell2oct(cell);
    if(vnei[ineiloc]==6){
      neicell=&(oct->cell[vcell[ineiloc]]);
    }
    else{
      if(oct->nei[ineiloc]->child!=NULL){
	neicell=&(oct->nei[ineiloc]->child->cell[vcell[ineiloc]]);
      }
      else{
	printf("big problem\n");
	abort();
      }
    }
  }


  /* oct=cell2oct(cell); */
  /* neioct=cell2oct(neicell); */
  /* dxcur=pow(0.5,oct->level); */

  // ============================ TRANSMISSIVE BOUNDARIES ====================
#ifdef TRANSXP
    if(ineiloc==1){
      if((oct->x+2.*dxcur)==1.){
	neicell=cell;
	face[0]=1;
	face[1]=1;
	face[2]=3;
	face[3]=3;
	face[4]=5;
	face[5]=5;
	face[6]=7;
	face[7]=7;
      }
    }
#endif


#ifdef TRANSYP
    if(ineiloc==3){
      if((oct->y+2.*dxcur)==1.){
	neicell=cell;
	face[0]=2;
	face[1]=3;
	face[2]=2;
	face[3]=3;
	face[4]=7;
	face[5]=6;
	face[6]=6;
	face[7]=7;
      }
    }
#endif

#ifdef TRANSZP
    if(ineiloc==5){
      if((oct->z+2.*dxcur)==1.){
	neicell=cell;
	face[0]=4;
	face[1]=5;
	face[2]=6;
	face[3]=7;
	face[4]=4;
	face[5]=5;
	face[6]=6;
	face[7]=7;
      }
    }
#endif


      
#ifdef TRANSXM
    if(ineiloc==0){
      if(oct->x==0.){
	neicell=cell;
	face[0]=0;
	face[1]=0;
	face[2]=2;
	face[3]=2;
	face[4]=4;
	face[5]=4;
	face[6]=6;
	face[7]=6;
      }
    }
#endif

#ifdef TRANSYM
    if(ineiloc==2){
      if(oct->y==0.){
	neicell=cell;
	face[0]=0;
	face[1]=1;
	face[2]=0;
	face[3]=1;
	face[4]=4;
	face[5]=5;
	face[6]=4;
	face[7]=5;
      }
    }
#endif

#ifdef TRANSZM
    if(ineiloc==4){
      if(oct->z==0.){
	neicell=cell;
	face[0]=0;
	face[1]=1;
	face[2]=2;
	face[3]=3;
	face[4]=0;
	face[5]=1;
	face[6]=2;
	face[7]=3;
      }
    }
#endif


  // ============================ END TRANSMISSIVE BOUNDARIES ====================


  if(neicell->child!=NULL){
    // optimal case
    //for(icell=0;icell<8;icell++) memcpy(&(stencil->oct[ioct].cell[icell].gdata),&(neicell->child->cell[face[icell]].gdata),sizeof(struct Gtype));
    //    for(icell=0;icell<8;icell++) memcpy(&(stencil->oct[ioct].cell[icell].gdata.p),&(neicell->child->cell[face[icell]].gdata.p),sizeof(REAL));
    for(icell=0;icell<8;icell++) stencil->oct[ioct].cell[icell].gdata.p=neicell->child->cell[face[icell]].gdata.p;
  }
  else{
    coarse2fine_grav(neicell,Wi);
    for(icell=0;icell<8;icell++){
      memcpy(&(stencil->oct[ioct].cell[icell].gdata),Wi+face[icell],sizeof(struct Gtype)); //
    }
  }




  // next order
  if(order==1){
    for(i=0;i<6;i++){
      if((i>>1)==(inei>>1)) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      if(visit[ioct2]) continue;
      visit[ioct2]=1;
      oldrecursive_neighbor_gather_oct_grav(ioct2, inei, i, -1, 2, neicell, stencil,visit);
    }
  }
  else if(order==2) {
    for(i=0;i<6;i++){
      if(((i>>1)==(inei>>1))||((i>>1)==(inei2>>1))) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      if(visit[ioct2]) continue;
      visit[ioct2]=1;
      oldrecursive_neighbor_gather_oct_grav(ioct2, inei, inei2, i, 3, neicell, stencil,visit);
    }
  }
}

void recursive_neighbor_gather_oct_grav(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct GGRID *stencil,char *visit){

  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;
  static int face[8]={0,1,2,3,4,5,6,7};
  REAL dxcur;

  struct Gtype Wi[8];
  struct OCT *oct;
  struct OCT *neioct;
  struct CELL *neicell;

  if(order==1){
    ineiloc=inei;
  }
  else if(order==2){
    ineiloc=inei2;
  }
  else if(order==3){
    ineiloc=inei3;

  }

  if(cell->child!=NULL){
    // the oct at the right level exists
    neicell=cell->child->nei[ineiloc];
  }
  else{
    getcellnei(cell->idx, vnei, vcell); // we get the neighbors
    oct=cell2oct(cell);
    if(vnei[ineiloc]==6){
      neicell=&(oct->cell[vcell[ineiloc]]);
    }
    else{
      if(oct->nei[ineiloc]->child!=NULL){
	neicell=&(oct->nei[ineiloc]->child->cell[vcell[ineiloc]]);
      }
      else{
	printf("big problem\n");
	abort();
      }
    }
  }


  /* oct=cell2oct(cell); */
  /* neioct=cell2oct(neicell); */
  /* dxcur=pow(0.5,oct->level); */

  // ============================ TRANSMISSIVE BOUNDARIES ====================
#ifdef TRANSXP
    if(ineiloc==1){
      if((oct->x+2.*dxcur)==1.){
	neicell=cell;
	face[0]=1;
	face[1]=1;
	face[2]=3;
	face[3]=3;
	face[4]=5;
	face[5]=5;
	face[6]=7;
	face[7]=7;
      }
    }
#endif


#ifdef TRANSYP
    if(ineiloc==3){
      if((oct->y+2.*dxcur)==1.){
	neicell=cell;
	face[0]=2;
	face[1]=3;
	face[2]=2;
	face[3]=3;
	face[4]=7;
	face[5]=6;
	face[6]=6;
	face[7]=7;
      }
    }
#endif

#ifdef TRANSZP
    if(ineiloc==5){
      if((oct->z+2.*dxcur)==1.){
	neicell=cell;
	face[0]=4;
	face[1]=5;
	face[2]=6;
	face[3]=7;
	face[4]=4;
	face[5]=5;
	face[6]=6;
	face[7]=7;
      }
    }
#endif


      
#ifdef TRANSXM
    if(ineiloc==0){
      if(oct->x==0.){
	neicell=cell;
	face[0]=0;
	face[1]=0;
	face[2]=2;
	face[3]=2;
	face[4]=4;
	face[5]=4;
	face[6]=6;
	face[7]=6;
      }
    }
#endif

#ifdef TRANSYM
    if(ineiloc==2){
      if(oct->y==0.){
	neicell=cell;
	face[0]=0;
	face[1]=1;
	face[2]=0;
	face[3]=1;
	face[4]=4;
	face[5]=5;
	face[6]=4;
	face[7]=5;
      }
    }
#endif

#ifdef TRANSZM
    if(ineiloc==4){
      if(oct->z==0.){
	neicell=cell;
	face[0]=0;
	face[1]=1;
	face[2]=2;
	face[3]=3;
	face[4]=0;
	face[5]=1;
	face[6]=2;
	face[7]=3;
      }
    }
#endif


  // ============================ END TRANSMISSIVE BOUNDARIES ====================


  if(neicell->child!=NULL){
    // optimal case
    //for(icell=0;icell<8;icell++) memcpy(&(stencil->oct[ioct].cell[icell].gdata),&(neicell->child->cell[face[icell]].gdata),sizeof(struct Gtype));
    //    for(icell=0;icell<8;icell++) memcpy(&(stencil->oct[ioct].cell[icell].gdata.p),&(neicell->child->cell[face[icell]].gdata.p),sizeof(REAL));
    for(icell=0;icell<8;icell++) stencil->oct[ioct].cell[icell].gdata.p=neicell->child->cell[face[icell]].gdata.p;
  }
  else{
    coarse2fine_grav(neicell,Wi);
    for(icell=0;icell<8;icell++){
      memcpy(&(stencil->oct[ioct].cell[icell].gdata),Wi+face[icell],sizeof(struct Gtype)); //
    }
  }




  // next order
  if(order==1){
    for(i=0;i<6;i++){
      if((i>>1)==(inei>>1)) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      if(visit[ioct2]) continue;
      visit[ioct2]=1;
      recursive_neighbor_gather_oct_grav(ioct2, inei, i, -1, 2, neicell, stencil,visit);
    }
  }
  else if(order==2) {
    for(i=0;i<6;i++){
      if(((i>>1)==(inei>>1))||((i>>1)==(inei2>>1))) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      if(visit[ioct2]) continue;
      visit[ioct2]=1;
      recursive_neighbor_gather_oct_grav(ioct2, inei, inei2, i, 3, neicell, stencil,visit);
    }
  }
}

//================================================================
//================================================================

struct OCT *oldgatherstencilgrav(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  
  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  char visit[27]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
  int ioct;

  //printf("let's gather\n");
  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	stencil[iread].oct[13].cell[icell].gdata.d=curoct->cell[icell].gdata.d;
	stencil[iread].oct[13].cell[icell].gdata.p=curoct->cell[icell].gdata.p;
      }
      visit[13]=1;

      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  visit[ioct]=1;
	  oldrecursive_neighbor_gather_oct_grav(ioct, inei, -1, -1, 1, cell, stencil+iread,visit);
	}
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  
  (*nread)=iread;
  return nextoct;
}

struct OCT *gatherstencilgrav(struct OCT *octstart, struct GGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  
  static int ix[6]={-1,1,0,0,0,0};
  static int iy[6]={0,0,-1,1,0,0};
  static int iz[6]={0,0,0,0,-1,1};
  char visit[27]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
  int ioct;

  //printf("let's gather\n");
  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	stencil[iread].oct[13].cell[icell].gdata.d=curoct->cell[icell].gdata.d;
	stencil[iread].oct[13].cell[icell].gdata.p=curoct->cell[icell].gdata.p;
      }
      visit[13]=1;

      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  visit[ioct]=1;
	  recursive_neighbor_gather_oct_grav(ioct, inei, -1, -1, 1, cell, stencil+iread,visit);
	}
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  
  (*nread)=iread;
  return nextoct;
}


//==============================================================================

 struct OCT *scatterstencilgrav(struct OCT *octstart, struct GGRID *stencil, int stride, struct CPUINFO *cpu)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int ioct[7]={12,14,10,16,4,22,13};
  
  nextoct=octstart;
  iread=0;

  //printf("let's scatter\n");
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(&(curoct->cell[icell].pnew),&(stencil[iread].pnew[icell]),sizeof(REAL));
	memcpy(&(curoct->cell[icell].res),&(stencil[iread].res[icell]),sizeof(REAL));
      }

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}

//==============================================================================
struct OCT *oldscatterstencilgrav(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int ioct[7]={12,14,10,16,4,22,13};
  
  nextoct=octstart;
  iread=0;

  //printf("let's scatter\n");
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(&(curoct->cell[icell].pnew),&(stencil[iread].New.cell[icell].pnew),sizeof(REAL));
	memcpy(&(curoct->cell[icell].res),&(stencil[iread].New.cell[icell].res),sizeof(REAL));
      }

      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}

//==============================================================================
//==============================================================================

struct OCT *scatterstencilforce(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int ioct[7]={12,14,10,16,4,22,13};
  
  nextoct=octstart;
  iread=0;

  //printf("let's scatter\n");
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;
      
      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	memcpy(curoct->cell[icell].f,stencil[iread].New.cell[icell].f,sizeof(REAL)*3);
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}




//============================================================================
int PoissonJacobi_single(struct GGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, int flag, REAL factdens){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int inei,icell;
  int i;
  REAL temp;
  REAL res;
  int vnei[6],vcell[6];
  int ioct[7]={12,14,10,16,4,22,13};

  struct Gtype *curcell;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
    for(i=0;i<nread;i++){ // we scan the octs
      
      
      temp=0.;
      res=0.;
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].gdata);
      
      // Computing the laplacian ===========================
 
      for(inei=0;inei<6;inei++){
 	temp+=stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].gdata.p;
 	/* if(level==6){ */
	/*   if(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].gdata.p!=0.) abort(); */
	/* } */
      }
 
      // setting up the residual
      res=temp;
      
      // we finish the laplacian
      temp=temp/6.0;
      temp=temp-dx*dx*curcell->d/6.0*factdens;

      // we finsih the residual
      res=res-6.0*curcell->p;
      res=res/(dx*dx)-factdens*curcell->d;

      // we store the new value of the potential
      stencil[i].pnew[icell]=temp;

      // we store the local residual
      if(flag) {
	stencil[i].res[icell]=factdens*curcell->d;
      }
      else{
	stencil[i].res[icell]=res;
      }

      // ready for the next cell
    }
    //ready for the next oct
  }
  return 0;
}


//============================================================================
int oldPoissonJacobi_single(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, int flag, REAL factdens){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int inei,icell;
  int i;
  REAL temp;
  REAL res;
  int vnei[6],vcell[6];
  int ioct[7]={12,14,10,16,4,22,13};

  struct Gtype *curcell;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
    for(i=0;i<nread;i++){ // we scan the octs
      
      
      temp=0.;
      res=0.;
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].gdata);
      
      // Computing the laplacian ===========================
 
      for(inei=0;inei<6;inei++){
 	temp+=stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].gdata.p;
 	/* if(level==6){ */
	/*   if(stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].gdata.p!=0.) abort(); */
	/* } */
      }
 
      // setting up the residual
      res=temp;
      
      // we finish the laplacian
      temp=temp/6.0;
      temp=temp-dx*dx*curcell->d/6.0*factdens;

      // we finsih the residual
      res=res-6.0*curcell->p;
      res=res/(dx*dx)-factdens*curcell->d;

      // we store the new value of the potential
      stencil[i].New.cell[icell].pnew=temp;

      // we store the local residual
      if(flag) {
	stencil[i].New.cell[icell].res=factdens*curcell->d;
      }
      else{
	stencil[i].New.cell[icell].res=res;
      }

      // ready for the next cell
    }
    //ready for the next oct
  }
  return 0;
}



 //============================================================================
REAL comp_residual(struct GGRID *stencil, int level, int curcpu, int nread,int stride, int flag){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int icell;
  int i;

  REAL residual=0.;
  REAL rloc;

  if(flag){
    for(i=0;i<nread;i++){ // we scan the octs
      for(icell=0;icell<8;icell++){ // we scan the cells
	rloc=pow(stencil[i].res[icell],2);
	residual+=rloc;
	// ready for the next cell
      }
      //ready for the next oct
    }
  }
  else{
    for(i=0;i<nread;i++){ // we scan the octs
      for(icell=0;icell<8;icell++){ // we scan the cells
	rloc=pow(stencil[i].res[icell],2);
	residual=(residual>rloc?residual:rloc);
	// ready for the next cell
      }
      //ready for the next oct
    }
  }

  return residual;
}
 
//-------------
REAL oldcomp_residual(struct HGRID *stencil, int level, int curcpu, int nread,int stride, int flag){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int icell;
  int i;

  REAL residual=0.;
  REAL rloc;

  if(flag){
    for(i=0;i<nread;i++){ // we scan the octs
      for(icell=0;icell<8;icell++){ // we scan the cells
	rloc=pow(stencil[i].New.cell[icell].res,2);
	residual+=rloc;
	// ready for the next cell
      }
      //ready for the next oct
    }
  }
  else{
    for(i=0;i<nread;i++){ // we scan the octs
      for(icell=0;icell<8;icell++){ // we scan the cells
	rloc=pow(stencil[i].New.cell[icell].res,2);
	residual=(residual>rloc?residual:rloc);
	// ready for the next cell
      }
      //ready for the next oct
    }
  }

  return residual;
}
 


//============================================================================
int Pot2Force(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL tsim){

  int inei,icell;
  int i;
  REAL temp;
  REAL res;
  int vnei[6],vcell[6];
  int ioct[7]={12,14,10,16,4,22,13};
  REAL floc[3];
  struct Gtype *curcell;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
    for(i=0;i<nread;i++){ // we scan the octs
      
      
      temp=0.;
      res=0.;
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].gdata);

      // computing the gradient

      floc[0]=0.5*(stencil[i].oct[ioct[vnei[1]]].cell[vcell[1]].gdata.p-stencil[i].oct[ioct[vnei[0]]].cell[vcell[0]].gdata.p)/dx*tsim;
      floc[1]=0.5*(stencil[i].oct[ioct[vnei[3]]].cell[vcell[3]].gdata.p-stencil[i].oct[ioct[vnei[2]]].cell[vcell[2]].gdata.p)/dx*tsim;
      floc[2]=0.5*(stencil[i].oct[ioct[vnei[5]]].cell[vcell[5]].gdata.p-stencil[i].oct[ioct[vnei[4]]].cell[vcell[4]].gdata.p)/dx*tsim;

      
      // store the force
      
      memcpy(stencil[i].New.cell[icell].f,floc,3*sizeof(REAL));

      // ready for the next cell
    }
    //ready for the next oct
  }
  return 0;
}

//=============================================================================

 REAL PoissonJacobi(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct GGRID *stencil, int stride, REAL tsim)
{
  REAL dxcur;
  int iter;
  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot;
  int nread;
  REAL fnorm,residual,residualold,dres;
  int icell;
  int nitmax;
  REAL factdens;
  REAL rloc;


  // Computing the factor of the density
  if(level>=param->lcoarse){
#ifndef TESTCOSMO
    factdens=4.0*M_PI;
#else
    //    factdens=6.0*tsim; WARNING JUST TESTING WITHOUT TSIM!!!
    factdens=6.0;
#endif
  }
  else{
    factdens=1.;
  }

  // Computing the max number for iteration

  if((level==param->mgridlmin)||(level>param->lcoarse)){
    nitmax=param->niter;
  }
  else{
    nitmax=param->nrelax;
  }

  dxcur=pow(0.5,level);
  
  double tall=0.,tcal=0.,tscat=0.,tgat=0.;
  double tglob=0.,tup=0.;

  double tstart,tstop,tt;
  
  for(iter=0;iter<nitmax;iter++){
    tstart=MPI_Wtime();

    // --------------- setting the first oct of the level
    nextoct=firstoct[level-1];
    nreadtot=0;

    double temps[10];

    // --------------- some inits for iterative solver
    if(iter==0){
      fnorm=0.;
      residual=0.;
    }
    else{
      residual=0.;
    }

    // Scanning the octs

    if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
      do {
	curoct=nextoct;
	nextoct=curoct->next; 
	// -------------  cleaning working arrays
	//memset(stencil,0,stride*sizeof(struct HGRID));
	temps[0]=MPI_Wtime();
	
	// ------------ gathering the stencil value values
	nextoct=gatherstencilgrav(curoct,stencil,stride,cpu, &nread);

	
	// ------------ solving the hydro
	temps[3]=MPI_Wtime();
	PoissonJacobi_single(stencil,level,cpu->rank,nread,stride,dxcur,(iter==0),factdens);
	temps[7]=MPI_Wtime();

	// ------------ computing the residuals

	rloc=comp_residual(stencil,level,cpu->rank,nread,stride,(iter==0));
	if(iter==0){
	  fnorm+=rloc;
	}
	else{
	  residual=(residual>rloc?residual:rloc);
	}
	
	// ------------ scatter back the data
	
	nextoct=scatterstencilgrav(curoct,stencil, nread, cpu);
	
	//nextoct=scatterstencilforce(curoct,stencil, nread, cpu);
	
	temps[9]=MPI_Wtime();

	tall+=temps[9]-temps[0];
	tcal+=temps[7]-temps[3];
	tscat+=temps[9]-temps[7];
	tgat+=temps[3]-temps[0];
	nreadtot+=nread;
      }while(nextoct!=NULL);
    }

    tt=MPI_Wtime();
    // at this stage an iteration has been completed : let's update the potential and compute the residual
    if(nreadtot>0){
      curoct=firstoct[level-1];
      if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
	nextoct=curoct;
	do{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
	  for(icell=0;icell<8;icell++){	
	    curoct->cell[icell].gdata.p=curoct->cell[icell].pnew;
	  }
	}while(nextoct!=NULL);
      }
    }
    tstop=MPI_Wtime();
    tup+=(tstop-tt);
    tglob+=(tstop-tstart);
    
    if(iter>0){
      dres=sqrt(residual);
      if((dres)<param->poissonacc){
	if(level>=param->lcoarse) break;
      }
    }
  }
  printf("CPU | level=%d iter=%d res=%e tgat=%e tcal=%e tscat=%e tall=%e tup=%e tglob=%e\n",level,iter,dres,tgat/iter,tcal/iter,tscat/iter,tall/iter,tup/iter, tglob/iter);
  return dres;
}

REAL OldPoissonJacobi(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL tsim)
{
  REAL dxcur;
  int iter;
  struct OCT *nextoct;
  struct OCT *curoct;
  int nreadtot;
  int nread;
  REAL fnorm,residual,residualold,dres;
  int icell;
  int nitmax;
  REAL factdens;
  REAL rloc;


  // Computing the factor of the density
  if(level>=param->lcoarse){
#ifndef TESTCOSMO
    factdens=4.0*M_PI;
#else
    //    factdens=6.0*tsim; WARNING JUST TESTING WITHOUT TSIM!!!
    factdens=6.0;
#endif
  }
  else{
    factdens=1.;
  }

  // Computing the max number for iteration

  if((level==param->mgridlmin)||(level>param->lcoarse)){
    nitmax=param->niter;
  }
  else{
    nitmax=param->nrelax;
  }

  dxcur=pow(0.5,level);
  
  double tall=0.,tcal=0.,tscat=0.,tgat=0.;
  double tglob=0.,tup=0.;

  double tstart,tstop,tt;
  
  for(iter=0;iter<nitmax;iter++){
    tstart=MPI_Wtime();

    // --------------- setting the first oct of the level
    nextoct=firstoct[level-1];
    nreadtot=0;

    double temps[10];

    // --------------- some inits for iterative solver
    if(iter==0){
      fnorm=0.;
      residual=0.;
    }
    else{
      residual=0.;
    }

    // Scanning the octs

    if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
      do {
	curoct=nextoct;
	nextoct=curoct->next; 
	// -------------  cleaning working arrays
	//memset(stencil,0,stride*sizeof(struct HGRID));
	temps[0]=MPI_Wtime();
	
	// ------------ gathering the stencil value values
	nextoct=oldgatherstencilgrav(curoct,stencil,stride,cpu, &nread);

	
	// ------------ solving the hydro
	temps[3]=MPI_Wtime();
	oldPoissonJacobi_single(stencil,level,cpu->rank,nread,stride,dxcur,(iter==0),factdens);
	temps[7]=MPI_Wtime();

	// ------------ computing the residuals

	rloc=oldcomp_residual(stencil,level,cpu->rank,nread,stride,(iter==0));
	if(iter==0){
	  fnorm+=rloc;
	}
	else{
	  printf("rloc=%e\n",sqrt(rloc));
	  residual=(residual>rloc?residual:rloc);
	}
	
	// ------------ scatter back the data
	
	nextoct=oldscatterstencilgrav(curoct,stencil, nread, cpu);
	
	//nextoct=scatterstencilforce(curoct,stencil, nread, cpu);
	
	temps[9]=MPI_Wtime();

	tall+=temps[9]-temps[0];
	tcal+=temps[7]-temps[3];
	tscat+=temps[9]-temps[7];
	tgat+=temps[3]-temps[0];
	nreadtot+=nread;
      }while(nextoct!=NULL);
    }

    tt=MPI_Wtime();
    // at this stage an iteration has been completed : let's update the potential and compute the residual
    if(nreadtot>0){
      curoct=firstoct[level-1];
      if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
	nextoct=curoct;
	do{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
	  for(icell=0;icell<8;icell++){	
	    curoct->cell[icell].gdata.p=curoct->cell[icell].pnew;
	  }
	}while(nextoct!=NULL);
      }
    }
    tstop=MPI_Wtime();
    tup+=(tstop-tt);
    tglob+=(tstop-tstart);
    
    if(iter>0){
      dres=sqrt(residual);
      if((dres)<param->poissonacc){
	if(level>=param->lcoarse) break;
      }
    }
  }
  printf("CPU | level=%d iter=%d res=%e tgat=%e tcal=%e tscat=%e tall=%e tup=%e tglob=%e\n",level,iter,dres,tgat/iter,tcal/iter,tscat/iter,tall/iter,tup/iter, tglob/iter);
  return dres;
}




//===============================================================================================

REAL PoissonMgrid(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct GGRID *stencil, int stride, REAL tsim)
{
  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  REAL dres;
  struct Gtype Wi[8];
  struct CELL* curcell;
  // pre-relaxation

#ifndef GPUAXL
  dres=PoissonJacobi(level,param,firstoct,cpu,stencil,stride,tsim);
#else
  dres=PoissonJacobiGPU(level,param,firstoct,cpu,stencil,stride,tsim);
#endif


  
  //  if(dres>param->poissonacc){
  // reduction
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
    do{ 
      curoct=nextoct;
      nextoct=curoct->next;
      curoct->parent->gdata.d=0.;
      curoct->parent->gdata.p=0.;
      for(icell=0;icell<8;icell++){
	curoct->parent->gdata.d+=curoct->cell[icell].res*0.125; // we average the residual and store it as the new density
      }
    }while(nextoct!=NULL);
  }

  // full relaxation at coarsest level or recursive call to mgrid
  
  if((level-1)==param->mgridlmin){
#ifndef GPUAXL
    PoissonJacobi(level-1,param,firstoct,cpu,stencil,stride,tsim);
#else
    PoissonJacobiGPU(level-1,param,firstoct,cpu,stencil,stride,tsim);
#endif
  }
  else{
    PoissonMgrid(level-1,param,firstoct,cpu,stencil,stride,tsim);
  }
  
  // prolongation + correction
  
  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
    do // sweeping level
	 {
	   curoct=nextoct;
	   nextoct=curoct->next;
	   curcell=curoct->parent;
	   coarse2fine_grav(curcell,Wi);
	   //for(icell=0;icell<8;icell++) memcpy(Wi+icell,&(curcell->gdata),sizeof(struct Gtype));
	   for(icell=0;icell<8;icell++) // looping over cells in oct
	     {
	       curoct->cell[icell].gdata.p-=Wi[icell].p; // we propagate the error and correct the evaluation
	     }
	 }while(nextoct!=NULL);
  }
  
  // post relaxation
#ifndef GPUAXL
  dres=PoissonJacobi(level,param,firstoct,cpu,stencil,stride,tsim);
#else
  dres=PoissonJacobiGPU(level,param,firstoct,cpu,stencil,stride,tsim);
#endif
    return dres;
}

//===================================================================================================================================

int PoissonForce(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL tsim){
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icomp,icell;
  struct PART *curp;
  struct PART *nexp;
  int nread,nreadtot;

  dxcur=pow(0.5,level);
  
  nextoct=firstoct[level-1];
  nreadtot=0;
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next; 
      // -------------  cleaning working arrays
      
      memset(stencil,0,stride*sizeof(struct HGRID));
      
      // ------------ gathering the stencil value values
      nextoct=oldgatherstencilgrav(curoct,stencil,stride,cpu, &nread);
	
	// ------------ Computing the potential gradient
      
      Pot2Force(stencil,level,cpu->rank,nread,stride,dxcur,tsim);
      
      // ------------ scatter back the data
      
      nextoct=scatterstencilforce(curoct,stencil, nread, cpu);
      nreadtot+=nread;
    }while(nextoct!=NULL);
  }

#ifdef WMPI
    mpi_exchange(cpu,sendbuffer,recvbuffer,5,1); // fx field
    mpi_exchange(cpu,sendbuffer,recvbuffer,6,1); // fy field
    mpi_exchange(cpu,sendbuffer,recvbuffer,7,1); // fz field
#endif
}

#endif

#ifdef WGRAV

//===================================================================================================================================
//===================================================================================================================================
int FillDens(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu){
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dx;
  int icomp,icell;
  struct PART *curp;
  struct PART *nexp;
  int nread;
  REAL locdens;
  REAL avgdens=0.;
  int nc=0;


  curoct=firstoct[level-1];
  if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
    nextoct=curoct;
    do{
      curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
      for(icell=0;icell<8;icell++){	
	locdens=0.;
#ifdef PIC
	locdens+=curoct->cell[icell].density;
#endif
#ifdef WHYDRO2
	locdens+=curoct->cell[icell].field.d;
#endif
	curoct->cell[icell].gdata.d=locdens;

#ifdef TESTCOSMO
	curoct->cell[icell].gdata.d-=1.;
#endif

	avgdens+=locdens;
	nc++;
      }
    }while(nextoct!=NULL);
  }


  return 0;
}

#endif

//===================================================================================================================================
//===================================================================================================================================
#ifdef WGRAV
int PoissonSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct GGRID *stencil, int stride, REAL aexp){

  REAL res;
  int igrid;
  struct Gtype Wi[8];
  struct CELL* curcell;
  int icell2;
  struct OCT* curoct;
  struct OCT* nextoct;
  int icell;
  if(cpu->rank==0) printf("Start Poisson Solver ");

#ifndef GPUAXL
  if(cpu->rank==0)  printf("on CPU\n");
#else
  if(cpu->rank==0)  printf("on GPU\n");
#endif
  //breakmpi();

  if((level==param->lcoarse)&&(param->lcoarse!=param->mgridlmin)){
    for(igrid=0;igrid<param->nvcycles;igrid++){ // V-Cycles
      printf("----------------------------------------\n");
      res=PoissonMgrid(level,param,firstoct,cpu,stencil,stride,aexp);
      if(res<param->poissonacc) break;
    }
  }
  else{
#ifndef GPUAXL
    PoissonJacobi(level,param,firstoct,cpu,stencil,stride,aexp);
#else
    PoissonJacobiGPU(level,param,firstoct,cpu,stencil,stride,aexp);
#endif

  }
	
  //once done we propagate the solution to level+1

  nextoct=firstoct[level-1];
  if(nextoct!=NULL){
    do // sweeping level
      {
	curoct=nextoct;
	nextoct=curoct->next;
	for(icell=0;icell<8;icell++) // looping over cells in oct
	  {
	    curcell=&(curoct->cell[icell]);
	    if(curcell->child!=NULL){

	      coarse2fine_grav(curcell,Wi);

	      for(icell2=0;icell2<8;icell2++){
		memcpy(&(curcell->child->cell[icell2].gdata.p),&(Wi[icell2].p),sizeof(REAL));
	      }

	    }
	  }
      }while(nextoct!=NULL);
  }
  return 0;
}
#endif
