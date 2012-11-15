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
 
void interpforce(REAL *f0, REAL *fp, REAL *Dx, REAL *Dy, REAL *Dz,REAL dx,REAL dy,REAL dz){
  
  int i;
  for(i=0;i<3;i++){
    fp[i]=f0[i]+Dx[i]*dx+Dy[i]*dy+Dz[i]*dz;
  }
 
}
 

//===============================================
void minmod2grav(struct Gtype *Um, struct Gtype *Up, struct Gtype *Ur){
  REAL r;
  REAL xi;
  REAL w=0.;

  Ur->d=(0.5*(1.+w)*Um->d+0.5*(1.-w)*Up->d);
  Ur->p=(0.5*(1.+w)*Um->p+0.5*(1.-w)*Up->p);
}


void minmod2grav_mix(struct Gtype *U1, struct Gtype *U2){
  REAL Dm=U1->d;
  REAL Dp=U2->d;
  REAL beta=1.; // 1 MINBEE 2 SUPERBEE

  if(Dp>0){
    U1->d=fmax(fmax(0.,fmin(beta*Dm,Dp)),fmin(Dm,beta*Dp));
    U2->d=U1->d;
  }
  else{
    U1->d=fmin(fmin(0.,fmax(beta*Dm,Dp)),fmax(Dm,beta*Dp));
    U2->d=U1->d;
  }



}

// ================== performs the difference between two Us

void diffUgrav(struct Gtype *U2, struct Gtype *U1, struct Gtype *UR){
  
  UR->d=U2->d- U1->d;
  UR->p=U2->p- U1->p;
}

void diffforce(REAL *f2, REAL *f1, REAL *fr){
  
  int i;
  for(i=0;i<3;i++){
    fr[i]=f2[i]-f1[i];
  }
 
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
    //memcpy(D+(2*dir+0),&Dm,sizeof(struct Gtype));
    //memcpy(D+(2*dir+1),&Dp,sizeof(struct Gtype));
    
    minmod2grav(&Dm,&Dp,D+dir);
}
  
  // Interpolation
  int ix,iy,iz;
  int icell;
  
for(iz=0;iz<2;iz++){
    for(iy=0;iy<2;iy++){
      for(ix=0;ix<2;ix++){
	icell=ix+iy*2+iz*4;
	
	//interpminmodgrav(&U0,&Up,D+ix,D+2+iy,D+4+iy,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation
	interpminmodgrav(&U0,&Up,D,D+1,D+2,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation
	memcpy(Wi+icell,&Up,sizeof(struct Gtype));
      }
    }
  }

}


// ==============================================================================
// ==============================================================================
void coarse2fine_gravlin(struct CELL *cell, struct Gtype *Wi){ 

  
  struct OCT * oct;
	  
  struct Gtype U0;
  struct Gtype Up;
  struct Gtype Um;
  struct Gtype D[6];
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
    
    diffUgrav(&U0,&Um,D+2*dir); 
    diffUgrav(&Up,&U0,D+2*dir+1); 

    /* memcpy(D+(2*dir+0),&Dm,sizeof(struct Gtype)); */
    /* memcpy(D+(2*dir+1),&Dp,sizeof(struct Gtype)); */
    
    minmod2grav_mix(D+2*dir,D+2*dir+1);
}
  
  // Interpolation
  int ix,iy,iz;
  int icell;
  
for(iz=0;iz<2;iz++){
    for(iy=0;iy<2;iy++){
      for(ix=0;ix<2;ix++){
	icell=ix+iy*2+iz*4;
	
	interpminmodgrav(&U0,&Up,D+ix,D+2+iy,D+4+iz,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation
	//interpminmodgrav(&U0,&Up,D,D+1,D+2,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation
	memcpy(Wi+icell,&Up,sizeof(struct Gtype));
      }
    }
  }

}


void coarse2fine_forcelin(struct CELL *cell, REAL *Wi){ 

  
  struct OCT * oct;
	  
  REAL U0[3];
  REAL Up[3];
  REAL Um[3];
  REAL *W;
  REAL D[18];
   /*  struct Gtype Up; */
  /* struct Gtype Um; */
  /* struct Gtype D[3]; */
  /* struct Gtype *W; */
  int inei2;
  int vcell[6],vnei[6];
  int dir;
  REAL dxcur;
  
  oct=cell2oct(cell);
  getcellnei(cell->idx, vnei, vcell); // we get the neighbors
  dxcur=pow(0.5,oct->level);
  
  W=cell->f;
  memcpy(U0,W,sizeof(REAL)*3);

  // Limited Slopes
  for(dir=0;dir<3;dir++){
    
    inei2=2*dir;
    if(vnei[inei2]==6){
      W=oct->cell[vcell[inei2]].f;
    }
    else{
      W=oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].f;
      
#ifdef TRANSXM
      if((oct->x==0.)&&(inei2==0)){
	W=cell->f;
      }
#endif
      
#ifdef TRANSYM
      if((oct->y==0.)&&(inei2==2)){
	W=cell->f;
      }
#endif
      
#ifdef TRANSZM
      if((oct->z==0.)&&(inei2==4)){
	W=cell->f;
      }
#endif
      
      
    }
    memcpy(Um,W,sizeof(REAL)*3);
    
    inei2=2*dir+1;
    if(vnei[inei2]==6){
      W=oct->cell[vcell[inei2]].f;
    }
    else{
      W=oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].f;
      
#ifdef TRANSXP
      if(((oct->x+2.*dxcur)==1.)&&(inei2==1)){
	W=cell->f;
      }
#endif
      
#ifdef TRANSYP
      if(((oct->y+2.*dxcur)==1.)&&(inei2==3)){
	W=cell->f;
      }
#endif
      
#ifdef TRANSZP
      if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
	W=cell->f;
      }
#endif
      
    }
    
    memcpy(Up,W,sizeof(REAL)*3);
    
    diffforce(U0,Um,D+2*dir*3); 
    diffforce(Up,U0,D+(2*dir+1)*3); 

 
}
  
  // Interpolation
  int ix,iy,iz;
  int icell;
  
for(iz=0;iz<2;iz++){
    for(iy=0;iy<2;iy++){
      for(ix=0;ix<2;ix++){
	icell=ix+iy*2+iz*4;
	
	interpforce(U0,Up,D+ix*3,D+6+iy*3,D+12+iz*3,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation
	//Up[2]=0.;Up[1]=0.;Up[0]=0.;
	memcpy(Wi+icell*3,Up,sizeof(REAL)*3);
      }
    }
  }

}

// ============================================================================================================
// ============================================================================================================
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
    for(icell=0;icell<8;icell++) stencil->oct[ioct].cell[icell].gdata.p=neicell->child->cell[face[icell]].gdata.p;
  }
  else{
    coarse2fine_gravlin(neicell,Wi);
    for(icell=0;icell<8;icell++){
      memcpy(&(stencil->oct[ioct].cell[icell].gdata),Wi+face[icell],sizeof(struct Gtype)); //
    }
  }



 
 
}

//================================================================
//================================================================


struct OCT *gatherstencilgrav(struct OCT *octstart, struct STENGRAV *gstencil, int stride, struct CPUINFO *cpu, int *nread)
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
  char visit[7]={0,0,0,0,0,0,0};
  int ioct;

  struct GGRID *stencil=gstencil->stencil;

  nextoct=octstart;
  if(nextoct!=NULL){
    do{ //sweeping levels
      curoct=nextoct;
      nextoct=curoct->next;

      if(curoct->cpu!=cpu->rank) continue;

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
	stencil[iread].oct[6].cell[icell].gdata.d=curoct->cell[icell].gdata.d;
	stencil[iread].oct[6].cell[icell].gdata.p=curoct->cell[icell].gdata.p;
      }
      visit[6]=1;

      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  //ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  ioct=inei;
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

struct OCT *scatterstencilgrav(struct OCT *octstart, struct STENGRAV *gstencil,int nread,int stride, struct CPUINFO *cpu)
{


  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;
  int level=octstart->level;
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
	curoct->cell[icell].pnew=gstencil->pnew[icell+iread*8];
	curoct->cell[icell].res =gstencil->res [icell+iread*8];
      }

#ifdef ONFLYRED
      if(level<=cpu->levelcoarse){
	curoct->parent->gdata.d=gstencil->resLR[iread];
	curoct->parent->gdata.p=0.;
      }
#endif  

      iread++;
    }while((nextoct!=NULL)&&(iread<nread));
  }
  return nextoct;
}

//==============================================================================
//==============================================================================
//==============================================================================

struct OCT *scatterstencilforce(struct OCT *octstart, struct STENGRAV *gstencil, int nread, int stride, struct CPUINFO *cpu,int dir)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  int ipos,iread;
  int icell;

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
	curoct->cell[icell].f[dir]=gstencil->pnew[icell+iread*8];
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<nread));
  }
  return nextoct;
}




//============================================================================
int PoissonJacobi_single(struct STENGRAV *gstencil, int level, int curcpu, int nread,int stride,REAL dx, int flag, REAL factdens){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int inei,icell;
  int i;
  REAL temp;
  REAL res;
  int vnei[6],vcell[6];
  //int ioct[7]={12,14,10,16,4,22,13};
  int ioct[7]={0,1,2,3,4,5,6};

  struct Gtype *curcell;
  struct GGRID *stencil=gstencil->stencil;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
    for(i=0;i<nread;i++){ // we scan the octs

#ifdef ONFLYRED
      if(icell==0) gstencil->resLR[i]=0.;
#endif
      
      temp=0.;
      res=0.;
      
      curcell=&(stencil[i].oct[ioct[6]].cell[icell].gdata);
      
      // Computing the laplacian ===========================
 
      for(inei=0;inei<6;inei++){
 	temp+=stencil[i].oct[ioct[vnei[inei]]].cell[vcell[inei]].gdata.p;
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
      gstencil->pnew[icell+i*8]=temp;
      // we store the local residual
      if(flag) {
	gstencil->res[i*8+icell]=factdens*curcell->d;
      }
      else{
	gstencil->res[icell+i*8]=res;
      }

#ifdef ONFLYRED
      // the low resolution residual
      gstencil->resLR[i]+=res*0.125;
#endif
      // ready for the next cell
    }

    //ready for the next oct
  }
  return 0;
}


//============================================================================



 //============================================================================
REAL comp_residual(struct STENGRAV *gstencil, int level, int curcpu, int nread,int stride, int flag){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int icell;
  int i;

  REAL residual=0.;
  REAL rloc;

  if(flag){
    for(i=0;i<nread;i++){ // we scan the octs
      for(icell=0;icell<8;icell++){ // we scan the cells
	rloc=pow(gstencil->res[icell+i*8],2);
	residual+=rloc;
	// ready for the next cell
      }
      //ready for the next oct
    }
  }
  else{
    for(i=0;i<nread;i++){ // we scan the octs
      for(icell=0;icell<8;icell++){ // we scan the cells
	rloc=pow(gstencil->res[icell+i*8],2);
	residual=(residual>rloc?residual:rloc);
      }
      //ready for the next oct
    }
  }

  return residual;
}
 
//-------------
 


//============================================================================
int Pot2Force(struct STENGRAV *gstencil, int level, int curcpu, int nread,int stride,REAL dx, REAL tsim, int dir){

  int inei,icell;
  int i;
  int vnei[6],vcell[6];
  //  int ioct[7]={12,14,10,16,4,22,13};
  int ioct[7]={0,1,2,3,4,5,6};
  REAL floc;
  struct GGRID *stencil=gstencil->stencil;

  for(icell=0;icell<8;icell++){ // we scan the cells
    getcellnei(icell, vnei, vcell); // we get the neighbors
    for(i=0;i<nread;i++){ // we scan the octs
      
      floc=0.5*(stencil[i].oct[ioct[vnei[(dir<<1)+1]]].cell[vcell[(dir<<1)+1]].gdata.p-stencil[i].oct[ioct[vnei[(dir<<1)]]].cell[vcell[(dir<<1)]].gdata.p)/dx*tsim;

      // store the force
      
      //gstencil->pnew[icell*nread+i]=floc;;
      gstencil->pnew[icell+i*8]=floc;;

      // ready for the next cell
    }
    //ready for the next oct
  }
  return 0;
}

//=============================================================================

REAL PoissonJacobi(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim)
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
  REAL res0;

  // Computing the factor of the density
  if(level>=param->lcoarse){
#ifndef TESTCOSMO
    factdens=4.0*M_PI;
#else
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
  //if(level==6) nitmax=10;

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
		
	
	// ------------ gathering the stencil value values
	nextoct=gatherstencilgrav(curoct,stencil,stride,cpu, &nread);
	
	// ------------ solving the hydro
	PoissonJacobi_single(stencil,level,cpu->rank,nread,stride,dxcur,(iter==0),factdens);

	// ------------ computing the residuals

	rloc=comp_residual(stencil,level,cpu->rank,nread,stride,(iter==0));

	if(iter==0){
	  fnorm+=rloc;
	}
	else{
	  //residual=(residual>rloc?residual:rloc);
	  residual+=rloc;
	}
	
	// ------------ scatter back the data
	
	nextoct=scatterstencilgrav(curoct,stencil,nread,stride, cpu);
	
	tall+=temps[9]-temps[0];
	tcal+=temps[7]-temps[3];
	tscat+=temps[9]-temps[7];
	tgat+=temps[3]-temps[0];
	nreadtot+=nread;
      }while(nextoct!=NULL);
    }

    if((iter==1)&&(level>param->lcoarse)) res0=residual;


    tt=MPI_Wtime();
    // at this stage an iteration has been completed : let's update the potential
    if(nreadtot>0){
      curoct=firstoct[level-1];
      if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
	nextoct=curoct;
	do{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
	  for(icell=0;icell<8;icell++){	
	    
	    REAL w;
	    if(level>param->lcoarse){
	      w=1.0;
	    }
	    else{
	      w=1.0;
	    }
	    curoct->cell[icell].gdata.p=curoct->cell[icell].pnew*w+(1.-w)*curoct->cell[icell].gdata.p;
 	  }
	}while(nextoct!=NULL);
      }
    }
    
    tstop=MPI_Wtime();
    tup+=(tstop-tt);
    tglob+=(tstop-tstart);
    
    if(iter>0){
      if(level<param->lcoarse){
	dres=sqrt(residual);
      }
      else{
	dres=sqrt(residual/fnorm);
      }
      if((dres)<param->poissonacc){
	if(level>=param->lcoarse) break;
      }
    }
  }

  if(level>param->lcoarse){
    printf("CPU | level=%d iter=%d res=%e res0=%e\n",level,iter,dres,sqrt(res0/fnorm));
  }
  else{
    printf("CPU | level=%d iter=%d res=%e\n",level,iter,dres);
  }

  return dres;
}





//===============================================================================================

REAL PoissonMgrid(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim)
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


#ifndef ONFLYRED
  // NOTE ON GPU the calculation is performed on the fly
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
#endif

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
	coarse2fine_gravlin(curcell,Wi);
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

int PoissonForce(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim){
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dxcur;
  int icomp,icell;
  struct PART *curp;
  struct PART *nexp;
  int nread,nreadtot;
  int idir;
  dxcur=pow(0.5,level);
  
  nextoct=firstoct[level-1];
  nreadtot=0;
  if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
    do {
      curoct=nextoct;
      nextoct=curoct->next; 
      
      // ------------ gathering the stencil value values
      nextoct=gatherstencilgrav(curoct,stencil,stride,cpu, &nread);
	
      
      for(idir=0;idir<3;idir++){
	// ------------ Computing the potential gradient
	Pot2Force(stencil,level,cpu->rank,nread,stride,dxcur,tsim,idir);
      
	// ------------ scatter back the data
	
	nextoct=scatterstencilforce(curoct,stencil, nread, stride, cpu,idir);
      }
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

  avgdens/=nc;

  


  return 0;
}

#endif

//===================================================================================================================================
//===================================================================================================================================
#ifdef WGRAV
int PoissonSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL aexp){

  REAL res;
  int igrid;
  struct Gtype Wi[8];
  struct CELL* curcell;
  int icell2;
  struct OCT* curoct;
  struct OCT* nextoct;
  int icell;
  double t[10];

 t[0]=MPI_Wtime();
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

	      coarse2fine_gravlin(curcell,Wi);
	      for(icell2=0;icell2<8;icell2++){
		//		Wi[icell2].p=0.;
		memcpy(&(curcell->child->cell[icell2].gdata.p),&(Wi[icell2].p),sizeof(REAL));
		//memcpy(&(curcell->child->cell[icell2].gdata.p),&(curcell->gdata.p),sizeof(REAL));
	      }

	    }
	  }
      }while(nextoct!=NULL);
  }

 t[9]=MPI_Wtime();
  
#ifndef GPUAXL
 printf("==== CPU POISSON TOTAL TIME =%e\n",t[9]-t[0]);
#else
 printf(" === GPU POISSON TOTAL TIME =%e\n",t[9]-t[0]);
#endif
  
  return 0;
}
#endif
