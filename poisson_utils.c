#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "oct.h"
#include <string.h>

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
      if((oct->nei[vnei[inei2]]->child->x-oct->x)>0.5){
		W=&(cell->gdata);
      }
#endif
      
#ifdef TRANSYM
      if((oct->nei[vnei[inei2]]->child->y-oct->y)>0.5){
	W=&(cell->gdata);
      }
#endif
      
#ifdef TRANSZM
      if((oct->nei[vnei[inei2]]->child->z-oct->z)>0.5){
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
      if((oct->nei[vnei[inei2]]->child->x-oct->x)<0.){
	W=&(cell->gdata);
      }
#endif
      
#ifdef TRANSYP
      if((oct->nei[vnei[inei2]]->child->y-oct->y)<0.){
	W=&(cell->gdata);
      }
#endif
      
#ifdef TRANSZP
      if((oct->nei[vnei[inei2]]->child->z-oct->z)<0.){
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
void recursive_neighbor_gather_oct_grav(int ioct, int inei, int inei2, int inei3, int order, struct CELL *cell, struct HGRID *stencil){

  int ix[6]={-1,1,0,0,0,0};
  int iy[6]={0,0,-1,1,0,0};
  int iz[6]={0,0,0,0,-1,1};
  int icell;
  int i;
  int ioct2;
  int vnei[6],vcell[6];
  int ineiloc;
  int face[8]={0,1,2,3,4,5,6,7};
  REAL dist;

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


  oct=cell2oct(cell);
  neioct=cell2oct(neicell);

  // ============================ TRANSMISSIVE BOUNDARIES ====================
#ifdef TRANSXP
    if(ineiloc==1){
      dist=neioct->x-oct->x;
      if(dist<0.){
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
      dist=neioct->y-oct->y;
      if(dist<0.){
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
      dist=neioct->z-oct->z;
      if(dist<0.){
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
      dist=neioct->x-oct->x;
      if(dist>0.5){
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
      dist=neioct->y-oct->y;
      if(dist>0.5){
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
      dist=neioct->z-oct->z;
      if(dist>0.5){
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
    for(icell=0;icell<8;icell++) memcpy(Wi+icell,&(neicell->child->cell[icell].gdata),sizeof(struct Gtype));
  }
  else{
    //coarse2fine_grav(neicell,Wi);
    for(icell=0;icell<8;icell++) memcpy(Wi+icell,&(neicell->gdata),sizeof(struct Gtype));
  }



  for(icell=0;icell<8;icell++){
    memcpy(&(stencil->oct[ioct].cell[icell].gdata),Wi+face[icell],sizeof(struct Gtype)); //
  }

  // next order
  if(order==1){
    for(i=0;i<6;i++){
      if((i/2)==(inei/2)) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      recursive_neighbor_gather_oct_grav(ioct2, inei, i, -1, 2, neicell, stencil);
    }
  }
  else if(order==2) {
    for(i=0;i<6;i++){
      if(((i/2)==(inei/2))||((i/2)==(inei2/2))) continue;
      ioct2=ioct+ix[i]+iy[i]*3+iz[i]*9; // oct position in stencil
      recursive_neighbor_gather_oct_grav(ioct2, inei, inei2, i, 3, neicell, stencil);
    }
  }
}

//================================================================
//================================================================

struct OCT *gatherstencilgrav(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread)
{
  struct OCT* nextoct;
  struct OCT* curoct;
  struct CELL *cell;

  int inei;
  int iread=0;
  int icell;
  //int ioct[7]={12,14,10,16,4,22,13};
  
  int ix[6]={-1,1,0,0,0,0};
  int iy[6]={0,0,-1,1,0,0};
  int iz[6]={0,0,0,0,-1,1};
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
	memcpy(&(stencil[iread].oct[13].cell[icell].gdata),&(curoct->cell[icell].gdata),sizeof(struct Gtype)); //
#ifdef TESTCOSMO
	stencil[iread].oct[13].cell[icell].gdata.d-=1.0;
#endif
      }

      //abort();
      cell=curoct->parent;
      
      //start recursive fill
      for(inei=0;inei<6;inei++)
	{
	  ioct=ix[inei]+iy[inei]*3+iz[inei]*9+13; // oct position in stencil
	  recursive_neighbor_gather_oct_grav(ioct, inei, -1, -1, 1, cell, stencil+iread);
	}
	  
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  (*nread)=iread;
  return nextoct;
}

//==============================================================================
//==============================================================================

struct OCT *scatterstencilgrav(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu)
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
	memcpy(&(curoct->cell[icell].pnew),&(stencil[iread].new.cell[icell].pnew),sizeof(REAL));
	memcpy(&(curoct->cell[icell].res),&(stencil[iread].new.cell[icell].res),sizeof(REAL));
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
	memcpy(curoct->cell[icell].f,stencil[iread].new.cell[icell].f,sizeof(REAL)*3);
      }
      iread++;
    }while((nextoct!=NULL)&&(iread<stride));
  }
  return nextoct;
}




//============================================================================
int PoissonJacobi_single(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, int flag, REAL factdens){

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
      }

      // setting up the residual
      res=temp;
      
      // we finish the laplacian
      temp=temp/6.0f;
      temp=temp-dx*dx*curcell->d/6.0f*factdens;

      // we finsih the residual
      res=res-6.0*curcell->p;
      res=res/(dx*dx)-factdens*curcell->d;

      // we store the new value of the potential
      stencil[i].new.cell[icell].pnew=temp;

      // we store the local residual
      if(flag) {
	stencil[i].new.cell[icell].res=factdens*curcell->d;
      }
      else{
	stencil[i].new.cell[icell].res=res;
      }

      // ready for the next cell
    }
    //ready for the next oct
  }
  return 0;
}



//============================================================================
int Pot2Force(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx){

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

      floc[0]=0.5*(stencil[i].oct[ioct[vnei[1]]].cell[vcell[1]].gdata.p-stencil[i].oct[ioct[vnei[0]]].cell[vcell[0]].gdata.p)/dx;
      floc[1]=0.5*(stencil[i].oct[ioct[vnei[3]]].cell[vcell[3]].gdata.p-stencil[i].oct[ioct[vnei[2]]].cell[vcell[2]].gdata.p)/dx;
      floc[2]=0.5*(stencil[i].oct[ioct[vnei[5]]].cell[vcell[5]].gdata.p-stencil[i].oct[ioct[vnei[4]]].cell[vcell[4]].gdata.p)/dx;

      if(floc[1]!=0.) abort();
      // store the force
      
      memcpy(stencil[i].new.cell[icell].f,floc,3*sizeof(REAL));

      // ready for the next cell
    }
    //ready for the next oct
  }
  return 0;
}

//=============================================================================

REAL PoissonJacobi(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL tsim)
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



  // Computing the factor of the density
  if(level>=param->lcoarse){
#ifndef TESTCOSMO
    factdens=4.0*M_PI;
#else
    factdens=6.0*tsim;
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
  
  for(iter=0;iter<nitmax;iter++){
    // --------------- setting the first oct of the level
    nextoct=firstoct[level-1];
    nreadtot=0;
    if((nextoct!=NULL)&&(cpu->noct[level-1]!=0)){
      do {
	curoct=nextoct;
	nextoct=curoct->next; 
	// -------------  cleaning working arrays
	
	memset(stencil,0,stride*sizeof(struct HGRID));
	
	// ------------ gathering the stencil value values
	nextoct=gatherstencilgrav(curoct,stencil,stride,cpu, &nread);
	
	// ------------ solving the hydro
	
	PoissonJacobi_single(stencil,level,cpu->rank,nread,stride,dxcur,(iter==0),factdens);
	      
	// ------------ scatter back the data
	
	nextoct=scatterstencilgrav(curoct,stencil, nread, cpu);
	nreadtot+=nread;
      }while(nextoct!=NULL);
    }
	  
    // at this stage an iteration has been completed : let's update the potential and compute the residual
    
    if(nreadtot>0){
      
      if(iter==0){
	fnorm=0.;
	residualold=0.;
      }
      else{
	residual=0.;
      }

      curoct=firstoct[level-1];
      if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
	nextoct=curoct;
	do{
	  curoct=nextoct;
	  nextoct=curoct->next;
	  if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
	  for(icell=0;icell<8;icell++){	
	    curoct->cell[icell].gdata.p=curoct->cell[icell].pnew;
	    if(iter==0){
	      fnorm+=pow(curoct->cell[icell].res,2);
	    }
	    else{
	      residual+=pow(curoct->cell[icell].res,2);
	    }
	  }
	}while(nextoct!=NULL);
      }
    }
    if(iter>0){
      dres=sqrt(fabs(residual-residualold)/fnorm);
      //dres=sqrt(fabs(residual-residualold)/fabs(residual+residualold)*0.5);
   
      residualold=residual;
      if(dres<param->poissonacc) break;
    }
  }
  printf("level=%d iter=%d res=%e fnorm=%e\n",level,iter,dres,fnorm);
  return dres;
}

//===============================================================================================

REAL PoissonMgrid(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL tsim)
{
  struct OCT *nextoct;
  struct OCT *curoct;
  int icell;
  REAL dres;
  struct Gtype Wi[8];
  struct CELL* curcell;
  // pre-relaxation

  dres=PoissonJacobi(level,param,firstoct,cpu,stencil,stride,tsim);

  
  if(dres>param->poissonacc){
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
       PoissonJacobi(level-1,param,firstoct,cpu,stencil,stride,tsim);
     }
     else{
       PoissonMgrid(level-1,param,firstoct,cpu,stencil,stride,tsim);
     }

     // prolongation + correction
     /* nextoct=firstoct[level-1]; */
     /* if(nextoct!=NULL){ */
     /*   do{  */
     /* 	 curoct=nextoct; */
     /* 	 nextoct=curoct->next; */
     /* 	 for(icell=0;icell<8;icell++) { */
     /* 	   curoct->cell[icell].gdata.p-=curoct->parent->gdata.p; // we propagate the error and correct the evaluation */
     /* 	 } */
     /*   }while(nextoct!=NULL); */
     /* } */

     nextoct=firstoct[level-1];
     if(nextoct!=NULL){
       do // sweeping level
	 {
	   curoct=nextoct;
	   nextoct=curoct->next;
	   curcell=curoct->parent;
	   coarse2fine_grav(curcell,Wi);
	   for(icell=0;icell<8;icell++) // looping over cells in oct
	     {
	       curoct->cell[icell].gdata.p-=Wi[icell].p; // we propagate the error and correct the evaluation
	     }
	 }while(nextoct!=NULL);
     }

     // post relaxation
     dres=PoissonJacobi(level,param,firstoct,cpu,stencil,stride,tsim);
  }

  return dres;
}

//===================================================================================================================================

int PoissonForce(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride){
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
      nextoct=gatherstencilgrav(curoct,stencil,stride,cpu, &nread);
	
	// ------------ Computing the potential gradient
      
      Pot2Force(stencil,level,cpu->rank,nread,stride,dxcur);
      
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
int FillDens(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu){
  struct OCT *nextoct;
  struct OCT *curoct;
  REAL dx;
  int icomp,icell;
  struct PART *curp;
  struct PART *nexp;
  int nread;
  REAL locdens;
  
  curoct=firstoct[level-1];
  if((curoct!=NULL)&&(cpu->noct[level-1]!=0)){
    nextoct=curoct;
    do{
      curoct=nextoct;
      nextoct=curoct->next;
      if(curoct->cpu!=cpu->rank) continue; // we don't update the boundary cells
      for(icell=0;icell<8;icell++){	
	locdens=0.f;
#ifdef PIC
	locdens+=curoct->cell[icell].density;
#endif
#ifdef WHYDRO2
	locdens+=curoct->cell[icell].field.d;
#endif
	curoct->cell[icell].gdata.d=locdens;
      }
    }while(nextoct!=NULL);
  }

  return 0;
}

#endif
