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


//=============================================================================
void update_pot_in_tree(struct OCT ** octList, int level,struct OCT ** firstoct,  struct CPUINFO *cpu, struct RUNPARAMS *param, REAL *distout, REAL *normpout){
  struct OCT *curoct;
  int icell, iOct;
  REAL dist=0.;
  REAL normp=0.;

  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    curoct=octList[iOct]; 

    for(icell=0;icell<8;icell++){	
      const REAL pnew  = curoct->GDATA_pnew[icell];
      const REAL pold  = curoct->GDATA_p   [icell];
      const REAL delta = pold-pnew;

      dist  += delta*delta;
      normp += pold *pold;
    }
    memcpy(curoct->GDATA_p,curoct->GDATA_pnew, 8*sizeof(REAL));
  }
  *normpout=normp;
  *distout=dist;
}
//================================================================
// ===============================================================
// ===============================================================

void interpminmodgrav(struct Gtype *U0, struct Gtype *Up, struct Gtype *Dx, struct Gtype *Dy, struct Gtype *Dz,REAL dx,REAL dy,REAL dz){
  
  Up->d =U0->d  +dx*Dx->d  +dy*Dy->d  +dz*Dz->d;
  Up->p =U0->p  +dx*Dx->p  +dy*Dy->p  +dz*Dz->p;
}

void interpminmodgrav_REAL(REAL *U0, REAL *Up, REAL Dx, REAL Dy, REAL Dz,REAL dx,REAL dy,REAL dz){
  *Up = *U0  + dx*Dx + dy*Dy +dz*Dz;
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

  
  struct OCT * oct=cell2oct(cell);;
	  
  struct Gtype U0;
  struct Gtype Up;
  struct Gtype Um;
  struct Gtype Dp,Dm;
  struct Gtype D[3];
  struct Gtype *W;
  int inei2;
  int dir;
  REAL dxcur=pow(0.5,oct->level);
    
  const int *vnei  = &NEIG_NEIP[cell->idx*6];
  const int *vcell = &NEIG_CELL[cell->idx*6];

  
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

//	interpminmodgrav(&U0,&Up,D,D+1,D+2,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation


//       interpminmodgrav_REAL(&U0.d,&Up.d, D[0].d, D[1].d, D[2].d,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation
//       interpminmodgrav_REAL(&U0.p,&Up.p, D[0].p, D[1].p, D[2].p,-0.25+ix*0.5,-0.25+iy*0.5,-0.25+iz*0.5); // Up contains the interpolation

	memcpy(Wi+icell,&Up,sizeof(struct Gtype));
      }
    }
  }

}


// ==============================================================================
// ==============================================================================
//void coarse2fine_gravlin(struct CELL *cell, struct Gtype *Wi){ 
void coarse2fine_gravlin(struct CELL *cell, REAL *Wid, REAL *Wip ){ 

  struct OCT * oct=cell2oct(cell);

  struct Gtype W;
  W.d = oct->GDATA_d[cell->idx];
  W.p = oct->GDATA_p[cell->idx];


  struct Gtype U0 =W;
  struct Gtype Up =W;
  struct Gtype Um =W;
  struct Gtype D[6]; 
  


  int inei2;
  int dir;
  REAL dxcur=pow(0.5,oct->level);;


  const int *vnei  = NEIGHBORS_NEIP[cell->idx];
  const int *vcell = NEIGHBORS_CELL[cell->idx];


  // Limited Slopes
  for(dir=0;dir<3;dir++){
    
    inei2=2*dir;
    if(vnei[inei2]==6){

      W.d =oct->GDATA_d[vcell[inei2]] ;
      W.p =oct->GDATA_p[vcell[inei2]] ;

    }
    else{
      if(oct->nei[vnei[inei2]]->child!=NULL){
	W.d = oct->nei[vnei[inei2]]->child->GDATA_d[vcell[inei2]] ;
	W.p = oct->nei[vnei[inei2]]->child->GDATA_p[vcell[inei2]] ;
      }
      else{
	W.d = cell2oct(oct->nei[vnei[inei2]]) -> GDATA_d[oct->nei[vnei[inei2]]->idx] ;
	W.p = cell2oct(oct->nei[vnei[inei2]]) -> GDATA_p[oct->nei[vnei[inei2]]->idx] ;
      }
      
#ifdef TRANSXM
      if((oct->x==0.)&&(inei2==0)){
	W.d = oct->GDATA_d[cell->idx];
	W.p = oct->GDATA_p[cell->idx];
      }
#endif
      
#ifdef TRANSYM
      if((oct->y==0.)&&(inei2==2)){
	W.d = oct->GDATA_d[cell->idx];
	W.p = oct->GDATA_p[cell->idx];      }
#endif
      
#ifdef TRANSZM
      if((oct->z==0.)&&(inei2==4)){
	W.d = oct->GDATA_d[cell->idx];
	W.p = oct->GDATA_p[cell->idx];
      }
#endif
    }

    Um=W;
    
    inei2=2*dir+1;
    if(vnei[inei2]==6){
//    W=&(oct->cell[vcell[inei2]].gdata);
     W.p=oct->GDATA_d[vcell[inei2]];
     W.d=oct->GDATA_p[vcell[inei2]];

    }
    else{
      if(oct->nei[vnei[inei2]]->child!=NULL){
//	W=&(oct->nei[vnei[inei2]]->child->cell[vcell[inei2]].gdata);
	W.d = cell2oct(oct->nei[vnei[inei2]]) -> GDATA_d[oct->nei[vnei[inei2]]->idx] ;
	W.p = cell2oct(oct->nei[vnei[inei2]]) -> GDATA_p[oct->nei[vnei[inei2]]->idx] ;
      }
      else{
//	W=&(oct->nei[vnei[inei2]]->gdata);
        W.p=oct->GDATA_d[vcell[inei2]];
        W.d=oct->GDATA_p[vcell[inei2]];
      }
      
#ifdef TRANSXP
      if(((oct->x+2.*dxcur)==1.)&&(inei2==1)){
	W.d = oct->GDATA_d[cell->idx];
	W.p = oct->GDATA_p[cell->idx];
      }
#endif
      
#ifdef TRANSYP
      if(((oct->y+2.*dxcur)==1.)&&(inei2==3)){
	W.d = oct->GDATA_d[cell->idx];
	W.p = oct->GDATA_p[cell->idx];
      }
#endif
      
#ifdef TRANSZP
      if(((oct->z+2.*dxcur)==1.)&&(inei2==5)){
	W.d = oct->GDATA_d[cell->idx];
	W.p = oct->GDATA_p[cell->idx];
      }
#endif
      
    }
    Up=W;


    diffUgrav(&U0,&Um,D+2*dir); 
    diffUgrav(&Up,&U0,D+2*dir+1); 
    
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
//	memcpy(Wi+icell,&Up,sizeof(struct Gtype));

	memcpy(Wid + icell,&Up.d,sizeof(REAL));
	memcpy(Wip + icell,&Up.p,sizeof(REAL));

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
  int dir;
  REAL dxcur;
  
  oct=cell2oct(cell);
    const int *vnei  = &NEIG_NEIP[cell->idx*6];
    const int *vcell = &NEIG_CELL[cell->idx*6];

  //  int vnei[6],vcell[6];   
    //getcellnei(cell->idx, vnei, vcell); // we get the neighbors

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
	memcpy(Wi+icell*3,Up,sizeof(REAL)*3);
      }
    }
  }

}

// ============================================================================================================
// ============================================================================================================

void recursive_neighbor_gather_oct_grav_trans(int ioct, int * face,  struct CELL *cell,  struct CELL *neicell, struct OCT* oct){
// ============================ TRANSMISSIVE BOUNDARIES ====================
#ifdef TRANSXP
	if(ioct==1){
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
	if(ioct==3){
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
	if(ioct==5){
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
	if(ioct==0){
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
	if(ioct==2){
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
	if(ioct==4){
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
}

//================================================================
//================================================================

//================================================================
//================================================================


//================================================================
void oct2Grid(struct GRID *stencil, struct OCT* oct, const int lMax, const int n, const int x, const int y, const int z, const int dx){

  int icell;

  const int x2 = x+dx;
  const int y2 = y+dx;
  const int z2 = z+dx;

  const int X[8] = { x,x2,x, x2,x, x2,x, x2 };
  const int Y[8] = { y,y, y2,y2,y, y, y2,y2 };
  const int Z[8] = { z,z, z, z, z2,z2,z2,z2 };

  if(oct->level<lMax) {
	  for(icell=0; icell<8; icell++){
		  struct CELL * curcell = &oct->cell[icell];

		  oct2Grid(stencil, curcell->child, lMax, n,X[icell],Y[icell],Z[icell],dx/2);	
    }  
  }else{	
	  for(icell=0; icell<8; icell++){
      const int cellID = Z[icell]*n*n+Y[icell]*n+X[icell];
      stencil->d[cellID] = oct->GDATA_d[icell];
			stencil->p[cellID] = oct->GDATA_p[icell]; 

    }
	}
}

//================================================================
void grid2Oct(struct GRID *stencil, struct OCT* oct, const int lMax, const int n, const int x, const int y, const int z, const int NBND, const int dx){
  int icell;

  const int x2 = x+dx;
  const int y2 = y+dx;
  const int z2 = z+dx;

  const int X[8] = { x,x2,x, x2,x, x2,x, x2 };
  const int Y[8] = { y,y, y2,y2,y, y, y2,y2 };
  const int Z[8] = { z,z, z, z, z2,z2,z2,z2 };

	if(oct->level<lMax) {
	  for(icell=0; icell<8; icell++){
		  struct CELL * curcell = &oct->cell[icell];
	    grid2Oct(stencil, curcell->child, lMax, n, X[icell],Y[icell],Z[icell], NBND, dx/2);	
    }	
  }else{
	  for(icell=0; icell<8; icell++){
      const int cellID = Z[icell]*n*n+Y[icell]*n+X[icell];
      oct->GDATA_pnew[icell] = stencil->pnew[cellID];
		  oct->GDATA_res [icell] = stencil->res [cellID];
    }    

#ifdef ONFLYRED
//  if(level<=cpu->levelcoarse){
	    struct CELL *parentCell = oct->parent; 
      struct OCT  *parentOct  = cell2oct(parentCell);

      const int N_LR      = (n-2*NBND)/2; 
      const int x_LR      = (x2-NBND)/2;
      const int y_LR      = (y2-NBND)/2;
      const int z_LR      = (z2-NBND)/2;

      const int cellID_LR = z_LR*N_LR*N_LR + y_LR*N_LR + x_LR;

      parentOct->GDATA_d[parentCell->idx] = stencil->resLR[cellID_LR];
	    parentOct->GDATA_p[parentCell->idx] = 0;
//}
#endif

	}
}

//================================================================
void setFace(struct GRID * stencil, struct OCT* oct, int lMax, int n, const int x, const int y, const int z, const int *neicell, const int dx){

  int icell;
  int cellID;

  const int x2 = x+dx;
  const int y2 = y+dx;
  const int z2 = z+dx;

  const int X[8] = { x,x2,x, x2,x, x2,x, x2 };
  const int Y[8] = { y,y, y2,y2,y, y, y2,y2 };
  const int Z[8] = { z,z, z, z, z2,z2,z2,z2 };

  if(oct->level<lMax ){
    for(icell=0; icell<4; icell++){
		  cellID = neicell[icell];
      struct CELL * curcell = &oct->cell[cellID];
    	setFace (stencil, curcell->child, lMax, n, X[cellID],Y[cellID],Z[cellID], neicell, dx/2);	
    }
  }else{
    for(icell=0; icell<8; icell++){
      cellID = Z[icell]*n*n+Y[icell]*n+X[icell];
      stencil->d[cellID] = oct->GDATA_d[icell];
		  stencil->p[cellID] = oct->GDATA_p[icell];
    }
  }
}

//================================================================

const static int NEI_CELL_FACE[24] = {	
	1,3,5,7,
	0,2,4,6,
 	2,3,6,7,
	0,1,4,5,
 	4,5,6,7,
	0,1,2,3			  
};

//================================================================
void setBound(struct GRID * stencil, struct OCT* oct, const int lMax, const int N, const int bound, const int dx){

  const int ncell = N - 2*bound;

  const int origine[18] = {
    2-ncell, bound,   bound,
    N-bound, bound,   bound,
    bound,   2-ncell, bound,
    bound,   N-bound, bound,
    bound,   bound,   2-ncell,
    bound,   bound,   N-bound
  };

  int ioct;
  for(ioct=0; ioct<6; ioct++){

		struct OCT* neioct = oct->nei[ioct]->child ;
    const int xc = origine[ioct*3 +0];    
    const int yc = origine[ioct*3 +1];    
    const int zc = origine[ioct*3 +2];    

    const int *neicell = NEI_CELL_FACE + ioct*4;

    setFace(stencil, neioct, lMax, N, xc,yc,zc, neicell, dx);
  }
}

//================================================================
void gatherstencilgrav_grid(struct GRID *grid, struct OCT* curoct, const int level, struct CPUINFO *cpu, struct RUNPARAMS *param, const int NCELL, const int NBND){

	const int N  = NCELL+2*NBND;
	const int dx = pow(2, param->gstride_grid); 
  oct2Grid(grid, curoct, level, N, NBND,NBND,NBND,dx);
  setBound(grid, curoct, level, N, NBND, dx);
}

//================================================================
void scatterstencilgrav_grid(struct GRID *grid, struct OCT* curoct, const int level, struct CPUINFO *cpu, struct RUNPARAMS *param, const int NCELL, const int NBND){
	const int dx = pow(2, param->gstride_grid); 
  grid2Oct(grid, curoct, level, NCELL+2*NBND, NBND,NBND,NBND, NBND,dx);
}

//============================================================================
int PoissonJacobi_single_grid(struct GRID *grid,int NCELL,int NBND, REAL dx, int flag, REAL factdens){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int i,j,k;
  REAL temp;
  REAL res;
  int idx, idxLR;

  //==============================================================
  // 2 DEF NCELL, idx,idxLR

  // we scan the cells of the grid


#ifdef ONFLYRED
  const int NCELL_LR = NCELL/2;
  memset(grid->resLR,0,(NCELL_LR*NCELL_LR*NCELL_LR)*sizeof(REAL));
#endif

  int Ntot = NCELL+2*NBND;
  for(k=0;k<NCELL;k++){
    for(j=0;j<NCELL;j++){
      for(i=0;i<NCELL;i++){
		idx=(i+NBND)+(j+NBND)*Ntot+(k+NBND)*Ntot*Ntot;

    temp=0.;
    res=0.;

    // Computing the laplacian ===========================
    temp+=grid->p[idx-1                            ]; // xm
    temp+=grid->p[idx+1                            ]; // xp
    temp+=grid->p[idx-(NCELL+2*NBND)               ]; // ym
    temp+=grid->p[idx+(NCELL+2*NBND)               ]; // yp
    temp+=grid->p[idx-(NCELL+2*NBND)*(NCELL+2*NBND)]; // zm
    temp+=grid->p[idx+(NCELL+2*NBND)*(NCELL+2*NBND)]; // zp

    // setting up the residual
    res=temp;

    // we finish the laplacian
    temp=temp/6.0;
    temp=temp-dx*dx*grid->d[idx]/6.0*factdens;

    // we finsih the residual
    res=res-6.0*grid->p[idx];
    res=res/(dx*dx)-factdens*grid->d[idx];

    // we store the new value of the potential

    grid->pnew[idx]=temp;

    // we store the local residual
    if(flag) {      grid->res[idx]=factdens*grid->d[idx];     }
    else{           grid->res[idx]=res;                       }

#ifdef ONFLYRED
    // the low resolution residual
    idxLR=i/2+j/2*NCELL_LR+k/2*NCELL_LR*NCELL_LR; // A OPTIMISER !!!!! 

    grid->resLR[idxLR]+=res*0.125;
#endif

      }
    }
  }
  return 0;
}

//============================================================================
REAL comp_residual_grid(struct GRID *grid, int NCELL, int NBND){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int i,j,k;
  int idx;
  REAL residual=0.;
  REAL rloc;

  for(k=0;k<NCELL;k++){
    for(j=0;j<NCELL;j++){
      for(i=0;i<NCELL;i++){
				idx=(i+NBND)+(j+NBND)*(NCELL+2*NBND)+(k+NBND)*(NCELL+2*NBND)*(NCELL+2*NBND);
		    rloc=grid->res[idx];
		    residual+=rloc*rloc;
      }
    }
  }

  return residual;
}


REAL PoissonJacobi_grid(struct OCT *** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim, struct GRID *grid)
{
  int iter, iOct, icell;
  int nread, nreadtot;
  int nitmax;
  int crit;

  REAL fnorm,residual,residualold,dres;
  REAL dxcur = pow(0.5,level);
  REAL factdens;
  REAL rloc;
  REAL res0=0.;
  REAL dist,normp,dresconv;

  double tall=0.,tcal=0.,tscat=0.,tgat=0., tglob=0.,tup=0.;
  double tstart,tstop,tt;
  double temps[10];


	const int NBND		=	2;
  const int NCELL		= pow(2,1+param->gstride_grid) ;
	const int N       = NCELL+2*NBND;

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

  if((level==param->mgridlmin)||(level>param->lcoarse))   	
    nitmax=param->niter;
  else
    nitmax=param->nrelax;
 

  for(iter=0;iter<nitmax;iter++){
    tstart=MPI_Wtime();

    // --------------- some inits for iterative solver
    if(iter==0){
      fnorm   =0.;
      residual=0.;
      nread   =0;
    }
    else{
      residual=0.;
    }

    // Scanning the octs
    const int newLevel = level-1-param->gstride_grid;
    struct OCT ** curOctList = octList[newLevel];

struct OCT  *curoct;


		for(iOct=0; iOct<cpu->locNoct[newLevel]; iOct++){
		  curoct = curOctList[iOct];

			temps[0]=MPI_Wtime();

			// ------------ gathering the stencil value 

			gatherstencilgrav_grid(grid, curoct, level, cpu,param, NCELL, NBND);

	    temps[2]=MPI_Wtime();

			// ------------ solving the hydro
			PoissonJacobi_single_grid(grid, NCELL, NBND, dxcur,  (iter==0), factdens);

			// ------------ computing the residuals

      rloc=comp_residual_grid(grid, NCELL, NBND);

			if(iter==0){
				fnorm+=rloc;
				//printf("rloc=%e\n",rloc);
			}
			else{
				//residual=(residual>rloc?residual:rloc);
				residual+=rloc;
			  //printf("iter=%d rloc=%e\n",iter,rloc);
			}
	
			// ------------ scatter back the data in the tree
			temps[4]=MPI_Wtime();

			scatterstencilgrav_grid(grid, curoct, level, cpu, param, NCELL, NBND);
	
  		temps[6]=MPI_Wtime();

			tall    += temps[6]-temps[0];
			tcal    += temps[4]-temps[2];
			tscat   += temps[6]-temps[4];
			tgat    += temps[2]-temps[0];
			nreadtot+= nread;
		}

    //printf("iter=%d nreadtot=%d\n",iter,nreadtot);

    tt=MPI_Wtime();
    // at this stage an iteration has been completed : let's update the potential in the tree

     update_pot_in_tree(octList[level-1], level,firstoct,cpu,param,&dist,&normp);

#ifdef WMPI
    // printf("iter=%d\n",iter);

    // reduction of relevant quantities

    if((iter<=param->niter)||(iter%1==0)){
      mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,(iter==0),level); // potential field exchange
      if(iter==0){
	      //if(level==7) printf("rank=%d fnorm=%e\n",cpu->rank,fnorm);
	      MPI_Allreduce(MPI_IN_PLACE,&fnorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }else{
	      MPI_Allreduce(MPI_IN_PLACE,&residual,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	      MPI_Allreduce(MPI_IN_PLACE,&dist,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	      MPI_Allreduce(MPI_IN_PLACE,&normp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
    }
#endif

    if((iter==1)&&(level>=param->lcoarse)) res0=residual;

    tstop=MPI_Wtime();
    tup+=(tstop-tt);
    tglob+=(tstop-tstart);
    
    if(iter>0){

      // here we test the convergence of the temporary solution
      dresconv=sqrt(dist/normp);

      // here we test the zero level of Poisson equation
      if(level<param->lcoarse){      	dres=sqrt(residual);            }
      else{                         	dres=sqrt(residual/fnorm);      }

      // we take the smallest
      crit=(dres<dresconv?0:1);
      dres=(dres<dresconv?dres:dresconv);

      if((dres)<param->poissonacc){
      	if(level>=param->lcoarse) break;
      }
    }
}
  if(cpu->rank==0){
    printf("CPU | tgat=%e tcal=%e tscat=%e\n",tgat,tcal,tscat);
    printf("CPU | level=%d iter=%d res=%e fnorm=%e",level,iter,dres,fnorm);
    if(level<=param->lcoarse) printf(" resraw=%e res0=%e crit=%d",sqrt(residual),sqrt(res0),crit);
    printf("\n");
  }
  return dres;
}


REAL PoissonMgrid_grid(struct OCT *** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim, struct GRID *grid)
{
  struct OCT *curoct;
  struct Gtype Wi[8];

  int icell, iOct;
  REAL dres;


  // pre-relaxation
#ifdef WMPI
  mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,1,level); // potential field exchange
#endif

#ifndef GPUAXL
  dres=PoissonJacobi_grid(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
#else
  if(level>=7){
    //dres=PoissonJacobiGPU(level,param,firstoct,cpu,stencil,stride,tsim);
    dres=PoissonJacobi_grid(octList, level,param,firstoct,cpu,stencil,stride,tsim, grid);
  }
  else{
    dres=PoissonJacobi_grid(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
  }
  //dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim);
#endif


#ifndef ONFLYRED
  // NOTE ON GPU the calculation is performed on the fly
  // reduction

  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    curoct=octList[level-1][iOct]; 

    curoct->parent->gdata.d=0.;
    curoct->parent->gdata.p=0.;

    REAL tmp=0;
    for(icell=0;icell<8;icell++){
     tmp += curoct->res[icell]*0.125; // we average the residual and store it as the new density
    }
     cell2oct(curoct->parent)->GDATA_p[curoct->parent->idx] = tmp;
  }
#endif

  // full relaxation at coarsest level or recursive call to mgrid
  
  if((level-1)==param->mgridlmin){
#ifdef WMPI
    mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,1,level-1); // potential field exchange
#endif

#ifndef GPUAXL
    PoissonJacobi_grid(octList, level-1,param,firstoct,cpu,stencil,stride,tsim,grid);
#else
    //PoissonJacobiGPU(level-1,param,firstoct,cpu,stencil,stride,tsim);
    PoissonJacobi_grid(octList, level-1,param,firstoct,cpu,stencil,stride,tsim,grid);
#endif
  }
  else{
    PoissonMgrid_grid(octList, level-1,param,firstoct,cpu,stencil,stride,tsim, grid);
  }
  
  // prolongation + correction
 
  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){

    curoct=octList[level-1][iOct]; 

//    REAL tmp[8];
//    coarse2fine_gravlin(curoct->parent, curoct->GDATA_d, tmp );

    for(icell=0;icell<8;icell++){ // looping over cells in oct 
        curoct->GDATA_p[icell] -=  cell2oct(curoct->parent)->GDATA_p[curoct->parent->idx] ;

//        curoct->GDATA_p[icell] -= tmp[icell];// we propagate the error and correct the evaluation
    }
  }


  // post relaxation
#ifdef WMPI		
  mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,1,level); // potential field exchange
#endif

#ifndef GPUAXL
  dres=PoissonJacobi_grid(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
#else
  dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim);

  if(level>=7){
    //dres=PoissonJacobiGPU(level,param,firstoct,cpu,stencil,stride,tsim);
    dres=PoissonJacobi_grid(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
  }
  else{
    dres=PoissonJacobi_grid(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
  }
#endif
  return dres;
}



//================================================================
//================================================================

//================================================================
//================================================================

//================================================================
//================================================================

void recursive_neighbor_gather_oct_grav(int ioct, struct CELL *cell, struct GGRID *stencil){

  struct CELL *neicell;
  static int face[8]={0,1,2,3,4,5,6,7};

  if(cell->child!=NULL){
    // the oct at the right level exists
    neicell=cell->child->nei[ioct];
  
    const int *vnei  = NEIGHBORS_NEIP[cell->idx];
    const int *vcell = NEIGHBORS_CELL[cell->idx];

    struct OCT  *oct=cell2oct(cell);

    if(vnei[ioct]==6){
      neicell=&(oct->cell[vcell[ioct]]);
    }
    else{
      if(oct->nei[ioct]->child!=NULL){
	      neicell=&(oct->nei[ioct]->child->cell[vcell[ioct]]);
	// ============================ TRANSMISSIVE BOUNDARIES ====================
        recursive_neighbor_gather_oct_grav_trans(ioct, &face[0], cell, neicell, oct);
      }
      else{
				printf("big problem\n");
				abort();
      }
    }
  }


  if(neicell->child!=NULL){   // optimal case

    memcpy( &stencil->p[ioct], &neicell->child->GDATA_p, 8*sizeof(REAL));

  }else{

    coarse2fine_gravlin(neicell, stencil->d[ioct], stencil->p[ioct]);

  }
}


//================================================================
//================================================================

void gatherstencilgrav(struct OCT*** octList, int iOct, struct GGRID *stencil, const int stride, struct CPUINFO *cpu, int *nread, int level){
  int iread;
  for(iread=0; iread<stride; iread++){
    struct OCT  * curoct = octList[level-1][iOct++];
    if(curoct==NULL) break;

    // filling the values in the central oct
    memcpy( &stencil[iread].d[6],  &curoct->GDATA_d, 8*sizeof(REAL));
    memcpy( &stencil[iread].p[6],  &curoct->GDATA_p, 8*sizeof(REAL));

    //start recursive fill
    int inei;
    for(inei=0;inei<6;inei++){
      recursive_neighbor_gather_oct_grav(inei, curoct->parent, &stencil[iread] );
    }
  }
  (*nread)=iread;
} 

void scatterstencilgrav(struct OCT *** octList, int iOct, struct OCT *octstart, struct STENGRAV *gstencil,int nread,int stride, struct CPUINFO *cpu,int level){

  int iread;
  for(iread=0; iread<stride; iread++){
		  struct OCT  *curoct = octList[level-1][iOct++];
      if(curoct==NULL) break;

      memcpy(&curoct->GDATA_pnew,&gstencil->pnew[iread*8], 8*sizeof(REAL));
      memcpy(&curoct->GDATA_res ,&gstencil->res [iread*8], 8*sizeof(REAL));
 

#ifdef ONFLYRED
      if(level<=cpu->levelcoarse){
				cell2oct(curoct->parent)->GDATA_d[curoct->parent->idx] = gstencil->resLR[iread];
				cell2oct(curoct->parent)->GDATA_p[curoct->parent->idx] = 0;
      }
#endif  

  }

  if(iread!=(nread)) abort();
}

//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================

void scatterstencilforce(struct OCT ** octList, int iOct, struct OCT *octstart, struct STENGRAV *gstencil, int nread, int stride, struct CPUINFO *cpu,int dir)
{
  struct OCT* curoct;
  int ipos,iread=0;
  int icell;
  int ilast;
  int vecpos;

#ifndef FASTGRAV
  ilast=nread;
#else
  ilast=stride;
#endif

  //printf("let's scatter\n");
  for(iread=0; iread<stride; iread++){
      curoct = octList[iOct++];
      if(curoct==NULL) break;

#ifdef FASTGRAV
      if(gstencil->valid[vecpos]!=1)	continue;
#endif

      // filling the values in the central oct
      for(icell=0;icell<8;icell++){
#ifndef FASTGRAV
	curoct->cell[icell].f[dir]=gstencil->pnew[icell+iread*8];
#else
	// NOTE FASTGRAV uses a different memory layout
	curoct->cell[icell].f[dir]=gstencil->pnew[icell*stride+vecpos];
#endif
      }
  }
}




//============================================================================
int PoissonJacobi_single(struct STENGRAV *gstencil, int level, int curcpu, int nread,int stride,REAL dx, int flag, REAL factdens){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  REAL temp;
  REAL res;
  REAL d, p;

  int inei,icell;
  int i;
  int *vnei, *vcell;
  int ioct[7]={0,1,2,3,4,5,6};
  int ilast;
  int ival=0;



#ifndef FASTGRAV
  struct GGRID *stencil=gstencil->stencil;
  ilast=nread;
#else
  ilast=stride;
#endif

      for(icell=0;icell<8;icell++){ // we scan the cells
        const int *vnei  = &NEIG_NEIP[icell*6];
        const int *vcell = &NEIG_CELL[icell*6];

//	    int vnei[6],vcell[6];   
  //  getcellnei(icell, vnei, vcell); // we get the neighbors

	//printf("Start openmp on thread %d on mpirank =%d icell=%d ival=%d nread=%d\n",omp_get_thread_num(),curcpu,icell,ival,nread);
	for(i=0;i<ilast;i++){ // we scan the octs

#ifdef FASTGRAV
      if(gstencil->valid[i]!=1){
	continue;
      }
#endif

#ifdef ONFLYRED
      if(icell==0) gstencil->resLR[i]=0.;
#endif
      
      temp=0.;
      res=0.;

#ifndef FASTGRAV      
        d = stencil[i].d[6][icell];
        p = stencil[i].p[6][icell];
#else
        d = gstencil->d[i+icell*stride];
        p = gstencil->p[i+icell*stride];
#endif      
      // Computing the laplacian ===========================
 
      for(inei=0;inei<6;inei++){
#ifndef FASTGRAV
 	temp+=stencil[i].p[ioct[vnei[inei]]][vcell[inei]];
#else
	int idxnei=gstencil->nei[i+vnei[inei]*stride];
	temp+=gstencil->p[idxnei+vcell[inei]*stride];
#endif
      }
 
      // setting up the residual
      res=temp;
      
      // we finish the laplacian
      temp=temp/6.0;
      temp=temp-dx*dx*d/6.0*factdens;
      
      //if(curcell->d>0) printf("temp=%e den=%e pot=%e\n",temp,curcell->d,curcell->p);

      // we finsih the residual
      res=res-6.0*p;
      res=res/(dx*dx)-factdens*d;

      // we store the new value of the potential
#ifndef FASTGRAV
      gstencil->pnew[icell+i*8]=temp;
#else
      gstencil->pnew[i+icell*stride]=temp;
#endif
      // we store the local residual
      if(flag) {
#ifndef FASTGRAV
	gstencil->res[i*8+icell]=factdens*d;
#else
	gstencil->res[i+icell*stride]=factdens*d;
#endif
     }
      else{
#ifndef FASTGRAV
	gstencil->res[icell+i*8]=res;
#else
	gstencil->res[icell*stride+i]=res;
#endif
      }

#ifdef ONFLYRED
      // the low resolution residual
      gstencil->resLR[i]+=res*0.125;
#endif
      // ready for the next cell
      if(icell==0) ival++;
      //printf("i=%d on thread %d on mpirank =%d \n",i,omp_get_thread_num(),curcpu);
      
      }
    
    //ready for the next oct
  }
  //printf("ival=%d nread=%d\n",ival,nread);

  return 0;
}


//============================================================================




 //============================================================================
REAL comp_residual(struct STENGRAV *gstencil, int level, int curcpu, int nread,int stride, int flag){

  // flag=1 means the residual contains the norm of the density
  // flag=0 means the resiual contains the actual residual of the Poisson Equation

  int icell;
  int i;
  int iread=0;

  REAL residual=0.;
  REAL rloc;
  int ilast;

  REAL powtmp ;

  ilast=nread;

  if(flag){
    for(i=0;i<ilast;i++){ // we scan the octs

#ifdef FASTGRAV
      if(gstencil->valid[i]!=1){
	continue;
      }
#endif


      for(icell=0;icell<8;icell++){ // we scan the cells
#ifndef FASTGRAV
	powtmp = gstencil->res[icell+i*8];      rloc = powtmp*powtmp;
#else
	powtmp = gstencil->res[icell*stride+i]; rloc = powtmp*powtmp;
#endif
	residual+=rloc;
	// ready for the next cell
      }
      //ready for the next oct
      iread++;
    }
  }
  else{
    for(i=0;i<ilast;i++){ // we scan the octs

#ifdef FASTGRAV
      if(gstencil->valid[i]!=1){
	continue;
      }
#endif
      for(icell=0;icell<8;icell++){ // we scan the cells
#ifndef FASTGRAV
	powtmp = gstencil->res[icell+i*8];rloc = powtmp*powtmp;
#else
	powtmp = gstencil->res[icell*stride+i];rloc = powtmp*powtmp;
#endif
	residual+=rloc;
	//residual=(residual>rloc?residual:rloc);
      }
      //ready for the next oct
      iread++;
    }
  }

  if(iread!=nread) abort();
  return residual;
}

//============================================================================
int Pot2Force(struct STENGRAV *gstencil, int level, int curcpu, int nread,int stride,REAL dx, REAL tsim, int dir){

  int inei,icell;
  int i;

  //  int ioct[7]={12,14,10,16,4,22,13};
  int ioct[7]={0,1,2,3,4,5,6};
  REAL floc;

  int ilast;

#ifndef FASTGRAV
  ilast=nread;
  struct GGRID *stencil=gstencil->stencil;
#else
  int idxR,idxL;
  int cR,cL;
  ilast=stride;
#endif



  for(icell=0;icell<8;icell++){ // we scan the cells

    const int *vnei  = &NEIG_NEIP[icell*6];
    const int *vcell = &NEIG_CELL[icell*6];

   // int vnei[6],vcell[6];   
   // getcellnei(icell, vnei, vcell); // we get the neighbors

    for(i=0;i<ilast;i++){ // we scan the octs
      
#ifndef FASTGRAV
//      floc=0.5*(stencil[i].oct[ioct[vnei[(dir<<1)+1]]].cell[vcell[(dir<<1)+1]].gdata.p - stencil[i].oct[ioct[vnei[(dir<<1)]]].cell[vcell[(dir<<1)]].gdata.p)/dx*tsim;

      floc=0.5*( stencil[i].p[ioct[vnei[(dir<<1)+1]]][vcell[(dir<<1)+1]] - stencil[i].p[ioct[vnei[(dir<<1)]]][vcell[(dir<<1)]] )/dx*tsim;

#else
      if(gstencil->valid[i]!=1){
	continue;
      }

      idxR=gstencil->nei[i+vnei[(dir<<1)+1]*stride];
      idxL=gstencil->nei[i+vnei[(dir<<1)  ]*stride];

      cR=vcell[(dir<<1)+1];
      cL=vcell[(dir<<1)  ];

      floc=0.5*(gstencil->p[idxR+cR*stride]-gstencil->p[idxL+cL*stride])/dx*tsim;
#endif

      // store the force
      
      //gstencil->pnew[icell*nread+i]=floc;;
#ifndef FASTGRAV
      gstencil->pnew[icell+i*8]=floc;;
#else
      gstencil->pnew[icell*stride+i]=floc;;
#endif

      // ready for the next cell
    }
    //ready for the next oct
  }
  return 0;
}


//============================================================================
void clean_vecpos(int level,struct OCT **octList, struct CPUINFO *cpu){
  int iOct;
  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      //octList[iOct]->vecpos=-1;
  }
}



//=============================================================================

REAL PoissonJacobi(struct OCT *** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim, struct GRID *grid)
{

  struct OCT *curoct;

  int iter, iOct, icell;
  int nread, nreadtot;
  int nitmax;
  int crit;

  REAL fnorm,residual,residualold,dres;
  REAL dxcur = pow(0.5,level);
  REAL factdens;
  REAL rloc;
  REAL res0=0.;
  REAL dist,normp,dresconv;

  double tall=0.,tcal=0.,tscat=0.,tgat=0., tglob=0.,tup=0.;
  double tstart,tstop,tt;
  double temps[10];

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

  if((level==param->mgridlmin)||(level>param->lcoarse))   	
    nitmax=param->niter;
  else
    nitmax=param->nrelax;
 
  for(iter=0;iter<nitmax;iter++){
    tstart=MPI_Wtime();

    // --------------- setting the first oct of the level

    nreadtot=0;

    // --------------- some inits for iterative solver
    if(iter==0){
      fnorm=0.;
      residual=0.;
      nread=0;
    }
    else{
      residual=0.;
    }

    // Scanning the octs

  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct+=stride){
    curoct=octList[level-1][iOct]; 
		
		temps[0]=MPI_Wtime();

		// ------------ gathering the stencil value 
		gatherstencilgrav(octList, iOct,stencil->stencil,stride,cpu, &nread, level);
	//  gatherstencilgrav_grid(grid, octList, level, 21, cpu,param);
	//  exit(0);
	
		temps[2]=MPI_Wtime();
		// ------------ solving the hydro
		PoissonJacobi_single(stencil,level,cpu->rank,nread,stride,dxcur,(iter==0),factdens);

		// ------------ computing the residuals

		rloc=comp_residual(stencil,level,cpu->rank,nread,stride,(iter==0));

		if(iter==0){
			fnorm+=rloc;
			//printf("rloc=%e\n",rloc);
		}
		else{
			//residual=(residual>rloc?residual:rloc);
			residual+=rloc;
			//printf("iter=%d rloc=%e\n",iter,rloc);
		}
	
		// ------------ scatter back the data in the tree
		temps[4]=MPI_Wtime();

		scatterstencilgrav(octList, iOct, curoct,stencil,nread,stride, cpu, level);

		temps[6]=MPI_Wtime();

		tall+=temps[6]-temps[0];
		tcal+=temps[4]-temps[2];
		tscat+=temps[6]-temps[4];
		tgat+=temps[2]-temps[0];
		nreadtot+=nread;
	}



    //printf("iter=%d nreadtot=%d\n",iter,nreadtot);


    tt=MPI_Wtime();
    // at this stage an iteration has been completed : let's update the potential in the tree

      if(nreadtot>0){
	update_pot_in_tree(octList[level-1], level,firstoct,cpu,param,&dist,&normp);
      }

#ifdef WMPI
    //printf("iter=%d\n",iter);

    // reduction of relevant quantities

    if((iter<=param->niter)||(iter%1==0)){
      mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,(iter==0),level); // potential field exchange
      if(iter==0){
	//if(level==7) printf("rank=%d fnorm=%e\n",cpu->rank,fnorm);
	MPI_Allreduce(MPI_IN_PLACE,&fnorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
      else{
	MPI_Allreduce(MPI_IN_PLACE,&residual,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE,&dist,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE,&normp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
    }
#endif

    if((iter==1)&&(level>=param->lcoarse)) res0=residual;

    tstop=MPI_Wtime();
    tup+=(tstop-tt);
    tglob+=(tstop-tstart);
    
    if(iter>0){

      // here we test the convergence of the temporary solution
      dresconv=sqrt(dist/normp);

      // here we test the zero level of Poisson equation
      if(level<param->lcoarse){
	dres=sqrt(residual);
      }
      else{
	dres=sqrt(residual/fnorm);
      }

      // we take the smallest
      crit=(dres<dresconv?0:1);
      dres=(dres<dresconv?dres:dresconv);

      if((dres)<param->poissonacc){
	if(level>=param->lcoarse) break;
      }
    }
  }

  if(cpu->rank==0) printf("CPU | tgat=%e tcal=%e tscat=%e\n",tgat,tcal,tscat);
#ifdef FASTGRAV
  if(cpu->noct[level-1]==nread){
    scatterstencilgrav(octList, iOct, curoct,stencil,nread,stride, cpu, level);
    if(nreadtot>0) update_pot_in_tree(octList, level,firstoct,cpu,param);
  }
#endif

  if(level>param->lcoarse){
    if(cpu->rank==0) printf("CPU | level=%d iter=%d res=%e fnorm=%e\n",level,iter,dres,fnorm);
  }
  else{
    if(cpu->rank==0) printf("CPU | level=%d iter=%d res=%e fnorm=%e resraw=%e res0=%e crit=%d\n",level,iter,dres,fnorm,sqrt(residual),sqrt(res0),crit);
  }
  return dres;
}



//=============================================================================




//===============================================================================================

REAL PoissonMgrid(struct OCT *** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim, struct GRID *grid)
{
  struct OCT *curoct;
  struct Gtype Wi[8];

  int icell, iOct;
  REAL dres;


  // pre-relaxation
#ifdef WMPI
  mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,1,level); // potential field exchange
#endif

#ifndef GPUAXL
  dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
#else
  if(level>=7){
    //dres=PoissonJacobiGPU(level,param,firstoct,cpu,stencil,stride,tsim);
    dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim, grid);
  }
  else{
    dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
  }
  //dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim);
#endif


#ifndef ONFLYRED
  // NOTE ON GPU the calculation is performed on the fly
  // reduction

  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    curoct=octList[level-1][iOct]; 

    curoct->parent->gdata.d=0.;
    curoct->parent->gdata.p=0.;

    for(icell=0;icell<8;icell++){
      cell2oct(curoct->parent)->GDATA_p[curoct->parent->idx] += curoct->res[icell]*0.125; // we average the residual and store it as the new density
//      curoct->parent->gdata.d += curoct->res[icell]*0.125; // we average the residual and store it as the new density
    }
  }
#endif

  // full relaxation at coarsest level or recursive call to mgrid
  
  if((level-1)==param->mgridlmin){
#ifdef WMPI
    mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,1,level-1); // potential field exchange
#endif

#ifndef GPUAXL
    PoissonJacobi(octList, level-1,param,firstoct,cpu,stencil,stride,tsim,grid);
#else
    //PoissonJacobiGPU(level-1,param,firstoct,cpu,stencil,stride,tsim);
    PoissonJacobi(octList, level-1,param,firstoct,cpu,stencil,stride,tsim,grid);
#endif
  }
  else{
    PoissonMgrid(octList, level-1,param,firstoct,cpu,stencil,stride,tsim, grid);
  }
  
  // prolongation + correction
 
  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){

    curoct=octList[level-1][iOct]; 

//    REAL tmp[8];
//    coarse2fine_gravlin(curoct->parent, curoct->GDATA_d, tmp );

    for(icell=0;icell<8;icell++){ // looping over cells in oct 
        curoct->GDATA_p[icell] -=  cell2oct(curoct->parent)->GDATA_p[curoct->parent->idx] ;

//        curoct->GDATA_p[icell] -= tmp[icell];// we propagate the error and correct the evaluation
    }
  }


  // post relaxation
#ifdef WMPI		
  mpi_exchange_level(cpu,cpu->sendbuffer,cpu->recvbuffer,2,1,level); // potential field exchange
#endif

#ifndef GPUAXL
  dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
#else
  //dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim);
  if(level>=7){
    //dres=PoissonJacobiGPU(level,param,firstoct,cpu,stencil,stride,tsim);
    dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
  }
  else{
    dres=PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,tsim,grid);
  }
#endif
  return dres;
}

//===================================================================================================================================

void PoissonForce(struct OCT *** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim){

  struct OCT *curoct;

  REAL dxcur;
  int icomp,icell, iOct;
  struct PART *curp;
  struct PART *nexp;
  int nread,nreadtot=0;
  int idir;
  dxcur=pow(0.5,level);
  
  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct+=stride){
      curoct = octList[level-1][iOct];
//      if(curoct==NULL) break;


#ifndef FASTGRAV
	// ------------ gathering the stencil value 
      gatherstencilgrav(octList, iOct,stencil->stencil,stride,cpu, &nread,level);
#else
      // ------------ some cleaning
      clean_vecpos(level,firstoct);
      
      //------------- gathering neighbour structure
      gatherstencilgrav_nei(curoct,stencil,stride,cpu, &nread);
      
      // ------------ gathering the stencil value 
      gatherstencilgrav(octList, iOct, firstoct[level-1],stencil,stride,cpu, &nread,level); //  the whole level must be parsed to recover the values

#endif
      
      
      for(idir=0;idir<3;idir++){
	// ------------ Computing the potential gradient
	Pot2Force(stencil,level,cpu->rank,nread,stride,dxcur,tsim,idir);
	
	// ------------ scatter back the data
	
	scatterstencilforce(octList[level-1], iOct, curoct,stencil, nread, stride, cpu,idir);
      }
      
      nreadtot+=nread;

  }

#ifdef WMPI
    mpi_exchange(cpu,cpu->sendbuffer,cpu->recvbuffer,5,1); // fx field
    mpi_exchange(cpu,cpu->sendbuffer,cpu->recvbuffer,6,1); // fy field
    mpi_exchange(cpu,cpu->sendbuffer,cpu->recvbuffer,7,1); // fz field
#endif
}

#endif

#ifdef WGRAV

//===================================================================================================================================
//===================================================================================================================================
int FillDens(struct OCT ** octList, int level,struct RUNPARAMS *param, struct CPUINFO *cpu){

  struct OCT *curoct;
  struct PART *curp;
  struct PART *nexp; 

  int icomp,icell,iOct;
  int nread;
  int nc=0;

  REAL dx;
  REAL locdens;
  REAL avgdens=0.;



  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
      curoct = octList[iOct];
      if(curoct==NULL) break;

      for(icell=0;icell<8;icell++){	
	locdens=0.;
#ifdef PIC
	locdens+=curoct->cell[icell].density;
#endif

#ifdef WHYDRO2
	locdens+=curoct->cell[icell].field.d;
#endif
	curoct->cell[icell].gdata.d=locdens;


	// WARNING HERE ASSUME PERIODIC BOUNDARY CONDITIONS AND REMOVE THE AVG DENSITY
	curoct->cell[icell].gdata.d-=1.;
	

	avgdens+=locdens;
	nc++;
      }
  }

  avgdens/=nc;

  


  return 0;
}

#endif

//===================================================================================================================================
//===================================================================================================================================
#ifdef WGRAV
int PoissonSolver(struct OCT *** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, struct GRID *gstencil_grid, int stride, REAL aexp){

  REAL res;
  int igrid;
  struct Gtype Wi[8];
  struct CELL* curcell;
  int icell2;
  struct OCT* curoct;

  int icell, iOct;
  double t[10];

  MPI_Barrier(cpu->comm);
 t[0]=MPI_Wtime();
  if(cpu->rank==0) printf("Start Poisson Solver ");

#ifndef GPUAXL
  if(cpu->rank==0)  printf("on CPU\n");
#else
  if(cpu->rank==0)  printf("on GPU\n");
#endif
  //breakmpi();

//=========== changement de structure de données==============
  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    curoct=octList[level-1][iOct]; 
    for(icell=0;icell<8;icell++){
      curoct->GDATA_d[icell] = curoct->cell[icell].gdata.d;
      curoct->GDATA_p[icell] = curoct->cell[icell].gdata.p;
    }
  }

  if((level==param->lcoarse)&&(param->lcoarse!=param->mgridlmin)){
    for(igrid=0;igrid<param->nvcycles;igrid++){ // V-Cycles
      if(cpu->rank==0) printf("----------------------------------------\n");
  //   res=PoissonMgrid(octList, level,param,firstoct,cpu,stencil,stride,aexp, gstencil_grid);
      res=PoissonMgrid_grid(octList, level,param,firstoct,cpu,stencil,stride,aexp, gstencil_grid);
      if(res<param->poissonacc) break;
    }
  }
  else{
#ifndef GPUAXL
    PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,aexp,gstencil_grid);
#else
    PoissonJacobi(octList, level,param,firstoct,cpu,stencil,stride,aexp,gstencil_grid);
    //PoissonJacobiGPU(level,param,firstoct,cpu,stencil,stride,aexp);
#endif

  }
	
  //once done we propagate the solution to level+1

  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    curoct=octList[level-1][iOct]; 
    for(icell=0;icell<8;icell++){ // looping over cells in oct   
      curcell=&(curoct->cell[icell]);
      if(curcell->child!=NULL){    
  
	coarse2fine_gravlin(curcell,curcell->child->GDATA_d,curcell->child->GDATA_p);
/*
	for(icell2=0;icell2<8;icell2++){
	//Wi[icell2].p=0.;
	  memcpy(&(curcell->child->cell[icell2].gdata.p),&(Wi[icell2].p),sizeof(REAL));
	//memcpy(&(curcell->child->cell[icell2].gdata.p),&(curcell->gdata.p),sizeof(REAL));
	}
*/
      }
    }
  }
  

//=========== retour a la structure de données standard
  for(iOct=0; iOct<cpu->locNoct[level-1]; iOct++){
    curoct=octList[level-1][iOct]; 
    for(icell=0;icell<8;icell++){
      curoct->cell[icell].gdata.d = curoct->GDATA_d[icell];
      curoct->cell[icell].gdata.p = curoct->GDATA_p[icell];
    }
  }



  MPI_Barrier(cpu->comm);
 t[9]=MPI_Wtime();
 if(cpu->rank==0){  
#ifndef GPUAXL
   printf("==== CPU POISSON TOTAL TIME =%e\n",t[9]-t[0]);
#else
   printf(" === GPU POISSON TOTAL TIME =%e\n",t[9]-t[0]);
#endif
 }
  return 0;
}
#endif
