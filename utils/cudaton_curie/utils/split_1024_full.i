#include "cudaton.i"

ll=exec("ls -d fields/snapaton*");

ncpux=4;
ncpuy=4;
ncpuz=4;
nbnd=16;

for(i=1;i<=numberof(ll);i++){

  
  
  //s=readsnap(ll(i),light=1);
  s=read1024(ll(i))
  ncx=dimsof(*s.dens)(2)/ncpux;
  ncy=dimsof(*s.dens)(3)/ncpuy;
  ncz=dimsof(*s.dens)(4)/ncpuz;

  // first we write the header and the source
  
  for(ic=0;ic<ncpux;ic++){
  for(jc=0;jc<ncpuy;jc++){
  for(kc=0;kc<ncpuz;kc++){
    rank=ic+jc*ncpux+kc*ncpux*ncpuy;
    
    newf=swrite(format=ll(i)+".p%05d",rank);
    newf;

    ff=open(newf,"wb");
    adress=0;
    ncx_bnd=int(ncx+2*nbnd);
    _write,ff,adress,ncx_bnd;adress+=sizeof(ncx_bnd);
    _write,ff,adress,s.nsource;adress+=sizeof(s.nsource);
    _write,ff,adress,s.time;adress+=sizeof(s.time);

    to=float((*s.temp)(ic*ncx+1:(ic+1)*ncx-1+1,jc*ncy+1:(jc+1)*ncy-1+1,kc*ncz+1:(kc+1)*ncz-1+1));
    tobnd=array(float,ncx+2*nbnd,ncy+2*nbnd,ncz+2*nbnd);
    tobnd(nbnd+1:nbnd+ncx,nbnd+1:nbnd+ncy,nbnd+1:nbnd+ncz)=to;

    so=float((*s.src)(ic*ncx+1:(ic+1)*ncx-1+1,jc*ncy+1:(jc+1)*ncy-1+1,kc*ncz+1:(kc+1)*ncz-1+1));
    sobnd=array(float,ncx+2*nbnd,ncy+2*nbnd,ncz+2*nbnd);
    sobnd(nbnd+1:nbnd+ncx,nbnd+1:nbnd+ncy,nbnd+1:nbnd+ncz)=so;

    ro=float((*s.dens)(ic*ncx+1:(ic+1)*ncx-1+1,jc*ncy+1:(jc+1)*ncy-1+1,kc*ncz+1:(kc+1)*ncz-1+1));
    robnd=array(float,ncx+2*nbnd,ncy+2*nbnd,ncz+2*nbnd);
    robnd(nbnd+1:nbnd+ncx,nbnd+1:nbnd+ncy,nbnd+1:nbnd+ncz)=ro;

    zerotab=float(robnd*0.);

    _write,ff,adress,zerotab;adress+=sizeof(zerotab);
    _write,ff,adress,zerotab;adress+=sizeof(zerotab);
    _write,ff,adress,zerotab;adress+=sizeof(zerotab);
    _write,ff,adress,zerotab;adress+=sizeof(zerotab);
    _write,ff,adress,zerotab;adress+=sizeof(zerotab);

    _write,ff,adress,tobnd;adress+=sizeof(tobnd);
    _write,ff,adress,sobnd;adress+=sizeof(sobnd);
    _write,ff,adress,robnd;adress+=sizeof(robnd);
    _write,ff,adress,s.aexp;adress+=sizeof(s.aexp);

    close,ff;

  }
  }
  }
  
 }