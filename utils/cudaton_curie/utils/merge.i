
istart=77;
istop=istart;

ncpux=2;
ncpuy=2;
ncpuz=2;
nbnd=16;

xavg=xavgm=aexp=[];
egy=xion=temp=src=den=[];
fx=fy=fz=[];

for(isnap=istart;isnap<=istop;isnap++){
  rootname=(swrite(format="/ccc/scratch/cont005/gch0003/aubertd/runclues/out512/snapaton.%05d",isnap));


  
  for(ic=0;ic<ncpux;ic++){
  for(jc=0;jc<ncpuy;jc++){
  for(kc=0;kc<ncpuz;kc++){
    rank=ic+jc*ncpux+kc*ncpux*ncpuy;
    swrite(format=rootname+".p%05d",rank);
    s=readsnap(swrite(format=rootname+".p%05d",rank));

    if(rank==0){
      ncxo=(s.ncells-2*nbnd);
      ncyo=(s.ncells-2*nbnd);
      nczo=(s.ncells-2*nbnd);

      ncx=(s.ncells-2*nbnd)*ncpux;
      ncy=(s.ncells-2*nbnd)*ncpuy;
      ncz=(s.ncells-2*nbnd)*ncpuz;

       src=array(float,ncx,ncy,ncz);
       den=array(float,ncx,ncy,ncz);
       xion=array(float,ncx,ncy,ncz);
       temp=array(float,ncx,ncy,ncz);
      egy=array(float,ncx,ncy,ncz);
      //fx=array(float,ncx,ncy,ncz);
      //fy=array(float,ncx,ncy,ncz);
      //fz=array(float,ncx,ncy,ncz);
      
      //      zion=array(float,ncx,ncy,ncz)+100.;
    }

    egy(ic*ncxo+1:(ic+1)*ncxo-1+1,jc*ncyo+1:(jc+1)*ncyo-1+1,kc*nczo+1:(kc+1)*nczo-1+1)=(*s.egy);
    //fx(ic*ncxo+1:(ic+1)*ncxo-1+1,jc*ncyo+1:(jc+1)*ncyo-1+1,kc*nczo+1:(kc+1)*nczo-1+1)=(*s.flx)(,,1);
    //fy(ic*ncxo+1:(ic+1)*ncxo-1+1,jc*ncyo+1:(jc+1)*ncyo-1+1,kc*nczo+1:(kc+1)*nczo-1+1)=(*s.flx)(,,2);
    //fz(ic*ncxo+1:(ic+1)*ncxo-1+1,jc*ncyo+1:(jc+1)*ncyo-1+1,kc*nczo+1:(kc+1)*nczo-1+1)=(*s.flx)(,,3);
     xion(ic*ncxo+1:(ic+1)*ncxo-1+1,jc*ncyo+1:(jc+1)*ncyo-1+1,kc*nczo+1:(kc+1)*nczo-1+1)=(*s.xion);
     src(ic*ncxo+1:(ic+1)*ncxo-1+1,jc*ncyo+1:(jc+1)*ncyo-1+1,kc*nczo+1:(kc+1)*nczo-1+1)=(*s.src);
     den(ic*ncxo+1:(ic+1)*ncxo-1+1,jc*ncyo+1:(jc+1)*ncyo-1+1,kc*nczo+1:(kc+1)*nczo-1+1)=(*s.dens);
     temp(ic*ncxo+1:(ic+1)*ncxo-1+1,jc*ncyo+1:(jc+1)*ncyo-1+1,kc*nczo+1:(kc+1)*nczo-1+1)=(*s.temp);
    
     INFO,(*s.xion);
     INFO,(*s.temp);
     
  }
  }
  }


  // www=where(den>0.01);
  // grow,xavg,xion(www)(avg);
  // grow,xavgm,(xion(www)*den(www))(sum)/(den(www)(sum));
  // grow,aexp,s.aexp;

 }
