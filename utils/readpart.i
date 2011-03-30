

func readpart(fname){
  //fname=exec("ls -d data/part.*")
  //x=y=[];for(i=1;i<=numberof(fname);i++){phase=readpart(fname(i));www=where(phase(7,)==1);grow,x,phase(1,www);grow,y,phase(2,www);}
  fp=open(fname,"rb");
  adress=0;
  npart=array(int);
  _read,fp,adress,npart;adress+=sizeof(npart);
  if(npart!=0){
  phase=array(float,7,npart);
  _read,fp,adress,phase;adress+=sizeof(phase);
  write,"found "+pr1(npart)+" particles";
  close,fp;
  return phase;
  }
  else{
    write,"No particles : empty file!";
  return 0;
  }
}



func mergepart(fname,ncpu){
  pf=[];
  for(i=0;i<ncpu;i++){
    fname2=swrite(format=fname+".p%05d",i);
    p=readpart(fname2);
    if(numberof(p)!=1) grow,pf,p;
  }
  return pf;
}


func part2cic(p,n){
dxn=1./n;

dens=array(double,n,n,n);

ic=int(p(1,)/dxn);
jc=int(p(2,)/dxn);
kc=int(p(3,)/dxn);

for(q=1;q<=dimsof(p)(0);q++){
  
  xc=ic(q)*dxn+dxn*0.5;
  yc=jc(q)*dxn+dxn*0.5;
  zc=kc(q)*dxn+dxn*0.5;

  // for particles with coordinates equals to 1
  if(ic(q)>n-1) continue;
  if(jc(q)>n-1) continue;
  if(kc(q)>n-1) continue;
  
  tx=abs(p(1,q)-xc)/dxn;
  ty=abs(p(2,q)-yc)/dxn;
  tz=abs(p(3,q)-zc)/dxn;
  
  dx=1.-tx;
  dy=1.-ty;

  dz=1.-tz;

  
  mx=1-(p(1,q)<xc)*2;
  my=1-(p(2,q)<yc)*2;
  mz=1-(p(3,q)<zc)*2;
  
  
  dens(ic(q)+1,jc(q)+1,kc(q)+1)+=dx*dy*dz;

  dens((ic(q)+1+mx)%n,(jc(q)+1)%n,(kc(q)+1)%n)+=tx*dy*dz;


  dens((ic(q)+1)%n,(jc(q)+1+my)%n,(kc(q)+1)%n)+=dx*ty*dz;
  dens((ic(q)+1)%n,(jc(q)+1)%n,(kc(q)+1+mz)%n)+=dx*dy*tz;

  dens((ic(q)+1+mx)%n,(jc(q)+1+my)%n,(kc(q)+1)%n)+=tx*ty*dz;
  dens((ic(q)+1)%n,(jc(q)+1+my)%n,(kc(q)+1+mz)%n)+=dx*ty*tz;
  dens((ic(q)+1+mx)%n,(jc(q)+1)%n,(kc(q)+1+mz)%n)+=tx*dy*tz;

  dens((ic(q)+1+mx)%n,(jc(q)+1+my)%n,(kc(q)+1+mz)%n)+=tx*ty*tz;


 }

 return dens;
}
