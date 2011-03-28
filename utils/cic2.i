



p=readpart("../data/part.00315.p00000");
n=256;
dxn=1./n;

dens=array(double,n,n,n);

ic=int(p(1,)/dxn);
jc=int(p(2,)/dxn);
kc=int(p(3,)/dxn);

for(q=1;q<=dimsof(p)(0);q++){
  
  xc=ic(q)*dxn+dxn*0.5;
  yc=jc(q)*dxn+dxn*0.5;
  zc=kc(q)*dxn+dxn*0.5;

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
