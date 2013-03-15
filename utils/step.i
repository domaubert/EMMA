

rep="data/";
imax=50;
lvl=6;
r=span(-0.5,0.5,pow(2,lvl)+1)(zcen);
// n=pow(2,lvl);
// x=r(,-:1:n)(,,-:1:n);
// y=r(-:1:n,)(,,-:1:n);
// zz=r(-:1:n,,)(-:1:n,,);
// RR=abs(x,y,zz);
// imid=pow(2,lvl-1);
// wp=where(r>=0.);
// rion=vion=tion=[];
// binr=span(r(dif)(1),0.7,256)

h0=73.;
om=0.27;
ov=0.73;


zr=array(float,imax+1,pow(2,lvl));
tr=array(float,imax+1,pow(2,lvl));
xr=array(float,imax+1,pow(2,lvl));
rr=r(-:1:imax+1,);
for(i=0;i<=imax;i+=1){
  
  x=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,706,a);
  xr(i+1,)=x(,16,16);
  zr(i+1,)=1./a-1.;
  rr(i+1,)=r*a;
  tr(i+1,)=(univAge(10000.,h0=h0,Omega_l=ov,Omega_m=om,silent=1)-univAge(1./a-1.,h0=h0,Omega_l=ov,Omega_m=om,silent=1))/(365*24*3600.*1e6);

  //hx=histo1d(RR,binr,wght=x)/(histo1d(RR,binr)+1e-6);grow,rion,interp(binr(zcen)(10:),hx(10:),0.5);
  //hx=histo1d(RR,binr,wght=x)/(4*pi*binr(zcen)^2*binr(dif))*(r(dif)(1))^3;grow,rion,interp(binr(zcen)(10:),hx(10:),0.5);
  //grow,rion,r(where((x(,imid,imid)<0.5)*(r>0.))(min));
  //grow,rion,interp(r(imid+1:),x(imid+1:,imid,imid),0.5);
  //grow,vion,numberof(where(x>0.01))/32.^3;
  //grow,tion,a;

 }

