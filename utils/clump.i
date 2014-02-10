
//rep="data_3g/";execut="./utils/oct2grid ";
//rep="data1g/";execut="./utils/oct2grid1 ";
rep="data/";execut="./utils/oct2grid ";

imax=2095;
imin=1;
lvl=7;
dxcur=1./pow(2,lvl+1);
r=span(0.,1.,pow(2,lvl)+1)(zcen);
n=pow(2,lvl);
xx=r(,-:1:n)(,,-:1:n);
y=r(-:1:n,)(,,-:1:n);
zz=r(-:1:n,,)(-:1:n,,);
RR=abs(x,y,zz);
imid=pow(2,lvl-1);
wp=where(r>=0.);
rion=vion=tion=[];
binr=span(r(dif)(1),1.,256);
// logxy,0,1;
// limits,0,1.05,1e-4,4e-3;
for(i=imin;i<=imax;i+=50){
  x=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,706,a,ncpu=64,execut=execut);
  // y=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,103,a,ncpu=32,execut="./utils/oct2grid ");
  // z=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,104,a,ncpu=32,execut="./utils/oct2grid ");
  // x=sqrt(x^2+y^2+z^2);
  //hx=histo1d(RR,binr,wght=x)/(histo1d(RR,binr)+1e-6);grow,rion,interp(binr(zcen)(15:),hx(15:),0.5);
  //hx=histo1d(RR,binr,wght=x)/(4*pi*binr(zcen)^2*binr(dif))*(r(dif)(1))^3;grow,rion,interp(binr(zcen)(10:),hx(10:),0.5);
  //grow,rion,r(where((x(,imid,imid)<0.5)*(r>0.))(min));
  grow,rion,interp(r(:110),x(:110,imid,imid),0.5);
  //  grow,vion,numberof(where(x>0.01))/(numberof(x)*1.0);
  grow,tion,a;
  //plg,hx*1e-6,binr(zcen),color="red";
 }

// PL,rion,tion,color="red",msize=.5;
// plg,rth,tion;
