
rep="data/";
imax=157;
lvl=7;
dxcur=1./pow(2,lvl+1);
r=span(-0.5,0.5,pow(2,lvl)+1)(zcen)-dxcur;
n=pow(2,lvl);
x=r(,-:1:n)(,,-:1:n);
y=r(-:1:n,)(,,-:1:n);
zz=r(-:1:n,,)(-:1:n,,);
RR=abs(x,y,zz);
imid=pow(2,lvl-1);
wp=where(r>=0.);
rion=vion=tion=[];
binr=span(r(dif)(1),0.7,256)
for(i=0;i<=imax;i+=5){
  
  x=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,706,a);
  hx=histo1d(RR,binr,wght=x)/(histo1d(RR,binr)+1e-6);grow,rion,interp(binr(zcen)(10:),hx(10:),0.1);
  //hx=histo1d(RR,binr,wght=x)/(4*pi*binr(zcen)^2*binr(dif))*(r(dif)(1))^3;grow,rion,interp(binr(zcen)(10:),hx(10:),0.5);
  //grow,rion,r(where((x(,imid,imid)<0.5)*(r>0.))(min));
  //grow,rion,interp(r(imid+1:),x(imid+1:,imid,imid),0.5);
  grow,vion,numberof(where(x>0.01))/(numberof(x)*1.0);
  grow,tion,a;
 }

zion=1./tion-1.;
