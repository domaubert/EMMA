require,"utils/readamr.i";

i=20;lvl=7;
rep="data/";
dxcur=1./32.;

e=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,701);
 fx=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,702,a);
xion=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,706,a);
l=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,0,a);
// t=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,707,a);
// s=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,705,a);
d=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,708,a);
// P=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,105,a);
d2=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,101,a);

//X=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,108,a);

r=span(-0.5,0.5,pow(2,lvl)+1)(zcen)-dxcur;
n=pow(2,lvl);
x=r(,-:1:n)(,,-:1:n);
y=r(-:1:n,)(,,-:1:n);
zz=r(-:1:n,,)(-:1:n,,);
RR=abs(x,y,zz);
binr=span(r(dif)(1),0.7,256);
hx=histo1d(RR(*),binr,wght=xion(*))/(histo1d(RR,binr)+1e-6);
hd=histo1d(RR(*),binr,wght=d(*))/(histo1d(RR,binr)+1e-6);

