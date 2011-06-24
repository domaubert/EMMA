#include "Chris/jgraph.i"

p=l=nd=[];
t=[];

for(i=1;i<=1080;i+=1){ 
  //  pname=swrite(format="data/l7_6/part.%05d",i);
  pname=swrite(format="../data/part.%05d",i);
  //  pp=readpart(pname,a);
  pp=mergepart(pname,2,a);
  s=sort(pp(7,));
  grow,p,[pp(,s)];
  grow,t,a;
  // ll=readcube(swrite(format="../data/lev3d.%05d.p00000",i));
  // grow,l,[ll(,,64)];
  // d=readcube(swrite(format="data/den3d.%05d.p00000",i));
  // ad=avg(d);
  // normd=4*pi*sqrt(((d-ad)^2)(*)(sum));
  // grow,nd,normd;
 }

epsilon=0.;
etot=0.5*(1.-epsilon)*abs(p(4,1,),p(5,1,),p(6,1,))^2+0.5*(epsilon)*abs(p(4,2,),p(5,2,),p(6,2,))^2-epsilon*(1-epsilon)/abs(p(1,1,)-p(1,2,),p(2,1,)-p(2,2,),p(3,1,)-p(3,2,));


//ws;
//JDrawCircle(0.5,0.5,0.2,number=100);
//plg,p(2,2,),p(1,2,),color="red";


r=abs(p(1,2,)-0.5,p(2,2,)-0.5,p(3,2,)-0.5);
r2=abs(p(1,2,)-p(1,1,),p(2,2,)-p(2,1,),p(3,2,)-p(3,1,));
v=abs(p(4,2,)-p(4,1,),p(5,2,)-p(5,1,),p(6,2,)-p(6,1,));
er=(p(1:3,2,)-p(1:3,1,))/(r(-:1:3,)+1e-6);
vr=((p(4:6,2,)-p(4:6,1,))*er)(sum,);
vt=sqrt(v^2-vr^2);

e=0.5*v^2-1./r2;
write,"deltae=",e(rms)/e(avg);
write,"deltae=",etot(rms)/etot(avg);
//ws;
//plg,p(2,2,),p(1,2,),color="red";
//JDrawCircle,0.5,0.5,0.1,number=50;
//plg,r;
//plg,r2,color="blue";
//plg,v,color="blue";
