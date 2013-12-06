#include "/home/daubert/Project/Quartz/utils/Zeldovich.i"
#include "/home/daubert/Project/Quartz/utils/readpart.i"

np=32;
x=span(0,1,np+1)(zcen);

pp=vp=[];
pt=vt=[];
aexp=[];

p=readpart("data/pstart.00000.p00000",a);
//p=mergepart(swrite(format="data/part.%05d",i),2,a);
p=p(,sort(p(7,)));
grow,pp,[p(1,:np)];
grow,vp,[p(4,:np)];
grow,aexp,a;
Zeldovich(a,0.5,xt,vxt,ng=np,omegam=0.999,omegav=0.001);
grow,pt,[xt];
grow,vt,[vxt];

for(i=0;i<=6;i+=1){
  i;
  p=readpart(swrite(format="./data/part.%05d.p00000",i),a);
  //p=mergepart(swrite(format="data/part.%05d",i),2,a);
  p=p(,sort(p(7,)));
  grow,pp,[p(1,:np)];
  grow,vp,[p(4,:np)];
  grow,aexp,a;
  Zeldovich(a,0.5,xt,vxt,ng=np,omegam=0.999,omegav=0.001);
  grow,pt,[xt];
  grow,vt,[vxt];
 }


deltax=sqrt(((pp-pt)^2)(sum,)/((pp-x)^2)(sum,));
deltav=sqrt(((vp-vt)^2)(sum,)/(vt^2)(sum,));
