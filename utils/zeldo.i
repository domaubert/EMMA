

x=span(0,1,65)(zcen);

pp=vp=[];
pt=vt=[];
aexp=[];
for(i=0;i<=180;i+=10){
  i;
  p=readpart(swrite(format="data/part.%05d.p00000",i),a);
  p=p(,sort(p(7,)));
  grow,pp,[p(1,:64)];
  grow,vp,[p(4,:64)];
  grow,aexp,a;
  Zeldovich(a,1.0,xt,vxt,ng=64);
  grow,pt,[xt];
  grow,vt,[vxt];
 }


deltax=sqrt(((pp-pt)^2)(sum,)/((pp-x)^2)(sum,));
deltav=sqrt(((vp-vt)^2)(sum,)/(vt^2)(sum,));
