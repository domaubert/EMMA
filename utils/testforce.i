#include "../utils/readamr.i"
#include "../utils/readpart.i"
  
p=readpart("data/part.00000.p00000");
fx=p(4,);
fy=p(5,);
fz=p(6,);
f=abs(fx,fy,fz);

x=p(1,)-0.5;
y=p(2,)-0.5;
z=p(3,)-0.5;
r=abs(x,y,z);
er=transpose([x,y,z]/(r+1e-8));
fr=er(1,)*fx+er(2,)*fy+er(3,)*fz;
ft=sqrt(f^2-fr^2);


