
#include "utils/readamr.i"

p=readcube("data/pot3d.00200.p00000");
d=readcube("data/den3d.00200.p00000",a);
INFO,p;

//pli,p(,,64);
dx=1./dimsof(p)(2);

lap=(roll(p,[1,0,0])+roll(p,[0,1,0])+roll(p,[0,0,1])+roll(p,[-1,0,0])+roll(p,[0,-1,0])+roll(p,[0,0,-1])-6.*p)/(dx^2);

d=(d-1.)*1.5*0.3/a

//lap=lap*a/(0.3*1.5)+1.;

