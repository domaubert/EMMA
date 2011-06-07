fx=readcube("data/fx.00003.p00000");
fy=readcube("data/fy.00003.p00000");
fz=readcube("data/fz.00003.p00000");
f=abs(fx,fy,fz);
r=span(-0.5,0.5,257)(zcen);
//plg,f(,128,128)*r^2,r,color="green";
