
#include "utils/read_grafic.i"

p=mergepart("data/part.00075",32,a);
L=100.;
h=0.67;
d=part2cic(p,256);
da=avg(d(*));
d=(d-da)/da;
//d=read_grafic("./utils/mpgrafic-0.2/ic_deltab",L,h);
//d=d(::2,::2,::2);
L/=h;



write,"==> Preparing FFT";
nx=dimsof(d)(0);
kf=fft_indgen(nx);
ka=array(0.,nx,nx,nx);



kx=kf(,-:1:nx)(,,-:1:nx);
ky=kf(-:1:nx,)(,,-:1:nx);
kz=kf(-:1:nx,)(-:1:nx,,);
ka=abs(kx,ky,kz)/L*2*pi;

write,"==> FFT dens";
nps=128;
fd=fft(d,[1,1,1]);
fd=double(fd*conj(fd))/nx^3;

write,"==> Computing PS";
kv=spanl(1e-3,1e3,nps+1)(zcen);
nk=(histo1d(ka(*),kv)+1e-10);
// P=histo1d(ka(*),kv,wght=fd(*))/(4*pi*kv(zcen)^2*kv(dif)/(2*pi/L)^3);
// P2=histo1d(ka(*),kv,wght=fd(*)^2)/(4*pi*kv(zcen)^2*kv(dif)/(2*pi/L)^3);
P=histo1d(ka(*),kv,wght=fd(*))/nk;
P2=histo1d(ka(*),kv,wght=fd(*)^2)/nk;
S=sqrt(P2-P^2);
EP=S/sqrt((4*pi*kv(zcen)^2*kv(dif)/(2*pi/L)^3));

//EP=EP*(L/nx)^3;
//P=P*(L/nx)^3; // Uncomment to comply with powspec4


// grafic data
vv=load("./utils/mpgrafic-0.2/power.dat");
dp=load("./utils/mpgrafic-0.2/dplus.dat");
dpt=load("./utils/mpgrafic-0.2/dplustab.dat");

gd=interp(dpt(,2),dpt(,1),a);

//ws;
PE=[P-3*EP,P+3*EP];
plshade,PE,kv(zcen),color=__rgb(,10);
plg,vv(,2)/h^3*(2*pi)^3*dp(1,3),vv(,1)*h,color="blue";
plg,vv(,2)/h^3*(2*pi)^3*gd^2,vv(,1)*h,color="green";
plg,P,kv(zcen),color="black",width=4;
// // powspec data
// pth=load("~/Project/BAO/test/powspec.dat");
// plg,pth(2:,3),pth(2:,1),color="green";

logxy,1,1;
