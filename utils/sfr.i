
#include "Bastien/cosmo.i"
#include "utils/readpart.i"

rep="./data";ncpu=16; srcint=5e15; tlife=20e6; // years
//rep="./dataref";ncpu=16; srcint=5e15; tlife=20e6; // years

nsnap=13;
lcoarse=8.;
bint=spanl(1e8,1e9,128);

delta=150.;
MPC=3.08e22; // m
H=67; // km/s/Mpc
Hor=H;
G=6.67e-11; // S.I.
mp=1.6e-27; // proton mass kg
omegam=0.3175;
omegab=0.049;
n=2^(3*lcoarse);
lbox=12./(H/100.)*MPC;// m
lorg=lbox/MPC;
msol=2e30;

// DM particle mass
H=H*1e3/MPC; // s^-1
rhoc=3.*H^2/(8.*pi*G);
mdm=((omegam-omegab)*rhoc*lbox^3)/n/msol;
munit=((omegam)*rhoc*lbox^3)/msol;
nH=((omegab)*rhoc*lbox^3)/mp;
munit;

// COmputing the SFR
s=mergepart(swrite(format=rep+"/star.%05d",nsnap),ncpu,a,star=1);
so=sort(s(11,));
s=s(,so);

//ms=indgen(dimsof(s)(0))*ms;//Solar masses
ts=s(11,); // Myrs
//plg,ms,ts;
//SFR=ms(dif)/ts(dif)/lorg^3;

SFR=histo1d(s(11,),bint,wght=s(8,))*munit/bint(dif)/lorg^3;

vt=(univAge(10000.,h0=Hor,Omega_m=omegam,lambda0=1.-omegam)-univAge(span(90,4,256),h0=Hor,Omega_m=omegam,lambda0=1.-omegam))/(3600.*24.*365);
zSFR=interp(span(90,4,256),vt,bint(zcen));

// emissivity of photons
//te=s(11,);
te=bint(zcen);
tlife=20e6;
srcint=5e15;
ne=te*0.;
for(i=1;i<=dimsof(s)(0);i++){
  heav=(te>=s(11,i))*(te<=(s(11,i)+tlife));
  ne+=integ(heav*(srcint*3600*24*365)*s(8,i)*munit*msol,te,te);
 }
ne/=nH;

// BOuwens
zB=[3.8,4.9,5.9,6.8,7.9,10.4];
SFRB=10.^[-1.44,-1.53,-1.84,-1.95,-2.20,-3.33]; // Ms/Mpc^3/yr
SFRB2=10.^[-1.06,-1.19,-1.59,-1.72,-2.05,-3.18]; // Ms/Mpc^3/yr
eSFRB=abs(SFRB*10.^[0.06,0.06,0.06,0.06,0.07,0.36]-SFRB);



// DISPLAY
ws,2;
plshade,[SFRB,SFRB2],zB,color=__rgb(,30);
PL,SFR(),zSFR(),color="green",msize=.3;
logxy,0,1;
limits,2,12,5e-5,1;
xytitles,"redshift z","SFR [MS/yr/Mpc^3^]";
ws,4;
plg,ne,zSFR;
plg,ne*0+1,zSFR,type=2;
logxy,0,1;
