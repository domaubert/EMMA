rep="../data";

nsnap=81;

lcoarse=6.;
delta=50.;
MPC=3.08e22; // m
H=67; // km/s/Mpc
Hor=H;
G=6.67e-11; // S.I.
omegam=0.3175;
omegab=0.045;
n=2^(3*lcoarse);
lbox=6./(H/100.)*MPC;// m
lorg=lbox/MPC;
msol=2e30;

// DM particle mass
H=H*1e3/MPC; // s^-1
rhoc=3.*H^2/(8.*pi*G);
mdm=((omegam-omegab)*rhoc*lbox^3)/n/msol;
munit=((omegam)*rhoc*lbox^3)/msol;
munit;

#if 1
s=mergepart(swrite(format=rep+"/star.%05d",nsnap),8,a,star=1);
so=sort(s(11,));
s=s(,so);

//ms=indgen(dimsof(s)(0))*ms;//Solar masses
ts=s(11,); // Myrs
//plg,ms,ts;
//SFR=ms(dif)/ts(dif)/lorg^3;

bint=span(1e8,1e9,16);
SFR=histo1d(s(11,),bint,wght=s(8,))*munit/bint(dif)/lorg^3;

vt=(univAge(10000.,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam)-univAge(span(90,4,256),h0=Hor,Omega_m=omegam,Omega_l=1.-omegam))/(3600.*24.*365);
zSFR=interp(span(90,4,256),vt,bint(zcen));
#endif

#if 0

ns=array(double,nsnap)*0.;
va=array(double,nsnap);
vt=array(double,nsnap);
for(i=1;i<=nsnap;i++){
  s=mergepart(swrite(format=rep+"/star.%05d",i),8,a,star=1);
  if(!is_void(s)){
    ns(i)=dimsof(s)(0)*ms;
  }
  va(i)=a;
  vt(i)=(univAge(10000.,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam)-univAge(1./a-1.,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam))/(3600.*24.*365);
 }

SFR=ns(dif)/vt(dif)/lorg^3;
zSFR=1./interp(va,vt,vt(zcen))-1.;

#endif

zB=[3.8,4.9,5.9,6.8,7.9,10.4];
SFRB=10.^[-1.44,-1.53,-1.84,-1.95,-2.20,-3.33]; // Ms/Mpc^3/yr
SFRB2=10.^[-1.06,-1.19,-1.59,-1.72,-2.05,-3.18]; // Ms/Mpc^3/yr
eSFRB=abs(SFRB*10.^[0.06,0.06,0.06,0.06,0.07,0.36]-SFRB);
