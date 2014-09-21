#include "Bastien/cosmo.i"

//rep="./data0_dyncool";ncpu=32;
rep="./data0_dyncoll_nosrc";ncpu=32;
//rep="./dataref";ncpu=16;

execut="./utils/oct2grid_cosmo ";
nsnap=19;
level=8;
nsamp=1;

MPC=3.08e22; // m
H=67; // km/s/Mpc
Hor=H;
G=6.67e-11; // S.I.
omegam=0.3175;
omegab=0.045;
n=128^3;
lbox=12./(H/100.)*MPC;// m
lorg=lbox/MPC;

H=H*1e3/MPC; // s^-1
rhoc=3.*H^2/(8.*pi*G);
rhobaryon=((omegab)*rhoc);
mH=(1.67262158e-27); // kg
nH=rhobaryon/mH; // atoms/m3


// stellar particle mass
ncell=pow(2,level);
//xion=array(float,ncell,ncell,ncell,nsnap);
//temp=array(float,ncell,ncell,ncell,nsnap);
xa=array(float,nsnap);
ea=array(double,nsnap);
ta=array(float,nsnap);

va=array(double,nsnap);
vt=array(double,nsnap);

for(i=1;i<=nsnap;i+=nsamp){
  xion=oct2cube(swrite(format=rep+"/grid.%05d",i),level,706,a,ncpu=ncpu,execut=execut);
  xa(i)=xion(avg,avg,avg,);
  ta(i)=oct2cube(swrite(format=rep+"/grid.%05d",i),level,707,a,ncpu=ncpu,execut=execut)(avg,avg,avg);
  ea(i)=(10.^oct2cube(swrite(format=rep+"/grid.%05d",i),level,701,a,ncpu=ncpu,execut=execut))(avg,avg,avg);
  va(i)=a;
  vt(i)=(univAge(10000.,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam)-univAge(1./a-1.,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam));

  
 }

www=where(va>0);
va=va(www);
vt=vt(www);
xa=xa(www);
ta=ta(www);
ea=ea(www);
vz=1./va-1.;

// ===== computing tau

// tab redshift
ze=array(float,2*numberof(va));
ze(:numberof(va)+1)=span(0,vz(0),numberof(va)+1);
ze(numberof(va)+2:2*numberof(va))=vz(1:-1)(::-1);
te=univAge(ze,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam);


// tab electron density
ne=array(float,2*numberof(va));
ne(:numberof(va))=nH;
ne(numberof(va)+1:2*numberof(va))=nH*xa(::-1);
ne=ne*(1+ze)^3; // physical density

//computing tau
clight=299792458.; // m
sigmat=6.65e-29; //m^2;
itau=te(dif)/ze(dif)*ne(zcen)*clight*sigmat;
tau=integ(itau,ze(zcen),ze(zcen));

window,0;
plg,1.-xa,vz;
logxy,0,1;
limits,4,12,1e-6,2;
xytitles,"redshift z","x_H";

xF=[5.5e-5,6.7e-5,6.5e-5,8e-5,1.1e-4,4e-4];
zF=[5.02,5.24,5.45,5.65,5.85,6.12];
xF1=1e-5*[4,4.5,4,4.9,7.8,12.];
xF2=1e-5*[7.1,9,9.5,13.,100000.,100000.];
plshade,[xF1,xF2],zF,color=__rgb(,30);
plg,1.-xa,vz;

window,1;
tauP=[array(0.097-0.038,numberof(tau)),array(0.097+0.038,numberof(tau))];
plshade,tauP,ze(zcen),color=__rgb(,30);
plg,tau,ze(zcen),color="red",width=5;
range,-0.01,0.16;
xytitles,"redshift z","!t";
