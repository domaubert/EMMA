#include "Bastien/cosmo.i"

//rep="./data0_dyncool";ncpu=32;
//rep="./data0_dyncoll_nosrc";ncpu=32;execut="./utils/oct2grid_cosmo ";
//rep="./dataref";ncpu=16;
//rep="./data_coarse_256_24MPC/";ncpu=256;execut="./utils/oct2grid "sbox=24.;
//rep="./data_12_noschaye/";ncpu=32;execut="./utils/oct2grid ";sbox=12.;
rep="./data/";ncpu=32;execut="./utils/alloct ";sbox=4.;
rep="./data_4_new_wsrc/";ncpu=32;execut="./utils/alloct ";sbox=4.;
rep="./data_4_new_wsrc_ministar_x100/";ncpu=32;execut="./utils/alloct ";sbox=4.;
//rep="./data_4_new_wsrc_ministar_x10/";ncpu=32;execut="./utils/alloct ";sbox=4.;
//rep="./data_4_new_wsrc_ministar_x1/";ncpu=32;execut="./utils/alloct ";sbox=4.;
//rep="./data_4_new_wsrc_ministar_x3/";ncpu=32;execut="./utils/alloct ";sbox=4.;
rep="./data_4_new_wsrc_ministar_x3_mono/";ncpu=32;execut="./utils/alloct ";sbox=4.;
//rep="./data_4_new_wsrc_ministar_x3_mono_noamr50/";ncpu=32;execut="./utils/alloct ";sbox=4.;
rep="./data16_TH035/";ncpu=256;execut="./utils/alloct ";sbox=16.;
rep="./data_8_TH_IDM00/";ncpu=256;execut="./utils/alloct ";sbox=16.;



#if 1
nsnap=22;
level=11;
nsamp=1;

MPC=3.08e22; // m
H=67; // km/s/Mpc
Hor=H;
G=6.67e-11; // S.I.
omegam=0.3175;
omegab=0.049;
n=128^3;
lbox=sbox/(H/100.)*MPC;// m
lorg=lbox/MPC;

H=H*1e3/MPC; // s^-1
rhoc=3.*H^2/(8.*pi*G);
rhobaryon=((omegab)*rhoc);
mH=(1.67262158e-27); // kg
nH=rhobaryon/mH; // atoms/m3


c=299792458.*0.1;
unitN=rhoc*omegam/mH;
sig=[3.00563459955e-22,5.04006266253e-23,7.89888816374e-24];


// stellar particle mass
ncell=pow(2,level);
//xion=array(float,ncell,ncell,ncell,nsnap);
//temp=array(float,ncell,ncell,ncell,nsnap);
xa=array(float,nsnap);
gamma=array(float,nsnap);
xam=array(float,nsnap);
dxa=array(float,nsnap);

va=array(double,nsnap);
vt=array(double,nsnap);


for(i=1;i<=nsnap;i+=nsamp){

#if 0
  xion=oct2cube(swrite(format=rep+"/grid.%05d",i),level,706,a,ncpu=ncpu,execut="utils/oct2grid ");
  dxion=oct2cube(swrite(format=rep+"/grid.%05d",i),level,1006,a,ncpu=ncpu,execut="utils/oct2grid ");
  xa(i)=xion(*)(avg);
  dxa(i)=dxion(*)(avg);
#endif

#if 1
  xion=alloct(swrite(format=rep+"/%05d/grid/grid.%05d",i,i),level,706,a,ncpu=ncpu,execut="utils/alloct ");
  xa(i)=(xion(5,)*pow(0.5,3*xion(4,)))(sum)/pow(0.5,xion(4,)*3)(sum);
#endif

#if 0
  ng=0.;
  gamma(i)=0.;
  for(ig=0;ig<3;ig++){
    nion=alloct(swrite(format=rep+"/grid.%05d",i),level,701+ig*10,a,ncpu=ncpu,execut="utils/alloct ");
    gamma(i)+=((nion(5,)*unitN*c*sig(ig+1)/a^3)*pow(0.5,3*nion(4,)))(sum)/pow(0.5,3*nion(4,))(sum);
  }
 #endif 


  //  dxion=alloct(swrite(format=rep+"/grid.%05d",i),level,1006,a,ncpu=ncpu,execut="utils/alloct ");
  //  d=alloct(swrite(format=rep+"/grid.%05d",i),level,101,a,ncpu=ncpu,execut="utils/alloct ");
  //  xam(i)=(d(5,)*xion(5,)*pow(0.5,3*xion(4,)))(sum)/(d(5,)*pow(0.5,xion(4,)*3))(sum);
  //  dxa(i)=(dxion(5,)*pow(0.5,3*dxion(4,)))(sum)/pow(0.5,dxion(4,)*3)(sum);

  

  
  //ta(i)=oct2cube(swrite(format=rep+"/grid.%05d",i),level,707,a,ncpu=ncpu,execut=execut)(avg,avg,avg);
  //ea(i)=(10.^oct2cube(swrite(format=rep+"/grid.%05d",i),level,701,a,ncpu=ncpu,execut=execut))(avg,avg,avg);
  va(i)=a;
  vt(i)=(univAge(10000.,h0=Hor,Omega_m=omegam,lambda0=1.-omegam)-univAge(1./a-1.,h0=Hor,Omega_m=omegam,lambda0=1.-omegam));

  
 }

www=where(va>0);
va=va(www);
vz=1./va-1.;
vt=vt(www);
xa=xa(www);
xam=xam(www);
dxa=dxa(www);

#if 0

//ta=ta(www);
//ea=ea(www);

// ===== computing tau

// tab redshift
ze=array(float,2*numberof(va));
ze(:numberof(va)+1)=span(0,vz(0),numberof(va)+1);
ze(numberof(va)+2:2*numberof(va))=vz(1:-1)(::-1);
te=univAge(ze,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam);


// tab electron density
davg=omegab/omegam;
ne=array(float,2*numberof(va));
ne(:numberof(va))=nH;
ne2=ne;
ne(numberof(va)+1:2*numberof(va))=nH*xa(::-1);
ne2(numberof(va)+1:2*numberof(va))=dxa(::-1)/davg*nH;
ne=ne*(1+ze)^3; // physical density
ne2=ne2*(1+ze)^3; // physical density

//computing tau
clight=299792458.; // m
sigmat=6.65e-29; //m^2;
itau=te(dif)/ze(dif)*ne(zcen)*clight*sigmat;
itau2=te(dif)/ze(dif)*ne2(zcen)*clight*sigmat;
tau=integ(itau,ze(zcen),ze(zcen));
tau2=integ(itau2,ze(zcen),ze(zcen));

#if 0

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

log10gamma_cal=[-13.1,-13.09,-13.51,-12.99,-12.94,-12.64,-12.37,-12.86,-12.71,-11.74,-12.32,-12.02,-11.99,-11.41,-11.73];
z_cal=[6.42,6.31,6.25,6.02,6.016,5.82,5.810,5.41,5.33,5.2,5.09,4.967,4.886,4.876,4.588];
e_cal=[0.53,0.46,0.64,0.42,0.55,0.43,0.40,0.46,0.48,1.15,0.69,0.34,0.54,0.58,0.54];
