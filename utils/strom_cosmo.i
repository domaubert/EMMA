#include "Bastien/cosmo.i"

H0=73.;
h0=H0*1e3/3.085677e22;
omegam=0.27;
omegav=0.73;

alpha=4.5;
 C=1.;
// ts=3e6*365*24*3600.;
// Ng=4000.;
// Ndot=Ng/(ts)*(alpha-1.)/alpha;

Ndot=5e52; //photons/sec

alphab=2.5918e-19;
n0=0.2;

zi=15.;
zf=0.;
ai=1./(1.+zi);
af=1./(1.+zf);
a=spanl(ai,af,128);
z=1./a-1.;


ti=univAge(10000,h0=H0,Omega_l=omegav,Omega_m=omegam)-univAge(zi,h0=H0,Omega_l=omegav,Omega_m=omegam);
t=univAge(10000,h0=H0,Omega_l=omegav,Omega_m=omegam)-univAge(z,h0=H0,Omega_l=omegav,Omega_m=omegam);

V=array(0.,numberof(z));
S=array(0.,numberof(z));

//Vmax=Ng/n0;


for(i=2;i<=numberof(z);i++){
  i;
  va=spanl(a(1),a(i),128);
  vf=sqrt(1./va^3+(1.-omegam)/omegam);
  floc=vf(0);
  vz=1./va-1.;
  vt=univAge(10000,h0=H0,Omega_l=omegav,Omega_m=omegam)-univAge(vz,h0=H0,Omega_l=omegav,Omega_m=omegam);
    
  F=-2./3.*C*alphab*n0/sqrt(omegam)/h0*(vf-floc);
  //F=-0.262*(vf-floc);
  //ndotloc=Ndot*(((vt-vt(1))<ts)+((vt-vt(1)+1e-18)/ts)^(-alpha)*((vt-vt(1))>=ts));
  ndotloc=Ndot;
  integrand=ndotloc/n0*exp(F);
  V(i)=integ(integrand,vt,vt(0));
 }

//V/=Vmax;
PARSEC=3.085677e16; // in m
unit_l=20.e6*PARSEC; // boxsize is unit_length

V/=(unit_l)^3;
R=(V/(4./3.*pi))^(1./3.);
