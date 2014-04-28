#include "Bastien/cosmo.i"

z=load("../HM2012/z.out");
z=z(1,);
dummy=load("../HM2012/e.out");
lam=dummy(,1)*1e-10; // wavelength in meters
em=dummy(,2:); // comoving emmisivity in erg/s/Hz/Mpc3
c=3e8; // speed of light (m/s);
nu=c/lam; // frequency (Hz);
hplanck=6.626e-34;
eev=hplanck*nu/1.6e-19;

www=where(eev>=13.6);
eionz=nionz=[];
eavg=[];
for(i=1;i<=numberof(z);i++){
  eion=integ(em(www,i),nu(www),nu(www)(0));
  nphot=integ(em(www,i)*1e-7/(hplanck*nu(www)),nu(www),nu(www)(0));
  write,"eion=",eion," from hnu=",nu(www)(1)," to hnu=", nu(www)(0);
  grow,eionz,abs(eion)*1e-7/(3.06e22)^3; // J/s/m3
  grow,eavg,abs(eion/nphot)*1e-7/1.6e-19;
  grow,nionz,abs(eion)*1e-7/(3.06e22)^3/(eavg(i)*1.6e-19); // n photons/s/m3
 }

#if 1
t=univAge(10000)-univAge(z);
st=sort(t);
ts=t(st);
nionzs=nionz(st); // number of ionising photons/m3/s
zs=z(st);
ncum=integ(nionzs,ts,ts);

#if 1
//====================================================================
H0=70e3/3.08e22; // s-1
om=0.3;
ov=0.7;


// redshift

coolz=xionz=[];
for(k=2;k<=numberof(zs);k+=10){
  k;
  ncum=integ(nionzs(k-1:k),ts(k-1:k),ts(k));
  zc=0.5*(zs(k)+zs(k-1));
  H=H0*sqrt(om*(1+zc)^3+ov);
  
  // temperature
  Temp=spanl(1e2,1e8,1024);
  // density
  navg=0.2; // average atom/m3
  nh=spanl(1e-3,1e10,1024)*navg*(1+zc)^3;

  cool=xion=[];
  for(j=1;j<=numberof(Temp);j++){
    temp=Temp(j);

    // nion
    ngamma=ncum*exp(-0*(nh/(200*navg))^1)*(1+zc)^3;
    
    // Photoionisation rate
    sigma=2.93055200728e-22; // cross-section (m2) for 20.27 eV photon
    cel=299792458.*0.01; // speed of light (m/s)

  
    //case A recombination rate
    lambda=2e0*157807e0/temp;
    alpha=1.269e-13*pow(lambda,1.503)/pow(1e0+pow(lambda/0.522,0.470),1.923)*1e-6; //m3/s

    // collisional ionisation rate
    T5=temp/1e5;
    beta=5.85e-11*sqrt(temp)/(1+sqrt(T5))*exp(-(157809e0/temp))*1e-6; // m3/s

    // solving equilibrium chem
    x=[];
    for(i=1;i<=numberof(nh);i++){
      //i;
      a=(alpha+beta)*nh(i);
      b=(cel*sigma*ngamma(i)+3*H-beta*nh(i));
      c=-cel*sigma*ngamma(i)-3*H;

      delta=b^2-4*a*c;
      q=-0.5*(b+sign(b)*sqrt(delta));
      x1=q/a;
      if(q==0){
        x2=x1;
      }
      else{
        x2=c/q;
      }
      if((x1>=0)*(x1<=1)) xs=x1;
      if((x2>=0)*(x2<=1)) xs=x2;
      //x1;x2;
      grow,x,xs;
    }

    x=min(1.,x*10.);
    
    // 1.-x;
    // plg,1.-x,nh,color="blue";
    // logxy,1,1;
    //plg,ncum/0.2,zs;


    //amp=1.2e-16;sig=1.;zavg=2;mz=1e-18;pz=1.2e-17;

    //ws;plg,amp/(sig*sqrt(2*pi))*exp(-(z-zavg)^2/(2.*sig^2))+mz*z+pz,z,color="red";plg,nionz,z;

    // === cooling
    CLUMPF=1.;
    nh2=nh*1e-6; // m-3  -> cm-3
    c1=exp(-157809.1e0/temp)*1.27e-21*sqrt(temp)/(1.+sqrt(temp/1e5))*x*(1.-x)*nh2*nh2*CLUMPF;


    // Case A Recombination Cooling

    c2=1.778e-29*temp*pow(2e0*157807e0/temp,1.965e0)/pow(1.+pow(2e0*157807e0/temp/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2*CLUMPF;
  

    // Collisional excitation cooling

    c3=exp(-118348e0/temp)*7.5e-19/(1+sqrt(temp/1e5))*x*(1.-x)*nh2*nh2*CLUMPF;


    // Bremmsstrahlung

    c4=1.42e-27*1.5e0*sqrt(temp)*x*x*nh2*nh2*CLUMPF;

    // Compton Cooling

    c5=1.017e-37*pow(2.727*(1+zc),4)*(temp-2.727*(1+zc))*nh2*x;
    //c5=(c5>0.?c5:0); // to be checked

    // Overall Cooling

    lambda=c1+c2+c3+c4+c5;// ! erg*cm-3*s-1
    lambda=lambda*1e-7/1.6e-19*1e6; // eV/m3/s
    // Heating

    heat=ngamma*cel*sigma*(20.24-13.6)*nh*(1.-x); // eV/m3/s

    // Net cooling
    Cool=heat-lambda;

    grow,cool,[Cool];
    grow,xion,[x];
  }

  grow,coolz,[cool];
  grow,xionz,[xion];

 }

T=Temp(-:1:numberof(nh),);
N=nh(,-:1:numberof(Temp));
