
#include "Bastien/cosmo.i"

func readinfo(dir){
  f=open(dir+"/param.info","r");
  ll=rdfile(f);
  close,f;

  vv=[];
  dummy=array(double);
  st=array(string);
  for(i=1;i<=numberof(ll);i++){
    if(strpart(ll(i),1:1)=="#") continue;
    a=sread(ll(i),format="%s %f",st,dummy);
    grow,vv,dummy;
  }

  return vv;
}


func monitor(dir,tpause=,zmax=,zmin=,clear=,col=,mstar=){
  if(is_void(zmax)) zmax=20;
  if(is_void(tpause)) tpause=30;
  if(is_void(clear)) clear=0;
  zB=[3.8,4.9,5.9,6.8,7.9,10.4];
  SFRB=10.^[-1.44,-1.53,-1.84,-1.95,-2.20,-3.33]; // Ms/Mpc^3/yr
  SFRB2=10.^[-1.06,-1.19,-1.59,-1.72,-2.05,-3.18]; // Ms/Mpc^3/yr
  
  tpause=int(tpause*1e3);
  vi=readinfo(dir);
  h0=vi(13);
  om=vi(10);
  ol=vi(11);
  if(is_void(mstar))mstar=vi(15);
  lbox=vi(2);
  mstar=om*(3.*(h0*1e3/3.08e22)^2)/8./pi/6.67e-11*(lbox*3.08e22/h0*100.)^3/2e30/128.^3/64.;//msol

  write,"########################################";
  
  vv=asciiRead(dir+"/param.avg");
  //  vv=load(dir+"/param.avg");
    
    z=vv(3,);
    a=vv(2,);
    lev=vv(5,);
    maxd=vv(6,);
    avgx=vv(7,);
    avgT=vv(8,);
    maxT=vv(9,);
    nstar=vv(10,);
    nSN=vv(11,);
    src=vv(12,);
    ts=(univAge(10000,h0=h0,Omega_m=om,lambda0=ol,silent=1)-univAge(vv(3,),h0=h0,Omega_m=om,lambda0=ol,silent=1))/(3600.*24*365.25);
    sfr=nstar(dif)*(ts(dif)>0)/(ts(dif)+1e-15)*mstar/(lbox/(h0/100.))^3;
    aa=write(format=" NSteps=%d z=%f [aexp=%f] === xHI=%e avg Temp=%e [K]\n",int(vv(1,0)),z(0),a(0),1.-avgx(0),avgT(0));    
    if(clear){
      WS,1,dpi=120;
      window,1,style="win42.gs";
    }

    plsys,1;
    plg,avgT,z,color=col;
    plg,maxT,z,color=col;
    limits,zmax,zmin;
    range,1e2,1e8;
    logxy,1,1;

    plsys,2;
    plg,SFRB,zB,color=__rgb(,30),type=2;
    plg,SFRB2,zB,color=__rgb(,30),type=2;
    plg,sfr,z(zcen),color=col;
    PL,sfr(0),z(zcen)(0),color=col;
    range,1e-5,1e1;
    limits,zmax,zmin;
    logxy,1,1;

    plsys,3;
    plg,lev,z,color=col;
    range,5,15;
    limits,zmax,zmin;
    logxy,1,0;

    plsys,4;
    plg,1.-avgx,z,color=col;
    plg,avgx,z,color=col;
    range,1e-6,2.;
    limits,zmax,zmin;
    logxy,1,1;

    plsys,5;
    plg,maxd,z,color=col;
    limits,zmax,zmin;
    range,1e0,1e10;
    logxy,1,1;

    plsys,6;
    plg,nstar+1,z,color=col;
    plg,nSN+1,z,color=col;
    range,1,1e8;
    limits,zmax,zmin;
    logxy,1,1;

    plsys,7;
    plg,src,z,color=col;
    limits,zmax,zmin;
    range,1e-4,1e2;
    logxy,1,1;


    redraw;
}
