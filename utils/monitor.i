
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


func monitor(dir,tpause=,zmax=,zmin=){
  if(is_void(zmax)) zmax=20;
  if(is_void(tpause)) tpause=30;

  
  tpause=int(tpause*1e3);
  vi=readinfo(dir);
  h0=vi(13);
  om=vi(10);
  ol=vi(11);
  mstar=vi(15);
  lbox=vi(2);


  write,"########################################";
  
  while(1){
    vv=asciiRead(dir+"/param.avg");
    
    z=vv(3,);
    a=vv(2,);
    lev=vv(5,);
    avgx=vv(7,);
    avgT=vv(8,);
    nstar=vv(10,);
    ts=(univAge(10000,h0=h0,Omega_m=om,Omega_l=ob,silent=1)-univAge(vv(3,),h0=h0,Omega_m=om,Omega_l=ol,silent=1))/(3600.*24*365.25);
    sfr=nstar(dif)/ts(dif)*mstar/(lbox/(h0/100.))^3;

    
    aa=write(format=" NSteps=%d z=%f [aexp=%f] === xHI=%e avg Temp=%e [K]\n",int(vv(1,0)),z(0),a(0),1.-avgx(0),avgT(0));    
    ws,1;
    window,1,style="win22.gs";

    plsys,1;
    plg,avgT,z;
    limits,zmax,zmin;
    range,1e2,1e7;
    logxy,1,1;

    plsys,2;
    PL,sfr,z(zcen);
    range,1e-5,1e1;
    limits,zmax,zmin;
    logxy,1,1;

    plsys,3;
    plg,lev,z;
    range,5,15;
    limits,zmax,zmin;
    logxy,1,1;

    plsys,4;
    plg,1.-avgx,z;
    range,1e-6,2.;
    limits,zmax,zmin;
    logxy,1,1;
    
    pause,tpause;
    
  }

}
