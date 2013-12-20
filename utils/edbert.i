
#include "Bastien/cosmo.i"
require,"utils/readamr.i";
require,"utils/readpart.i";

th=span(1e-3,2*pi,4096);
lambda=sin(th/2.)^2*((th-sin(th))/pi)^(-8./9.);
Vth=lambda*sin(th)*(th-sin(th))/(1.-cos(th)+1e-10)^2;



#if 1
//i=1;plh,log10(hd(,i)+1e-10),log10(binr(zcen)/rt(i));plg,log10(D),log10(lb);plh,log10(hdb(,i)+1e-10),log10(binr(zcen)/rt(i)),color="green";plg,log10(DB*0.01),log10(lbb),color="green";plh,vd(,i)/vt(i),log10(binr(zcen)/rt(i)),color="blue";plg,VB,log10(lbb),color="blue";
//wkl;

imin=100;
imax=100;
level=12;
ncpu=4;
n=2^level;
R=0.05;
deltai=0.2;
astart=1e-3;
om=1.0;
ob=0.01;
deltaib=deltai*ob;
field=1;

rita=R/deltai*astart; // initial turnaround radius 

// r =(span(-0.5,0.5,(n+1))^2)(zcen)(,-:1:n)(,,-:1:n);
// r+=(span(-0.5,0.5,(n+1))^2)(zcen)(-:1:n,)(,,-:1:n);
// r+=(span(-0.5,0.5,(n+1))^2)(zcen)(-:1:n,)(-:1:n,,);
// r=sqrt(r);
//ws;
binr=spanl(0.01,1e1,32);

hd=array(float,numberof(binr)-1,imax-imin+1);
hl=array(float,numberof(binr)-1,imax-imin+1);
hdb=array(float,numberof(binr)-1,imax-imin+1);
nb=array(float,numberof(binr)-1,imax-imin+1);
ehdb=array(float,numberof(binr)-1,imax-imin+1);
ed=array(float,numberof(binr)-1,imax-imin+1);
vd=array(float,numberof(binr)-1,imax-imin+1);
evd=array(float,numberof(binr)-1,imax-imin+1);
evd=array(float,numberof(binr)-1,imax-imin+1);
ts=array(float,imax-imin+1);
a=array(float,imax-imin+1);
rt=array(float,imax-imin+1);
mt=array(float,imax-imin+1);
rtb=array(float,imax-imin+1);
mtb=array(float,imax-imin+1);
vt=array(float,imax-imin+1);
rhoH=array(float,imax-imin+1);

ti=(univAge(10000.,h0=1.,k=0.,Omega_l=0.,Omega_m=1.,silent=1)-univAge(1./astart-1.,h0=1.,k=0.,Omega_l=0.,Omega_m=1.,silent=1)); 
tita=0.75*pi/deltai^1.5*ti;


for(i=imin;i<=imax;i++){
  i;
  //dg=oct2cube(swrite(format="~/Project/Quartz/data/grid.%05d",i),level,101,t,ncpu=ncpu,execut="~/Project/Quartz/utils/oct2grid ");
  l=oct2cell(swrite(format="~/Project/Quartz/data/grid.%05d",i),level,0,aexp,ncpu=4,execut="~/Project/Quartz/utils/oct2cell ");
  p=mergepart(swrite(format="~/Project/Quartz/data/part.%05d",i),ncpu,t);
  www=where(p(8,)==max(p(8,)));
  cdm=p(:3,www)(,avg);
  vcdm=p(4:6,www)(,avg);

  tau=((univAge(10000.,h0=1.,k=0.,Omega_l=0.,Omega_m=1.,silent=1)-univAge(1./aexp-1.,h0=1.,k=0.,Omega_l=0.,Omega_m=1.,silent=1)))/ti;
  tsurts=((univAge(10000.,h0=1.,k=0.,Omega_l=0.,Omega_m=1.,silent=1)-univAge(1./aexp-1.,h0=1.,k=0.,Omega_l=0.,Omega_m=1.,silent=1)))/2.*(1e3/3.0856e22); 
  tsurtita=((univAge(10000.,h0=1.,k=0.,Omega_l=0.,Omega_m=1.,silent=1)-univAge(1./aexp-1.,h0=1.,k=0.,Omega_l=0.,Omega_m=1.,silent=1)))/tita;

  //rt(i-imin+1)=(0.75*pi)^(-(8./9.))*(deltai)^(1./3.)*R*tau^(8./9.)*astart/aexp;
  rt(i-imin+1)=rita*(tsurtita)^(8./9.)/aexp;
  
  rtb(i-imin+1)=(0.75*pi)^(-(8./9.))*(deltaib)^(1./3.)*R*tau^(8./9.)*astart/aexp;
  mt(i-imin+1)=4./3.*pi*rt(i-imin+1)^3*(om-ob);
  mtb(i-imin+1)=4./3.*pi*rt(i-imin+1)^3*(ob);
  vt(i-imin+1)=3.*aexp^0.5*rt(i-imin+1); // here we assume om=1.
  ts(i-imin+1)=t;

  if(field){
    l=oct2cell(swrite(format="~/Project/Quartz/data/grid.%05d",i),level,0,aexp,ncpu=4,execut="~/Project/Quartz/utils/oct2cell ");
    db=oct2cell(swrite(format="~/Project/Quartz/data/grid.%05d",i),level,101,aexp,ncpu=4,execut="~/Project/Quartz/utils/oct2cell ");
    vx=oct2cell(swrite(format="~/Project/Quartz/data/grid.%05d",i),level,102,aexp,ncpu=4,execut="~/Project/Quartz/utils/oct2cell ");
    vy=oct2cell(swrite(format="~/Project/Quartz/data/grid.%05d",i),level,103,aexp,ncpu=4,execut="~/Project/Quartz/utils/oct2cell ");
    vz=oct2cell(swrite(format="~/Project/Quartz/data/grid.%05d",i),level,104,aexp,ncpu=4,execut="~/Project/Quartz/utils/oct2cell ");
    a(i-imin+1)=aexp;
    
    vx=vx(5,)-vcdm(1)+2.*(db(1,)-cdm(1))*aexp^0.5;
    vy=vy(5,)-vcdm(2)+2.*(db(2,)-cdm(2))*aexp^0.5;
    vz=vz(5,)-vcdm(3)+2.*(db(3,)-cdm(3))*aexp^0.5; // total velocity (hubble flow included)

    vb=(vx*(db(1,)-cdm(1))+vy*(db(2,)-cdm(2))+vz*(db(3,)-cdm(3)))/abs(db(1,)-cdm(1),db(2,)-cdm(2),db(3,)-cdm(3)); // radial velocity
    vb=vb/vt(i-imin+1);

    
    rb=abs(db(1,)-cdm(1),db(2,)-cdm(2),db(3,)-cdm(3))/rt(i-imin+1);
    //vb=(((db(1,)-0.5)*vx(5,)+(db(2,)-0.5)*vy(5,)+(db(3,)-0.5)*vz(5,))/rt(i-imin+1)/rb)/vt(i-imin+1)+2./3.*rb; // divided bby aexp because of supercomoving description
    

    nb(,i-imin+1)=histo1d(rb,binr)+1e-10;
    hl(,i-imin+1)=histo1d(rb,binr,wght=l(5,))/nb(,i-imin+1);
    hdb(,i-imin+1)=histo1d(rb,binr,wght=db(5,)*pow(0.5,3*db(4,))/mtb(i-imin+1))/(3*binr(zcen)^2*binr(dif));
    vd(,i-imin+1)=histo1d(rb,binr,wght=vb*db(5,)*pow(0.5,3*db(4,)))/(histo1d(rb,binr,wght=db(5,)*pow(0.5,3*db(4,)))+1e-10);
    evd(,i-imin+1)=histo1d(rb,binr,wght=vb^2*db(5,)*pow(0.5,3*db(4,)))/(histo1d(rb,binr,wght=db(5,)*pow(0.5,3*db(4,)))+1e-10);
  }

  rp=abs(p(1,)-cdm(1),p(2,)-cdm(2),p(3,)-cdm(3))/rt(i-imin+1);
  hd(,i-imin+1)=histo1d(rp,binr,wght=p(8,)/mt(i-imin+1))/(3*binr(zcen)^2*binr(dif));



  // fma;
  // limits,-2.5,2.5,-2.5,2.5;
  // pl,(p(2,)-0.5)/rt(i-imin+1),(p(1,)-0.5)/rt(i-imin+1);
  // pause(100);
  
 }



lb=indgen(50)*0.02;
//M=[0.561,0.909,1.20,1.45,1.67,1.88,2.06,2.24,2.45,2.56,2.70,2.92,3.00,3.09,3.19,3.31,3.45,3.68,3.85,3.88,3.91,3.95,3.98,4.02,4.06,4.1,4.14,4.18,4.23,4.27,4.32,4.37,4.42,4.47,4.53,4.58,4.64,4.70,4.76,4.82,4.89,4.95,5.02,5.09,5.16,5.24,5.31,5.39,5.47,5.55];
D=[1.58e4,4.41e3,1.09e3,5.19e2,3.12e2,2.53e2,1.70e2,9.69e1,5.42e1,5.16e1,5.57e1,2.23e1,2.10e1,2.02e1,2.0e1,2.08e1,2.41e1,4.75e1,3.6,3.38,3.19,3.03,2.88,2.74,2.62,2.51,2.41,2.32,2.24,2.16,2.09,2.02,1.96,1.91,1.86,1.81,1.77,1.72,1.68,1.65,1.61,1.58,1.55,1.53,1.50,1.47,1.45,1.43,1.40,1.39];

// db=mb(dif)/lb(dif)/3./lb(zcen)^2;
// lb=lb(zcen);


lbb=[0.017,0.035,0.052,0.069,0.087,0.104,0.122,0.139,0.156,0.174,0.191,0.208,0.226,0.243,0.260,0.278,0.295,0.312,0.330,0.347,0.34701];
DB=[2.22e4,4.71e3,1.84e3,9.40e2,5.54e2,3.59e2,2.47e2,1.79e2,1.34e2,1.03e2,8.11e1,6.51e1,5.3e1,4.36e1,3.62e1,3.04e1,2.57e1,2.19e1,1.88e1,1.61e1,4.02];
VB=-[1.32e-3,1.45e-3,2.64e-3,4.47e-3,6.72e-3,9.39e-3,1.25e-2,1.61e-2,2.03e-2,2.49e-2,3.05e-2,3.65e-2,4.31e-2,5.13e-2,6.07e-2,7.07e-2,8.18e-2,9.44e-2,1.09e-1,1.27e-1,1.43];
www=where(lb>lbb(0));
grow,DB,D(www(2:));
grow,VB,interp(Vth,lambda,lb(www(2:)));
grow,lbb,lb(www(2:));


// mbb=[0,0.505,0.825,1.10,1.35,1.58,1.79,1.98,2.17,2.34,2.51,2.67,2.82,2.96,3.10,3.23,3.35,3.47,3.59,3.70,3.8];
// dbb=mbb(dif)/lbb(dif)/3./lbb(zcen)^2;
// lbb=lbb(zcen);
//ws;

 // wk;
 // window,style="win22.gs",dpi=120;
if(imin==imax){
  i=imin;
  rfact=rt(i-imin+1);
  // window,0;
  // ws1,dpi=30;
  plsys,1;
  //PL,log10(nb(,i-imin+1)+1e-10),log10(binr(zcen)*rfact);
  plg,hl(,i-imin+1),log10(binr(zcen)*rfact);
  limits,-3.,0.8,5,13;
  
  // window,1;
  // ws1,dpi=30;
  plsys,2;
  PL,log10(hd(,i-imin+1)+1e-10),log10(binr(zcen)*rfact),msize=.2;
  plg,log10(D),log10(lb*rfact),width=3;
  limits,-3.,0.8,-2,5;
  
  // window,2;
  // ws1,dpi=30;
  plsys,3;
  PL,log10(hdb(,i-imin+1)+1e-10),log10(binr(zcen)*rfact),msize=.2;
  //pleb,log10(hdb(,i-imin+1)+1e-10),log10(binr(zcen)*rfact/rt(i)),color="green",dy=log10(ehdb(,i-imin+1)+1e-10);
  plg,log10(DB),log10(lbb*rfact),color="black",width=3;
  limits,-3.,0.8,-2,5;

  // window,3;
  // ws1,dpi=30;
  plsys,4;
  PL,vd(,i-imin+1),log10(binr(zcen)*rfact),color="blue",msize=.2  ;
  //pleb,vd(,i-imin+1),log10(binr(zcen)*rfact),dy=sqrt(evd(,i-imin+1)-vd(,i-imin+1)^2)*1.;
  plg,VB,log10(lbb*rfact),color="black",width=3;
  //plg,VB,lbb,color="blue";
  limits,-3.,0.8,-2,2;
  //  limits,-1,0;
 }
              
redraw;
