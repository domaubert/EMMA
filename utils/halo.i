#include "~/qtest/utils/hop/readhop.i"

dur_src=[];
//rep="./data_4_new_wsrc_ministar/";istart=32;istop=32;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;
//rep="./data_4_new_wsrc_ministar_x1/";istart=37;istop=37;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;
//rep="~/qtest/data_4_new_wsrc_ministar_x3/";istart=37;istop=37;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;dur_src=30.;flux_src=3e16;
//rep="./data_4_new_wsrc_ministar_x3_mono/";istart=37;istop=37;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;
//rep="./data_4_new_wsrc_ministar_x3_mono_vb/";istart=11;istop=37;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;
//rep="./data_4_new_wsrc_ministar_x10/";istart=37;istop=37;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;
//rep="./data_4_new_wsrc_ministar_x100/";istart=35;istop=35;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;
//rep="~/qtest/data_3MYR/";istart=12;istop=12;ncpu=32;lmax=15;sbox=12.;lcoarse=7.;dur_src=3.;flux_src=3e16;
//rep="~/qtest/data/";istart=36;istop=36;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;dur_src=5.;flux_src=1e16;
rep="~/qtest/data/";istart=36;istop=36;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;dur_src=5.;flux_src=1e16;
rep="~/qtest/data_8_256/";istart=14;istop=14;ncpu=32;lmax=15;sbox=8.;lcoarse=8.;dur_src=10.;flux_src=2e16;
//rep="~/qtest/data_4_new_wsrc_ministar_x1/";istart=37;istop=37;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;dur_src=20.;flux_src=5e15;
//rep="~/qtest/data_4_new_wsrc_ministar_x3/";istart=37;istop=37;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;dur_src=20.;flux_src=1.5e16;
//rep="~/qtest/data_4_new_wsrc_ministar_x10/";istart=37;istop=37;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;dur_src=20.;flux_src=5e16;
//rep="~/qtest/data_4_new_wsrc_ministar_x100/";istart=35;istop=35;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;dur_src=20.;flux_src=5e17;

restore,openb("~/qtest/utils/mag/f1600_SB99_5Myrs.yor"),f1600,age;
f1600*=flux_src;

#if 1
HSFR=HM200=HM200S=HR=HS=HMHOP=HMHOPS=HPOND=vz=HSTOT=HFB=HFB2=HL=[];
tnew=dur_src*1e6;



for(isnap=istart;isnap<=istop;isnap+=1){
#if 1
  d=alloct(swrite(format=rep+"/grid.%05d",isnap),13,101,a,ncpu=ncpu,execut="~/qtest/utils/alloct ");
  xd=d(:3,);
  zz=1./a-1.;
  zz;

  tsnap= univAge(1000,h0=67.,Omega_l=0.6825,Omega_m=0.3175)-univAge(zz,h0=67.,Omega_l=0.6825,Omega_m=0.3175);
  tsnap/=(3600.*24*365*1e6);
  p=mergepart(swrite(format=rep+"part.%05d",isnap),ncpu,a);
  s=mergepart(swrite(format=rep+"star.%05d",isnap),ncpu,a,star=1);
  den=read_dens(swrite(format=rep+"hop.%05d.den",isnap));
  tag=read_tag(swrite(format=rep+"reg.%05d.tag",isnap));
#endif
  MPC=3.08e22; // m
  H=67; // km/s/Mpc
  h=H/100.;
  lbox=sbox/(H/100.)*MPC;// m
  lorg=lbox/MPC;
  G=6.67e-11; // S.I.
  mp=1.6e-27; // proton mass kg
  msol=2e30;
  omegam=0.3175;
  omegab=0.049;
  bint=spanl(1e8,1e9,32);

  dxmin=3*(1./2^lmax);
  nhalo=max(tag);
  H=H*1e3/MPC; // s^-1
  rhoc=3.*H^2/(8.*pi*G);
  n=dimsof(p)(0);
  mdm=((omegam-omegab)*rhoc*lbox^3)/n/msol;
  munit=((omegam)*rhoc*lbox^3)/msol;

  write,"nhalo=",nhalo;
  ms=sfrnew=m200=mage=sfr=mhop=ph=m200d=m200d2=m2002=l200d=mag1600=[];
  newest=tsnap;
  for(ih=0;ih<nhalo;ih+=1){
    www=where(tag==ih);
    if(ih%10==0) write,ih,numberof(www);

    grow,mhop,p(8,www)(sum);

    mxden=den(www)(mxx);
    cdm=p(,www)(,mxden);
    r=abs(p(1,)-cdm(1),p(2,)-cdm(2),p(3,)-cdm(3));
    rd=abs(xd(1,)-cdm(1),xd(2,)-cdm(2),xd(3,)-cdm(3));

    
    if(numberof(s)>0){
      rs=abs(s(1,)-cdm(1),s(2,)-cdm(2),s(3,)-cdm(3));
    }
    else{
      rs=[];
    }

    www=where(r<0.01);
    sr=sort(r(www));
    rin=r(www)(sr);
  
    rhoa=dimsof(p)(0)*1.0;
    rhorin=indgen(numberof(rin))/(4./3.*pi*(rin+1e-15)^3);
    w200=where(rhorin>(rhoa*200.));
    r200=rin(w200)(0);

    if(numberof(s)>0){
      wsin=where(rs<r200);
    }
    else{
      wsin=[];
    }



    //r200;

    grow,ph,[cdm];
    wpin=where(r<r200);
    wpin2=where(r<1.0*r200);
    wpind=where(rd<r200);
    wpind2=where(rd<1.0*r200);

    if(numberof(wpind)>0){
      mdens=(d(5,wpind)*pow(0.5,3.*d(4,wpind)))(sum);
      levin=d(4,wpind)(avg);
      if(levin==0) error;
    }
    else{
      mdens=0.;
      levin=0.;
    }

    
    if(numberof(wpind2)>0){
      mdens2=((d(5,wpind2)>500*omegab/omegam)*d(5,wpind2)*pow(0.5,3.*d(4,wpind2)))(sum);
    }
    else{
      mdens2=0.;
    }


    
    if(numberof(wsin)>0){

      flux_in=(s(8,wsin)*munit/1e6*interp(f1600,age,(newest*1e6-s(11,wsin))/1e6))(sum);
      M1600loc=-2.5*log10(flux_in)-48.6;
      grow,mag1600,M1600loc;
      grow,ms,s(8,wsin)(sum);
      grow,mage,median(s(11,wsin));
      grow,sfr,[histo1d(s(11,wsin),bint,wght=s(8,wsin))*munit/bint(dif)/lorg^3];
      wnew=where((newest*1e6-s(11,wsin))<tnew);
      if(numberof(wnew)>0){
        grow,sfrnew,s(8,wsin)(wnew)(sum)*munit/tnew/lorg^3;
      }
      else{
        grow,sfrnew,[0.];
      }
      
    }
    else{
      grow,mag1600,[0.];
      grow,ms,[0.];
      grow,sfrnew,[0.];
      grow,mage,[0.];
      grow,sfr,[bint(zcen)*0.];
    }
  
    grow,m200,p(8,wpin)(sum);
    grow,m2002,p(8,wpin2)(sum);
    grow,m200d,mdens;
    grow,m200d2,mdens2;
    grow,l200d,levin;
    // binr=spanl(dxmin,r200,24);
    // rho=histo1d(r,binr)/(4*pi*binr(zcen)^2)/binr(dif);
    // erho=sqrt(histo1d(r,binr))/(4*pi*binr(zcen)^2)/binr(dif);
    // plg,rho,binr(zcen);
  }


  nm=32;
  binM=spanl(10*mdm,1e11,nm);

  hsfr=hsfrnew=[]; // sfr per halo mass bin
  for(i=1;i<=dimsof(sfr)(2);i++){
    grow,hsfr,[histo1d(m200*munit,binM,wght=sfr(i,))/(histo1d(m200*munit,binM)+1e-15)];
  }


  www=where(ms>0); // looking for luminous haloes
  if(numberof(www)>0){
    hm200s=histo1d(m200(www)*munit*h,binM)/binM(dif)/sbox^3; // luminous haloes mass histogram
    hmhops=histo1d(mhop(www)*munit*h,binM)/binM(dif)/sbox^3; // luminous haloes mass histogram
  }
  else{
    h200s=binM(dif)*0.;
    hmhops=binM(dif)*0.;
  }
  
  hm200=histo1d(m200*munit*h,binM)/binM(dif)/sbox^3; // halo mass histogram
  hmhop=histo1d(mhop*munit*h,binM)/binM(dif)/sbox^3; // halo mass histogram
  
  hr=histo1d(mhop*munit*h,binM,wght=ms/mhop)/(histo1d(mhop*munit*h,binM)+1e-15); // halo ms/ms mass histogram
  hsnew=histo1d(mhop*munit*h,binM,wght=sfrnew)/(histo1d(mhop*munit*h,binM)+1e-15); // halo ms/ms mass histogram
  hsnewtot=histo1d(mhop*munit*h,binM,wght=sfrnew); // halo ms/ms mass histogram
  hpond=histo1d(mhop*munit*h,binM,wght=ms*munit*h)/(histo1d(mhop*munit*h,binM)+1e-15);

  w0=where(m200d>0);
  if(numberof(w0)>0){
    hfb=histo1d(mhop(w0)*munit*h,binM,wght=m200d(w0)/(m200(w0)+m200d(w0)))/(histo1d(mhop(w0)*munit*h,binM)+1e-15); // halo ms/ms mass histogram
    hl=histo1d(mhop(w0)*munit*h,binM,wght=l200d(w0))/(histo1d(mhop(w0)*munit*h,binM)+1e-15); // halo ms/ms mass histogram
  }
  else{
    hfb=binM(zcen)*0;
    hl=binM(zcen)*0;
  }
  
    w0=where(m200d>0);
  
    if(numberof(w0)>0){
      hfb2=histo1d(mhop(w0)*munit*h,binM,wght=m200d2(w0)/(m2002(w0)+m200d2(w0)))/(histo1d(mhop(w0)*munit*h,binM)+1e-15); // halo ms/ms mass histogram
    }
    else{
      hfb2=binM(zcen)*0;
    }
         

  
  grow,HSFR,[hsfr];
  grow,HM200,[hm200];
  grow,HM200S,[hm200s];
  grow,HMHOP,[hmhop];
  grow,HMHOPS,[hmhops];
  grow,HR,[hr];
  grow,HS,[hsnew];
  grow,HSTOT,[hsnewtot];
  grow,HPOND,[hpond];
  grow,HFB,[hfb];
  grow,HFB2,[hfb2];
  grow,HL,[hl];
  grow,vz,zz;
 }
#endif

#if 0

nh=numberof(sfrnew);
posh=ph(:3,);
bind=spanl(0.01,0.7,32);

HD=HDLUM=CL=CLL=[];
for(i=1;i<=nh;i++){
  if(i%10==0) write,i;
  p0=posh(,i)(,-:1:nh);
  pp=(posh-p0)-((posh-p0)>0.5)+((posh-p0)<-0.5);
  dist=abs(pp(1,),pp(2,),pp(3,));
  grow,HD,[histo1d(dist,bind)/(4.*pi*bind(zcen)^2*bind(dif))];
  grow,HDLUM,[histo1d(dist,bind,wght=(sfrnew>0))/(4.*pi*bind(zcen)^2*bind(dif))];


  dist(where(dist==0))=1e6;
  closest=dist(sort(dist))(1);
  dist(where(sfrnew==0))=1e6;
  closestL=dist(sort(dist))(1);
  
  grow,CL,closest;
  grow,CLL,closestL;
 }






