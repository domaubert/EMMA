#include "utils/hop/readhop.i"

//rep="data_coarse_256_24MPC_alt_th/";isnap=24;ncpu=256;lmax=12;sbox=24.;lcoarse=8.;
//rep="./data_cosmo_corse_c0_ok_src0/";isnap=19;ncpu=32;lmax=12;sbox=12.;lcoarse=7.;
//rep="./data_cosmo_corse_c0_ok/";isnap=19;ncpu=32;lmax=12;sbox=12.;lcoarse=7.;
//rep="./data_cosmo_corse_c1_ok/";isnap=19;ncpu=32;lmax=12;sbox=12.;lcoarse=7.;
//rep="./data_4_noschaye/";istart=31;istop=31;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;
//rep="./data/";istart=16;istop=16;ncpu=32;lmax=15;sbox=12.;lcoarse=7.;
//rep="./data_12_noschaye_nosrc/";istart=21;istop=21;ncpu=32;lmax=15;sbox=12.;lcoarse=7.;
rep="./data_12_noschaye/";istart=15;istop=15;ncpu=32;lmax=15;sbox=12.;lcoarse=7.;
rep="./data_4_schaye/";istart=27;istop=27;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;
rep="./data/";istart=27;istop=27;ncpu=32;lmax=15;sbox=4.;lcoarse=7.;

#if 1
HSFR=HM200=HM200S=HR=HS=HMHOP=HMHOPS=HPOND=[];
tnew=2e7;



for(isnap=istart;isnap<=istop;isnap+=4){
#if 1
  p=mergepart(swrite(format=rep+"part.%05d",isnap),ncpu,a);
  s=mergepart(swrite(format=rep+"star.%05d",isnap),ncpu,a,star=1);
  den=read_dens(swrite(format=rep+"hop.%05d.den",isnap));
  tag=read_tag(swrite(format=rep+"reg.%05d.tag",isnap));
#endif
  z=1./a-1.;
  MPC=3.08e22; // m
  H=67; // km/s/Mpc
  h=H/100.;
  lbox=sbox/(H/100.)*MPC;// m
  lorg=lbox/MPC;
  G=6.67e-11; // S.I.
  mp=1.6e-27; // proton mass kg
  omegam=0.3175;
  omegab=0.049;
  bint=spanl(1e8,1e9,32);

  dxmin=3*(1./2^lmax);
  nhalo=max(tag);
  H=H*1e3/MPC; // s^-1
  rhoc=3.*H^2/(8.*pi*G);
  mdm=((omegam-omegab)*rhoc*lbox^3)/n/msol;
  munit=((omegam)*rhoc*lbox^3)/msol;

  write,"nhalo=",nhalo;
  ms=sfrnew=m200=mage=sfr=mhop=ph=[];
  newest=s(11,max);
  for(ih=0;ih<nhalo;ih+=1){
    www=where(tag==ih);
    if(ih%10==0) write,ih,numberof(www);

    grow,mhop,p(8,www)(sum);

    mxden=den(www)(mxx);
    cdm=p(,www)(,mxden);
    r=abs(p(1,)-cdm(1),p(2,)-cdm(2),p(3,)-cdm(3));
    rs=abs(s(1,)-cdm(1),s(2,)-cdm(2),s(3,)-cdm(3));
    www=where(r<0.01);
    sr=sort(r(www));
    rin=r(www)(sr);
  
    rhoa=dimsof(p)(0)*1.0;
    rhorin=indgen(numberof(rin))/(4./3.*pi*(rin+1e-15)^3);
    w200=where(rhorin>(rhoa*200.));
    r200=rin(w200)(0);
    //r200;

    grow,ph,[cdm];
    wsin=where(rs<r200);
    wpin=where(r<r200);
    
    if(numberof(wsin)>0){
      grow,ms,s(8,wsin)(sum);
      grow,mage,median(s(11,wsin));
      grow,sfr,[histo1d(s(11,wsin),bint,wght=s(8,wsin))*munit/bint(dif)/lorg^3];
      wnew=where(newest-s(11,wsin)<tnew);
      if(numberof(wnew)>0){
        grow,sfrnew,s(8,wsin)(wnew)(sum)*munit/tnew/lorg^3;
      }
      else{
        grow,sfrnew,[0.];
      }
      
    }
    else{
      grow,ms,[0.];
      grow,sfrnew,[0.];
      grow,mage,[0.];
      grow,sfr,bint(zcen)*0.;
    }
  
    grow,m200,p(8,wpin)(sum);
    // binr=spanl(dxmin,r200,24);
    // rho=histo1d(r,binr)/(4*pi*binr(zcen)^2)/binr(dif);
    // erho=sqrt(histo1d(r,binr))/(4*pi*binr(zcen)^2)/binr(dif);
    // plg,rho,binr(zcen);
  }


  nm=32;
  binM=spanl(1e7,5e12,nm);

  hsfr=hsfrnew=[]; // sfr per halo mass bin
  for(i=1;i<=dimsof(sfr)(2);i++){
    grow,hsfr,[histo1d(m200*munit,binM,wght=sfr(i,))/(histo1d(m200*munit,binM)+1e-15)];
  }


  www=where(ms>0); // looking for luminous haloes
  hm200=histo1d(m200*munit*h,binM)/binM(dif)/sbox^3; // halo mass histogram
  hm200s=histo1d(m200(www)*munit*h,binM)/binM(dif)/sbox^3; // luminous haloes mass histogram
  hmhop=histo1d(mhop*munit*h,binM)/binM(dif)/sbox^3; // halo mass histogram
  hmhops=histo1d(mhop(www)*munit*h,binM)/binM(dif)/sbox^3; // luminous haloes mass histogram
  hr=histo1d(mhop*munit*h,binM,wght=ms/mhop)/(histo1d(mhop*munit*h,binM)+1e-15); // halo ms/ms mass histogram
  hsnew=histo1d(mhop*munit*h,binM,wght=sfrnew)/(histo1d(mhop*munit*h,binM)+1e-15); // halo ms/ms mass histogram
  hpond=histo1d(mhop*munit*h,binM,wght=ms*munit*h)/(histo1d(mhop*munit*h,binM)+1e-15);
  
  
  grow,HSFR,[hsfr];
  grow,HM200,[hm200];
  grow,HM200S,[hm200s];
  grow,HMHOP,[hmhop];
  grow,HMHOPS,[hmhops];
  grow,HR,[hr];
  grow,HS,[hsnew];
  grow,HPOND,[hpond];
 }
#endif

#if 1

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






