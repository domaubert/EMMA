


//rep="./data/";ncpu=32;execut="./utils/oct2grid ";sbox=12.;isnap=19;
//rep="./data_4_noschaye/";ncpu=32;execut="./utils/oct2grid ";sbox=4.;isnap=32;
rep="~/qtest/data/";ncpu=32;execut="~/qtest/utils/alloct ";sbox=12.5;isnap=13;
//rep="./data_4_new_wsrc/";ncpu=32;execut="./utils/oct2grid ";sbox=4.;isnap=31;
//rep="./data_4_new_wsrc_ministar_x100/";ncpu=32;execut="./utils/oct2grid ";sbox=4.;isnap=30;


level=13;

bind=spanl(1e-2,1e6,64);
bint=spanl(1e1,1e6,64);
binx=spanl(1e-7,1e0,64);

t=alloct(swrite(format=rep+"grid.%05d",isnap),15,707,a,ncpu=32,execut=execut)(5,);
d=alloct(swrite(format=rep+"grid.%05d",isnap),15,101,a,ncpu=32,execut=execut);
x=alloct(swrite(format=rep+"grid.%05d",isnap),15,706,a,ncpu=32,execut=execut)(5,);
p=alloct(swrite(format=rep+"grid.%05d",isnap),15,105,a,ncpu=32,execut=execut)(5,);
l=d(4,);
d=d(5,);

# if 1
hphase=histo2d([d(*),t(*)],bind,bint)/(bind(dif)(,-:1:numberof(bint)-1))/(bint(dif)(-:1:numberof(bind)-1,));
hphasem=histo2d([d(*),t(*)],bind,bint,wght=d(*)/pow(2.,3*l(*)))/(bind(dif)(,-:1:numberof(bint)-1))/(bint(dif)(-:1:numberof(bind)-1,));
hphasexm=histo2d([d(*),1.-x(*)],bind,binx,wght=d(*)/pow(2.,3*l(*)))/(bind(dif)(,-:1:numberof(bint)-1))/(binx(dif)(-:1:numberof(bind)-1,));
hphasex=histo2d([d(*),1.-x(*)],bind,binx)/(bind(dif)(,-:1:numberof(bint)-1))/(binx(dif)(-:1:numberof(bind)-1,));


#if 0
hphaseUV=histo2d([d(*),n(*)],bind,bind)/(bind(dif)(,-:1:numberof(bind)-1))/(bind(dif)(-:1:numberof(bind)-1,));
hphasemUV=histo2d([d(*),n(*)],bind,bind,wght=d(*))/(bind(dif)(,-:1:numberof(bind)-1))/(bind(dif)(-:1:numberof(bind)-1,));
#endif

window,4,style="win12.gs";
plsys,1;
plfc,log10(hphasem+1e-10),bint(zcen)(-:1:numberof(bind)-1,),bind(zcen)(,-:1:numberof(bint)-1);
logxy,1,1;
plsys,2;
plfc,log10(hphase+1e-10),bint(zcen)(-:1:numberof(bind)-1,),bind(zcen)(,-:1:numberof(bint)-1);
//plfc,log10(hphasemUV+1e-10),bind(zcen)(-:1:numberof(bind)-1,),bind(zcen)(,-:1:numberof(bind)-1);
logxy,1,1;

window,3,style="win12.gs";
plsys,1;
plfc,log10(hphasexm+1e-10),binx(zcen)(-:1:numberof(bind)-1,),bind(zcen)(,-:1:numberof(binx)-1);
logxy,1,1;
plsys,2;
plfc,log10(hphasex+1e-10),binx(zcen)(-:1:numberof(bind)-1,),bind(zcen)(,-:1:numberof(binx)-1);
//plfc,log10(hphasemUV+1e-10),bind(zcen)(-:1:numberof(bind)-1,),bind(zcen)(,-:1:numberof(bind)-1);
logxy,1,1;

