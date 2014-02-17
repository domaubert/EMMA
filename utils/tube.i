

vv=load("utils/test1.dat");
plg,vv(,2),vv(,1);
#if 0
imax=22;
imin=1;
dt=t=[];
for(i=imin;i<=imax;i++){
  d=oct2cube(swrite(format="data/grid.%05d",i),9,101,a,ncpu=64,execut="utils/oct2grid ",zmin=0,zmax=0.05);
  grow,dt,d(,1,1);
  grow,t,a;
 }
