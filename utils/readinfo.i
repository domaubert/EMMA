struct AMRINFO{
  int lcoarse;
  int lmax;
  int lmap;
  int nmax;
  float thresh;
  pointer t,dt,res;
};

func readinfo(fname){
  ff=open(fname,"r");
  lcoarse=lmax=lmap=nmax=array(int);
  thresh=array(float);
  n=read(ff,format="levelcoarse=%d",lcoarse);
  n=read(ff,format="levelmax=%d",lmax);
  n=read(ff,format="levelmap=%d",lmap);
  n=read(ff,format="nstepmax=%d",nmax);
  n=read(ff,format="amr threshold=%e",thresh);

  idx=array(int);
  t=array(float,nmax);
  dt=array(float,nmax);
  res=array(float,nmax);

  for(i=1;i<=nmax;i++){
    n=read(ff,format="istep=%d t=%e dt=%e coarseres=%e",idx,t(i),dt(i),res(i));
  }

  close,ff;

  s=AMRINFO();
  s.lcoarse=lcoarse;
  s.lmax=lmax;
  s.lmap=lmap;
  s.thresh=thresh;
  s.nmax=nmax;
  s.t=&t;
  s.dt=&dt;
  s.res=&res;

  return s;
}
