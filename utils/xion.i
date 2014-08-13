rep="../data";

nsnap=19;
level=7;

MPC=3.08e22; // m
H=67; // km/s/Mpc
Hor=H;
G=6.67e-11; // S.I.
omegam=0.3175;
omegab=0.045;
n=128^3;
lbox=12./(H/100.)*MPC;// m
lorg=lbox/MPC;


// stellar particle mass
ncell=pow(2,level);
xion=array(float,ncell,ncell,ncell,nsnap);
temp=array(float,ncell,ncell,ncell,nsnap);
va=array(double,nsnap);
vt=array(double,nsnap);

for(i=1;i<=nsnap;i++){
  xion(,,,i)=oct2cube(swrite(format=rep+"/grid.%05d",i),level,706,a,ncpu=8,execut="../utils/oct2grid ")
  temp(,,,i)=oct2cube(swrite(format=rep+"/grid.%05d",i),level,707,a,ncpu=8,execut="../utils/oct2grid ")
  va(i)=a;
  vt(i)=(univAge(10000.,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam)-univAge(1./a-1.,h0=Hor,Omega_m=omegam,Omega_l=1.-omegam))/(3600.*24.*365);
 }
