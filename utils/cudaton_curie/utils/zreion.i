
z=array(-1.,128,128,128);

for(i=0;i<=80;i++){
  i;
  s=readsnap(swrite(format="runs2/out/snap.%05d",i));
  www=where((z<0)*((*s.xion)>0.5)*((*s.dens)>2e-2));
  if(numberof(www)>0){
    z(www)=1./(s.aexp)-1.;
  }
 }

