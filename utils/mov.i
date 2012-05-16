for(i=1;i<=357;i++){
  d=mergefield(swrite(format="data/denhyd.%05d",i),"data/cpu3d.00000",4,7,a);
  l=mergefield(swrite(format="data/lev3d.%05d",i),"data/cpu3d.00000",4,7,a);
  fma;
  plotamr(l(,64,));
  hcps(swrite(format="data/snap2.%05d",i));
  system("convert -density 256x256 "+swrite(format="data/snap2.%05d.ps",i)+ " "+swrite(format="data/snap2.%05d.jpg",i));
 }
