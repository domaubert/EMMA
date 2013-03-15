require,"utils/readamr.i";

z=t=[];
for(i=0;i<=10;i++){

  lvl=6;
  rep="data/";


 d1=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,1,a);
 d=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,101,a);
 p=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,105,a);
 l=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,0,a);
 tt=oct2cube(swrite(format=rep+"grid.%05d",i),lvl,707,a);

 grow,t,tt(*)(avg);
 grow,z,1./a-1.;
 }
