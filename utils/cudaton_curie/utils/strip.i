require,"utils/cudaton.i";

ll=exec("ls -d fields/snapaton*");
for(i=1;i<=numberof(ll);i++){
  newf=swrite(format="fields/snaplight_%03d",i);
  newf;
  s=readsnap(ll(i),keep=1);

  ff=open(newf,"wb");
  adress=0;
  _write,ff,adress,s.ncells;adress+=sizeof(s.ncells);
  _write,ff,adress,s.nsource;adress+=sizeof(s.nsource);
  _write,ff,adress,s.time;adress+=sizeof(s.time);

  dummy=*s.src;_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=*s.dens;_write,ff,adress,dummy;adress+=sizeof(dummy);
  _write,ff,adress,s.aexp;adress+=sizeof(s.aexp);

  close,ff;
 }
