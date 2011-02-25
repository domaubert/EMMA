

func readpart(fname){
  //fname=exec("ls -d data/part.*")
  //x=y=[];for(i=1;i<=numberof(fname);i++){phase=readpart(fname(i));www=where(phase(7,)==1);grow,x,phase(1,www);grow,y,phase(2,www);}
  fp=open(fname,"rb");
  adress=0;
  npart=array(int);
  _read,fp,adress,npart;adress+=sizeof(npart);
  if(npart!=0){
  phase=array(float,7,npart);
  _read,fp,adress,phase;adress+=sizeof(phase);
  write,"found "+pr1(npart)+" particles";
  close,fp;
  return phase;
  }
  else{
    write,"No particles : empty file!";
  return 0;
  }
}



func mergepart(fname,ncpu){
  pf=[];
  for(i=0;i<ncpu;i++){
    fname2=swrite(format=fname+".p%05d",i);
    p=readpart(fname2);
    if(numberof(p)!=1) grow,pf,p;
  }
  return pf;
}

