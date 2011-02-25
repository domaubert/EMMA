

// fp=open("dens.dat","rb");
// adress=0;
// dens=array(float,128^3);
// _read,fp,adress,dens;
// close,fp;
//dens=reform(dens,[3,128,128,128]);

func readmap(fname){
  fp=open(fname,"rb");
  adress=0;
  nmap=array(int);
  _read,fp,adress,nmap;adress+=sizeof(nmap);
  map=array(float,nmap^2);
  _read,fp,adress,map;
  close,fp;

  map=reform(map,[2,nmap,nmap]);
  return map;
}

func readcube(fname){
  fp=open(fname,"rb");
  adress=0;
  nmap=array(int);
  _read,fp,adress,nmap;adress+=sizeof(nmap);
  map=array(float,nmap^3);
  _read,fp,adress,map;
  close,fp;

  map=reform(map,[3,nmap,nmap,nmap]);
  return map;
}

func mapcpu(fname,ncpu){

  for(i=0;i<ncpu;i++)
    {
      fname2=swrite(format=fname+".p%05d",i);
      c=readcube(fname2);
      if(i==0) cglob=c*0.;
      www=where(c==i);
      cglob(www)=i;
    }
  
  return cglob;
}

func mergefield(fname,cname,ncpu,level){

  field=array(float,2^level,2^level,2^level);
  c=mapcpu(cname,ncpu);

  for(i=0;i<ncpu;i++)
    {
      fname2=swrite(format=fname+".p%05d",i);
      f=readcube(fname2);
      www=where(c==i);
      field(www)=f(www);
    }
  
    return field;

}
