

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

func plotamr(levmap,color=,width=)
{
  nx=dimsof(levmap)(2);
  ny=dimsof(levmap)(3);

  lcur=max(levmap);
  dx=1;
  imod=1;
  while(lcur>=1){
    lcur;
    www=where2(levmap==lcur);
    if(numberof(www)!=0){
      for(i=1;i<dimsof(www)(0);i++){
        xmin=www(1,i)-1;
        ymin=www(2,i)-1;
        if((xmin%dx==0)*(ymin%dx==0)){
          plbox,xmin,xmin+dx,ymin,ymin+dx,color=color,width=width;
        }
      }
    }
    lcur=lcur-1;
    dx*=2;
  }
  
}

