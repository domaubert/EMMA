

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

func readcube(fname,&time){
  fp=open(fname,"rb");
  adress=0;
  nmap=array(int);
  time=array(float);
  _read,fp,adress,nmap;adress+=sizeof(nmap);
  _read,fp,adress,time;adress+=sizeof(time);
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

func mergefield(fname,cname,ncpu,level,&time){

  field=array(float,2^level,2^level,2^level);
  c=mapcpu(cname,ncpu);
  time=array(float);
  for(i=0;i<ncpu;i++)
    {
      fname2=swrite(format=fname+".p%05d",i);
      f=readcube(fname2,time);
      www=where(c==i);
      field(www)=f(www);
    }
  
    return field;

}

func plotamr(levmap,color=,width=,lmin=)
{
  if(is_void(lmin)) lmin=1;
  nx=dimsof(levmap)(2);
  ny=dimsof(levmap)(3);

  lmap=int(log(nx)/log(2));
  lcur=max(levmap);
  dx=2^(lmap-lcur);
  imod=1;
  while(lcur>=lmin){
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

func ext_amr1D(levmap,field,color=,lmin=)
{
  if(is_void(lmin)) lmin=min(levmap);
  nx=dimsof(levmap)(2);
  //  ny=dimsof(levmap)(3);

  lmap=int(log(nx)/log(2));
  lcur=max(levmap);
  dx=2^(lmap-lcur);
  imod=1;
  xi=span(0,1,nx+1)(zcen);
  y=x=lv=[];
  while(lcur>=lmin){
    lcur;
    www=where(levmap==lcur);
    if(numberof(www)!=0){
      xloc=xi(www)(::int(pow(2,max(levmap)-lcur)));
      yloc=field(www)(::int(pow(2,max(levmap)-lcur)));
      lvloc=levmap(www)(::int(pow(2,max(levmap)-lcur)));
      grow,x,xloc;
      grow,y,yloc;
      grow,lv,lvloc;
    }
    lcur=lcur-1;
    dx*=2;
  }
  s=sort(x);
  return [x(s),y(s),lv(s)];
  
}


func oct2cube(fname,lvl,field,&time,ncpu=,execut=){
  if(is_void(execut)) execut="~/Project/Quartz/utils/oct2grid ";
  if(is_void(ncpu)) ncpu=1;
  time=array(double);
  commande=execut+fname+" "+pr1(lvl)+" "+pr1(field)+" "+fname+".f"+pr1(field)+" "+pr1(ncpu);
  system(commande);
  cube=readcube(swrite(format=fname+".f"+pr1(field)+".p%05d",ncpu-1),time);
  return cube;
                                                                                             
}


