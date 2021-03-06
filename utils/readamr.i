

// fp=open("dens.dat","rb");
// adress=0;
// dens=array(float,128^3);
// _read,fp,adress,dens;
// close,fp;
//dens=reform(dens,[3,128,128,128]);


func readfield(fname,&time){
  fp=open(fname,"rb");
  adress=0;
  ncells=array(int,3);
  time=array(float);
  bounds=array(float,6);
  _read,fp,adress,ncells;adress+=sizeof(ncells);
  _read,fp,adress,time;adress+=sizeof(time);
  _read,fp,adress,bounds;adress+=sizeof(bounds);
  dummy=array(float,ncells(1),ncells(2));
  field=array(float,ncells(1),ncells(2),ncells(3));
  for(i=1;i<=ncells(3);i++){
    _read,fp,adress,dummy;adress+=sizeof(dummy);
    field(,,i)=reform(dummy,[2,ncells(1),ncells(2)]);
  }
  close,fp;
  return field;
}

func readavg(fname,&time){
  fp=open(fname,"rb");
  adress=0;
  nc=array(int);
  time=array(float);
  _read,fp,adress,nc;adress+=sizeof(nc);
  _read,fp,adress,time;adress+=sizeof(time);
  nc;
  time;
  data=array(float,6);
  _read,fp,adress,data;
  close,fp;
  return data;
}


func readmap(fname,&time){
  fp=open(fname,"rb");
  adress=0;
  nmapx=array(int);
  nmapy=array(int);
  aexp=array(float);

  _read,fp,adress,nmapx;adress+=sizeof(nmapx);
  _read,fp,adress,nmapy;adress+=sizeof(nmapy);
  _read,fp,adress,aexp;adress+=sizeof(aexp);

  ntot=nmapx*nmapy;
  res=array(float,ntot*4);
  
  _read,fp,adress,res;adress+=sizeof(res);
  
  close,fp;  

  res=reform(res,[3,nmapx,nmapy,4]);
  return res;

}


func readcube(fname,&time){
  fp=open(fname,"rb");
  adress=0;
  nmapx=array(int);
  nmapy=array(int);
  nmapz=array(int);
  time=array(float);
  _read,fp,adress,nmapx;adress+=sizeof(nmapx);
  _read,fp,adress,nmapy;adress+=sizeof(nmapy);
  _read,fp,adress,nmapz;adress+=sizeof(nmapz);
  _read,fp,adress,time;adress+=sizeof(time);
  nmapx;
  nmapy;
  nmapz;
  
  map=array(float,nmapx,nmapy,nmapz);
  dummy=array(float,nmapx*nmapy);
  for(i=1;i<=nmapz;i++){
    _read,fp,adress,dummy;adress+=sizeof(dummy);
    map(,,i)=reform(dummy,[2,nmapx,nmapy]);
  }
  close,fp;

  //map=reform(map,[3,nmapx,nmapy,nmapz]);
  return map;
}


func readcell(fname,&time){
  fp=open(fname,"rb");
  adress=0;
  nc=array(int);
  time=array(float);
  _read,fp,adress,nc;adress+=sizeof(nc);
  _read,fp,adress,time;adress+=sizeof(time);
  nc;
  time;
  map=array(float,nc*5);
  _read,fp,adress,map;
  close,fp;

  map=reform(map,[2,5,nc]);
  time;
  return map;
}


func savesplit(fname,field,time){
  fp=open(fname,"wb");
  adress=0;
  nslice=int(dimsof(field)(0));
  time=float(time);
  _write,fp,adress,nslice;adress+=sizeof(nslice);
  _write,fp,adress,time;adress+=sizeof(time);
  for(i=0;i<nslice/256;i++){
    i;
    ff=float(field(,,i*256+1:(i+1)*256));
    _write,fp,adress,ff;adress+=sizeof(ff);
  }
  close,fp;
}

func readsplit(fname,&time){
  fp=open(fname,"rb");
  adress=0;
  nslice=array(int);
  time=array(float);
  _read,fp,adress,nslice;adress+=sizeof(nslice);
  _read,fp,adress,time;adress+=sizeof(time);
  write,"nslice=",nslice;
  field=array(float,nslice,nslice,nslice);
  ff=array(float,nslice,nslice,256);
  for(i=0;i<nslice/256;i++){
    i;
    _read,fp,adress,ff;adress+=sizeof(ff);
    field(,,i*256+1:(i+1)*256)=ff;
  }
  close,fp;
  return field;
}

// ====================================================================
// ====================================================================
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

func plotamr(levmap,color=,width=,lmin=,lbox=)
{
  if(is_void(lmin)) lmin=1;
  if(is_void(lbox)) lbox=1.;
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
          plbox,(xmin/lbox),(xmin+dx)/lbox,(ymin/lbox),(ymin+dx)/lbox,color=color,width=width;
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
  xi=span(-0.5,0.5,nx+1)(zcen);
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


func oct2cube(fname,lvl,field,&time,ncpu=,execut=,zmin=,zmax=,xmin=,xmax=,ymin=,ymax=,mono=,proj=,fout=,cen=,dcen=){
  if(is_void(execut)) execut="../utils/oct2grid ";
  if(is_void(ncpu)) ncpu=1;
  if(is_void(zmin)) zmin=0.;
  if(is_void(zmax)) zmax=1.;
  if(is_void(xmin)) xmin=0.;
  if(is_void(xmax)) xmax=1.;
  if(is_void(ymin)) ymin=0.;
  if(is_void(ymax)) ymax=1.;
  if(is_void(mono)) mono=-1;
  if(is_void(proj)) proj=0;
  if(is_void(fout)) fout=fname+".f"+pr1(field);
  if(!is_void(cen)){
    cen=cen(,1);
    if(is_void(dcen)) dcen=array(0.05,3);
    if(dimsof(dcen)(1)==0) dcen=array(dcen,3);
    xmin=cen(1)-dcen(1);
    xmax=cen(1)+dcen(1);
    ymin=cen(2)-dcen(2);
    ymax=cen(2)+dcen(2);
    zmin=cen(3)-dcen(3);
    zmax=cen(3)+dcen(3);
  }

  time=array(double);
  commande=execut+" "+fname+" "+pr1(lvl)+" "+pr1(field)+" "+fout+" "+pr1(ncpu)+" 0 "+pr1(mono)+" "+pr1(xmin)+" "+pr1(xmax)+" "+pr1(ymin)+" "+pr1(ymax)+" "+pr1(zmin)+" "+pr1(zmax)+" "+pr1(proj);
  //commande;
  system(commande);
  cube=readcube(fname+".f"+pr1(field),time);
  return cube;
}


func alloct(fname,lvl,field,&time,ncpu=,execut=,zmin=,zmax=,mono=){
  if(is_void(execut)) execut="../utils/alloct ";
  if(is_void(ncpu)) ncpu=1;
  if(is_void(zmin)) zmin=0.;
  if(is_void(zmax)) zmax=1.;
  if(is_void(mono)) mono=-1;
  time=array(double);
  commande=execut+" "+fname+" "+pr1(lvl)+" "+pr1(field)+" "+fname+".f"+pr1(field)+" "+pr1(ncpu)+" 0 "+pr1(mono);//+" "+pr1(zmin)+" "+pr1(zmax);
  commande;
  system(commande);
  cube=readcell(fname+".f"+pr1(field),time);
  return cube;
}



func gensilo(fname,lvl,field,&time,ncpu=,execut=,fnameout=,zmin=,zmax=,xmin=,xmax=,ymin=,ymax=,mono=){
  if(is_void(execut)) execut="~/Project/Quartz/utils/oct2grid ";
  if(is_void(ncpu)) ncpu=1;
  if(is_void(fnameout)) fnameout=fname;
  if(is_void(zmin)) zmin=0.;
  if(is_void(zmax)) zmax=1.;
  if(is_void(xmin)) xmin=0.;
  if(is_void(xmax)) xmax=1.;
  if(is_void(ymin)) ymin=0.;
  if(is_void(ymax)) ymax=1.;
  if(is_void(mono)) mono=-1;
  time=array(double);
  //  commande=execut+fname+" "+pr1(lvl)+" "+pr1(field)+" "+fnameout+" "+pr1(ncpu)+" 1";
  commande=execut+" "+fname+" "+pr1(lvl)+" "+pr1(field)+" "+fname+".f"+pr1(field)+" "+pr1(ncpu)+" 1 "+pr1(mono)+" "+pr1(xmin)+" "+pr1(xmax)+" "+pr1(ymin)+" "+pr1(ymax)+" "+pr1(zmin)+" "+pr1(zmax);
  commande;
  system(commande);
}


