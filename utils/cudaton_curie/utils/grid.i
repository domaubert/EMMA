func dumpgrid(fname,x,t,d,a)
{
  n=int(dimsof(x)(0));
  adress=0;
  ff=open(fname,"wb");
  _write,ff,adress,a;adress+=sizeof(a);
  _write,ff,adress,n;adress+=sizeof(n);
  for(i=1;i<=(n/256);i++) {dummy=x(,,(i-1)*256+1:i*256);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  for(i=1;i<=(n/256);i++) {dummy=t(,,(i-1)*256+1:i*256);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  for(i=1;i<=(n/256);i++) {dummy=d(,,(i-1)*256+1:i*256);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  //  _write,ff,adress,t;adress+=sizeof(t);
  //  _write,ff,adress,d;adress+=sizeof(d);
  
//   for(i=1;i<=n;i++) {dummy=t(,,i);_write,ff,adress,dummy;adress+=sizeof(dummy);}
//   for(i=1;i<=n;i++) {dummy=d(,,i);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  close,ff;
  return 0;
}

func dumpgrid1(fname,x,a)
{
  n=int(dimsof(x)(0));
  adress=0;
  ff=open(fname,"wb");
  _write,ff,adress,a;adress+=sizeof(a);
  _write,ff,adress,n;adress+=sizeof(n);
  for(i=1;i<=(n/256);i++) {dummy=x(,,(i-1)*256+1:i*256);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  //  for(i=1;i<=(n/256);i++) {dummy=t(,,(i-1)*256+1:i*256);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  //  for(i=1;i<=(n/256);i++) {dummy=d(,,(i-1)*256+1:i*256);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  //  _write,ff,adress,t;adress+=sizeof(t);
  //  _write,ff,adress,d;adress+=sizeof(d);
  
//   for(i=1;i<=n;i++) {dummy=t(,,i);_write,ff,adress,dummy;adress+=sizeof(dummy);}
//   for(i=1;i<=n;i++) {dummy=d(,,i);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  close,ff;
  return 0;
}

func readgrid(fname,&x,&t,&d,&a)
{
  adress=0;
  ff=open(fname,"rb");
  a=array(float);
  n=array(int);
  _read,ff,adress,n;adress+=sizeof(n);
  _read,ff,adress,a;adress+=sizeof(a);
  dummy=array(float,n,n,n/(n/256));
  x=array(float,n,n,n);
  t=array(float,n,n,n);
  d=array(float,n,n,n);
  for(i=1;i<=(n/256);i++) {i;_read,ff,adress,dummy;adress+=sizeof(dummy);x(,,(i-1)*256+1:i*256)=dummy;}
    for(i=1;i<=(n/256);i++) {i;_read,ff,adress,dummy;adress+=sizeof(dummy);t(,,(i-1)*256+1:i*256)=dummy;}
    for(i=1;i<=(n/256);i++) {i;_read,ff,adress,dummy;adress+=sizeof(dummy);d(,,(i-1)*256+1:i*256)=dummy;}
//   for(i=1;i<=(n/256);i++) {_read,ff,adress,dummy;adress+=sizeof(dummy);}
  //  _write,ff,adress,t;adress+=sizeof(t);
  //  _write,ff,adress,d;adress+=sizeof(d);
  
//   for(i=1;i<=n;i++) {dummy=t(,,i);_write,ff,adress,dummy;adress+=sizeof(dummy);}
//   for(i=1;i<=n;i++) {dummy=d(,,i);_write,ff,adress,dummy;adress+=sizeof(dummy);}
  close,ff;
  return 0;
}

func readgridd(fname,&d,&a)
{
  adress=0;
  ff=open(fname,"rb");
  a=array(float);
  n=array(int);
  _read,ff,adress,a;adress+=sizeof(a);
  _read,ff,adress,n;adress+=sizeof(n);
  dummy=array(float,n,n,n/(n/256));
  d=array(float,n,n,n);
  for(i=1;i<=(n/256);i++) {adress+=sizeof(dummy);}
  for(i=1;i<=(n/256);i++) {adress+=sizeof(dummy);}
  for(i=1;i<=(n/256);i++) {_read,ff,adress,dummy;adress+=sizeof(dummy);d(,,(i-1)*256+1:i*256)=dummy;}
  close,ff;
  return 0;
}

func readgrid1(fname,&x,&a)
{
  adress=0;
  ff=open(fname,"rb");
  a=array(float);
  n=array(int);
  _read,ff,adress,a;adress+=sizeof(a);
  _read,ff,adress,n;adress+=sizeof(n);
  dummy=array(float,n,n,n/(n/256));
  x=array(float,n,n,n);
  for(i=1;i<=(n/256);i++) {_read,ff,adress,dummy;adress+=sizeof(dummy);x(,,(i-1)*256+1:i*256)=dummy;}
  close,ff;
  return 0;
}


func readg(fname,&x,&t,&a)
{
  adress=0;
  ff=open(fname,"rb");
  a=array(float);
  n=array(int);
  _read,ff,adress,n;adress+=sizeof(n);
  _read,ff,adress,a;adress+=sizeof(a);
  x=array(float,n,n,n);
  t=array(float,n,n,n);
  _read,ff,adress,x;adress+=sizeof(x);
  _read,ff,adress,t;
  x=x(17:n-16,17:n-16,17:n-16);
  t=t(17:n-16,17:n-16,17:n-16);
  close,ff;

}
