
func read_grafic(filename,&L,&h,nslice=,bigendian=,decendian=)
{

  
  np1=np2=np3=array(int);
  dx=x1o=x2o=x3o=astart=omegam=omegav=h0=array(float);

  adress=0;
  ff=open(filename,"rb");

  if(bigendian) sun_primitives,ff;
  if(decendian) dec_primitives,ff;
  
  expsize=3*sizeof(int)+8*sizeof(float);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,np1;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,np2;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,np3;adress+=sizeof(dummy);
  dummy=array(float);_read,ff,adress,dx;adress+=sizeof(dummy);
  dummy=array(float);_read,ff,adress,x1o;adress+=sizeof(dummy);
  dummy=array(float);_read,ff,adress,x2o;adress+=sizeof(dummy);
  dummy=array(float);_read,ff,adress,x3o;adress+=sizeof(dummy);
  dummy=array(float);_read,ff,adress,astart;adress+=sizeof(dummy);
  dummy=array(float);_read,ff,adress,omegam;adress+=sizeof(dummy);
  dummy=array(float);_read,ff,adress,omegav;adress+=sizeof(dummy);
  dummy=array(float);_read,ff,adress,h0;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

  if(dummy!=expsize) error,"I/O header error";
  if(!is_void(nslice)) np3=nslice;
  
  write,"*********";
  write,"Cosmology";
  write,"*********";
  write,"omegam=",omegam;
  write,"omegav=",omegav;
  write,"h0=",h0;
  write,"astart=",astart;
  write,"zstart=",1./astart-1.;
  write,"dx=",dx;
  field=array(float,np1,np2,np3);
  
  for(i=1;i<=np3;i++)
    {
      dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
      dummy=array(float,np1,np2);_read,ff,adress,dummy;adress+=sizeof(dummy);
      field(,,i)=dummy;
      dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
    }


  close,ff;

  L=dx*np1*h0/100.;
  h=h0/100.;
  return field;
  
}


func grafic2stef(filename,dd)
{
  ff=open(filename,"wb");

  n=array(int(dimsof(dd)(2)),3);
  adress=0;
  dummy=int(sizeof(n));_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=n;_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=int(sizeof(n));_write,ff,adress,dummy;adress+=sizeof(dummy);

  dummy=int(sizeof(dd(*)));_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=dd(*);_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=int(sizeof(dd(*)));_write,ff,adress,dummy;adress+=sizeof(dummy);
  close,ff;

}

func write_stef(filename,vec)
{
  vec=float(vec);
  
  ff=open(filename,"wb");
  adress=0;
  n=int(dimsof(vec)(2:4));
  
  dummy=int(sizeof(n));_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=n(1);_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=n(2);_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=n(3);_write,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=int(sizeof(n));_write,ff,adress,dummy;adress+=sizeof(dummy);

  dummy=int(sizeof(vec));_write,ff,adress,dummy;adress+=sizeof(dummy);
  _write,ff,adress,vec;adress+=sizeof(vec);
  dummy=int(sizeof(vec));_write,ff,adress,dummy;adress+=sizeof(dummy);

  close,ff;

  
}



func read_noise(filename,nslice=)
{
  np1=np2=np3=iseed=array(int);
  dx=x1o=x2o=x3o=astart=omegam=omegav=h0=array(float);

  adress=0;
  ff=open(filename,"rb");

  expsize=4*sizeof(int);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,np1;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,np2;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,np3;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,iseed;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

  //  if(dummy!=expsize) error,"I/O header error";
  if(!is_void(nslice)) np3=nslice;
  
  write,"*********";
  write,"Cosmology";
  write,"*********";
  write,"omegam=",omegam;
  write,"omegav=",omegav;
  write,"h0=",h0;
  write,"astart=",astart;

  field=array(float,np1,np2,np3);
  
  for(i=1;i<=np3;i++)
    {
      dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
      dummy=array(float,np1,np2);_read,ff,adress,dummy;adress+=sizeof(dummy);
      field(,,i)=dummy;
      dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
    }


  close,ff;

  return field;
  
}

