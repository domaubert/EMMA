//#include "ascii2ls.i"

func read_size(sizefile)
{
  ff=open(sizefile,"r");
  npart=array(int);
  npartingroup=array(int);
  ngroup=array(int);
  read,ff,npart;
  read,ff,npartingroup;
  read,ff,ngroup;

  write,"npart,npart in groups, ngroups =",npart,npartingroup,ngroup;

  vsize=array(int,2,ngroup);
  read,ff,vsize;
  close,ff;

  return vsize(2,);
    
  
}

func read_tag(tagfile,verbose=)
{
  if(is_void(verbose)) verbose=0;
  ff=open(tagfile,"rb");

  if(verbose) write,"Reading "+tagfile;
  adress=0;
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  if(verbose)  write,"Npart in simulation=",dummy;
  Npart=dummy;
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  if(verbose)  write,"N halos found in simulation=",dummy;
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

  if(verbose) write,"Tagging Particles";
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  tagarray=array(int,Npart);_read,ff,adress,tagarray;
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

  close,ff;
  return tagarray;
}

func read_dens(tagfile,verbose=)
{
  if(is_void(verbose)) verbose=0;
  ff=open(tagfile,"rb");

  if(verbose) write,"Reading "+tagfile;
  adress=0;
  dummy=array(int);_read,ff,adress,dummy;adress+=4;
  if(verbose)  write,"Npart in simulation=",dummy;
  Npart=dummy;

  if(verbose) write,"Getting Densities";
  tagarray=array(float,Npart);_read,ff,adress,tagarray;

  close,ff;
  return tagarray;
}



func read_pos(posfile)
{
  ff=open(posfile,"rb");
  adress=0;
  n=array(int);
  _read,ff,adress,n;write,n;adress+=sizeof(n);
  pos=array(float,3*n);
  _read,ff,adress,pos;
  close,ff;
  return reform(pos,[2,3,n]);
}
   
func read_group_pos(groupposfile)
{
  ff=open(groupposfile,"rb");
  adress=0;
  n=array(int);
  _read,ff,adress,n;write,n;adress+=sizeof(n);
  pos=array(float,3*n);
  _read,ff,adress,pos;
  close,ff;
  return reform(pos,[2,3,n]);
}
