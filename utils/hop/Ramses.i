
/*
  The Ramses.i Library contains routines which deals with Ramses
  snapshots structures
*/


struct amr_part{
  double pos(3);
  float vel(3);
  float mass;
  int id;
}


struct amr_header{

  // Run Parameters
  
  int ncpu;
  int ndim;
  int n_loc(3);
  int nlevelmax;
  int ngridmax;
  int nstepcoarse;

  // Physical Parameters

  float scale;
  float t;
  float aexp;
  float hexp;
  float omegam;
  float omegav;
  float omegak;
  float omegab;
  float scale_l;
  float scale_d;
  float scale_t;
}

struct amr_grid{
  int level;
  float x(3);
  int index;
  int sons(8);
}

struct amr_phys{
  int level;
  float dens(8);
  float vel(8,3);
  float pressure(8);
}

struct amr_dens{
  int level;
  float dens(8);
}

struct amr{

  amr_grid grid;
  amr_phys gaz;
}

struct amr_light{

  amr_grid grid;
  amr_dens gaz;
}


func Read_AMR(dirname,&header,&apart,levelmax=,verbose=,megaverbose=,gaz=,part=,ncpu=,trig=,ngridmax=,npartmax=)
{

  /* DOCUMENT

  Read_AMR(dirname,levelmax=,verbose=,resort=,megaverbose=)

  reads the amr structure contained in DIRNAME up to LEVELMAX
  if RESORT==1 the routine sort the Tree structure relative to the indexes
  
  */

  if(is_void(verbose)) verbose=0;
  if(is_void(megaverbose)) megaverbose=0;
  if(megaverbose) verbose=1;
  if(is_void(resort)) resort=0;
  if(is_void(gaz)) gaz=0;
  if(is_void(part)) part=0;
  if(is_void(ngridmax)) ngridmax=100000;
  if(is_void(npart)) npartmax=512^3;
  
  ll=exec("ls -d "+dirname+"/amr_*");
  lg=exec("ls -d "+dirname+"/hydro_*");
  lp=exec("ls -d "+dirname+"/part_*");
  if(verbose) write," Found "+pr1(numberof(ll))+" files";
  if(is_void(ncpu))  ncpu=numberof(ll);
  if(is_void(ncpumin)) ncpumin=1;
  if(is_void(ncpumax)) ncpumax=ncpu;
  
  grid=[];
  phys=[];
  apart=[];
  
  //  write,"***************************************";
  //  write," HAVE YOU CHECKED TRIGUEDINA STYLE ????";
  //  write,"***************************************";
  
  gridt=Read_GRID_local(ll(1),header,verbose=megaverbose,levelmax=levelmax,trig=trig);

  if(gaz)
    {
      grid=array(amr_grid,ngridmax);
      phys=array(amr_phys,ngridmax);
    }
  if(part)
    {
      apart=array(amr_part,npartmax);
    }
  offgrid=1;
  offpart=1;

  for(i=ncpumin;i<=ncpumax;i++)
    {
      if(megaverbose) write,"Reading "+ll(i);
      if(megaverbose) write,"getting grid";

      if(gaz)
        {
          gridt=Read_GRID_local(ll(i),header,verbose=megaverbose,levelmax=levelmax,trig=trig);
          grid(offgrid:offgrid+numberof(gridt)-1)=gridt;

          if(megaverbose) write,"getting gaz";
          physt=Read_PHYS_local(lg(i),verbose=megaverbose,levelmax=levelmax);
          phys(offgrid:offgrid+numberof(physt)-1)=physt;

        }

      if(part)
        {
          if(megaverbose) write,"getting particles";
          partt=Read_PART_local(lp(i),verbose=megaverbose);
          apart(offpart:offpart+numberof(partt)-1)=partt;
        }

      offgrid=offgrid+numberof(gridt);
      offpart=offpart+numberof(partt);
    }

  if(!is_void(grid))
    {
      grid=grid(:offgrid);
      phys=phys(:offgrid);
      res=array(amr,numberof(grid));
      res.grid=grid;
      if(gaz) res.gaz=phys;

    }

  if(!is_void(apart)) apart=apart(:offpart-1);

  
  grid=[];
  phys=[];

  return res;
}

func Read_AMR_LIGHT(dirname,&header,&apart,levelmax=,verbose=,megaverbose=,gaz=,part=,ncpu=,trig=,ngridmax=,npartmax=,ncpumin=,ncpumax=)
{

  /* DOCUMENT

  Read_AMR(dirname,levelmax=,verbose=,resort=,megaverbose=)

  reads the amr structure contained in DIRNAME up to LEVELMAX
  if RESORT==1 the routine sort the Tree structure relative to the indexes
  
  */

  if(is_void(verbose)) verbose=0;
  if(is_void(megaverbose)) megaverbose=0;
  if(megaverbose) verbose=1;
  if(is_void(resort)) resort=0;
  if(is_void(gaz)) gaz=0;
  if(is_void(part)) part=0;
  if(is_void(ngridmax)) ngridmax=1000000;
  if(is_void(npart)) npartmax=512^3;
  if(is_void(trig)) trig=0;
  
  ll=exec("ls -d "+dirname+"/amr_*");
  lg=exec("ls -d "+dirname+"/hydro_*");
  lp=exec("ls -d "+dirname+"/part_*");
  if(verbose) write," Found "+pr1(numberof(ll))+" files";
  if(is_void(ncpu))  ncpu=numberof(ll);
  if(is_void(ncpumin)) ncpumin=1;
  if(is_void(ncpumax)) ncpumax=ncpu;

  
  grid=[];
  phys=[];
  apart=[];
  
  //  write,"***************************************";
  //  write," HAVE YOU CHECKED TRIGUEDINA STYLE ????";
  //  write,"***************************************";
  
  gridt=Read_GRID_local(ll(1),header,verbose=megaverbose,levelmax=levelmax,trig=trig);

  if(gaz)
    {
      grid=array(amr_grid,ngridmax);
      phys=array(amr_dens,ngridmax);
    }
  if(part)
    {
      apart=array(amr_part,npartmax);
    }
  offgrid=1;
  offpart=1;

  for(i=ncpumin;i<=ncpumax;i++)
    {
      if(megaverbose) write,"Reading "+ll(i);
      if(megaverbose) write,"getting grid";

      if(gaz)
        {
          gridt=Read_GRID_local(ll(i),header,verbose=megaverbose,levelmax=levelmax,trig=trig);
          grid(offgrid:offgrid+numberof(gridt)-1)=gridt;

          if(megaverbose) write,"getting gaz";
          physt=Read_LIGHT_local(lg(i),verbose=megaverbose,levelmax=levelmax);
          phys(offgrid:offgrid+numberof(physt)-1)=physt;

        }

      if(part)
        {
          if(megaverbose) write,"getting particles";
          partt=Read_PART_local(lp(i),verbose=megaverbose);
          apart(offpart:offpart+numberof(partt)-1)=partt;
        }

      offgrid=offgrid+numberof(gridt);
      offpart=offpart+numberof(partt);
    }

  if(!is_void(grid))
    {
      grid=grid(:offgrid);
      phys=phys(:offgrid);
      res=array(amr_light,numberof(grid));
      res.grid=grid;
      if(gaz) res.gaz=phys;

    }

  if(!is_void(apart)) apart=apart(:offpart-1);

  
  grid=[];
  phys=[];

  return res;
}


func Read_PART_local(filename,verbose=)
{

  /* DOCUMENT

  pp=array(float,3,512^3+1000);istart=1;for(i=1;i<=64;i++) {ff=swrite(format="/scratch/cont003/p568da/cosmo400/512/output_00001/part_00001.out%05d",i);ff;aa=Read_PART_local(ff);n=dimsof(aa.pos)(0);pp(,istart:istart+n-1)=aa.pos;istart+=n;}
   */
    if(is_void(verbose)) verbose=0;


  // *****************************
  // Checking and Opening the file
  // *****************************

  if(open(filename,"rb",1))
    {
      ff=open(filename,"rb");
    }   
  else
    {
      error,filename+" does not exist !";
    }
  
  // *****************************
  // Check For Byte Swap
  // *****************************
  

  check_bs=0;
  adress=0;

  sun_primitives,ff;
  check_bs=CheckSwap(ff,4,adress);
  
  if(!check_bs)
    {
      dec_primitives,ff;
      check_bs=CheckSwap(ff,4,adress);
    }

  if(!check_bs) error,"Havent found the right Byte Swapping";
  
  adress=0;

  ncpu=ndim=npart=seed=nstar_tot=nsink=array(int);
  mstar_tot=mstar_lost=array(double);
  
  ncpu=readf(ff,ncpu,adress);
  ndim=readf(ff,ndim,adress);
  npart=readf(ff,npart,adress);
  
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
  dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);

//   dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);dummy;
//   dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);
//   dummy=array(int);_read,ff,adress,dummy;adress+=sizeof(dummy);dummy;

  //   nstar_tot=readf(ff,nstar_tot,adress);
//   mstar_tot=readf(ff,mstar_tot,adress);
//   mstar_lost=readf(ff,mstar_lost,adress);
//   nsink=readf(ff,nsink,adress);

  
  if(verbose) write,"Found "+pr1(npart)+" particles";

  res=array(amr_part,npart);

  for(i=1;i<=ndim;i++)
    {
      res.pos(i,)=readf(ff,res.pos(i,),adress);
    }
  
  for(i=1;i<=ndim;i++)
    {
      res.vel(i,)=readf(ff,res.vel(i,),adress);
    }
  res.mass=readf(ff,res.mass,adress);
  res.id=readf(ff,res.id,adress);


  close,ff;
  return res;
}


func Read_PHYS_local(filename,levelmax=,verbose=)
{
  if(is_void(verbose)) verbose=0;


  // *****************************
  // Checking and Opening the file
  // *****************************

  if(open(filename,"rb",1))
    {
      ff=open(filename,"rb");
    }   
  else
    {
      error,filename+" does not exist !";
    }
  
  // *****************************
  // Check For Byte Swap
  // *****************************
  

  check_bs=0;
  adress=0;

  sun_primitives,ff;
  check_bs=CheckSwap(ff,4,adress);
  
  if(!check_bs)
    {
      dec_primitives,ff;
      check_bs=CheckSwap(ff,4,adress);
    }

  if(!check_bs) error,"Havent found the right Byte Swapping";

  // ******************************
  // Read the small physical header
  // ******************************

  ncpu=nvar=ndim=nlevelmax=array(int);
  gamma=array(float);
  adress=0;

  ncpu=readf(ff,ncpu,adress);
  nvar=readf(ff,nvar,adress);
  ndim=readf(ff,ndim,adress);
  nlevelmax=readf(ff,nlevelmax,adress);
  gamma=readf(ff,gamma,adress);
  
  if(is_void(levelmax)) levelmax=nlevelmax;

  Vphys=[]; //Will contain an array of phys structures

  ilevel=array(int);
  ncache=array(int);

  
  for(i=1;i<=levelmax;i++) // Loop over refinments levels;
    {
      ilevel=readf(ff,ilevel,adress);
      ncache=readf(ff,ncache,adress);

      if(verbose) write,"Parsing level "+pr1(i);
      if(verbose) write,"Ngrids (curr. level)= "+pr1(ncache);

      if(ncache>0)    // Reading the physical quantities
        {
          physt=array(amr_phys,ncache);
          physt.level=array(ilevel,ncache);
          
          for(icell=1;icell<=2^ndim;icell++)
            {

              // Density
              physt.dens(icell,)=readf(ff,physt.dens(icell,),adress);
              
              //Velocity
              for(idim=1;idim<=ndim;idim++)
                {
                  physt.vel(icell,idim,)=readf(ff,physt.vel(icell,idim,),adress);
                }
              
              //Pressure
              physt.pressure(icell,)=readf(ff,physt.pressure(icell,),adress);


              //skipping other variables

              for(iskip=ndim+3;iskip<=nvar;iskip++)
                {
                  dummy=array(float,ncache);
                  dummy=readf(ff,dummy,adress);
                }
              
            }   

          // grow vector
          grow,Vphys,physt;
        }
      
    }

  close,ff;

  return Vphys;
}

func Read_LIGHT_local(filename,levelmax=,verbose=)
{
  if(is_void(verbose)) verbose=0;


  // *****************************
  // Checking and Opening the file
  // *****************************

  if(open(filename,"rb",1))
    {
      ff=open(filename,"rb");
    }   
  else
    {
      error,filename+" does not exist !";
    }
  
  // *****************************
  // Check For Byte Swap
  // *****************************
  

  check_bs=0;
  adress=0;

  sun_primitives,ff;
  check_bs=CheckSwap(ff,4,adress);
  
  if(!check_bs)
    {
      dec_primitives,ff;
      check_bs=CheckSwap(ff,4,adress);
    }

  if(!check_bs) error,"Havent found the right Byte Swapping";

  // ******************************
  // Read the small physical header
  // ******************************

  ncpu=nvar=ndim=nlevelmax=array(int);
  gamma=array(float);
  adress=0;

  ncpu=readf(ff,ncpu,adress);
  nvar=readf(ff,nvar,adress);
  ndim=readf(ff,ndim,adress);
  nlevelmax=readf(ff,nlevelmax,adress);
  gamma=readf(ff,gamma,adress);
  
  if(is_void(levelmax)) levelmax=nlevelmax;

  Vphys=[]; //Will contain an array of phys structures

  ilevel=array(int);
  ncache=array(int);

  
  for(i=1;i<=levelmax;i++) // Loop over refinments levels;
    {
      ilevel=readf(ff,ilevel,adress);
      ncache=readf(ff,ncache,adress);

      if(verbose) write,"Parsing level "+pr1(i);
      if(verbose) write,"Ngrids (curr. level)= "+pr1(ncache);

      if(ncache>0)    // Reading the physical quantities
        {
          physt=array(amr_dens,ncache);
          physt.level=array(ilevel,ncache);
          
          for(icell=1;icell<=2^ndim;icell++)
            {

              // Density
              physt.dens(icell,)=readf(ff,physt.dens(icell,),adress);
              
              //Velocity SKIP
              for(idim=1;idim<=ndim;idim++)
                {
                  dummy=array(float,ncache);
                  dummy=readf(ff,dummy,adress);
                }
              
              //Pressure SKIP
              dummy=array(float,ncache);
              dummy=readf(ff,dummy,adress);

              //skipping other variables

              for(iskip=ndim+3;iskip<=nvar;iskip++)
                {
                  dummy=array(float,ncache);
                  dummy=readf(ff,dummy,adress);
                }
              
            }   

          // grow vector
          grow,Vphys,physt;
        }
      
    }

  close,ff;

  return Vphys;
}



func Read_GRID_local(filename,&hdr,levelmax=,verbose=,trig=)
{
  /* DOCUMENT

  Reads the "local" AMR Structure
  
  
  */

  if(is_void(verbose)) verbose=0;

  // Reads the header
  hdr=Read_header(filename,adress,ff,trig=trig);
  if(is_void(levelmax)) levelmax=hdr.nlevelmax;
  

  
  ilevel=array(int);
  ncache=array(int);
  Vncache=array(0,levelmax);
  
  Vgrid=[]; //Will contain an array of grid structures

  for(i=1;i<=levelmax;i++) // Loop over refinments levels;
    {
      ilevel=readf(ff,ilevel,adress);
      ncache=readf(ff,ncache,adress);
      Vncache(i)=ncache;
        
      if(verbose) write,"Parsing level "+pr1(i);
      if(verbose) write,"Ngrids (curr. level)= "+pr1(ncache);

      if(ncache>0)    // Reading the Grid Positions
        {
          gridt=array(amr_grid,ncache);
          gridt.level=array(ilevel,ncache);
            
          for(idim=1;idim<=hdr.ndim;idim++)
            {
              gridt.x(idim,)=readf(ff,gridt.x(idim,),adress);
            }
          grow,Vgrid,gridt;
        }
    }

  // Tree Structure for Grids

  offset=0;
  for(i=1;i<=levelmax;i++) // Loop over refinments levels;
    {
      if(Vncache(i)>0)
        {
          // Looking for indexes
          Vgrid(offset+1:offset+Vncache(i)).index=readf(ff,Vgrid(offset+1:offset+Vncache(i)).index,adress);
          for(k=1;k<=2^hdr.ndim;k++)// Looking for Sons
            {
              Vgrid(offset+1:offset+Vncache(i)).sons(k)=readf(ff,Vgrid(offset+1:offset+Vncache(i)).sons(k),adress);
            }
          offset+=Vncache(i);
        }
    }
  

  close,ff;
  
  return Vgrid;
  
}




func Read_header(filename,&adress,&ff,trig=)
{
  /* DOCUMENT

  Read_header(filename,&adress)
     
  Reads and return the header of contained in FILENAME.
  ADRESS contains the current pointer adress for further reading.
  
   */

  if(is_void(trig))
    {
      trig=1;
      write,"*****************************";
      write," WARNING TRIGUEDINA STYLE";
      write,"*****************************";

    }
  

  // *****************************
  // Checking and Opening the file
  // *****************************

  if(open(filename,"rb",1))
    {
      ff=open(filename,"rb");
    }   
  else
    {
      error,filename+" does not exist !";
    }

  // *****************************
  // Check For Byte Swap
  // *****************************


  check_bs=0;
  adress=0;

  sun_primitives,ff;
  check_bs=CheckSwap(ff,4,adress);

  if(!check_bs)
    {
      dec_primitives,ff;
      check_bs=CheckSwap(ff,4,adress);
    }

  if(!check_bs) error,"Havent found the right Byte Swapping";
  // *****************************
  // Reading the Header
  // *****************************
  
  hdr=amr_header();
  adress=0;
  hdr.ncpu=readf(ff,hdr.ncpu,adress);
  hdr.ndim=readf(ff,hdr.ndim,adress);
  hdr.n_loc=readf(ff,hdr.n_loc,adress);
  hdr.nlevelmax=readf(ff,hdr.nlevelmax,adress);
  hdr.ngridmax=readf(ff,hdr.ngridmax,adress);
  hdr.nstepcoarse=readf(ff,hdr.nstepcoarse,adress);
  hdr.scale=readf(ff,hdr.scale,adress);
  hdr.t=readf(ff,hdr.t,adress,closer=0);
  hdr.aexp=readf(ff,hdr.aexp,adress,openr=0,closer=0);
  hdr.hexp=readf(ff,hdr.hexp,adress,openr=0);
  hdr.omegam=readf(ff,hdr.omegam,adress,closer=0);
  hdr.omegav=readf(ff,hdr.omegav,adress,openr=0,closer=0);
  hdr.omegak=readf(ff,hdr.omegak,adress,openr=0,closer=0);
  hdr.omegab=readf(ff,hdr.omegab,adress,openr=0);

  if(!trig)
    {
      hdr.scale_l=readf(ff,hdr.scale_l,adress,closer=0);
      hdr.scale_d=readf(ff,hdr.scale_d,adress,openr=0,closer=0);  //PROBABLY WRONG
      hdr.scale_t=readf(ff,hdr.scale_t,adress,openr=0);
    }
  
  return hdr;
}






func readf(file,x,&address,openr=,closer=)
{
  /* DOCUMENT

  readf(file,x,&address,openr=,closer=)

  Read a fortran unformatted record.
  It looks the object X in FILE at adress ADDRESS and returns it.
  By default, this routines looks for a full reccord that contains only X.
  If several objects coexists in the same record, one should decide if the
  record has to be opened or closed using the OPENR/CLOSER=1/0 options.

  e.g.

  ff=open("toto.dat","rb");
  adress=0;
  hdr.nstepcoarse=readf(ff,hdr.nstepcoarse,adress);
  hdr.scale=readf(ff,hdr.scale,adress);
  hdr.t=readf(ff,hdr.t,adress,closer=0);
  hdr.aexp=readf(ff,hdr.aexp,adress,openr=0,closer=0);
  hdr.hexp=readf(ff,hdr.hexp,adress,openr=0);

  
   */
  if(is_void(openr)) openr=1;
  if(is_void(closer)) closer=1;
  
  if(openr){ dummy=array(int);_read,file,address,dummy;address+=sizeof(dummy);}
  _read,file,address,x;address+=sizeof(x);
  if(closer){  dummy=array(int);_read,file,address,dummy;address+=sizeof(dummy);}
  
  return x;
}

func CheckSwap(ff,n,adress)
{
  /* DOCUMENT

  CheckSwap(ff,n,adress)

  return 1 if the integer located at ADRESS in the file FF is equal to N
  
   */

  dummy=array(int);_read,ff,adress,dummy;
  return (dummy==n);
}


func DoubleImage(field,n=)
{
  if(is_void(n)) n=1;
  if(n==0)
    {
      return field;
    }
  else
    {
      for(pp=1;pp<=n;pp++)
        {
          ni=dimsof(field)(2);
          nj=dimsof(field)(3);
          res=array(structof(field),2*ni,2*nj);
          ii=indgen(ni)(,-:1:nj)(*);
          jj=indgen(nj)(-:1:ni,)(*);
          res(*)(2*ii-1+(2*jj-2)*(2*ni))=field(*)(ii+(jj-1)*ni);
          res(*)(2*ii+(2*jj-2)*(2*ni))=field(*)(ii+(jj-1)*ni);
          res(*)(2*ii-1+(2*jj-1)*(2*ni))=field(*)(ii+(jj-1)*ni);
          res(*)(2*ii+(2*jj-1)*(2*ni))=field(*)(ii+(jj-1)*ni);
          field=res;
        }

      return res;
    }

}


func GetInterp(amesh,levelmin=,levelmax=,massp=,maxd=,iz=)
{

  if(is_void(dposc)) dposc=[.1,.1,.1];
  if(is_void(posc)) posc=[.5,.5,.5];
  if(is_void(massp)) massp=1;
  if(is_void(maxd)) maxd=0;


  ngridmax=2^(levelmax-1)*2;
  res=array(float,ngridmax,ngridmax);
  if(is_void(iz)) iz=indgen(2^levelmin);
      
  write,"final grid will be "+pr1(ngridmax)+"x"+pr1(ngridmax);

  ni=1;
  for(i=levelmax;i>=levelmin;i--)
  //for(i=levelmin;i<=levelmax;i++)
    {
      ngrid=2^i;
      rest=array(float(0),ngrid,ngrid);

      dx=1./ngrid;
      www=where((amesh.grid.level==i)*(amesh.grid.sons(sum,)==0));
      if(numberof(www)==0)
        {
          write,"nolevel found";
          ni=ni*2;
          continue;
        }
      
      index_z=int(amesh(www).grid.x(3,)/(dx));
      amt=amesh(www);
      for(iz=1;iz<=ngrid;iz++)
        {
          //          rest2=array(float(0),ngrid,ngrid);

          wz=where(index_z==iz);
          if(numberof(wz)==0) continue;
          amtz=amt(wz);
          
          index_x=int(amtz.grid.x(1,)/(dx));
          index_y=int(amtz.grid.x(2,)/(dx));

          rest(index_x+(index_y-1)*ngrid)+=(amtz.gaz.dens([1,3],))(sum,);
          rest(index_x+1+(index_y-1)*ngrid)+=(amtz.gaz.dens([2,4],))(sum,);
          rest(index_x+1+index_y*ngrid)+=(amtz.gaz.dens([6,8],))(sum,);
          rest(index_x+index_y*ngrid)+=(amtz.gaz.dens([5,7],))(sum,);
        }
      rest=DoubleImage(rest,n=levelmax-i);

      //www0=where(rest!=0);
      //res(www0)=rest(www0);

      res=max(res,rest);

      ni=ni*2;

    }

  return res;
}

func GetDens(dir,&hdr,levelmin=,levelmax=,massp=,maxd=,ncpu=)
{

  if(is_void(dposc)) dposc=[.1,.1,.1];
  if(is_void(posc)) posc=[.5,.5,.5];
  if(is_void(massp)) massp=1;
  if(is_void(maxd)) maxd=0;
  

  ngridmax=2^(levelmax-1)*2;
  res=resf=array(float,ngridmax,ngridmax);
  if(is_void(iz)) iz=indgen(2^levelmin);
      
  write,"final grid will be "+pr1(ngridmax)+"x"+pr1(ngridmax);

  ll=exec("ls -d "+dir+"/amr_*");
  INFO,ll;
  if(is_void(ncpu)) ncpu=numberof(ll);
  

  for(i=levelmax;i>=levelmin;i--)
    //for(i=levelmin;i<=levelmax;i++)
    {
      ngrid=2^i;
      rest=array(float(0),ngrid,ngrid);
      
      dx=1./ngrid;
      
      for(ifile=1;ifile<=ncpu;ifile++)
        {
          write,"ifile=",ifile;
          amesh=Read_AMR_LIGHT(dir,hdr,apart,gaz=1,trig=0,part=0,ncpumin=ifile,ncpumax=ifile);
          
          www=where((amesh.grid.level==i)*(amesh.grid.sons(sum,)==0));
          if(numberof(www)==0)
            {
              write,"nolevel found";
              continue;
            }
          
          index_z=int(amesh(www).grid.x(3,)/(dx));
          amt=amesh(www);
          for(iz=1;iz<=ngrid;iz++)
            {
              //          rest2=array(float(0),ngrid,ngrid);
              
              wz=where(index_z==iz);
              if(numberof(wz)==0) continue;
              amtz=amt(wz);
              
              index_x=int(amtz.grid.x(1,)/(dx));
              index_y=int(amtz.grid.x(2,)/(dx));
              
              rest(index_x+(index_y-1)*ngrid)+=(amtz.gaz.dens([1,3],))(sum,);
              rest(index_x+1+(index_y-1)*ngrid)+=(amtz.gaz.dens([2,4],))(sum,);
              rest(index_x+1+index_y*ngrid)+=(amtz.gaz.dens([6,8],))(sum,);
              rest(index_x+index_y*ngrid)+=(amtz.gaz.dens([5,7],))(sum,);
            }
          //www0=where(rest!=0);
          //res(www0)=rest(www0);
        }
      
      rest=DoubleImage(rest,n=levelmax-i);
      res=max(res,rest);
    }

  return res;
}


func PlotAMR(amesh,levelmin=,levelmax=,posc=,dpos=,fill=,dens=,velx=,vely=,velz=,vel=,pressure=,logz=,ologz=)
{
  if(is_void(levelmin)) levelmin=2;
  if(is_void(levelmax)) levelmax=max(amesh.grid.level);
  if(is_void(posc)) posc=[.5,.5,.5];
  if(is_void(dpos)) dpos=[.5,.5,.05];
  if(is_void(logz)) logz=0;
  if(is_void(ologz)) ologz=1e-6;

  opt=[1,0,0]; //cell display
  if(fill)      opt=[0,1,0,0,0,0,0,0];
  if(dens)      opt=[0,1,1,0,0,0,0,0];
  if(velx)      opt=[0,1,0,1,0,0,0,0];
  if(vely)      opt=[0,1,0,0,1,0,0,0];
  if(velz)      opt=[0,1,0,0,0,1,0,0];
  if(vel)       opt=[0,1,0,0,0,0,1,0];
  if(pressure)  opt=[0,1,0,0,0,0,0,1];
  
  
  ii=indgen(numberof(amesh.grid));
  www=where((amesh.grid.level>=levelmin)*(amesh.grid.level<=levelmax));
  amesh=amesh(www);
  ii=ii(www);
  www=where((abs(amesh.grid.x(3,)-posc(3))<dpos(3))*(abs(amesh.grid.x(2,)-posc(2))<dpos(2))*(abs(amesh.grid.x(1,)-posc(1))<dpos(1)));
  //www=where((abs(amesh.grid.x(1,)-posc(1))<1./(2^(amesh.grid.level)))*(abs(amesh.grid.x(2,)-posc(2))<1./(2^(amesh.grid.level)))*(abs(amesh.grid.x(3,)-posc(3))<1./(2^(amesh.grid.level))));
  amesh=amesh(www);
  ii=ii(www);

  sz=sort(amesh.grid.x(3,));
  amesh=amesh(sz(::-1));

  if(numberof(www)==0) error,"The selected grid is empty in the current range";
  
  cc=Colors(levelmax-levelmin+1);
  xgf=ygf=array(float,numberof(www)*4*4);
  zgf=array(float,numberof(www)*4);

  off=0;
  offz=0;

  //pos matrix

  //  mx=transpose([[-.5,0,0,-.5],[0.,0.5,0.5,0.],[-.5,0,0,-.5],[0.,0.5,0.5,0.],[-.5,0,0,-.5],[0.,0.5,0.5,0.],[-.5,0,0,-.5],[0.,0.5,0.5,0.]]);
  //  my=transpose([[-.5,-.5,0,0],[-.5,-.5,0,0],[0,0,.5,.5],[0,0,.5,.5],[-.5,-.5,0,0],[-.5,-.5,0,0],[0,0,.5,.5],[0,0,.5,.5]]);

  //  mx=float([-.5,0,0,-.5,0.,0.5,0.5,0.,-.5,0,0,-.5,0.,0.5,0.5,0.,-.5,0,0,-.5,0.,0.5,0.5,0.,-.5,0,0,-.5,0.,0.5,0.5,0.]);
  //  my=float([-.5,-.5,0,0,-.5,-.5,0,0,0,0,.5,.5,0,0,.5,.5,-.5,-.5,0,0,-.5,-.5,0,0,0,0,.5,.5,0,0,.5,.5]);

  mx=float([-.5,0,0,-.5,0.,0.5,0.5,0.,-.5,0,0,-.5,0.,0.5,0.5,0]);
  my=float([-.5,-.5,0,0,-.5,-.5,0,0,0,0,.5,.5,0,0,.5,.5]);

  
  www2=indgen(numberof(amesh));
  
  
  for(i=1;i<=numberof(www2);i++)
    {
      a=1./(2^(amesh.grid(www2(i)).level-1));//cell size
      
      if(opt(2))
        {
          
          xg=mx*a+amesh.grid(www2(i)).x(1);
          yg=my*a+amesh.grid(www2(i)).x(2);
          
          xgf(off+1:off+16)=xg;
          ygf(off+1:off+16)=yg;
          
          
          if(opt(3))
            {
              zgf(offz+1:offz+4)=amesh(www2(i)).gaz.dens(:4);
            }
          else if(opt(4))
            {
                  zgf(offz+1:offz+4)=amesh(www2(i)).gaz.vel(,1)(:4);
            }
          else if(opt(5))
            {
              zgf(offz+1:offz+4)=amesh(www2(i)).gaz.vel(,2)(:4);
            }
          else if(opt(6))
            {
              zgf(offz+1:offz+4)=amesh(www2(i)).gaz.vel(,3)(:4);
            }
          else if(opt(7))
            {
              zgf(offz+1:offz+4)=abs(amesh(www2(i)).gaz.vel(,3)(:4),amesh(www2(i)).gaz.vel(,1)(:4),amesh(www2(i)).gaz.vel(,2)(:4));
            }
          else if(opt(8))
            {
              zgf(offz+1:offz+4)=amesh(www2(i)).gaz.pressure(:4);
            }
          else
            {
              zgf(offz+1:offz+4)=j/10.;
            }
          
          off=off+16;
          offz+=4;
          
          
        }
      
      if(opt(1))
        {
          xg=span(-.5,.5,3)(,-:1:3)*a+amesh.grid(www2(i)).x(1);
          yg=span(-.5,.5,3)(-:1:3,)*a+amesh.grid(www2(i)).x(2);
        }
      
      if(opt(1)) plm,yg,xg,boundary=0,color=cc(amesh.grid.level(www2(i))-levelmin+1);
    }
  //   }

  write,"Done";
  if(opt(2))
    {
      
     if(logz)
       {
         plfp,log(ologz+zgf),ygf,xgf,array(4,numberof(zgf));
       }
     else
       {
         plfp,zgf,ygf,xgf,array(4,numberof(zgf));
       }
    }
  
  
}

func PlotAMRfast(amesh,levelmin=,levelmax=,posc=,dpos=,fill=,dens=,velx=,vely=,velz=,vel=,pressure=,logz=,ologz=)
{
  if(is_void(levelmin)) levelmin=2;
  if(is_void(levelmax)) levelmax=max(amesh.grid.level);
  if(is_void(posc)) posc=[.5,.5,.5];
  if(is_void(dpos)) dpos=[.5,.5,.05];
  if(is_void(logz)) logz=0;
  if(is_void(ologz)) ologz=1e-6;

  opt=[1,0,0]; //cell display
  if(fill)      opt=[0,1,0,0,0,0,0,0];
  if(dens)      opt=[0,1,1,0,0,0,0,0];
  if(velx)      opt=[0,1,0,1,0,0,0,0];
  if(vely)      opt=[0,1,0,0,1,0,0,0];
  if(velz)      opt=[0,1,0,0,0,1,0,0];
  if(vel)       opt=[0,1,0,0,0,0,1,0];
  if(pressure)  opt=[0,1,0,0,0,0,0,1];
  
  
  www=where((amesh.grid.level>=levelmin)*(amesh.grid.level<=levelmax));
  amesh=amesh(www);
  www=where((abs(amesh.grid.x(3,)-posc(3))<dpos(3))*(abs(amesh.grid.x(2,)-posc(2))<dpos(2))*(abs(amesh.grid.x(1,)-posc(1))<dpos(1)));
  amesh=amesh(www);
  

  
  sz=sort(amesh.grid.x(3,));
  amesh=amesh(sz(::-1));

  if(numberof(www)==0) error,"The selected grid is empty in the current range";
  
  cc=Colors(levelmax-levelmin+1);
  xgf=ygf=array(float,numberof(www)*4*4);
  zgf=array(float,numberof(www)*4);

  off=0;
  offz=0;

  //pos matrix


  mx=float([-.5,0,0,-.5,0.,0.5,0.5,0.,-.5,0,0,-.5,0.,0.5,0.5,0]);
  my=float([-.5,-.5,0,0,-.5,-.5,0,0,0,0,.5,.5,0,0,.5,.5]);

  a=1./(2^(amesh.grid.level-1));//cell size
      
  if(opt(2))
    {
          
      xg=mx*a+amesh.grid().x(1);
      yg=my*a+amesh.grid().x(2);
          
      xgf=xg;
      ygf=yg;
          
      
      if(opt(3))
        {
          zgf=amesh().gaz.dens(:4);
        }
      else if(opt(4))
        {
          zgf=amesh().gaz.vel(,1)(:4);
        }
      else if(opt(5))
        {
          zgf(offz+1:offz+4)=amesh(www2(i)).gaz.vel(,2)(:4);
        }
      else if(opt(6))
        {
          zgf(offz+1:offz+4)=amesh(www2(i)).gaz.vel(,3)(:4);
        }
      else if(opt(7))
        {
          zgf(offz+1:offz+4)=abs(amesh(www2(i)).gaz.vel(,3)(:4),amesh(www2(i)).gaz.vel(,1)(:4),amesh(www2(i)).gaz.vel(,2)(:4));
        }
      else if(opt(8))
        {
          zgf(offz+1:offz+4)=amesh(www2(i)).gaz.pressure(:4);
        }
      else
        {
          zgf(offz+1:offz+4)=j/10.;
        }
    }          
      
      
  if(opt(1))
    {
      xg=span(-.5,.5,3)(,-:1:3)*a+amesh.grid(www2(i)).x(1);
      yg=span(-.5,.5,3)(-:1:3,)*a+amesh.grid(www2(i)).x(2);
    }
      
  if(opt(1)) plm,yg,xg,boundary=0,color=cc(amesh.grid.level(www2(i))-levelmin+1);
    
  if(opt(2))
    {
      
     if(logz)
       {
         plfp,log(ologz+zgf),ygf,xgf,array(4,numberof(zgf));
       }
     else
       {
         plfp,zgf,ygf,xgf,array(4,numberof(zgf));
       }
    }
  
  
}



func __PlotAMR(amesh,levelmin=,levelmax=,posc=,dpos=,fill=,dens=,velx=,vely=,velz=,vel=,pressure=,logz=,ologz=)
{
  if(is_void(levelmin)) levelmin=2;
  if(is_void(levelmax)) levelmax=max(amesh.grid.level);
  if(is_void(posc)) posc=[.5,.5,.5];
  if(is_void(dpos)) dpos=[.5,.5,.05];
  if(is_void(logz)) logz=0;
  if(is_void(ologz)) ologz=1e-6;

  opt=[1,0,0]; //cell display
  if(fill)      opt=[0,1,0,0,0,0,0,0];
  if(dens)      opt=[0,1,1,0,0,0,0,0];
  if(velx)      opt=[0,1,0,1,0,0,0,0];
  if(vely)      opt=[0,1,0,0,1,0,0,0];
  if(velz)      opt=[0,1,0,0,0,1,0,0];
  if(vel)       opt=[0,1,0,0,0,0,1,0];
  if(pressure)  opt=[0,1,0,0,0,0,0,1];
  
  
  ii=indgen(numberof(amesh.grid));
  www=where((amesh.grid.level>=levelmin)*(amesh.grid.level<=levelmax));
  amesh=amesh(www);
  ii=ii(www);
  www=where((abs(amesh.grid.x(3,)-posc(3))<dpos(3))*(abs(amesh.grid.x(2,)-posc(2))<dpos(2))*(abs(amesh.grid.x(1,)-posc(1))<dpos(1)));
  amesh=amesh(www);
  ii=ii(www);

  if(numberof(www)==0) error,"The selected grid is empty in the current range";
  
  cc=Colors(levelmax-levelmin+1);
  xgf=ygf=array(float,numberof(www)*4*4);
  zgf=array(float,numberof(www)*4);

  off=0;
  offz=0;


  mx=float([-.5,0,0,-.5,0.,0.5,0.5,0.,-.5,0,0,-.5,0.,0.5,0.5,0]);
  my=float([-.5,-.5,0,0,-.5,-.5,0,0,0,0,.5,.5,0,0,.5,.5]);


  res=array(float,2^levelmax,2^levelmax); // will contain the final image
  rest=res; // Temporary image
  
  for(j=levelmin;j<=levelmax;j++)
    {

      www2=where(amesh.grid.level==j);
      aoct=2^levelmax/(2^(j-1)); //OCT size in terms of levelmax CELL size
      poct=amesh(www2).grid.x*2^levelmax/(2^(j-1)); // OCT position in terms of levelmax CELL size
    }
  
  for(i=1;i<=numberof(www2);i++)
    {
      a=1./(2^(amesh.grid(www2(i)).level-1));//cell size
      
      if(opt(2))
        {
          if(anyof(amesh.grid(www2(i)).sons!=0)) continue;
          
          xg=mx*a+amesh.grid(www2(i)).x(1);
          yg=my*a+amesh.grid(www2(i)).x(2);
          
          xgf(off+1:off+16)=xg;
          ygf(off+1:off+16)=yg;
          
          
          if(opt(3))
            {
              zgf(offz+1:offz+4)=amesh(www2(i)).gaz.dens(:4);
            }
          else if(opt(4))
            {
                  zgf(offz+1:offz+4)=amesh(www2(i)).gaz.vel(,1)(:4);
            }
          else if(opt(5))
            {
              zgf(offz+1:offz+4)=amesh(www2(i)).gaz.vel(,2)(:4);
            }
          else if(opt(6))
            {
              zgf(offz+1:offz+4)=amesh(www2(i)).gaz.vel(,3)(:4);
            }
          else if(opt(7))
            {
              zgf(offz+1:offz+4)=abs(amesh(www2(i)).gaz.vel(,3)(:4),amesh(www2(i)).gaz.vel(,1)(:4),amesh(www2(i)).gaz.vel(,2)(:4));
            }
          else if(opt(8))
            {
              zgf(offz+1:offz+4)=amesh(www2(i)).gaz.pressure(:4);
            }
          else
            {
              zgf(offz+1:offz+4)=j/10.;
            }
          
          off=off+16;
          offz+=4;
          
          
        }
      
      if(opt(1))
        {
          xg=span(-.5,.5,3)(,-:1:3)*a+amesh.grid(www2(i)).x(1);
          yg=span(-.5,.5,3)(-:1:3,)*a+amesh.grid(www2(i)).x(2);
        }
      
      if(opt(1)) plm,yg,xg,boundary=0,color=cc(amesh.grid.level(www2(i))-levelmin+1);
    }
  //   }

  write,"Done";
  if(opt(2))
    {
      
     if(logz)
       {
         plfp,log(ologz+zgf),ygf,xgf,array(4,numberof(zgf));
       }
     else
       {
         plfp,zgf,ygf,xgf,array(4,numberof(zgf));
       }
    }
  
  
}







func PlotOct(grid,coord=,color=,width=,type=)
{
  if(is_void(coord)) coord=1;

  a=1./(2^(grid.level));

  x=array(float,2,4);
  x(,1)=[grid.x(coord)+a,grid.x(coord+1)+a];
  x(,2)=[grid.x(coord)-a,grid.x(coord+1)+a];
  x(,3)=[grid.x(coord)-a,grid.x(coord+1)-a];
  x(,4)=[grid.x(coord)+a,grid.x(coord+1)-a];

  for(i=0;i<=3;i++)
    {
      plg,[x(2,i),x(2,i+1)],[x(1,i),x(1,i+1)],color=color,width=width,type=type;
    }
  
  return 0;
}
