

func readdens(fname,nprocx,nprocy,nprocz,&dens,&pot)
{

  nfile=nprocx*nprocy*nprocz;
  
  for(ifile=0;ifile<nfile;ifile++)
    {
      fich=swrite(format=fname+".p%05d",ifile);
      fich;
      fbuf=open(swrite(format=fname+".p%05d",ifile),"rb");
      adress=0;
      nx=array(int);_read,fbuf,adress,nx;adress+=sizeof(nx);nx;
      ny=array(int);_read,fbuf,adress,ny;adress+=sizeof(ny);ny;
      nz=array(int);_read,fbuf,adress,nz;adress+=sizeof(nz);nz;
      nb=array(int);_read,fbuf,adress,nb;adress+=sizeof(nb);nb;

      k=ifile/(nprocx*nprocy);
      j=(ifile-k*(nprocx*nprocy))/nprocx;
      i=ifile-k*(nprocx*nprocy)-j*nprocx ;
      
      nxnb=nx-2*nb;
      nynb=ny-2*nb;
      nznb=nz-2*nb;

      if(ifile==0)
        {
          dens=array(float,nxnb*nprocx,nynb*nprocy,nznb*nprocz);
          pot=array(float,nxnb*nprocx,nynb*nprocy,nznb*nprocz);
        }

      temp=array(float,nz,ny,nx);
      _read,fbuf,adress,temp;adress+=sizeof(temp);
      dens(k*nznb+1:(k+1)*nznb,j*nynb+1:(j+1)*nynb,i*nxnb+1:(i+1)*nxnb)=temp(nb+1:nb+nznb,nb+1:nb+nynb,nb+1:nb+nxnb);

      
      _read,fbuf,adress,temp;adress+=sizeof(temp);
      pot(k*nznb+1:(k+1)*nxnb,j*nynb+1:(j+1)*nynb,i*nxnb+1:(i+1)*nxnb)=temp(nb+1:nb+nznb,nb+1:nb+nynb,nb+1:nb+nxnb);

      close,fbuf;
    }

  
}
