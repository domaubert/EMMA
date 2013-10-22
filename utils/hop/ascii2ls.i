//#include "sfit.i"

func myread(file,n,&h,nlines=,ncol=,s=){
  /* DOCUMENT
     reads ascii files, where the text header is n lines long, stores the header into h.
     if nlines is not given:
     computes file length (number of data lines) from
     cat -n thomas04.data | tail -n 1 | cut -c -10   (for instance)
     so be careful if the number of lines is VERY LARGE
     if ncol is not given:
     also it computes the number of columns by :
     ncol=numberof(split2words(line(n+1),sep=" "));
     s=1 silent mode
     

  */

  // Limited header reading only to the cases where n>0
  
  local tab;
  
  //if(is_void(nlines)) nlines=str2int(exec("cat -n "+file+" | tail -n 1 | cut -c -10 ")(1))-n;
  if(is_void(nlines)) nlines=str2int(exec("cat "+file+" | wc -l")(1))-n;
  f=open(file,"r");
  if(n>0) h=rdline(f,n);
  a=rdline(f,1);
  if(is_void(ncol)) ncol=numberof(split2words(a,sep=" "));

  if(s!=1) write,ncol,nlines;
  
  close,f;
  f=open(file,"r");
  if (n>0) h=rdline(f,n);
  tab=array(0.,ncol,nlines);
  read,f,tab;
  return tab;
};
  

func stextread(file,n,&h,nlines=,ncol=,s=,dsep=){
  /* DOCUMENT
     raw ascii reading routine: returns columns with text inside
     separator in text file is " " by default (dsep=)
     stextread means "slow" text read, indeed it is not expected to be fast
     NOTE if ncol=1, it just returns the whole file read line by line (can be convenient)
  */

  local tab;

  if(is_void(dsep)) dsep=" ";
  
  //if(is_void(nlines)) nlines=str2int(exec("cat -n "+file+" | tail -n 1 | cut -c -10 ")(1))-n;
  if(is_void(nlines)) nlines=str2int(exec("cat "+file+" | wc -l")(1))-n;
  f=open(file,"r");
  if(n>0) h=rdline(f,n);
  a=rdline(f,1);
  if(is_void(ncol)) ncol=numberof(split2words(a,sep=" "));

  if(s!=1) write,ncol,nlines;
  
  close,f;
  f=open(file,"r");
  if (n>0) h=rdline(f,n);
  tab=array(string,ncol,nlines);
  for(i=1;i<=nlines;i++){
    dummy=rdline(f);
    //dummy;
    tab(,i)=split2words(dummy,sep=dsep)(:ncol);
  };
  return tab;
};



