
nsnap=20;


X=Y=DE=R=FX=FY=F=[];

for(ip=0;ip<=4;ip++){
  x=y=z=egy=t=fx=fy=[];
  x2=y2=[];
  tp=array(double);
  for(i=0;i<=nsnap;i++){
    fname=swrite(format="data/part.%05d.p00000",i);
    p=readpart(fname,tp);
    www=where(p(7,)==ip);
    grow,x,p(1,www);
  
    grow,y,p(2,www);
    grow,egy,p(4,www);
    grow,fx,p(5,www);
    grow,fy,p(6,www);
    grow,t,tp;
  }


  r=abs(x-0.5,y-0.5);
  de=egy/egy(avg)-1.;
  grow,X,[x];
  grow,Y,[y];
  grow,FX,[fx];
  grow,FY,[fy];
  grow,F,[abs(fx,fy)];
  grow,DE,[de];
  grow,R,[r];

 }
