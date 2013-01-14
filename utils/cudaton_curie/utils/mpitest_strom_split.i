

#include "utils/cudaton.i"


//ll=exec("ls -d runs/strom128comov.0*.p*");
nc=128;
nsplit=4;
bnd=16;
x=array(float,nsplit*nc,nsplit*nc,nsplit*nc);
t=array(float,nsplit*nc,nsplit*nc,nsplit*nc);
d=array(float,nsplit*nc,nsplit*nc,nsplit*nc);
//fname="/home/cont003/teyssier/dom/prod/runs/strom/1024/out/strom1024.%05d.p%05d";
fname="runs/out/snap.%05d.p%05d";

//ll=exec("ls -d "+"/home/cont003/teyssier/dom/prod/runs/strom/1024/out/strom1024.0*.p00000");
ll=exec("ls -d "+"runs/out/snap.0*.p00000");
//ll=ll(:24);
xion=eng=flx=array(float,nsplit*nc,nsplit*nc,numberof(ll));
//zion=iion=array(float,nsplit*nc,nsplit*nc,nsplit*nc);
xavg=tavg=xrms=trms=tt=aa=array(float,numberof(ll));
va=vr=array(float,numberof(ll));

//
for(isnap=1;isnap<=numberof(ll);isnap++)
  {
  isnap;
  for(k=0;k<=nsplit-1;k++)
    {
      for(j=0;j<=nsplit-1;j++)
        {
          for(i=0;i<=nsplit-1;i++)
            {
              rank=i+j*nsplit+k*nsplit*nsplit;
              s=readsnappave(swrite(format=fname,isnap,rank),nbnd=bnd,npave=1,getflux=1);
              swrite(format=fname,isnap,rank);
              x(i*(nc)+1:(i+1)*(nc),j*(nc)+1:(j+1)*(nc),k*(nc)+1:(k+1)*(nc))=(*s.xion)(,,);
              t(i*(nc)+1:(i+1)*(nc),j*(nc)+1:(j+1)*(nc),k*(nc)+1:(k+1)*(nc))=(*s.temp)(,,);
              d(i*(nc)+1:(i+1)*(nc),j*(nc)+1:(j+1)*(nc),k*(nc)+1:(k+1)*(nc))=(*s.egy)(,,);
            }
        }
    }
  INFO,x;
  
  xion(,,isnap)=x(,,nc*nsplit/2);
  eng(,,isnap)=t(,,nc*nsplit/2);
  flx(,,isnap)=d(,,nc*nsplit/2);
  //  vr(isnap)=interp(span(0,1.,nc*nsplit/2),x(nc+1:,nc,nc),0.5);
  //  va(isnap)=s.aexp;

  // 
  //   binx=(indgen(512)-1)/(2.*nc);
  //   spos=255/(2.*nc);
  //     vr(isnap)=interp((binx-spos)(257:),x(257:,256,256),0.5);
  
  
//   xavg(isnap)=x(*)(avg);
//   xrms(isnap)=x(*)(rms);
//   tavg(isnap)=t(*)(avg);
//   trms(isnap)=t(*)(rms);
//   tt(isnap)=s.time;
//   aa(isnap)=s.aexp;
//   www=where((x(,,)(*)>0.5)*(zion(*)==0));
//   if(numberof(www)>0) {
//     zion(*)(www)=s.aexp;
//     iion(*)(www)=isnap;
//   }

}





#if 0
pli,f(,,18);
plg,[0,2*nc],[nc,nc];
plg,[nc,nc],[0,2*nc];

#if 0
pli,f(,,128);
