

#include "utils/cudaton.i"


//ll=exec("ls -d runs/strom128comov.0*.p*");
nc=512;
nsplitx=4;
nsplity=4;
nsplitz=2;
bnd=16;
x=array(float,nc,nc,nc);
t=array(float,nc,nc,nc);
//d=array(float,nsplit*nc,nsplit*nc,nsplit*nc);
//fname="/home/cont003/teyssier/dom/prod/runs/strom/1024_pave/out/strom1024_pave.%05d.p%05d";
fname="/home/cont003/teyssier/dom/cudatonGC/trunk/arch3/runs/out/snap.%05d.p%05d";
//ll=exec("ls -d "+"/home/cont003/teyssier/dom/prod/runs/strom/1024_pave/out/strom1024_pave.0*.p00000");
ll=exec("ls -d "+"/home/cont003/teyssier/dom/cudatonGC/trunk/arch3/runs/out/snap.0*.p00000");
//ll=ll(:24);
xion=temp=array(float,nc,nc,numberof(ll));
//zion=iion=array(float,nsplit*nc,nsplit*nc,nsplit*nc);
xavg=tavg=xrms=trms=tt=aa=array(float,numberof(ll));
va=vr=array(float,numberof(ll));

//for(isnap=2;isnap<=numberof(ll);isnap++)
for(isnap=1;isnap<=numberof(ll);isnap++)
{
  isnap;
  for(k=0;k<=nsplitz-1;k++)
    {
      for(j=0;j<=nsplity-1;j++)
        {
          for(i=0;i<=nsplitx-1;i++)
            {
              rank=i+j*nsplitx+k*nsplitx*nsplity;
              s=readsnappave(swrite(format=fname,isnap,rank),nbnd=bnd,npave=nsplitx/nsplitz,getflux=1);
              swrite(format=fname,isnap,rank);
              INFO,(*s.xion);
              x(i*(nc/nsplitx)+1:(i+1)*(nc/nsplitx),j*(nc/nsplity)+1:(j+1)*(nc/nsplity),k*(nc/nsplitz)+1:(k+1)*(nc/nsplitz))=(*s.xion)(,,);
              t(i*(nc/nsplitx)+1:(i+1)*(nc/nsplitx),j*(nc/nsplity)+1:(j+1)*(nc/nsplity),k*(nc/nsplitz)+1:(k+1)*(nc/nsplitz))=(*s.temp)(,,);
              //              d(i*(nc)+1:(i+1)*(nc),j*(nc)+1:(j+1)*(nc),k*(nc)+1:(k+1)*(nc))=(*s.dens)(,,);
            }
        }
    }

  
  xion(,,isnap)=x(,,nc/2);
  temp(,,isnap)=t(,,nc/2);
  

  // va(isnap)=s.aexp;
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
