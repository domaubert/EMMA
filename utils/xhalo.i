
rep="data_4_schaye/";nmax=9;dcen=0.005;lmax=15;snap=36;
WS,10;
WS,11;
WS,12;
window,10,style="win55.gs";
window,11,style="win55.gs";
window,12,style="win55.gs";
mapx=mapy=mapz=[];

i=0;
for(ih=1;ih<=25;ih+=1){
  i++;
  cen=ph(:3,ih);
  dz=oct2cube(swrite(format=rep+"grid.%05d",snap),lmax,101,a,ncpu=32,execut="utils/oct2grid ",cen=cen,dcen=dcen);
  p=mergepart(swrite(format=rep+"part.%05d",snap),32,a,cen=cen,dcen=dcen);
  s=mergepart(swrite(format=rep+"star.%05d",snap),32,a,star=1,cen=cen,dcen=dcen);
  www=where2(dz==max(dz))(,1);
  // grow,mapz,[dz(,,www(3))];
  // grow,mapy,[dz(,www(2),)];
  // grow,mapx,[dz(www(1),,)];
  grow,mapz,[dz(,,avg)];
  grow,mapy,[dz(,avg,)];
  grow,mapx,[dz(avg,,)];

  window,10;
  plsys,i;
  pli,log10(mapz(,,i)+1e-2),cen(1)-dcen,cen(2)-dcen,cen(1)+dcen,cen(2)+dcen,cmin=-2,cmax=8;
  ol=limits();
  // window,11;
  // plsys,i;
  //pl,p(2,),p(1,),color="black";
  if(numberof(s)(0)>0) pl,s(2,),s(1,),color="red";
  //  limits,ol(1),ol(2),ol(3),ol(4);
 }
