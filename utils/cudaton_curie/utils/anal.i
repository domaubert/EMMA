
ll=exec("ls -d runs2/out/snap_c01.00.*");
istop=numberof(ll);
aglob=xglob=tglob=[];
for(i=1;i<istop;i++){
  i;
  dummy=readg(ll(i),x,t,a);
  grow,aglob,a;
  grow,xglob,avg(x);
  grow,tglob,avg(t);
 }
