
vv=load("hpc/energystat.txt");
//vv=load("hpc/elog/e7_7h.txt");
lmax=10;
lc=7;
ks=10;

vv2=[];
ns=numberof(vv(,0));
for(i=1;i<=ns;i++){
  tag=1;
  for(l=lmax;l>=lc;l--){
    if(vv(i-(lmax-l),0)!=l){
      tag=0;
      break;
    }
  }
  if(tag){
    grow,vv2,[vv(i,)];
  }
 }
vv=transpose(vv2);
info,vv;
//error();

ns=numberof(vv(,1));
de=ae=[];
for(i=ks+1;i<=ns;i++){
  i;
  t=vv(i,4);
  t0=vv(i-ks,4);
  u=vv(i,5);
  u0=vv(i-ks,5);
  ei=vv(i,7);
  ei0=vv(i-ks,7);
  ea=(vv(i-ks:i,4)+vv(i-ks:i,5)+vv(i-ks:i,7))(avg);
  rhs=integ(vv(i-ks:i,5)/vv(i-ks:i,1),vv(i-ks:i,1),vv(i,1));
  //grow,de,(t+u+ei-(t0+u0+ei0)-rhs)/(t+u+ei+t0+u0+ei0)*2;
  // grow,de,(t+u+ei-(t0+u0+ei0)-rhs)/ea;
  grow,de,(t+u+ei-(t0+u0+ei0)-rhs)/(u-u0);
  // grow,de,(t+u+ei-(t0+u0+ei0)-rhs)/rhs;
  grow,ae,vv(i-ks:i,1)(avg);
 }
