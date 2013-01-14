
ll=exec("ls -d ~/Project/cudaton_cosmoromain/runs2/out/teyssier/c*");
WS,1,dpi=150;
//window,style="win22.gs";
palette,"idl-33.gp";
c=["red","blue","black","green"];
for(i=1;i<=numberof(ll);i++){

  ll2=exec("ls -d "+ll(i)+"/sn*");
  
  //  plsys,i;
  isnap=3;
  for(j=isnap;j<=isnap;j++){
    readg(ll2(j),x,t,a);
    a;
    //pli,log(1-x(,,64)+1e-6),cmin=-6,cmax=1;
    nc=contour(yc,xc,0.5,x(,,64),span(0,128,129)(zcen)(-:1:128,),span(0,128,129)(zcen)(,-:1:128));
    PL,yc,xc,color=c(i),msize=.1;
  }
 }
limits,0,128,0,128;
