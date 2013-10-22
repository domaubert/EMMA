// SOME ROUTINES FOR PLOTTING

func splitplot2(y,x,cut=,color=,width=,type=,noerase=,rg=,nmax=)
/* DOCUMENT 
     splitplot(y,x,cut=,color=,width=,type=,noerase=)
   will split over numberof(y)/cut subplot the curve y(x)
   (at most 7x4 plots)
   SEE ALSO:
 */
{
  if(is_void(cut)) cut=500;
  if(is_void(nmax)) nmax=20;
  nn=numberof(y);
 if(is_void(x)) x=indgen(nn);
   qq=min(long((nn/cut/2))+1,nmax);
  if(nn>20*4*cut) write,"WARNING: the plot is truncated\n";
  for(j=0;j<=qq-1;j++)
    {
	if(j*2*cut==nn) continue;
     if(is_void(noerase)) window,j,style="struc.gs"; else   window,j;
         for(i=0;i<=min(2,nn/cut-j*2)-1;i++)
        {
          plsys,2-i;
      if((i+j*2)*cut==min(((i+j*2)+1)*cut,nn)) continue;
    	  plh,y((i+j*2)*cut+1:min(((i+j*2)+1)*cut,nn)),
      x((i+j*2)*cut+1:min(((i+j*2)+1)*cut,nn)),color=color,width=width,type=type;
  if(!is_void(rg)) range,rg(1),rg(2)   ;

  }
	 redraw; 
    }
    
}


func wsp(win,dpi=,width=,height=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100) in 2 plsys per page mode
   include logxy+animate

   DElADAUBE
   
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=75;
  if (is_void(width)) width=long(18.*dpi)+1;
  if (is_void(height)) height=long(6.54*dpi)+1;
 if (is_void(win) && (win= current_window())<0) win= 0;
 // window, win;   fma;  animate, 0; 
window, win,dpi=dpi,width=width,height=height,style="struc.gs";  limits;  fma;
}

func wspJ(ne){
/* DOCUMENT
   opens a spectrum plotting window
*/
//  window,ne,style="nst.gs",width=800,height=500;
  window,ne,style=usdir+"/Yorick/STECKMAP/Gist/nst.gs",width=800,height=500;
};

func wspJ2(ne){
/* DOCUMENT
  a longer, thinner version of wspJ, probably better for publications and visu
  looks nice but legends are missing on the left x>1 in system 0
  
*/
//  window,ne,style="nst.gs",width=800,height=500;
  window,ne,dpi=45,style=usdir+"/Yorick/STECKMAP/Gist/nst2.gs",width=1000,height=500;
};

func wspJ3(ne){
/* DOCUMENT
  a longer, thinner version of wspJ, probably better for publications and visu
*/
//  window,ne,style="nst.gs",width=800,height=500;
  window,ne,dpi=90,style=usdir+"/Yorick/STECKMAP/Gist/nst3.gs",width=800,height=200;
};




func plJ(a,b,color,fun=,type=){

  /* DOCUMENT
     nice plots for spectra , see wspJ
  */
  if(is_void(color)) color="black";
  if(is_void(fun)) fun="plh";
  if(is_void(type)) type=1;

  
  
  plsys,2;

  if(fun=="plg") plg,a(:int(numberof(a)/2.)),b(:int(numberof(a)/2.)),color=color,type=type;
  if(fun=="plh") plh,a(:int(numberof(a)/2.)),b(:int(numberof(a)/2.)),color=color,type=type;
  plsys,1;
  if(fun=="plh") plh,a(int(numberof(a)/2.):),b(int(numberof(a)/2.):),color=color,type=type;
  if(fun=="plg") plg,a(int(numberof(a)/2.):),b(int(numberof(a)/2.):),color=color,type=type;
  
};

func plJ2(a,b,wins,wine,color,init=,fun=,cmin=,cmax=){

  /* DOCUMENT
     nice plots for spectra , see wspJ
     in n wspJ windows
     wins is the number of first window
     wine is the number of last window
  */

  if(init==1) {for(win=wins;win<=wine;win++){ws,win;wspJ;};};

  na=numberof(a);
  nw=wine-wins;
  nia=int(na/(2*(nw+1)));
  for(i=0;i<=nw;i++){
    window,(i+wins);
    plsys,2;
    if(fun=="plg") plg,a(2*nia*i+1:2*nia*i+nia),b(2*nia*i+1:2*nia*i+nia),color=color;
    if(fun!="plg") plh,a(2*nia*i+1:2*nia*i+nia),b(2*nia*i+1:2*nia*i+nia),color=color;
    if(!is_void(cmin)) range,cmin,cmax;
    plsys,1;
    if(fun=="plg") plg,a(2*nia*i+nia+1:2*nia*i+2*nia),b(2*nia*i+nia+1:2*nia*i+2*nia),color=color;
    if(fun!="plg") plh,a(2*nia*i+nia+1:2*nia*i+2*nia),b(2*nia*i+nia+1:2*nia*i+2*nia),color=color;
    if(!is_void(cmin)) range,cmin,cmax;

  };

  for(i=0;i<=nw;i++){window,(i+wins);replotAll;};
};


func spepsi(wins,wine,name,&l){
  /* DOCUMENT
     save from window wins to wine
     returns the filenames in l;
  */
  l=[];
  for(win=wins;win<=wine;win++){window,win;epsi,name+"-"+pr1(win),nodisp=1,norescale=1;grow,l,name+"-"+pr1(win)+".eps";};
};

func speps(wins,wine,name,&l){
  /* DOCUMENT
     save from window wins to wine
     returns the filenames in l;
  */
  l=[];
  for(win=wins;win<=wine;win++){window,win;eps,name+"-"+pr1(win)+".ps",nodisp=1;grow,l,name+"-"+pr1(win)+".ps";};
};

func sphcpq(wins,wine,name,&l){
  /* DOCUMENT
     save from window wins to wine
     returns the filenames in l;
  */
  l=[];
  for(win=wins;win<=wine;win++){window,win;hcpq,name+"-"+pr1(win);grow,l,name+"-"+pr1(win)+".ps";};
};




func myexportfig(ll,fname=,pp=,ps=,nodisp=)
/* DOCUMENT exportfig(ll,fname=,pp=,ps=) creates a latex file (default crap.tex)
   with a table containing the figures in the list ll
   EXAMPLE
   #include "Eric/system.i"
   ll =exec("ls clusters/*.jpg");
   exportfig(ll(1:21));
   or
   exportfig(ll(1:21),pp=3);      
   KEYWORDS
   fname = tex file name
   pp = number of fig per line (default 4 )
   SEE ALSO: dumppdf,eps,pdf,jpg
*/
{
if (is_void(fname))  fname="crap.tex";
  ff =open(fname,"w")

    ll=ll;
  //  ll=ll(::-1);
    
    nn= numberof(ll);
  //  for(i=1;i<=nn;i++) ll(i)=_filterdot(ll(i));
  if (is_void(pp))  pp=4;
 if(pp==3){ qq=9; q1=0.40; q2=0.25;}
 if(pp==4){ qq=24; q1=1.; q2=0.2;}
 if(pp==2){ qq=6; q1=0.6; q2=0.32;}
 
write,ff,"\\documentclass{article}";
write,ff,"\\usepackage[T1]{fontenc}";
write,ff,"\\usepackage[latin1]{inputenc}";
write,ff,"\\usepackage{graphics}";
write,ff,"\\textwidth=16.5cm";
write,ff,"\\textheight=27.5cm"
write,ff,"\\topmargin=-3.cm";
write,ff,"\\evensidemargin=-2.5cm";
write,ff,"\\oddsidemargin=0cm";    
write,ff,"\\begin{document}";

 for(k=0;k<=nn/qq;k++)
   {//write,ff,"{\\centering \\begin{tabular}{c c c c}";
     write,ff,"{\\centering "
     //     write,ff,"{";
   
     tt=((k<nn/qq)?qq:(nn%qq));
 for(i=1;i<=tt ;i++)
   { 
     if(is_void(ps)) tt1=strpart(ll(i+k*qq),1:-4); else tt1=ll(i+k*qq);
  write,ff,
            format="\\resizebox*{%2f\\columnwidth}{%2f\\textheight}{\\rotatebox{-90}{\\includegraphics{%s}}} ",
    //format="\\resizebox{%2f\\textheight}{%2f\\columnwidth}{\\rotatebox{-90}{\\includegraphics{%s}}} ",
                       q1,q2,tt1;
      if (i % pp ) {write,ff," ";} else if (i<nn) {  write,ff,"  ";} else {write,ff,"";} 
    }

 write,ff,"\\par}";
   }
 write,ff,"\\vspace{0.3cm}";
write,ff,"\\end{document}";
close,ff;

 col="cp "+fname+" ./crap.tex"; 
 tra=exec(col);
 col="latex crap.tex";
 tra=exec(col);
 //dvina=strreplace(fname,"tex","dvi");
 psna=strreplace(fname,"tex","eps");
 //col="cp crap.dvi "+pr1(dvina);
 //exec(col);
 // col="dvips -o "+pr1(psna)+" "+pr1(dvina);
 //exec(col);
 //col="gv "+psna;
 //exec(col);
 
//$ latex crap.tex
 
$ dvips -o crap.ps crap.dvi
  if(is_void(nodisp)) exec("gv -upsidedown crap.ps");     
$pwd

  col="cp -f crap.ps "+psna;
 exec(col); 
  
 return fname;
};






func markLick(namel,u,names,yi=,ys=,color=,wname=,nv=,ylabel=,height=){
  /* DOCUMENT
     plots shades in bands contained in u, names containing the names of the indices having the bands in u
     see invmxpaper.i for an example
     upload lick4paper     first
  */
  
  local ind;
  if(is_void(nv)) nv=2;
  if(is_void(ylabel)) ylabel=0.3;
  if(is_void(yi)) yi=0.;
  if(is_void(ys)) ys=2.;
//  if(is_void(color)) color=__gray80;
  if(is_void(color)) color=char([200,200,200]);

  
  if(is_void(wname)) wname=[];
  if(numberof(namel)>1) {
    markLick,namel(1),u,names,yi=yi,ys=ys,color=color,wname=wname,ylabel=ylabel,height=height;
    markLick,namel(2:),u,names,yi=yi,ys=ys,color=color,wname=wname,ylabel=ylabel,height=height;
    return;
  };
  ind=where(names==namel);
  if (is_void(dimsof(ind))) return "bad name";
  plshade,[array(yi,nv),array(ys,nv)],span(u(1,ind),u(2,ind),nv),color=color,edge=[];
  if (!is_void(wname)) plt1,namel,u(2,ind),ylabel,orient=1,tosys=1,height=height;
  return [];
};

func mlabel(s,p1,p2,c=,t=,w=,h=,np=,dx=,dy=,auto=,yshift=,xshift=){
  /* DOCUMENT
     plots labels in the rectangle defined by p1 (left bottom corner) p2 (top right corner), equally spaced, with colors specified in array c, types in array t, label in array s
     in this routine the p1 and p2 must be specified in the same coordinates as the plot (system=1)
     s the string array containg the labels.
     c can be char colors (char([200,200,200]) for example since 15/07/2008
     
     FIX ME for slog!!
     
     see also PLmlabel
  */

  local nl,np,xl,yl;
  slog=0;
  if(is_void(sys)) sys=1;
  if(is_void(np)) np=20;
  nl=numberof(s);
  if(is_void(yshift)) yshift=0.04;
  if(is_void(xshift)) xshift=0.1;

  if(auto==1){
    if(is_void(p1)) p1=0.1;
    if(is_void(p2)) p2=[0.9-yshift*(nl-1),0.9];
    bla=limits();
    xlab1=p1*(bla(2)-bla(1))+bla(1);
    xlab2=(p1+xshift)*(bla(2)-bla(1))+bla(1);
    ylab1=p2(1)*(bla(4)-bla(3))+bla(3);
    ylab2=p2(2)*(bla(4)-bla(3))+bla(3);
    p1=[xlab1,ylab1];
    p2=[xlab2,ylab2];
  };
  
  if(is_void(c)) c="black";
  if(dimsof(c)(0)!=nl) {
      c=c(-:1:nl);
      write,"assign same colors everywhere";
  };
  if(is_void(h)) h=1;
  if(numberof(h)!=nl) h=h(-:1:nl);
  if(is_void(t)) t=1;
  if(numberof(t)!=nl) t=t(-:1:nl);
  if(is_void(w)) w=1;
  if (numberof(w)!=nl) w=w(-:1:nl);
  if(is_void(dx)) dx=(p2(1)-p1(1))/15.;
  if(is_void(dy)) dy=-(p2(2)-p1(2))/float(nl)/5.;
  //  if(ayto==1) 

  
  xa=span(p1(1),p2(1),np);
  ya=span(p2(2),p1(2),nl)(,-:1:np);
  
  for(i=1;i<=nl;i++){
      dimsof(c)(1);
//      c(,i); // debug only
      plh,ya(i,),xa,type=t(i),color=(dimsof(c)(1)<2?c(i):c(,i)),width=w(i); // changed to handle char values for colors (example char([200,200,200])
    //    plt,s(i),xa(0)+dx,ya(i,1)+dy,tosys=1;
    add_txt,s(i),pos=[xa(0)+dx,((slog==0)?ya(i,1)+dy:ya(i,1)*(1.+dy))],sys=1,height=h(i);

    
  };
  if (deb==1) error;
};

func mmlabel(s,c=,t=,w=,np=,dx=,dy=){
  /* DOCUMENT
     plots labels using function mlabel with same arguments except coordinates of the box are taken by mouse
  */
  a=mouse(1,1,"draw a box");
  //write,a(1:4);
  a([1,3])=a([1,3])(sort(a([1,3])));
  a([2,4])=a([2,4])(sort(a([2,4])));
  mlabel(s,a([1,2]),a([3,4]),c=c,t=t,w=w,np=np,dx=dx,dy=dy);
  return a(1:4);
};
  
func mPLmlabel(s,c=,t=,w=,np=,dx=,dy=,inc=,m=,ms=){
  /* DOCUMENT
     plots labels using function PLmlabel with same arguments except coordinates of the box are taken by mouse
  */
  a=mouse(1,1,"draw a box");
  a([1,3])=a([1,3])(sort(a([1,3])));
  a([2,4])=a([2,4])(sort(a([2,4])));
  //write,a(1:4);
  PLmlabel(s,a([1,2]),a([3,4]),c=c,t=t,w=w,np=np,dx=dx,dy=dy,inc=inc,m=m,ms=ms);
  return a(1:4);
};

func PLmlabel(s,p1,p2,c=,t=,w=,np=,dx=,dy=,inc=,m=,ms=,h=,slog=,auto=,yshift=,stretch=){
  /* DOCUMENT
     plots labels for symbols (PL style) in the rectangle defined by p1 (left bottom corner) p2 (top right corner), equally spaced, with colors specified in array c, types in array t, label defined in array s for PL, printed height in array h
     Fine-tune position with dx,dy
     slog=1 -> positions are logarithmic, useful when logxy,1,1
     see also mlabel
     stretch is useful when using auto=1: stretch = 1 will give 
  */

  local nl,np,xl,yl;
  nl=numberof(s);
  if(is_void(yshift)) yshift=0.04;
  if(is_void(auto)) auto=0;
  if(is_void(stretch)) stretch=1.;
  if(auto==1){
    if(is_void(p1)) p1=0.1;
    if(is_void(p2)) p2=[0.9-yshift*(nl-1),0.9];
    bla=limits();
    xlab=p1*(bla(2)-bla(1))+bla(1);
    ylab1=p2(1)*(bla(4)-bla(3))+bla(3);
    ylab2=p2(2)*(bla(4)-bla(3))*stretch+bla(3);
    p1=[xlab,ylab1];
    p2=[xlab,ylab2];
  };
  
  if(is_void(np)) np=20;
  if(is_void(c)) c="black";
  if(numberof(c)!=nl) c=c(-:1:nl);
  if(is_void(h)) h=1;
  if(numberof(h)!=nl) h=h(-:1:nl);
  if(is_void(t)) t=1;
  if(numberof(t)!=nl) t=t(-:1:nl);
  if(is_void(w)) w=1;
  if (numberof(w)!=nl) w=w(-:1:nl);
  if (is_void(inc)) inc="white";
  if (!is_void(inc)) { if ((numberof(inc)!=nl)) inc=inc(-:1:nl);};
  if(is_void(m)) m=1;
  if (numberof(m)!=nl) m=m(-:1:nl);
  if(is_void(ms)) ms=1;
  if (numberof(ms)!=nl) ms=ms(-:1:nl);
  if(is_void(dx)) dx=(p2(1)-p1(1))/15.;
  if(is_void(dy)) dy=-(p2(2)-p1(2))/float(nl-1)/3.;
  if(is_void(slog)) slog=0;
  
  xa=p1(1);
  ya=(slog==0)?span(p2(2),p1(2),nl):spanl(p2(2),p1(2),nl);
  
  
  for(i=1;i<=nl;i++){
    PL,ya(i,),xa,type=t(i),color=c(i),width=w(i),incolor=inc(i),msize=ms(i),marker=m(i);
    //    write,inc(i);
    //plt,s(i),xa(0)+dx,ya(i,1)+dy,tosys=1;
    //write,xa;
    //write,i;
    add_txt,s(i),pos=[xa(0)+dx,((slog==0)?ya(i,1)+dy:ya(i,1)*(1.+dy))],sys=1,height=h(i);
    
  };
  if (deb==1) error;
};



func xyleg(s1,s2,s3,cs1=,cs2=,cs3=,dx1=,dy1=,dx2=,dy2=,dx3=,dy3=){
  /* DOCUMENT
     default position is  0.37,0.37;  0.13,0.62;
     default text is log(age[yr])
  */

  
  if(is_void(s1)) s1="log(age[yr])";
  if(is_void(cs1)) cs1=[0.37,0.37];
  if(is_void(cs2)) cs2=[0.13,0.62];
  if(is_void(cs3)) cs3=[0.64,0.62];
  if(is_void(dx1)) dx1=0.;
  if(is_void(dx2)) dx2=0.;
  if(is_void(dx3)) dx3=0.;
  if(is_void(dy1)) dy1=0.;
  if(is_void(dy2)) dy2=0.;
  if(is_void(dy3)) dy3=0.;

  
  
  plt,s1,cs1(1)+dx1,cs1(2)+dy1,tosys=0;
  if (!is_void(s2)) plt,s2,cs2(1)+dx2,cs2(2)+dy2,tosys=0,orient=1;
  if(!is_void(s3)) plt,s3,cs3(1)+dx3,cs3(2)+dy3,tosys=0,orient=1;
};



func mypler(y,x, dy=,dx=,marker=, width=, color=,incolor=,msize=,
          mean=,skew=,kurto=,med=,yweight=,xweight=,save=,noerror=,line=,logx=,logy=,type=,trimmed=,errorweight=
          ,shading=,ticks=)
  /* DOCUMENT

     pler modified to account for ticks. No big deal
     
   pler(y,x, marker=, width=, color=,incolor=,
   msize=,mean=,skew=,kurto=,med=)
   plot without fuss with error bars x and y
   can be either arrays or sorted lists
   mean= divides by sqrt(number of data points in bin)
   skew= computes skewness as well
   kurto= computes kurtosis as well
   med= computes median and quartile instead of mean and rms
   trimmed= trims 100*trimmed points on both side of the distribution
   yweight= can be a list a vector or an array of dimsof(y)
   xweight= itou
   noerror= doesn't plot error bar
   line= uses lines to join points
   logx= plots in log  without logxy
   logy=  itou
   note color=-1   and/or marker will pick random colors resp. marker
   SEE ALSO: splitstat,median,quartile
 */
{
  
  if(trimmed)
    {
      yy=xx=[];
      wwy=wwx=[];
      if(typeof(y)=="list")
          for(i=1;i<=_len(y);i++)
       {       yy=_cat(yy,trim(_car(y,i),wy,lcut=trimmed,cutedges=1));
       if(typeof(yweight)=="list") wwy=_cat(wwy,trim(_car(yweight,i)(wy),lcut=trimmed,cutedges=1));}
      else
      if(typeof(x)=="list")
          for(i=1;i<=_len(x);i++)
       {     xx=_cat(xx,trim(_car(x,i),wx,lcut=trimmed,cutedges=1)); 
       if(typeof(xweight)=="list") wwx=_cat(wwx,trim(_car(xweight,i)(wx),lcut=trimmed,cutedges=1));}
      else
      for(i=1;i<=dimsof(y)(2);i++)
        {
          yy=_cat(yy,trim(y(i,),wy,lcut=trimmed,cutedges=1)); 
        if(is_array(yweight)) wwy=_cat(wwy,yweight(wy));
        if (numberof(x))    if (dimsof(x)(1)==2)
          {
            xx=_cat(xx,trim(x(i,),wx,lcut=trimmed,cutedges=1));
            if(is_array(xweight)) wwx=_cat(wwx,xweight(wx));
          }
        }
      
      if (numberof(x)) if (dimsof(x)(1)==2)  x=xx;
      ss=pler(yy,x,marker=marker, width=width, color=color,incolor=incolor,msize=msize,
        mean=mean,skew=skew,kurto=kurto,med=med,yweight=wwy,xweight=wwx,save=save,noerror=noerror,
        line=line,logx=logx,logy=logy,type=type,errorweight=errorweight);
      
     if(save) return ss;
     return;
    }

  if( is_void(dy) && is_void(dx))
    {

      if (typeof(y)=="list")
        {
          if(typeof(yweight)=="list")
            {
              
              yw=yweight;  sw=[]; nn=_len(yw); for(i=1;i<=nn;i++){grow,sw,(1.*_nxt(yw))(sum); }
              w1=where(!sw); if(is_array(w1)) sw(w1)=1.;
              yy=y; yw=yweight; 
            y1=[]; for(i=1;i<=nn;i++){ grow,y1,(1./sw(i)*_nxt(yw)*_nxt(yy))(sum);}
             yy=y;yw=yweight;
            ye=[]; for(i=1;i<=nn;i++){ grow,ye,sqrt((1./sw(i)*_nxt(yw)*(_nxt(yy)-y1(i,-))^2)(sum));}
            if(!is_void(mean)) {
              yw=yweight;  yn=[]; for(i=1;i<=nn;i++){ grow,yn,numberof(_nxt(yw));}
              ye = ye/sqrt(yn); }
            }
          else
            {
          if(is_void(med)) y1=lst2arr(_map(AVG,y)); else {  y1=lst2arr(_map( median,y));};
          if(is_void(med)) ye=lst2arr(_map(RMS,y)); else ye=lst2arr(_map(quartile,y));
          yn=lst2arr(_map(numberof,y));
          if(!is_void(skew))     ys=lst2arr(_map(skewness,y));
          if(!is_void(kurto))     yk=lst2arr(_map(kurtosis,y));
          if(!is_void(mean)) ye = ye/sqrt(yn);
            }
        }
      else
        {
          if(!is_void(yweight))
            {
              if (dimsof(yweight)(1) == 1)
                {
                  yweight *= 1/sum(double(yweight));
                  if(is_void(med))
                    {
                      y1=(y*yweight(-,))(,sum);
                      ye=sqrt((yweight(-,)*((y-y1(,-))^2))(,sum));
                      if(!is_void(mean)) ye = sqrt(((yweight(-,))^2*((y-y1(,-))^2))(,sum));
                    }
                }
              else
                {
                  yweight *=1. /double(yweight(,sum));
                  y1=(y*yweight)(,sum);
                  ye=sqrt((yweight*((y-y1(,-))^2))(,sum));
                  if(!is_void(mean)) ye = sqrt((yweight^2*((y-y1(,-))^2))(,sum));
                }
            }
          else
            {
              if(is_void(med))
                {
                  ye=y(,rms);    y1=y(,avg);
                  if(!is_void(skew)) ys=lst2arr(_map(skewness,arr2lst(transpose(y))));
                  if(!is_void(kurto)) yk=lst2arr(_map(kurtosis,arr2lst(transpose(y))));

                }
              else
                {
                  y1=median(y,2);ye=quartile(transpose(y));
                }
              if(!is_void(mean)) ye = ye/sqrt(dimsof(y)(0));
            }
        }
      if (is_void(x)) x=1.*indgen(dimsof(y1)(0));

      if (typeof(x)=="list")
        {
          if(typeof(xweight)=="list")
            {
              xw=xweight;  sw=[]; nn=_len(xw); for(i=1;i<=nn;i++){ grow,sw,(1.*_nxt(xw))(sum);}
             w1=where(!sw); if(is_array(w1)) sw(w1)=1.;
             xx=x; xw=xweight;
            x1=[]; for(i=1;i<=nn;i++){ grow,x1,(1./sw(i)*_nxt(xw)*_nxt(xx))(sum);}
            xx=x;xw=xweight;
            xe=[]; for(i=1;i<=nn;i++){ grow,xe,sqrt((1./sw(i)*_nxt(xw)*(_nxt(xx)-x1(i,-))^2)(sum));}
            if(!is_void(mean)) {
              xw=xweight;  xn=[]; for(i=1;i<=nn;i++){ grow,xn,numberof(_nxt(xw));}
              xe = xe/sqrt(xn); }
            }
          else
            {
              if(is_void(med))       x1=lst2arr(_map(AVG,x)); else x1=lst2arr(_map(median,x));
              if(is_void(med))       xe=lst2arr(_map(RMS,x)); else xe=lst2arr(_map(quartile,x));
              xn=lst2arr(_map(numberof,x));
              if(!is_void(skew))     xs=lst2arr(_map(skewness,x));
              if(!is_void(kurto))     xk=lst2arr(_map(kurtosis,x));
            if(!is_void(mean)) xe = xe/sqrt(xn);
            }
        }
      else  if (dimsof(x)(1) > 1) 
        {
          if(!is_void(xweight)){
            if (dimsof(xweight)(1) == 1)
              {
                xweight *= 1/sum(double(xweight));
                if(is_void(med))
                  {
                    x1=(x*xweight(-,))(,sum);
                    xe=sqrt((xweight(-,)*((x-x1(,-))^2))(,sum));
                    if(!is_void(mean)) xe = sqrt(((xweight(-,))^2*((x-x1(,-))^2))(,sum));
                  }
              }
            else
              {
                xweight *=1. /double(xweight(,sum));
                x1=(x*xweight)(,sum);
                xe=sqrt((xweight*((x-x1(,-))^2))(,sum));
                if(!is_void(mean)) xe = sqrt((xweight^2*((x-x1(,-))^2))(,sum));
              }
          }
          else
            {
              if(is_void(med))
                {
                  x1=x(,avg);  xe=x(,rms);
                   if(!is_void(skew)) xs=lst2arr(_map(skewness,arr2lst(transpose(x))));
                  if(!is_void(kurto)) xk=lst2arr(_map(kurtosis,arr2lst(transpose(x))));
                }
              else
                {
                  x1=median(x,2);
                  xe=quartile(transpose(x));
                }
              if(!is_void(mean)) xe =xe/sqrt(dimsof(x)(0)); 
            }
        }
      else {
        x1=x; xe=0.;
      }
    } else{
      if(is_void(x)) x=indgen(numberof(y));
      if(!is_void(dx))  xe=dx; else xe=x*0;
      if(!is_void(dy))  ye=dy;  else ye=y*0;
      y1=y; x1=x;
    }
  if(errorweight) {ye *= errorweight; xe *= errorweight;}
  
  if (is_void(color)) color=incolor;
if(numberof(color)==1) if(color==-1) { color=Colors(8)(1+long(random()*6)); incolor=color;}

  if(!is_void(logx)) {
    ww=where((x1>0)*(x1>xe)*(x1+xe>0));
    if(is_array(xe>0))
      {xlo=log10((x1-xe)(ww));
      xhi=log10((x1+xe)(ww));} else
        {xlo=xhi=log10(x1);}
    x1=log10(x1(ww));
    y1=y1(ww);
  }
  if(!is_void(logy)) {
    ww=where((y1>0)*(y1>ye)*(y1+ye>0));
    if(is_array(xe>0) ){
      ylo=log10((y1-ye)(ww));
      yhi=log10((y1+ye)(ww));
    } else {ylo=yhi=log10(y1);}
    y1=log10(y1(ww));
    x1=x1(ww);
  }

  if(is_void(xlo)){ xlo=x1-xe; xhi=x1+xe;}
  if(is_void(ylo)){ ylo=y1-ye; yhi=y1+ye;}
  if(is_void(noerror)&is_void(shading)) 
    {   plp,y1,x1,xlo=xlo,xhi=xhi,ylo=ylo,yhi=yhi,ticks=ticks,color=color,width=width,type=1,symbol=0;
  PL,y1,x1,marker=marker,incolor=incolor,width=width,color=color,msize=msize;}
else
  if(is_void(noerror)&!is_void(shading))  plshade,[ylo,yhi],x1,color=color,edge=1;
 else
   PL,y1,x1,marker=marker,incolor=incolor,width=width,color=color,msize=msize;
  if(!is_void(line)) plg,y1,x1,width=width,color=color,type=type;
  if(!is_void(skew))
    {
      PL,y1+ye*ys,x1,marker=5,incolor=incolor,width=width,color=color,msize=0.4;
      if(numberof(xs))  PL,y1,x1+xe*xs,marker=5,incolor=incolor,width=width,color=color,msize=0.4;
      if(!is_void(save)& is_void(kurto)) return [x1,xe,y1,ye,ys];
    }
  if(!is_void(kurto))
    {
      PL,y1+yk*ye,x1,marker=6,incolor=incolor,width=width,color=color,msize=0.3;
      PL,y1-yk*ye,x1,marker=6,incolor=incolor,width=width,color=color,msize=0.3;
      if(numberof(xk))
        {
          PL,y1,x1+xe*xk,marker=6,incolor=incolor,width=width,color=color,msize=0.4;
          PL,y1,x1-xe*xk,marker=6,incolor=incolor,width=width,color=color,msize=0.4;
        }
      if(!is_void(save) & is_void(skew)) return [x1,xe,y1,ye,yk];
      if(!is_void(save)) return [x1,xe,y1,ye,ys,yk];
    }
  if(!is_void(save) & !is_void(noerror)) return [x1,y1];
  if(!is_void(save)) return [x1,xe,y1,ye];
}

func blank(filename){
  /* DOCUMENT
     prints a blank 
     => filename has no suffix, is eps
     marche pas
  */
  local a;
  extern com;
  filename=filename+".eps";
  com="cp /home6/ocvirk/perso/modeles/galaxie/POP/SDSS/blank.eps "+filename;
  a=exec(com);
};

func myPL(y,x,marker=,width=,color=,incolor=,msize=,line=,type=){
  /* DOCUMENT
     modofication of PL
     msize can be a vector of same size as y
     
  */
  local nx;
  nx=numberof(x);
  for(i=1;i<=nx;i++){
    PL,y(i),x(i),marker=marker,width=width,color=color,incolor=incolor,msize=msize(i),line=line,type=type;
  };
};


func plcc(y,x,z,zbins=,w=,l=){
  /* DOCUMENT
     plots y versus x, codes z with color in 4 bins=> zbins is only a 3-array, monotoniclally increasing, w:width, l=1 => labels
     PUT IT IN plot2ls.i
  */
  if (is_void(w)) w=1;
  if(is_void(zbins)) zbins=[z(min),z(avg),z(max)];
  //if(is_void(zbins)) zbins=[(z(min)+z(max))/2.-z(rms),(z(min)+z(max))/2.,(z(min)+z(max))/2.+z(rms)];

  
  i1=where(z<=zbins(1));
  i2=where((z>=zbins(1))&(z<=zbins(2)));
  i3=where((z>=zbins(2))&(z<=zbins(3)));
  i4=where(z>=zbins(3));
  
  PL,y(i4),x(i4),marker=4,color="red",width=w;
  PL,y(i3),x(i3),marker=3,color="black",width=w;
  PL,y(i2),x(i2),marker=2,color="blue",width=w;
  PL,y(i1),x(i1),marker=1,color="green",width=w;


  

  if (!is_void(l)) {
  plt,("k<"+pr1(zbins(1))),0.5,0.82,color="green";
  plt,pr1(zbins(1))+"<k<"+pr1(zbins(2)),0.5,0.805,color="blue";
  plt,pr1(zbins(2))+"<k<"+pr1(zbins(3)),0.5,0.79,color="black";
  plt,pr1(zbins(3))+"<k",0.5,0.775,color="red";
  };
};

func mytitle(title,dx=,dy=){
  if(is_void(dx)) dx=0.22;
  if(is_void(dy)) dy=0.875;
  plt,title,0.22,0.875,tosys=0;
  return;
};
  


func hcpq(file){
  /*DOCUMENT
    wrapper allowing printing using hcp.
    does hcp_file(file);hcp;hcp_finish;
  */
  hcp_file(file);hcp;hcp_finish;
};


func pfits(file){
  /*DOCUMENT
    specific to the file format decided for comparing fits with Mina
    loads the fits file and plots in current display
  */

  local a;
  a=fits_read(file,h);
  ws;
  plh,a(,2),a(,1);
  plh,a(,3),a(,1),color="red";
  plh,a(,4),a(,1),color="green";
  plh,a(,5),a(,1),color="blue";
  plh,a(,6),a(,1);
  if(!is_void(fits_get(h,"SSPmodel"))) pltitle,fits_get(h,"SSPmodel")+"   fit";
  xyleg,"Angstrom",[];
  return;
};


func ylims(py1,py2){
  /* DOCUMENT
     vertically rescales the display to make some room for labelling the plotted curves
     with PLmlabel
  */

  if (is_void(py2)) py2=1.5;

  local yli,xlims;
  
  yli=limits();
  xlims=yli([1,2]);
  yli=yli([3,4]);
  if(!is_void(py1)) yli(1)*=py1;
  if(!is_void(py2)) yli(2)*=py2;
  limits,xlims(1),xlims(2),yli(1),yli(2);
};
     

func plminmax(y,x,type=,width=,color=){
  /* DOCUMENT
     plots vertical error bars without horizontal ticks corresponding to the min and max of each bin.
  */
  nb=dimsof(y)(2);
  if(is_void(x)) x=indgen(nb);
  for(i=1;i<=nb;i++){
    plh,[min(y(i,)),max(y(i,))],array(x(i),2),type=type,width=width,color=color;
  };
};


func plSFR(file,bab,plsad=,log=,sadrange1=,sadrange2=,plmass=){
  /* DOCUMENT
     plots the sfr and MC in file.
     (Should be used with am absolute flux basis, having bab as binning (such as bab (zcen)=b.ages and also 1./bab = 1/deltat 
     should remedy this by making the input of bRbasis directly M/L / deltat NO THAT DOESNT WORK< CONVERGES BADLY!!) SO IN THE FILE WE HAVE EITHER THE SAD AND TRANSFORM IT INTO SFR OR WE PUT THE SFR IN THE RESULT FILE DIRECTLY
     sadrange1-2 keywords are used as the range of the plot
     plsad if the SAD is to be plotted. In that case no need for bab, a dummy bab is generated.
     NB: Its easy to switch between results in file is SAD or SFR by setting MsL to its real value or to an array of 1. The same holds for bab.
     USES MsLinterp, so MsL should have a value, for example MSL(base,filter)
     NB also: SFRs are normalized to 1
  */
  if(is_void(N)) N=1;
  upload,file,s=1;
  nab=numberof(ages);
  if(plsad==1) bab=log10(indgen(numberof(bab)));
  if(plmass==1) bab=log10(indgen(numberof(bab)));
  sads=gres(:2*nab,);
  sads(:nab,)=sads(:nab,)^2;
  sfrs=SAD2SFR(sads(:nab,),sads(nab+1:2*nab,),10^bab);
  masses=SAD2MASS(sads(:nab,),sads(nab+1:2*nab,));

  
  gres=sfrs;
  if(plsad==1) gres=sads;
  if(plmass==1) gres=masses;
  
  //  deltat=((10^bab)(dif))(,-:1:nMC);
  ws;
  //  gres(:nab,)/=sqrt(deltat);
  if(N==1) gres(:nab,)=gres(:nab,)/((sqrt(((gres^2)(:nab,))(sum,)))(-:1:nab,));
  //  gres(:nab,)/=sqrt(deltat);
  plh,median((gres)(:nab,),2),log10(ages)+6,width=3;
  plh,(gres)(:nab,1),log10(ages)+6,width=3,color="red";
  plminmax,(gres)(:nab,),log10(ages)+6;
  if(is_void(plsad)&(is_void(plmass))){
    xyleg,"log(age[yr])","normalized SFR (Msol/yr)";
    plt,"SFR",0.218,0.8,tosys=0;
  };
  if(plsad==1) {
    xyleg,"log(age[yr])","flux fractions";
    plt,"SAD",0.218,0.8,tosys=0;
  };
  if(plmass==1) {
    xyleg,"log(age[yr])","Msun";
    plt,"MASS",0.218,0.8,tosys=0;
  };
     
  pltitle,dfile;
  ylims;
  if(log==1) logxy,0,1;
  if(log==1) range,1.e-10,1.e-1;
  if(!is_void(sadrange1)) range,sadrange1,sadrange2;
  //  if(log!=1) ylims;
};

func plAMR(file){

  /* DOCUMENT
     warning, will bug for resfiles without MC,
  */
  
  upload,file,s=1;
  nab=numberof(ages);
  ws;
  //  gres(:nab,)=gres(:nab,)/((sqrt(((gres^2)(:nab,))(sum,)))(-:1:nab,));
  plh,log10(Zrescalem1((median((gres)(nab+1:2*nab,),2)))/0.02),log10(ages)+6,width=3;
  plminmax,log10(Zrescalem1((gres)(nab+1:2*nab,))/0.02),log10(ages)+6;
  xyleg,"log(age[yr])","log10(Z/Z_sun_)";
  pltitle,dfile;
  plt,"AMR",0.218,0.8,tosys=0;
  if(dimsof(gres)(3)>1) plt,"avg!c!^2_ ="+pr1(_ki(2:)(avg)),0.218,0.78,tosys=0;
  plt,"!c!^2__1="+pr1(_ki(1)),0.218,0.76,tosys=0;
  
  range,-2.,0.5;

};

func prSAD_AMR(filelist,b,dir=,yp=,log=,plsad=,sadrange1=,sadrange2=,DEBUG=){
  /* DOCUMENT
     plots and prints ad psnups SAD/SFR and AMR
     yp= is the first argument of ypsnup: has to be >=2
  */

  local nl;
  //  ws,1;
  nl=numberof(filelist);
  if(yp==1) yp=2;
  
  //  if(is_void(dir)) dir="./PS/";
  //  exec("mkdir  "+dir);
  //  pslistSFR=strreplace(filelist,".res*","-SFR.ps");
  //  pslistamr=strreplace(filelist,".res*","-AMR.ps");

  pslistSFR=filelist+"-SFR.ps";
  pslistAMR=filelist+"-AMR.ps";
  pslistMASS=filelist+"-MASS.ps";
  pslistSAD=filelist+"-SAD.ps";
  
  for(i=1;i<=nl;i++){

    if(DEBUG==1){  // HORRIBLE DEBUG !!
      write,filelist(i);
      if((i==1)|(i==2)|(i==5)|(i==4)) b=bc03;
      if((i==3)|(i==6)) b=miles;
     };   // END OF HORRIBLE DEBUG
    
    write,b.filename;
    bab=b.bab;
    //    MsL=MSL(b,"flat_wide");  // required by MsLinterp
    MsL=b.MsLratio;
    _m=b.met;                // required also by MsLinterp
    if(plsad==1) MsL=MsL*0.+1.;
    
    plSFR(filelist(i),bab,log=log,plsad=[],sadrange1=sadrange1,sadrange2=sadrange2);
    eps,pslistSFR(i),nodisp=1;
    plAMR(filelist(i));
    eps,pslistAMR(i),nodisp=1;
    plSFR(filelist(i),bab,log=log,plsad=1,sadrange1=sadrange1,sadrange2=sadrange2);
    eps,pslistSAD(i),nodisp=1;
    plSFR(filelist(i),bab,log=log,plsad=[],plmass=1,sadrange1=sadrange1,sadrange2=sadrange2);
    eps,pslistMASS(i),nodisp=1;
  };

  if(!is_void(yp)) {
    st=usdir+"/bin/ypsnup -"+pr1(yp)+" \" "+strglue(pslistAMR,sep="  ") +"   "+ strglue(pslistSFR,sep="  ")+" "+strglue(pslistSAD,sep="  ")+" "+ strglue(pslistMASS,sep="  ")+" \" ";
      prac=exec(st);
  };
      
  return pslistSFR;
  
};
    
func plbigsp(filelist,n,fname=,nodisp=){
  /* DOCUMENT
     plots big spectra 
  */
  nl=numberof(filelist);
  for(i=1;i<=nl;i++){
    upload,filelist(i);
    wspJ,1;
    plJ2,d,x0,1,1+n,init=1;
    plJ2,pmodel1,x0,1,1+n,"red";
    //plJ2,W*1.e-2,x0,1,4,"blue";
    if(!is_void(npec)) plJ2,npec,x0,1,1+n,"green";
    window,1;
    pltitle,filelist(i);
    //spepsi(1,3,ll(i)+".bsp",l);
    speps(1,1+n,filelist(i),l);
    ba=myexportfig(l,pp=,ps=1,nodisp=nodisp,fname=(is_void(fname)?[]:fname(i)));
  };

  
  //  plt,"avg!c!^2_ ="+pr1(_ki(avg)),0.218,0.78,tosys=0;
  //  plt,"!c!^2__1="+pr1(_ki(min)),0.218,0.76,tosys=0;
  //  range,-2.,0.5;
};


func mycolorbar(levs){
  /* DOCUMENT
     make color_bar so that colors reflect values in levs rather than index
  */

  local n,colors,levs;
  //  n =min(15,dimsof(z)(3));
  n=numberof(levs);
  colors= bytscl(span(1,n+1,n+1),cmin=0.5 ,cmax=n+1.5);
  color_bar,levs,colors,vert=1,font="timesI",height=12,format="%3.2f";
};

func color_bar(levs, colors, vert=, labs=, adjust=, ecolor=, height=,
               vport=, format=, font=,edges=)
/* DOCUMENT
   MODIF of original color_bar in Eric/plot.i, adding edges option
   see color_bar
         or color_bar, levs, colors
     Draw a color bar below the current coordinate system.  If LEVS is
     not specified uses plfc_levs (set by previous call to plfc).  If
     COLORS is specified, it should have one more value than LEVS,
     otherwise equally spaced colors are chosen, or plfc_colors if
     plfc_levs was used.  With the vert=1 keyword the color bar appears
     to the left of the current coordinate system (vert=0 is default).
     By default, color_bar will attempt to label some of the color
     interfaces.  With the labs= keyword, you can force the labelling
     algorithm as follows: labs=0 supresses all labels, labs=n forces
     a label at every nth interface, labs=[i,n] forces a label at every
     nth interface starting from interface i (0<=i<=numberof(LEVS)).

     You can specify the viewport coordinates by keyword
     VPORT=[xmin,xmax,ymin,ymax]; by default the colorbar is drawn next to
     the current viewport. You can use the ADJUST keyword to move the bar
     closer to (adjust<0) or further from (adjust>0) the viewport.

     You can specify the string format for labels with keyword FORMAT (default
     "%g"), the font type with keyword FONT (default "helvetica") and the font
     height with keyword HEIGHT (default 14 points).
     
   SEE ALSO: plfc
 */
{
  if (is_void(levs)) {
    if (is_void(plfc_levs)) error, "no levels specified";
    levs= plfc_levs;
    n= numberof(levs)+1;
    if (is_void(colors)) colors= plfc_colors;
  } else {
    n= numberof(levs)+1;
    if (is_void(colors)) colors= bytscl(span(1,n,n),cmin=0.5,cmax=n+0.5);
  }
  if (n != numberof(colors))
    error, "numberof(colors) must be one more than numberof(levs)";
  
  if (is_void(vport)) vport= viewport();    
  if (is_void(adjust)) adjust= 0.0;
  dx= dy= 0.0;
  if (vert) {
    x= (vport(2)+adjust+[0.022,0.042])(-:1:n+1,);
    dx= 0.005;
    y= span(vport(3),vport(4),n+1)(,-:1:2);
  } else {
    y= (vport(3)-adjust-[0.045,0.065])(-:1:n+1,);
    dy= -0.005;
    x= span(vport(1),vport(2),n+1)(,-:1:2);
  }
  sys= plsys(0);
  plf,[colors],y,x,edges=edges,ecolor=ecolor, legend="";
  plsys, sys;
  
  if (is_void(labs) || labs(0)>0) {
    if (numberof(levs)>1) {
      dz= levs(dif);
      if (numberof(dz)!=numberof(levs)-1 ||
          anyof((dz>0.0)!=(dz(1)>0.0)) || !dz(1))
        error, "levs must be monotone 1D";
      levs= levs(1:0);
      levs= grow([2*levs(1)-levs(2)],levs,[2*levs(0)-levs(-1)]);
    } else {
      levs= double(levs(1));
      if (!levs) levs= [-1.0,levs,1.0];
      else levs= [0.0,levs,2*levs];
    }
    if (numberof(labs)<2) {
      if (is_void(labs)) labs= (n-1)/4 + 1;
      orig= where(levs<1.0e-9*max(levs(dif)));
      if (numberof(orig)==1) labs= [orig(1)%labs,labs];
      else labs= [(n%labs)/2,labs];
    }
    list= where(indgen(0:n)%labs(2)==labs(1));
    x= x(list,);
    y= y(list,);
    if (is_void(format)) format= "%g";
    labs= swrite(format=format,levs(list));
    plsys, 0;
    pldj, x(,2),y(,2),x(,2)+dx,y(,2)+dy, legend="";
    plsys, sys;
    if (is_void(font)) font= "helvetica";
    plt1, labs,x(,2)+2*dx,y(,2)+2*dy, justify=(vert?"LH":"CT"),
      height=height,
      font=font;
  }
}




func ccolorbar(xmin,xmax,nlabels,edges=,nc=,vert=){
    /* DOCUMENT
       tentative sensible implementation of a continuous color bar with publishing quality, i.e. a continous color gradient, no frame by default (edges=0 by default), and reasonable labeling (integers favored, will work best if xmin and xmax are integers or reasonable as well)
    */

    if(is_void(edges)) edges=0;
    if(is_void(vert)) vert=1;
    if(is_void(nc)) nc=50; // number of colors in the color bar (200->quite continuous already)

    color_bar,span(xmin,xmax,nc*(nlabels-1)+1),labs=[1,nc],edges=edges,vert=vert;
};

    


func pldif(x1,y1,x2,y2){
  /* DOCUMENT
     plots lines/vectors between (x1,y1) and (x2,y2)
  */

  nx=numberof(x1);
  for(i=1;i<=nx;i++){
    xs=[x1(i),x2(i)];
    ys=[y1(i),y2(i)];

    plg,ys,xs;
  };
};
    


func cossmooth(a,level){
    /* DOCUMENT
       cosmetic smoothing:
       use for badly sampled data:
       actually spline representation sould do the same
    */
};



func Logxy(a,b)
/* DOCUMENT
     alternative to logxy with different ticks
   SEE ALSO:
 */
{
  logxy,a,b;
  get_style,land,sys,legs,clegs;
 lens=sys.ticks.vert.tickLen;
  if(b)
    {
      sys.ticks.vert.tickLen(4:5)=0;
      sys.ticks.vert.logAdjMinor=1.1;
    }
  if(a)
    {
      sys.ticks.horiz.tickLen(4:5)=0;
      sys.ticks.horiz.logAdjMinor=1.1;
    }
 set_style,land,sys,legs,clegs;
} 


func wsnodisp(i,file)
  /* DOCUMENT
     creates a non-window for plotting directly into ps
     Useful in batch mode
     but must be used with hcp
  */
{
  window,i,display="",hcp=file;
};



func WSL_nodisp(win,file)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100) in portrait mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  if (is_void(win) && (win= current_window())<0) win= 0;
 if (is_void(dpi))  dpi=100;
 //  window,win,wait=0; winkill;

  width=long(10.54*dpi)+1;
  height=long(7.04*dpi)+1;
  
  window, win,dpi=dpi,style=sdir+"Gist/large.gs",height=height,width=width,display="",hcp=file;
//limits; 
//  fma;  animate, 0; 
//  window, win,wait=0;  limits;  fma; logxy,0,0;
}



