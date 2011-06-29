

// s=readout("out00025.dat");ss=sort(*s.id);PL,(*s.vel)(1,ss)(1000:1300),(*s.pos)(1,ss)(1000:1300),msize=.3;Zeldovich(s.time,2.,x,vx,ng=128);plg,vx,x;


//============================================================
func dplus(a,omegam,omegav)
{
  eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav);
  ai=spanl(a/1e6,a,1024);
  integrand=ddplus(ai,omegam,omegav);
  
  res=eta/a*integ(integrand,ai,ai(0));
  return res;
}

//============================================================
func ddplus(a,omegam,omegav)
{
  eta=sqrt(omegam/a+omegav*a*a+1.-omegam-omegav);
  res=2.5/(eta^3);
  return res;
}

//============================================================
func fomega(a,omegam,omegav)
{
  omegak=1.0-omegam-omegav;
  eta=sqrt(omegam/a+omegav*a*a+omegak);
  res=(2.5/dplus(a,omegam,omegav)-1.5*omegam/a-omegak)/(eta*eta);
  return res;
}

//============================================================
func dladt(a,omegam,omegav)
{
  eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav);
  res=a*eta;
  return res;
}

//============================================================

func flahav(a,omegam,omegav)
{
  z=1./a-1.;
  X=1.+omegam*z+omegav*(1./(1+z)^2-1);
  va=spanl(1e-6,a,1024);
  aX=1.+omegam*(1./va-1)+omegav*(va^2-1);
  res=(omegav*a^2-omegam*0.5/a)/X-1+a*X^(-1.5)/integ(aX^(-1.5),va,va(0));
  return res;
  
    
}

func Hubblea(a,omegam,omegav)
{
  return sqrt(omegam/a^3+omegav+(1-omegav-omegam)/a^2);
}


func Zeldovich(aini,across,&x,&vx,ng=,omegam=,omegav=)
{
  if(is_void(ng))ng=64;
  if(is_void(omegam)) omegam=0.3;
  if(is_void(omegav))  omegav=0.7;
  amp=1./(dplus(across,omegam,omegav)*2*pi/ng);
  amp;
  vfact=dplus(aini,omegam,omegav);
  x=(indgen(ng)-0.5)/ng+vfact*amp*sin((indgen(ng)-0.5)*2*pi/ng)/ng;
  
  vfact=fomega(aini,omegam,omegav)*Hubblea(aini,omegam,omegav)*dplus(aini,omegam,omegav);
  vfact;
  vx=aini^2*amp*sin(2*pi*(indgen(ng)-0.5)/ng)*vfact/ng;
  return 0;
}

