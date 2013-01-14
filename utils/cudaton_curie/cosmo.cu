

float a2t(float a, float omegav, float h0)
{
  return 2./(3.*h0)*asinh(sqrt(a*a*a*omegav/(1-omegav)))/sqrt(omegav);
}

float t2a(float t, float omegav, float h0)
{
  return powf((1.-omegav)/omegav*powf(sinh(3*h0*t*sqrtf(omegav)/2.),2),1./3.);
}

#ifdef COSMO
#ifndef FLAT_COSMO

float a2tgen(float a, float omegam, float omegav, float Hubble0)
{
  float h0=(Hubble0/100.)*3.08568e22/1e3; // H0 s-1 ==> km/s/Mpc
  float omegar=4.165e-5/(h0*h0);
  float omegak=1.-(omegam+omegav+omegar);
  float integ;
  float h=a/NINTEG;
  float acur;
  int i;
  float t=0.;
  for(i=0;i<NINTEG;i++)
    {
      acur=h*(i+0.5);
      integ=1./sqrt(omegak+omegam/acur+omegav*acur*acur+omegar/acur/acur);
      t+=integ;
    }

  t*=h;
  
  t/=Hubble0;

  return t;
}

#endif
#endif
