

nflux_emma=4e16; // photons/sec/kg
dur_emma=5.; //Myrs


#if 1

// --- > SB99 Spectrum Input
vv=asciiRead("~/qtest/utils/mag/sb99_spec_z=0.001.dat");
lam=vv(1,); // angstrom
h=6.626e-34; //SI
c=299792458.; //SI
egy=h*c/(lam*1e-10)/1e-7; // ergs
flux=pow(10.,vv(2:,)); //ergs/s/angstrom
age=grow(grow(span(1,20.,20),span(30,100,8)),span(200,900,8));

nage=numberof(age);

//--> ionizing photons flux

nflux=array(0.,nage);
for(i=1;i<=nage;i++){
 
  nflux(i)=integ(flux(i,)/egy*(lam<912.),lam,lam(0));
 }

nflux/=(1e6*2e30);// photons/sec/kg

// including ages <1Myr
nfluxbis=grow(nflux(1),nflux);
agebis=grow(0,age);

// sb99 average emission over the same duration as emma
nflux_avg_dur_emma=integ(nfluxbis,agebis,dur_emma)/dur_emma;
nflux_avg_dur_emma;
//nflux_avg_dur_emma=integ(nfluxbis,agebis,dur_emma);

// we get the fudge factor to be applied to spectra
//fact_flux=nflux_emma/nflux_avg_dur_emma;
fact_flux=1./nflux_avg_dur_emma;

//Normalised fluxes for magnitudes calculations

flambda=flux*fact_flux;// ergs/angstrom/sec
flambda*=1e10; // ergs/m/sec
fnu=flambda*pow(lam*1e-10,2)(-:1:36,)/c; //ergs/sec/Hz
// putting everything at 10 pc
pc10=3.0856e16*10.*1e2; //cm
fnu/=(4.*pi*pc10*pc10); //ergs/sec/Hz/cm2

f1600=interp(fnu,lam,1600.,2); //ergs/sec/Hz/cm2 at 1600 Angstrom for 1e6 Ms particle

save,createb(swrite(format="~/qtest/utils/mag/f1600_SB99_%dMyrs.yor",int(dur_emma))),f1600,age;

M1600=-2.5*log10(f1600/1e6)-48.6;
f5500=interp(fnu,lam,5500,2); // at 5500 Angstrom for V mag
f5500Vega=3640e-23;// vega flux ergs/s/cm2/Hz
MV=-2.5*log10(f5500/f5500Vega);

