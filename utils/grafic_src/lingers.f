C Developed by Edmund Bertschinger, Chung-Pei Ma and Paul Bode
C              Department of Physics
C              MIT, Room 6-207
C              Cambridge, MA 02139
C              edbert@mit.edu
C
C  For conditions of distribution and use, see the accompanying 
C  LICENSE file.
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	program linger
c  This program integrates the linearized equations of general relativity,
c  the Boltzmann equations and the fluid equations for scalar perturbations
c  of a spatially flat Friedmann-Robertson-Walker universe in synchronous
c  gauge.  The time variable is conformal time dtau=dt/a(t) and the spatial
c  dependence is Fourier transformed.  ak is the comoving wavevector;
c  comoving distances are x=r/a(t), with a(t)=1 today.  The units of both
c  length and time are Mpc.
c
c  Presently restricted to a flat (K=0) universe.
c
	implicit double precision (a-h,o-z)
c
	common /lingerinc/ omegab,omegac,omegav,omegan,h0,tcmb,yhe,
     &                     nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmax,nqmax,iq0,iq1,iq2
	common /genparm/ grhom,grhog,grhor,adotrad
	common /initcase/ initfl
c
c  Numerical integration parameters.
c  Dependent integration variables: a,ahdot,eta,delta_c,theta_c,delta_b,
c  theta_b, 2*(lmax+1) photon moments (total intensity plus polarization),
c  (lmax+1) massless neutrino moments, nqmax*(lmaxnu+1) massive neutrino
c  moments.
	parameter (lmax0=10000,lmaxnu=50,nqmax0=128,qmax=25.0d0)
	parameter (nvar0=7+3*(lmax0+1)+nqmax0*(lmaxnu+1))
	dimension y(nvar0)
	parameter (nkmax=50000)
	integer kdone(nkmax)
c
	common /starter/ taumax,dk,dlk,akmin,nk,lmax1,ifulbol,kdone
c
c  Dverk parameters.
	parameter (tol=1.0d-6,nw=nvar0)
c  To select dverk options, call dverk with ind=2.  The first 9 entries
c  of c(i) must be set to zero except for those options to be set.  A useful
c  option is c(7) = maximum number of function evaluations.
	dimension c(24),w(nw,9)
	external derivs
c
	call setup
c
c  Loop over wavenumbers.
	do 30 ik=1,nk
	if (kdone(ik).eq.1) go to 30
	if (ifulbol.eq.1) then
	  ak=ik*dk
	else
	  ak=akmin*exp((ik-1)*dlk)
	end if
c  Begin when wave is far outside horizon.
c  Conformal time (in Mpc) in the radiation era, for photons plus 3 species
c  of relativistic neutrinos.
	taustart=0.001d0/ak
c  Make sure to start early in the radiation era.
	taustart=min(taustart,0.1d0)
c  Start when massive neutrinos are strongly relativistic.
	if (amnu.ne.0.0d0) then
	  arel=1.d-3/amnu
	  taustart=min(taustart,arel/adotrad)
	end if
	ak2=ak*ak
	lmax=1.5*ak*taumax+10.0
	lmax=min(lmax,lmax1)
	iq0=11+3*lmax
	iq1=iq0+nqmax
	iq2=iq1+nqmax
	nvar=7+3*(lmax+1)+nqmax*(lmaxnu+1)
c
	call initial(y,taustart)
	tau=taustart
c
	ind=1
	nout=10
c  Begin timestep loop.
	do 10 itime=1,1000000
	taumax1=taumax*exp(-0.3*nout)
	a=y(1)
	a2=a*a
c  Timestep constraints: Courant time dtau1, expansion time dtau2.
c  Note that dverk takes many small steps per call, so it is not
c  necessary to satisfy stringent timestep constraints here.
	dtau1=10.0d0/ak
	dtau2=1.0d0*tau
	dtau=min(dtau1,dtau2)
	tauend=min(tau+dtau,taumax1)
c  Advance one timestep.
	call dverk(nvar,derivs,tau,y,tauend,tol,ind,c,nw,w)
	if (ind.lt.0) then
	  write(*,*) 'Dverk error: ind,ak=',ind,ak
	  write(*,*) 'c(1-24)=',c
	  write(*,*) 'y(1-10)=',(y(i),i=1,10)
	  go to 30
	end if
c  Evaluate metric and massive neutrino perturbations for output.
	a=y(1)
	a2=a*a
	ahdot=y(2)
	eta=y(3)
	deltac=y(4)
	thetac=y(5)
	deltab=y(6)
	thetab=y(7)
	deltag=y(8)
	thetag=y(9)
	shearg=y(10)/2.0d0
	deltar=y(10+2*lmax)
	thetar=y(11+2*lmax)
	shearr=y(12+2*lmax)/2.0d0
	call nu1(a,rhonu,pnu)
c  8*pi*G*rho*a**2.
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
	call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
	deltan=drhonu/rhonu
	thetan=ak*fnu/(rhonu+pnu)
	shearn=shearnu/(rhonu+pnu)
c  8*pi*G*rho*delta*a**2.
	dgrho=grhom*(omegac*deltac+omegab*deltab)/a
     2    +(grhog*deltag+grhor*(nnur*deltar+nnunr*drhonu))/a2
c  Add a seed if desired.
	if (initfl.eq.4) dgrho=dgrho+grhom/a
c  8*pi*G*(rho+P)*theta*a**2.
	dgtheta=grhom*(omegac*thetac+omegab*thetab)/a
     2    +4.0d0/3.0d0*(grhog*thetag+nnur*grhor*thetar)/a2
     3    +nnunr*grhor*ak*fnu/a2
c  8*pi*G*(rho+P)*sigma*a**2.
	dgshear=4.0d0/3.0d0*(grhog*shearg+nnur*grhor*shearr)/a2
     2    +nnunr*grhor*shearnu/a2
	adotoa=sqrt(grho/3.0d0)
c  Check energy constraint equation.
	econ=(adotoa*ahdot/a-2.0d0*ak2*eta-dgrho)/grho
	if (tau.eq.taumax1) then
c  Output data.
c  Compute gauge transformation variable etatophi=phi-eta to convert
c  from synchronous gauge variable eta to conformal Newtonian gauge
c  variable phi.  Output in place of thetac since the latter is 0 in
c  synchronous gauge.
	  etadot=0.5d0*dgtheta/ak2
	  etatophi=-0.5d0*adotoa/ak2*(ahdot/a+6.0d0*etadot)
c  This is inelegant, but it's the only way that works for all
c  the various preprocessors.
c  Append statement for Dec:
	  open(10,file='linger.dat',status='old',access='append')
	  write(10,'(i7,1x,19(1pe11.4,1x))') ik,ak,a,tau,ahdot,eta,
     &      deltac,deltab,deltag,deltar,deltan,etatophi,thetab,thetag,
     &      thetar,thetan,shearg,shearr,shearn,econ
	  close(10)
c
	  if (nout.eq.0) then
c  Output CMB data.
	    open(11,file='lingerg.dat',status='old',
     &                           access='append',form='unformatted')
	    write(11) ik,ak,tau,lmax
c  Correct variables to Delta_l.
	    write(11) 0.25d0*y(8),y(9)/(3.d0*ak),(0.25d0*y(8+l),
     &                   l=2,lmax)
	    write(11) (0.25d0*y(9+lmax+l),l=0,lmax)
	    close(11)
	    go to 30
	  end if
	  nout=nout-1
	end if
10	continue
30	continue
c
	stop
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine setup
c
	implicit double precision (a-h,o-z)
c
	common /lingerinc/ omegab,omegac,omegav,omegan,h0,tcmb,yhe,
     &                     nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmax,nqmax,iq0,iq1,iq2
	common /genparm/ grhom,grhog,grhor,adotrad
	common /initcase/ initfl
c
c  Numerical integration parameters.
c  Dependent integration variables: a,ahdot,eta,delta_c,theta_c,delta_b,
c  theta_b,2*(lmax+1) photon moments (total intensity plus polarization),
c  (lmax+1) massless neutrino moments, nqmax*(lmaxnu+1) massive neutrino
c  moments.
	parameter (lmax0=10000,lmaxnu=50,nqmax0=128,qmax=25.0d0)
	parameter (nkmax=50000)
	integer kdone(nkmax)
	common /starter/ taumax,dk,dlk,akmin,nk,lmax1,ifulbol,kdone
c
c  const=7*pi**4/120.
	parameter (const=5.68219698d0,zeta3=1.20205690d0)
c
	parameter (tol=1.0d-6)
c
c  Internal variables.
	dimension denl(lmax0),dlfdlq(nqmax0)
	common /store/ denl,dlfdlq
	save /store/
c
	external dtauda
c
c  Read initial parameters.
1	write(*,*)
     2 'Enter Omega_b, Omega_c, Omega_v, Omega_nu (e.g. .05 .30 0.65 0)'
	read(*,*) omegab,omegac,omegav,omegan
	omega=omegab+omegac+omegav+omegan
	omegak=1.0d0-omega
	if (abs(omegak).gt.1.0d-3) then
	  write(*,*) 'Linger currently works only for a flat universe'
	  write(*,*) '  You request Omega_curvature=',omegak
	  write(*,*) '  Illegal choice.  Try again:'
	  go to 1
	end if
2	write(*,*)
     2    'Enter H0, Tcmb, Y_He, N_nu(massive) (e.g. 65 2.726 0.24 0)'
	read(*,*) h0,tcmb,yhe,nnunr
	if (h0.lt.25.d0.or.h0.gt.100.d0) then
	  write(*,*)
     2      '  Warning: H0 has units of km/s/Mpc.  Your value is weird.'
	end if
	if (tcmb.lt.2.7d0.or.tcmb.gt.2.8d0) then
	  write(*,*)
     2      '  Warning: Tcmb has units of K.  Your value is weird.'
	end if
	if (yhe.lt.0.2d0.or.yhe.gt.0.3d0) then
	  write(*,*)
     2      '  Warning: Y_He is the Helium fraction of baryons.',
     3      '  Your value is weird.'
	end if
	if (nnunr.lt.0.or.nnunr.gt.3) then
	  write(*,*)
     2      'Warning: N_nu(massive) should be 0, 1, 2, or 3'
	  write(*,*) '  Illegal choice.  Try again:'
	  go to 2
	end if
	nnur=3-nnunr
3	continue
c	write(*,*)	
c     2    'Enter 1 for full Boltzmann equation for CMB',
c     3    ' (lmax<=10000, linear k)'
c	write(*,*)
c     2    '   or 0 for matter transfer functions only (lmax=100, log k)'
c	read(*,*) ifulbol
	write(*,*) 'Computer matter transfer functions only, lmax=100'
	ifulbol=0
	if (ifulbol.ne.0.and.ifulbol.ne.1) then
	  write(*,*) '  Illegal choice.  Try again:'
	  go to 3
	end if
	if (ifulbol.eq.1) then
	  write(*,*) 'Enter kmax (1/Mpc), nk, zend'
	  read(*,*) akmax,nk,zend
	  dk=akmax/nk
	  akmin=dk
	  lmax1=lmax0
	else
	  write(*,*) 'Enter kmin (1/Mpc), kmax (1/Mpc), nk, zend ',
     2 '(e.g. 1.e-4 50 100 0)'
	  read(*,*) akmin,akmax,nk,zend
	  if (nk.gt.1) then
	    dlk=log(akmax/akmin)/(nk-1)
	  else
	    dlk=0.0d0
	  end if
	  lmax1=100
	end if
	write(*,*) 'Initial conditions cases:'
	write(*,*) '    1 for isentropic (adiabatic) fluctuations,'
	write(*,*) '    2 for cdm entropy/isocurvature fluctuations, or'
	write(*,*) '    3 for baryon entropy/isocurvature',
     2                    ' fluctuations, or'
	write(*,*) '    4 for seed/isocurvature fluctuations'
	write(*,*) 'Enter 1, 2, 3, or 4 now'
	read(*,*) initfl
	if (initfl.lt.1.or.initfl.gt.4) then
	  write(*,*) 'Illegal choice!'
	  stop
	end if
c
	do l=1,lmax0
	  denl(l)=1.0d0/(2*l+1)
	end do
c  Evaluate general cosmology constants.
c  grho gives the contribution to the expansion rate from: (g) photons,
c  (r) one flavor of relativistic neutrino (2 degrees of freedom),
c  (m) nonrelativistic matter (for Omega=1).  grho is actually
c  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
c  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
c  (Used only to set the initial conformal time.)
	grhom=3.3379d-11*h0*h0
	grhog=1.4952d-13*tcmb**4
	grhor=3.3957d-14*tcmb**4
c  adotrad gives the relation a(tau) in the radiation era:
	adotrad=2.8948d-7*tcmb*tcmb
c
c  Initialize neutrino mass and storage.
	if (nnunr.eq.0.or.omegan.eq.0.0) then
	  amnu=0.0d0
	  nqmax=0
	else
c  amnu=m_nu*c**2/(k_B*T_nu0).
c  The 1.0e-18 is to prevent certain compilers (e.g. Sun's) from
c  thinking that a divide by zero will occur.
	  amnu=const/(1.5d0*zeta3)*grhom/grhor*omegan/(nnunr+1.0e-18)
	  nqmax=nqmax0
	  dq=qmax/nqmax
	    do i=1,nqmax
	    q=i*dq
	    expq=exp(-q)
	    dlfdlq(i)=-q/(1.0d0+expq)
	  end do
	end if
c
c  Stopping time.
	amax=1.0d0/(1.0d0+zend)
	taumax=rombint(dtauda,0.0d0,amax,tol)
c
	taumin=0.001d0/akmax
	taumin=min(taumin,0.1d0)
	if (amnu.ne.0.0d0) then
c  Initialize massive neutrinos.
	  call initnu1(amnu)
c  Prepare taumin for inithermo.
	  arel=1.d-3/amnu
	  taumin=min(taumin,arel/adotrad)
	end if
c
c  Initialize baryon temperature and ionization fractions vs. time.
	call inithermo(taumin,taumax)
c
	do i=1,nkmax
	  kdone(i)=-1
	end do
c
c  Check whether linger.dat exists.
	open(10,err=6,file='linger.dat',status='unknown')
	rewind 10
	read(10,*,end=6) omegab1,omegac1,omegav1,omegan1
	read(10,*,end=6) h01,tcmb1,yhe1,nnur1,nnunr1,initfl1
c  If it does, check that it has the correct parameters.
	e1=abs(omegab1-omegab)
	e2=abs(omegac1-omegac)
	e3=abs(omegav1-omegav)
	e4=abs(omegan1-omegan)
	e5=abs(h01-h0)
	e6=abs(tcmb1-tcmb)
	e7=abs(yhe1-yhe)
	err=e1+e2+e3+e4+e5+e6+e7
	if (err.gt.1.0d-4.or.nnur1.ne.nnur.or.nnunr1.ne.nnunr.or.
     &      initfl1.ne.initfl) then
	  write(*,*) 'linger.dat file does not match input parameters'
	  write(*,*) 'input: ',omegab,omegac,omegav,omegan,h0,
     &                tcmb,yhe,nnur,nnunr,initfl
	  write(*,*) 'output: ',omegab1,omegac1,omegav1,omegan1,h01,
     &                tcmb1,yhe1,nnur1,nnunr1,initfl1
	  close(10)
	  stop
	end if
	close(10)
c  Set flag iflin=1 when linger.dat exists and is satisfactory.
	iflin=1
	go to 7
6	close(10)
c  Set flag iflin=0 when linger.dat does not exist or is flawed.
	iflin=0
c  Check whether lingerg.dat exists.
7	open(11,err=8,file='lingerg.dat',status='unknown',
     &       form='unformatted')
	rewind 11
	read(11,end=8) omegab1,omegac1,omegav1,omegan1,h01,tcmb1,yhe1,
     &                 nnur1,nnunr1,initfl1
c  If it does, check that it has the correct parameters.
	e1=abs(omegab1-omegab)
	e2=abs(omegac1-omegac)
	e3=abs(omegav1-omegav)
	e4=abs(omegan1-omegan)
	e5=abs(h01-h0)
	e6=abs(tcmb1-tcmb)
	e7=abs(yhe1-yhe)
	err=e1+e2+e3+e4+e5+e6+e7
	if (err.gt.1.0d-4.or.nnur1.ne.nnur.or.nnunr1.ne.nnunr.or.
     &      initfl1.ne.initfl) then
	  write(*,*) 'lingerg.dat file does not match input parameters'
	  write(*,*) 'input: ',omegab,omegac,omegav,omegan,h0,
     &                tcmb,yhe,nnur,nnunr,initfl
	  write(*,*) 'output: ',omegab1,omegac1,omegav1,omegan1,h01,
     &                tcmb1,yhe1,nnur1,nnunr1,initfl1
	  close(11)
	  stop
	end if
	read(11,end=8) i
	read(11,end=8) dummin,dummax,nkdum,dumz,dumtau
c  Check that ifulbol, akmin, akmax, nk, zend, taumax are correct.
	if ((i.ne.ifulbol).or.(dummin.ne.akmin).or.(dummax.ne.akmax).or.
     2    (nkdum.ne.nk).or.(dumz.ne.zend).or.(dumtau.ne.taumax)) then
	  write(*,*) 'lingerg.dat file does not match input parameters'
	  write(*,*) 'input: ',ifulbol,akmin,akmax,nk,zend,taumax
	  write(*,*) 'output: ',i,dummin,dummax,nkdum,dumz,dumtau
	  close(11)
	  stop
	end if
c  Find which k have already been done (from lingerg.dat only).
	do i=1,nkmax
	  read(11,end=8) ik,ak,tau,lmax
	  read(11,end=8) (deltat,l=0,lmax)
	  read(11,end=8) (deltap,l=0,lmax)
	  kdone(ik)=1
	end do
	if (iflin.eq.0) then
	  write(*,*) 'Warning: lingerg.dat exists and is correct, but',
     &               ' linger.dat'
	  write(*,*) '         does not exist or has a bad header.'
	  write(*,*) '    Did you misplace linger.dat?'
	  write(*,*) 'If not, please move/remove linger.dat and try',
     &               ' again.'
	  close(11)
	  stop
	end if
8	close(11)
	if (kdone(1).eq.-1) then
c  Data files were empty or non-existent.  Write their headers.
	  open(10,file='linger.dat',status='unknown')
	  rewind 10
	  write(10,'(4(1pe12.6,1x))') omegab,omegac,omegav,omegan
	  write(10,'(3(1pe12.6,1x),3(i2,2x))') h0,tcmb,yhe,nnur,nnunr,
     &           initfl
	  close(10)
	  open(11,file='lingerg.dat',status='unknown',
     &           form='unformatted')
	  rewind 11
	  write(11) omegab,omegac,omegav,omegan,h0,tcmb,yhe,nnur,nnunr,
     &           initfl
	  write(11) ifulbol
	  write(11) akmin,akmax,nk,zend,taumax
	  close(11)
	end if
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine derivs(n,x,y,yprime)
c  Evaluate the time derivatives of the perturbations.
c
	implicit double precision (a-h,o-z)
	dimension y(n),yprime(n)
c
	common /lingerinc/ omegab,omegac,omegav,omegan,h0,tcmb,yhe,
     &                     nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmax,nqmax,iq0,iq1,iq2
	common /genparm/ grhom,grhog,grhor,adotrad
	common /initcase/ initfl
c
	parameter (lmax0=10000,lmaxnu=50,nqmax0=128,qmax=25.0d0)
c  Temperatures at which to switch off tightly-coupled photon-baryon
c  approximation.
	parameter (tsw1=2.0d4,tsw2=2.0d5)
c  Internal variables.
	dimension denl(lmax0),dlfdlq(nqmax0),akv(nqmax0),akov(nqmax0)
	common /store/ denl,dlfdlq
c
	tau=x
	a=y(1)
	ahdot=y(2)
	eta=y(3)
c  CDM.
	deltac=y(4)
	thetac=y(5)
c  Baryons.
	deltab=y(6)
	thetab=y(7)
c  Photons.
	deltag=y(8)
	thetag=y(9)
	shearg=y(10)/2.0d0
c  Polarization term.
	polter=y(10)+y(9+lmax)+y(11+lmax)
c  Massless neutrinos.
	deltar=y(10+2*lmax)
	thetar=y(11+2*lmax)
	shearr=y(12+2*lmax)/2.0d0
c
	a2=a*a
	a3=a2*a
	call thermo(tau,tempb,cs2,xe)
c  akthom is the Thomson opacity coefficient.
	akthom=2.3048d-9*(1-yhe)*omegab*h0*h0
c  Thomson opacity a*n_e*sigma_T (1/Mpc).
	opac=xe*akthom/a2
c  Photon mass density over baryon mass density.
	photbar=grhog/(grhom*omegab*a)
	pb43=4.0d0/3.0d0*photbar
c  Compute expansion rate.
	call nu1(a,rhonu,pnu)
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
	adotoa=sqrt(grho/3.0d0)
	yprime(1)=adotoa*a
	gpres=((grhog+grhor*nnur)/3.0d0+grhor*nnunr*pnu)/a2
     2       -grhom*omegav*a2
c  Evaluate metric and massive neutrino perturbations.
	call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
	deltan=drhonu/rhonu
	thetan=ak*fnu/(rhonu+pnu)
c  8*pi*G*delta_rho*a**2 and 8*pi*G*delta_P*a**2.
	dgrho=grhom*(omegac*deltac+omegab*deltab)/a
     2    +(grhog*deltag+grhor*(nnur*deltar+nnunr*drhonu))/a2
	dgpres=(grhog*deltag+grhor*nnur*deltar)/a2/3.0d0
     2    +grhor*nnunr*dpnu/a2
c  Add a seed if desired.
	if (initfl.eq.4) dgrho=dgrho+grhom/a
	dahdotdtau=-(dgrho+3.0d0*dgpres)*a
	yprime(2)=dahdotdtau
c  Force energy conservation.
	hdot=(2.0d0*ak2*eta+dgrho)/adotoa
c  8*pi*G*(rho+P)*theta*a**2.
	dgtheta=grhom*(omegac*thetac+omegab*thetab)/a
     2    +4.0d0/3.0d0*(grhog*thetag+nnur*grhor*thetar)/a2
     3    +nnunr*grhor*ak*fnu/a2
	etadot=0.5d0*dgtheta/ak2
	yprime(3)=etadot
c  8*pi*G*(rho+P)*sigma*a**2.
	dgshear=4.0d0/3.0d0*(grhog*shearg+nnur*grhor*shearr)/a2
     2    +nnunr*grhor*shearnu/a2
c  CDM equations of motion.
	deltacdot=-thetac-0.5d0*hdot
	yprime(4)=deltacdot
	thetacdot=-adotoa*thetac
	yprime(5)=thetacdot
c  Baryon equations of motion.
	deltabdot=-thetab-0.5d0*hdot
	yprime(6)=deltabdot
c  Need photon perturbation for first-order correction to tightly-coupled
c  baryon-photon approximation.
	deltagdot=4.0d0/3.0d0*(-thetag-0.5d0*hdot)
	drag=opac*(thetag-thetab)
	if (tempb.lt.tsw1) then
c  Treat baryons and photons as uncoupled.
	  thetabdot=-adotoa*thetab+ak2*cs2*deltab+pb43*drag
	else
c  Treat baryons and photons as tightly coupled.
c  Zeroth-order approximation to baryon velocity.
	  thetabdot=(-adotoa*thetab+ak2*cs2*deltab
     2  +ak2*pb43*(0.25d0*deltag-shearg))/(1.0d0+pb43)
c  (\ddot a)/a.
	  adotdota=0.5d0*(adotoa*adotoa-gpres)
c  First-order approximation to baryon-photon slip, thetabdot-thetagdot.
	  slip=2.0d0*pb43/(1.0d0+pb43)*adotoa*(thetab-thetag)
     2     +1.0d0/opac*(-adotdota*thetab-adotoa*ak2*0.5d0*deltag
     3     +ak2*(cs2*deltabdot-0.25d0*deltagdot))/(1.0d0+pb43)
c  First-order approximation to baryon velocity.
	  thetabdot=thetabdot+pb43/(1.0d0+pb43)*slip
	end if
	yprime(7)=thetabdot
c  Photon total intensity and polarization equations of motion.
	yprime(8)=deltagdot
	thetagdot=(-thetabdot-adotoa*thetab+ak2*cs2*deltab)/pb43
     2   +ak2*(0.25d0*deltag-shearg)
	yprime(9)=thetagdot
	if (tempb.lt.tsw2) then
c  Treat baryons and photons as uncoupled.
	  yprime(10)=8.0d0/15.0d0*thetag-0.6d0*ak*y(11)-opac*y(10)
     2        +4.0d0/15.0d0*hdot+8.0d0/5.0d0*etadot+0.1d0*opac*polter
c  Polarization equations for l = 0, 1, 2.
	  yprime(9+lmax)=-ak*y(10+lmax)-opac*y(9+lmax)+0.5d0*opac*polter
	  yprime(10+lmax)=ak/3.0d0*(y(9+lmax)-2.0d0*y(11+lmax))
     2             -opac*y(10+lmax)
	  yprime(11+lmax)=ak*(0.4d0*y(10+lmax)-0.6d0*y(12+lmax))
     2             -opac*y(11+lmax)+0.1d0*opac*polter
	    do 10 l=3,lmax-1
	    yprime(8+l)=ak*denl(l)*(l*y(7+l)-(l+1)*y(9+l))-opac*y(8+l)
	    yprime(9+lmax+l)=ak*denl(l)*(l*y(8+lmax+l)-(l+1)*
     2                 y(10+lmax+l))-opac*y(9+lmax+l)
10	  continue
	else
c  Treat baryons and photons as tightly coupled (with negligible polarization).
	  yprime(10)=0.0d0
	  yprime(9+lmax)=0.0d0
	  yprime(10+lmax)=0.0d0
	  yprime(11+lmax)=0.0d0
	    do 15 l=3,lmax-1
	    yprime(8+l)=0.0d0
	    yprime(9+lmax+l)=0.0d0
15	  continue
	end if
c  Truncate moment expansion.
c	yprime(8+lmax)=ak*lmax*y(7+lmax)/(2*lmax+1)-opac*y(8+lmax)
c	yprime(9+2*lmax)=ak*lmax*y(8+2*lmax)/(2*lmax+1)-opac*y(9+2*lmax)
	yprime(8+lmax)=ak*y(7+lmax)-(lmax+1)/tau*y(8+lmax)
     2                -opac*y(8+lmax)
	yprime(9+2*lmax)=ak*y(8+2*lmax)-(lmax+1)/tau*y(9+2*lmax)
     2                -opac*y(9+2*lmax)
c  Massless neutrino equations of motion.
	deltardot=4.0d0/3.0d0*(-thetar-0.5d0*hdot)
	yprime(10+2*lmax)=deltardot
	thetardot=ak2*(0.25d0*deltar-shearr)
	yprime(11+2*lmax)=thetardot
	yprime(12+2*lmax)=8.0d0/15.0d0*thetar-0.6d0*ak*y(13+2*lmax)
     2             +4.0d0/15.0d0*hdot+8.0d0/5.0d0*etadot
	  do 20 l=3,lmax-1
	  yprime(10+2*lmax+l)=ak*denl(l)*(l*y(9+2*lmax+l)
     2                    -(l+1)*y(11+2*lmax+l))
20	continue
c  Truncate moment expansion.
c	yprime(10+3*lmax)=ak*lmax*y(9+3*lmax)/(2*lmax+1)
	yprime(10+3*lmax)=ak*y(9+3*lmax)-(lmax+1)/tau*y(10+3*lmax)
c  Massive neutrino equations of motion.
	if (nqmax.eq.0) return
	dq=qmax/nqmax
	  do i=1,nqmax
	  q=i*dq
	  aq=a*amnu/q
	  v=1.0d0/sqrt(1.0d0+aq*aq)
	  akv(i)=ak*v
	  akov(i)=ak/v
	end do
c  l = 0, 1, 2,lmaxnu.
	  do 30 i=1,nqmax
	  ind=iq0+i-1
	  yprime(ind)=-akv(i)*y(ind+nqmax)+hdot*dlfdlq(i)/6.0d0
	  ind=iq1+i-1
	  yprime(ind)=akv(i)*(y(ind-nqmax)-2*y(ind+nqmax))/3
	  ind=iq2+i-1
	  yprime(ind)=akv(i)*(2*y(ind-nqmax)-3*y(ind+nqmax))/5
     2           -(hdot/15.0d0+2.0d0/5.0d0*etadot)*dlfdlq(i)
	  ind=10+3*lmax+i+lmaxnu*nqmax
c  Truncate moment expansion.
c	  yprime(ind)=akv*lmaxnu*y(ind-nqmax)/(2*lmaxnu+1)
	  yprime(ind)=akv(i)*y(ind-nqmax)-(lmaxnu+1)/tau*y(ind)
30	continue
	  do 50 l=3,lmaxnu-1
	    do 40 i=1,nqmax
	    ind=10+3*lmax+i+l*nqmax
	    yprime(ind)=akv(i)*denl(l)*(l*y(ind-nqmax)-(l+1)*y(ind+nqmax))
40	  continue
50	continue
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine thermo(tau,tempb,cs2b,xeb)
c  Compute unperturbed baryon temperature, sound speed squared,
c  and ionization fraction by interpolating pre-computed tables.
c
	implicit double precision (a-h,o-z)
	parameter (nthermo=10000)
	dimension tb(nthermo),cs2(nthermo),xe(nthermo)
	dimension dtb(nthermo),dcs2(nthermo),dxe(nthermo)
	common /thermod/ taumin,dlntau,tb,cs2,xe,dtb,dcs2,dxe
c
	d=log(tau/taumin)/dlntau+1.0d0
	i=int(d)
	d=d-i
	if (i.lt.1) then
c  Linear interpolation if out of bounds (should not occur).
	  tempb=tb(1)+(d+i-1)*dtb(1)
	  cs2b=cs2(1)+(d+i-1)*dcs2(1)
	  xeb=xe(1)+(d+i-1)*dxe(1)
	else if (i.ge.nthermo) then
	  tempb=tb(nthermo)+(d+i-nthermo)*dtb(nthermo)
	  cs2b=cs2(nthermo)+(d+i-nthermo)*dcs2(nthermo)
	  xeb=xe(nthermo)+(d+i-nthermo)*dxe(nthermo)
	else
c  Cubic spline interpolation.
	  tempb=tb(i)+d*(dtb(i)+d*(3.0d0*(tb(i+1)-tb(i))-2.0d0*dtb(i)
     2         -dtb(i+1)+d*(dtb(i)+dtb(i+1)+2.0d0*(tb(i)-tb(i+1)))))
	  cs2b=cs2(i)+d*(dcs2(i)+d*(3.0d0*(cs2(i+1)-cs2(i))
     2         -2.0d0*dcs2(i)-dcs2(i+1)+d*(dcs2(i)+dcs2(i+1)
     2         +2.0d0*(cs2(i)-cs2(i+1)))))
	  xeb=xe(i)+d*(dxe(i)+d*(3.0d0*(xe(i+1)-xe(i))-2.0d0*dxe(i)
     2         -dxe(i+1)+d*(dxe(i)+dxe(i+1)+2.0d0*(xe(i)-xe(i+1)))))
	end if
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine inithermo(taumin,taumax)
c  Compute and save unperturbed baryon temperature and ionization fraction
c  as a function of time.  With nthermo=10000, xe(tau) has a relative accuracy
c  (numerical integration precision) better than 1.e-5.
c
	implicit double precision (a-h,o-z)
c
	common /lingerinc/ omegab,omegac,omegav,omegan,h0,tcmb,yhe,
     &                     nnur,nnunr
	common /genparm/ grhom,grhog,grhor,adotrad
c
	parameter (barssc0=9.1820d-14)
c
	parameter (nthermo=10000)
	dimension tb(nthermo),cs2(nthermo),xe(nthermo)
	dimension dtb(nthermo),dcs2(nthermo),dxe(nthermo)
	common /thermod/ tauminn,dlntau,tb,cs2,xe,dtb,dcs2,dxe
	save /thermod/
c
	thomc0=5.0577d-8*tcmb**4
	tauminn=taumin
	dlntau=log(taumax/taumin)/(nthermo-1)
c
c  Initial conditions: assume radiation-dominated universe.
	tau0=taumin
	adot0=adotrad
	a0=adotrad*taumin
c  Assume that any entropy generation occurs before taumin.
c  This gives wrong temperature before pair annihilation, but
c  the error is harmless.
	tb(1)=tcmb/a0
	xe0=1.0d0
	x1=0.0d0
	x2=1.0d0
	xe(1)=xe0+0.25d0*yhe/(1.0d0-yhe)*(x1+2*x2)
	barssc=barssc0*(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(1))
	cs2(1)=4.0d0/3.0d0*barssc*tb(1)
c
	  do 10 i=2,nthermo
	  tau=taumin*exp((i-1)*dlntau)
	  dtau=tau-tau0
c  Integrate Friedmann equation using inverse trapezoidal rule.
	  a=a0+adot0*dtau
	  a2=a*a
	  call nu1(a,rhonu,pnu)
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
	  adot=sqrt(grho/3.0d0)*a
	  a=a0+2.0d0*dtau/(1.0d0/adot0+1.0d0/adot)
c  Baryon temperature evolution: adiabatic except for Thomson cooling.
c  Use  quadrature solution.
	  tg0=tcmb/a0
	  ahalf=0.5d0*(a0+a)
	  adothalf=0.5d0*(adot0+adot)
c  fe=number of free electrons divided by total number of free baryon
c  particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
c  more accuracy is required (unlikely) then this can be iterated with
c  the solution of the ionization equation.
	  fe=(1.d0-yhe)*xe(i-1)/(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(i-1))
	  thomc=thomc0*fe/adothalf/ahalf**3
	  etc=exp(-thomc*(a-a0))
	  a2t=a0*a0*(tb(i-1)-tg0)*etc-tcmb/thomc*(1.d0-etc)
	  tb(i)=tcmb/a+a2t/(a*a)
c  Integrate ionization equation.
	  tbhalf=0.5d0*(tb(i-1)+tb(i))
	  call ionize(tbhalf,ahalf,adothalf,dtau,xe0)
	  call ionhe(tb(i),a,xe0,x1,x2)
	  xe(i)=xe0+0.25d0*yhe/(1.0d0-yhe)*(x1+2*x2)
c  Baryon sound speed squared (over c**2).
	  adotoa=adot/a
	  dtbdla=-2.0d0*tb(i)-thomc*adothalf/adot*(a*tb(i)-tcmb)
	  barssc=barssc0*(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(i))
	  cs2(i)=barssc*tb(i)*(1-dtbdla/tb(i)/3.0d0)
c
	  a0=a
	  tau0=tau
	  adot0=adot
10	continue
c
	call splini
	call splder(xe,dxe,nthermo)
	call splder(cs2,dcs2,nthermo)
	call splder(tb,dtb,nthermo)
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ionize(tempb,a,adot,dtau,xe)
c  Integrate the ionization fraction xe for hydrogen semi-implicitly
c  from tau to tau+dtau, treating tempb, a, and adot as constants
c  during this interval.
c
	implicit double precision (a-h,o-z)
c
	common /lingerinc/ omegab,omegac,omegav,omegan,h0,tcmb,yhe,
     &                     nnur,nnunr
	common /genparm/ grhom,grhog,grhor,adotrad
c  Ionization temperature and coefficient.
	parameter (tion=1.5789d5,beta0=43.082)
c  Two-photon decay rate (in 1/Mpc).
	parameter (dec2g=8.468d14)
c  Switch for fully implicit (switch=1.0) or semi-implicit (switch=0.5) scheme.
	parameter (switch=0.5d0)
c
c  Recombination coefficient (in sqrt(K)/Mpc).
	alpha0=2.3866d-6*(1-yhe)*omegab*h0*h0
c  Coefficient for correction of radiative decay (dimensionless).
	crec=8.0138d-26*(1-yhe)*omegab*h0*h0
c  Recombination and ionization rates.
	phi2=0.448d0*log(tion/tempb)
	phi2=max(phi2,0.0d0)
	alpha=alpha0/sqrt(tempb)*phi2/(a*a*a)
	beta=tempb*phi2*exp(beta0-tion/tempb)
c  Peebles' correction factor.
	if (tempb.lt.200.0d0) then
	  cpeebles=1.0d0
	else
	  cp1=crec*dec2g*(1.0d0-xe)/(a*adot)
	  cp2=crec*tempb*phi2*exp(beta0-0.25d0*tion/tempb)*
     2        (1.0d0-xe)/(a*adot)
	  cpeebles=(1.0d0+cp1)/(1.0d0+cp1+cp2)
	end if
c  Integrate dxe=bb*(1-xe)-aa*xe*xe by averaging rhs at current tau
c  (fraction 1-switch) and future tau (fraction switch).
	aa=a*dtau*alpha*cpeebles
	bb=a*dtau*beta*cpeebles
	b1=1.0d0+switch*bb
	bbxe=bb+xe-(1.0d0-switch)*(bb*xe+aa*xe*xe)
	rat=switch*aa*bbxe/(b1*b1)
c  Prevent roundoff error.
	if (rat.lt.1.0d-6) then
	  xe=bbxe/b1*(1.0d0-rat)
	else
	  xe=b1/(2.0d0*switch*aa)*(sqrt(4.0d0*rat+1.0d0)-1.0d0)
	end if
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ionhe(tempb,a,x0,x1,x2)
c  Compute the helium ionization fractions using the Saha equation.
c  x0 is the hydrogen ionization fraction n(H+)/n(H) (input),
c  x1 is the helium first ionization fraction n(He+)/n(He)
c    (input and output), and
c  x2 is the helium second ionization fraction n(He++)/n(He)
c    (input and output).
c
	implicit double precision (a-h,o-z)
c
	common /lingerinc/ omegab,omegac,omegav,omegan,h0,tcmb,yhe,
     &                     nnur,nnunr
	common /genparm/ grhom,grhog,grhor,adotrad
c  Ionization temperatures.
	parameter (tion1=2.855e5,tion2=6.313e5)
c
c  Constant for electron partition function per baryon.
	b0=2.150e24/((1.0d0-yhe)*omegab*h0*h0)
c  Electron partition function per baryon.
	b=b0*a*a*a*tempb*sqrt(tempb)
c  Dimensionless right-hand sides in Saha equations.
	r1=4.0d0*b*exp(-tion1/tempb)
	r2=b*exp(-tion2/tempb)
c
c  Solve coupled equations iteratively.
	c=0.25d0*yhe/(1.0d0-yhe)
	err=1.0d0
	niter=0
10	niter=niter+1
	if (err.lt.1.0d-12) return
	  xe=x0+c*(x1+2*x2)
	  x2new=r1*r2/(r1*r2+xe*r1+xe*xe)
	  x1=xe*r1/(r1*r2+xe*r1+xe*xe)
	  err=abs(x2new-x2)
	  x2=x2new
	  go to 10
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine nu1(a,rhonu,pnu)
c  Compute massive neutrino density and pressure in units of the mean
c  density of one flavor of massless neutrinos.  Use cubic splines to
c  interpolate from a table.
c
	implicit double precision (a-h,o-z)
c
	common /cosmoparm/ ak,ak2,amnu,lmax,nqmax,iq0,iq1,iq2
	parameter (nrhopn=10000)
	dimension r1(nrhopn),p1(nrhopn)
	dimension dr1(nrhopn),dp1(nrhopn)
	common /nu1d/ amin,dlna,r1,p1,dr1,dp1
c
	if (amnu.eq.0.0d0) then
	  rhonu=1.0d0
	  pnu=1.0d0/3.0d0
	  return
	end if
c
	d=log(a/amin)/dlna+1.0d0
	i=int(d)
	d=d-i
	if (i.lt.1) then
c  Use linear interpolation, bounded by results for massless neutrinos.
	  rhonu=r1(1)+(d-1)*dr1(1)
	  pnu=p1(1)+(d-1)*dp1(1)
	  rhonu=min(exp(rhonu),1.0d0)
	  pnu=min(exp(pnu),0.3333333333d0)
	else if (i.gt.nrhopn) then
c  This should not happen, unless the user evolves to z<0!
	  rhonu=r1(nrhopn)+(d+i-nrhopn)*dr1(nrhopn)
	  pnu=p1(nrhopn)+(d+i-nrhopn)*dp1(nrhopn)
	  rhonu=exp(rhonu)
	  pnu=exp(pnu)
	else
c  Cubic spline interpolation.
	  rhonu=r1(i)+d*(dr1(i)+d*(3.0d0*(r1(i+1)-r1(i))-2.0d0*dr1(i)
     2         -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2.0d0*(r1(i)-r1(i+1)))))
	  pnu=p1(i)+d*(dp1(i)+d*(3.0d0*(p1(i+1)-p1(i))-2.0d0*dp1(i)
     2         -dp1(i+1)+d*(dp1(i)+dp1(i+1)+2.0d0*(p1(i)-p1(i+1)))))
	  rhonu=exp(rhonu)
	  pnu=exp(pnu)
	end if

c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine initnu1(amnu)
c  Initialize interpolation tables for massive neutrinos.
c  Use cubic splines interpolation of log rhonu and pnu vs. log a.
	implicit double precision (a-h,o-z)
c
	parameter (nrhopn=10000)
	dimension r1(nrhopn),p1(nrhopn)
	dimension dr1(nrhopn),dp1(nrhopn)
	common /nu1d/ amin,dlna,r1,p1,dr1,dp1
	save /nu1d/
c
	amin=1.0d-9
	dlna=-log(amin)/(nrhopn-1)
c
c  Check whether correct interpolation table already exists on disk.
	open(12,err=10,file='nu1.dat',status='old',form='unformatted')
	read(12,end=10) amnu1
	err=abs(amnu1/amnu-1.0d0)
	if (err.gt.1.0e-8) then
	  go to 10
	else
	  read(12,end=10) (r1(j),p1(j),j=1,nrhopn)
	  close(12)
c  Yes, interpolation table was okay.
	  go to 20
	end if

c  No, correct interpolation tables aren't on disk: compute and save.
10	close(12)
	do i=1,nrhopn
	  a=amin*exp((i-1)*dlna)
	  call ninu1(a,rhonu,pnu)
	  r1(i)=log(rhonu)
	  p1(i)=log(pnu)
	end do

20	call splini
	call splder(r1,dr1,nrhopn)
	call splder(p1,dp1,nrhopn)

	open(12,file='nu1.dat',status='unknown',form='unformatted')
	rewind 12
	write(12) amnu
	write(12) (r1(j),p1(j),j=1,nrhopn)
	close(12)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ninu1(a,rhonu,pnu)
c  Compute the mean density and pressure of one flavor of massive neutrinos,
c  in units of the mean density of one flavor of massless neutrinos.
c  Relative accuracy is better than 1.e-10 for all a.
c
	implicit double precision (a-h,o-z)
c
	common /cosmoparm/ ak,ak2,amnu,lmax,nqmax,iq0,iq1,iq2
c
	parameter (qmax=30.0d0,nq=1000,nq1=nq+1)
	dimension dum1(nq1),dum2(nq1)
c  const=7*pi**4/120.
	parameter (const=5.68219698d0)
c
	if (amnu.eq.0.0d0) then
	  rhonu=1.0d0
	  pnu=1.0d0/3.0d0
	  return
	end if
c
c  q is the comoving momentum in units of k_B*T_nu0/c.
c  Integrate up to qmax and then use asymptotic expansion for remainder.
	dq=qmax/nq
	dum1(1)=0.0d0
	dum2(1)=0.0d0
	  do 10 i=1,nq
	  q=i*dq
	  aq=a*amnu/q
	  v=1.0d0/sqrt(1.0d0+aq*aq)
	  qdn=dq*q*q*q/(exp(q)+1.0d0)
	  dum1(i+1)=qdn/v
	  dum2(i+1)=qdn*v
10	continue
	call splint(dum1,rhonu,nq1)
	call splint(dum2,pnu,nq1)
c  Apply asymptotic corrrection for q>qmax and normalize by relativistic
c  energy density.
	rhonu=(rhonu+dum1(nq1)/dq)/const
	pnu=(pnu+dum2(nq1)/dq)/const/3.0d0
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine nu2(a,drhonu,fnu,dpnu,shearnu,psi0,psi1,psi2)
c  Compute the perturbations of density, energy flux, pressure, and
c  shear stress of one flavor of massive neutrinos, in units of the mean
c  density of one flavor of massless neutrinos, by integrating over momentum.
c
	implicit double precision (a-h,o-z)
c
	common /cosmoparm/ ak,ak2,amnu,lmax,nqmax,iq0,iq1,iq2
c
	parameter (lmax0=10000,lmaxnu=50,nqmax0=128,qmax=25.0d0)
	dimension psi0(nqmax0),psi1(nqmax0),psi2(nqmax0)
c  const=7*pi**4/120.
	parameter (const=5.68219698d0)
c
c  MAXJ=3 stops Richardson extrapolation at three levels, so that the
c  integration scheme is 8th-order Newton-Cotes.
	parameter (MAXITER=7,MAXJ=3)
	dimension g0(4),gmax(4),g(4,MAXJ+1)
c
	if (nqmax.eq.0) then
	  drhonu=0.0
	  fnu=0.0
	  dpnu=0.0
	  shearnu=0.0
	  return
	end if
c
c  Make sure that nqmax=2**MAXITER.
	if (nqmax.ne.2**MAXITER) then
	  write(*,*) 'Must set MAXITER=log_2(nqmax) in nu2'
	  stop
	end if
c  q is the comoving momentum in units of k_B*T_nu0/c.
c  Integrate up to qmax and then use asymptotic expansion for remainder.
	h=0.5d0*qmax
	q=qmax
	aq=a*amnu/q
	v=1.0d0/sqrt(1.0d0+aq*aq)
	qdn=h*q*q*q/(exp(q)+1.0d0)
	gmax(1)=qdn*psi0(nqmax)/v
	gmax(2)=qdn*psi0(nqmax)*v
	gmax(3)=qdn*psi1(nqmax)
	gmax(4)=qdn*psi2(nqmax)*v
	g(1,1)=gmax(1)
	g(2,1)=gmax(2)
	g(3,1)=gmax(3)
	g(4,1)=gmax(4)
	gf1=gmax(1)
	gf2=gmax(2)
	gf3=gmax(3)
	gf4=gmax(4)
	nint=1
	nn=nqmax/2
c
	do 40 iter=1,MAXITER
c  Calculate next trapezoidal rule approximation to integral.
	  g0(1)=0.0d0
	  g0(2)=0.0d0
	  g0(3)=0.0d0
	  g0(4)=0.0d0
	    do 20 k=1,nint
	    q=(k+k-1)*h
	    iq=(k+k-1)*nn
	    aq=a*amnu/q
	    v=1.0d0/sqrt(1.0d0+aq*aq)
	    qdn=h*q*q*q/(exp(q)+1.0d0)
	    g0(1)=g0(1)+qdn*psi0(iq)/v
	    g0(2)=g0(2)+qdn*psi0(iq)*v
	    g0(3)=g0(3)+qdn*psi1(iq)
	    g0(4)=g0(4)+qdn*psi2(iq)*v
20	  continue
	  g0(1)=0.5d0*g(1,1)+g0(1)
	  g0(2)=0.5d0*g(2,1)+g0(2)
	  g0(3)=0.5d0*g(3,1)+g0(3)
	  g0(4)=0.5d0*g(4,1)+g0(4)
	  h=0.5d0*h
	  nint=nint+nint
	  nn=nn/2
	  jmax=min(iter,MAXJ)
	  fourj=1.0d0
	    do 30 j=1,jmax
c  Use Richardson extrapolation.
	    fourj=4.0d0*fourj
	    g1=g0(1)+(g0(1)-g(1,j))/(fourj-1.0d0)
	    g2=g0(2)+(g0(2)-g(2,j))/(fourj-1.0d0)
	    g3=g0(3)+(g0(3)-g(3,j))/(fourj-1.0d0)
	    g4=g0(4)+(g0(4)-g(4,j))/(fourj-1.0d0)
	    g(1,j)=g0(1)
	    g(2,j)=g0(2)
	    g(3,j)=g0(3)
	    g(4,j)=g0(4)
	    g0(1)=g1
	    g0(2)=g2
	    g0(3)=g3
	    g0(4)=g4
30	  continue
	  err1=1.0d0-gmax(1)/g0(1)
	  err2=1.0d0-gmax(2)/g0(2)
	  err3=1.0d0-gmax(3)/g0(3)
	  err4=1.0d0-gmax(4)/g0(4)
	  gmax(1)=g0(1)
	  gmax(2)=g0(2)
	  gmax(3)=g0(3)
	  gmax(4)=g0(4)
	  g(1,jmax+1)=g0(1)
	  g(2,jmax+1)=g0(2)
	  g(3,jmax+1)=g0(3)
	  g(4,jmax+1)=g0(4)
40	continue
	drhonu=(g0(1)+gf1*2.0d0/qmax)/const
	dpnu=(g0(2)+gf2*2.0d0/qmax)/const/3.0d0
	fnu=(g0(3)+gf3*2.0d0/qmax)/const
	shearnu=(g0(4)+gf4*2.0d0/qmax)/const*2.0d0/3.0d0
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine initial(y,tau)
c  Initial conditions.
	implicit double precision (a-h,o-z)
c
	common /lingerinc/ omegab,omegac,omegav,omegan,h0,tcmb,yhe,
     &                     nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmax,nqmax,iq0,iq1,iq2
	common /genparm/ grhom,grhog,grhor,adotrad
	common /initcase/ initfl
c
	parameter (lmax0=10000,lmaxnu=50,nqmax0=128,qmax=25.0d0)
	parameter (nvar0=7+3*(lmax0+1)+nqmax0*(lmaxnu+1))
	dimension y(nvar0)
c
	a=tau*adotrad
	a2=a*a
	call nu1(a,rhonu,pnu)
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
	adotoa=sqrt(grho/3.0d0)
	gpres=((grhog+grhor*nnur)/3.0d0+grhor*nnunr*pnu)/a2
     2       -grhom*omegav*a2
	s=grho+gpres
	fracnu=grhor*4.0d0/3.0d0*(nnur+nnunr)/a2/s
c  Use yrad=rho_matter/rho_rad to correct initial conditions for
c  matter+radiation.
	yrad=grhom*(omegac+omegab)*a/(grhog+grhor*(nnur+nnunr*rhonu))
c
c  Choose one of the following four cases for initial conditions, or
c  add your own.  Comment out the other cases.
c
	if (initfl.eq.1) then
c-------------------------------------------------------------------------------
c  First case.
c  Isentropic ("adiabatic") initial conditions.
	psi=-1.0d0
	C=(15.0d0+4.0d0*fracnu)/20.0d0*psi
	akt2=(ak*tau)**2
	h=C*akt2*(1.0d0-0.2d0*yrad)
	eta=2.0d0*C-(5.0d0+4.0d0*fracnu)/6.0d0/(15.0d0+4.0d0*fracnu)*
     2              C*akt2*(1.0d0-yrad/3.0d0)
	f1=(23.0d0+4.0d0*fracnu)/(15.0d0+4.0d0*fracnu)
	deltac=-0.5d0*h
	deltag=-2.0d0/3.0d0*h*(1.0d0-akt2/36.0d0)
	deltab=0.75d0*deltag
	deltar=-2.0d0/3.0d0*h*(1.0d0-akt2/36.0d0*f1)
	thetac=0.0d0
	thetag=-C/18.0d0*akt2*akt2/tau
	thetab=thetag
	thetar=f1*thetag
	shearr=4.0d0/15.0d0*ak2/s*psi*(1.0d0+7.0d0/36.0d0*yrad)
	ahdot=2.0d0*C*ak2*tau*a*(1.0d0-0.3d0*yrad)
c-------------------------------------------------------------------------------
	else if (initfl.eq.2) then
c  Second case.
c  Isocurvature CDM initial conditions: perturb only CDM as a --> 0.
	delta0=1.0d0
	h=delta0*yrad*(1.0d0/(1.0d0+omegab/omegac)-0.5d0*yrad)
	deltac=delta0-0.5d0*h
c  Compensate perturbation with everything else.
	deltag=-2.0d0/3.0d0*h
	deltab=-0.75d0*deltag
	deltar=deltag
	thetac=0.0d0
	thetag=-h/12.0d0*ak2*tau
	thetab=thetag
	thetar=thetag
	shearr=0.0d0
	ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
	eta=-h/6.0d0
c-------------------------------------------------------------------------------
	else if (initfl.eq.3) then
c  Third case.
c  Isocurvature baryon initial conditions: perturb only baryons as a --> 0.
	delta0=1.0d0
	h=delta0*yrad*(1.0d0/(1.0d0+omegac/omegab)-0.5d0*yrad)
	deltab=delta0-0.5d0*h
c  Compensate perturbation with everything else.
	deltac=-0.5d0*h
	deltag=-2.0d0/3.0d0*h
	deltar=deltag
	thetac=0.0d0
	thetag=-h/12.0d0*ak2*tau
	thetab=thetag
	thetar=thetag
	shearr=0.0d0
	ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
	eta=-h/6.0d0
c-------------------------------------------------------------------------------
	else if (initfl.eq.4) then
c  Fourth case.
c  Isocurvature seed initial conditions: everything is unperturned as a --> 0.
	delta0=1.0d0
	h=delta0*yrad*(1.0d0/(1.0d0+omegac/omegab)-0.5d0*yrad)
c  Compensate perturbation with everything else.
	deltab=-0.5d0*h
	deltac=-0.5d0*h
	deltag=-2.0d0/3.0d0*h
	deltar=deltag
	thetac=0.0d0
	thetag=-h/12.0d0*ak2*tau
	thetab=thetag
	thetar=thetag
	shearr=0.0d0
	ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
	eta=-h/6.0d0
c-------------------------------------------------------------------------------
	else
	  write(*,*) 'initfl must equal 1-4! initfl=',initfl
	  stop
	end if
c
	deltan=deltar
	thetan=thetar
c
	y(1)=a
	y(2)=ahdot
	y(3)=eta
c  CDM.
	y(4)=deltac
	y(5)=thetac
c  Baryons.
	y(6)=deltab
	y(7)=thetab
c  Photons (total intensity and polarization).
	y(8)=deltag
	y(9)=thetag
	shearg=0.0d0
	y(9+lmax)=0.0d0
	y(10+lmax)=0.0d0
	  do 10 l=2,lmax
	  y(8+l)=0.0d0
	  y(9+lmax+l)=0.0d0
10	continue
c  Massless neutrinos.
	y(10+2*lmax)=deltar
	y(11+2*lmax)=thetar
	y(12+2*lmax)=shearr*2.0d0
	  do 20 l=3,lmax
	  y(10+2*lmax+l)=0.0d0
20	continue
c  Massive neutrinos.
	if (nqmax.eq.0) go to 50
	dq=qmax/nqmax
	  do 40 i=1,nqmax
	  q=i*dq
	  aq=a*amnu/q
	  v=1.0d0/sqrt(1.0d0+aq*aq)
	  akv=ak*v
	  expq=exp(-q)
	  dlfdlq=-q/(1.0d0+expq)
	  y(iq0+i-1)=-0.25d0*dlfdlq*deltan
c  Divide by v to get first-order correction for neutrino mass.
	  y(iq1+i-1)=-dlfdlq*thetan/akv/3.0d0
	  y(iq2+i-1)=-0.5d0*dlfdlq*shearr
	    do 30 l=3,lmaxnu
	    ind=10+3*lmax+i+l*nqmax
	    y(ind)=0.0d0
30	  continue
40	continue
c  Check energy constraint equation.
50	call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
	deltan=drhonu/rhonu
	thetan=ak*fnu/(rhonu+pnu)
	shearn=shearnu/(rhonu+pnu)
	dgrho=grhom*(omegac*deltac+omegab*deltab)/a
     2    +(grhog*deltag+grhor*(nnur*deltar+nnunr*drhonu))/a2
c  Add a seed if desired.
	if (initfl.eq.4) dgrho=dgrho+grhom/a
	econ=(adotoa*ahdot/a-2.0d0*ak2*eta-dgrho)/grho
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function dtauda(a)
	implicit double precision (a-h,o-z)
c
	common /lingerinc/ omegab,omegac,omegav,omegan,h0,tcmb,yhe,
     &                     nnur,nnunr
	common /genparm/ grhom,grhog,grhor,adotrad
c
c  Compute expansion rate.
	call ninu1(a,rhonu,pnu)
c  8*pi*G*rho*a**4.
	grho2=grhom*(omegac+omegab)*a+(grhog+grhor*(nnur+nnunr*rhonu))
     2       +grhom*omegav*a**4
	dtauda=sqrt(3.0d0/grho2)
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine splder(y,dy,n)
c  Splder fits a cubic spline to y and returns the first derivatives at
c  the grid points in dy.  Dy is equivalent to a 4th-order Pade
c  difference formula for dy/di.
c
	implicit double precision (a-h,o-z)
	dimension f(20001),g(20001),y(n),dy(n)
	common /spline/ g
c
	n1=n-1
	if (n1.gt.20000) write(*,*) 'Spline array overflow!!! n1=',
     2      n1,'>20000'
c  Quartic fit to dy/di at boundaries, assuming d3y/di3=0.
	f(1)=(-10.0d0*y(1)+15.0d0*y(2)-6.0d0*y(3)+y(4))/6.0d0
	f(n)=(10.0d0*y(n)-15.0d0*y(n1)+6.0d0*y(n-2)-y(n-3))/6.0d0
c  Solve the tridiagonal system
c  dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1,
c  with dy(1)=f(1), dy(n)=f(n).
	  do 10 i=2,n1
	  f(i)=g(i)*(3.0d0*(y(i+1)-y(i-1))-f(i-1))
10	continue
	dy(n)=f(n)
	  do 20 i=n1,1,-1
	  dy(i)=f(i)-g(i)*dy(i+1)
20	continue
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine splini
c  Splini must be called before splder to initialize array g in common.
c
	implicit double precision (a-h,o-z)
	dimension g(20001)
	common /spline/ g
	save /spline/
c
	g(1)=0.0d0
	  do 10 i=2,20001
	  g(i)=1.0d0/(4.0d0-g(i-1))
10	continue
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine splint(y,z,n)
c  Splint integrates a cubic spline, providing the ouput value
c  z = integral from 1 to n of s(i)di, where s(i) is the spline fit
c  to y(i).
c
	implicit double precision (a-h,o-z)
	dimension y(n)
c
	n1=n-1
c  Cubic fit to dy/di at boundaries.
c	dy1=(-11.0d0*y(1)+18.0d0*y(2)-9.0d0*y(3)+2.0d0*y(4))/6.0d0
	dy1=0.0d0
	dyn=(11.0d0*y(n)-18.0d0*y(n1)+9.0d0*y(n-2)-2.0d0*y(n-3))/6.0d0
c
	z=0.5d0*(y(1)+y(n))+(dy1-dyn)/12.0d0
	  do 10 i=2,n1
	  z=z+y(i)
10	continue
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function rombint(f,a,b,tol)
c  Rombint returns the integral from a to b of f(x)dx using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).  The function
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.
c
	parameter (MAXITER=16,MAXJ=5)
	implicit double precision (a-h,o-z)
	dimension g(MAXJ+1)
	external f
c
	h=0.5d0*(b-a)
	gmax=h*(f(a)+f(b))
	g(1)=gmax
	nint=1
	error=1.0d20
	i=0
10	  i=i+1
	  if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol))
     2      go to 40
c  Calculate next trapezoidal rule approximation to integral.
	  g0=0.0d0
	    do 20 k=1,nint
	    g0=g0+f(a+(k+k-1)*h)
20	  continue
	  g0=0.5d0*g(1)+h*g0
	  h=0.5d0*h
	  nint=nint+nint
	  jmax=min(i,MAXJ)
	  fourj=1.0d0
	    do 30 j=1,jmax
c  Use Richardson extrapolation.
	    fourj=4.0d0*fourj
	    g1=g0+(g0-g(j))/(fourj-1.0d0)
	    g(j)=g0
	    g0=g1
30	  continue
	  if (abs(g0).gt.tol) then
	    error=1.0d0-gmax/g0
	  else
	    error=gmax
	  end if
	  gmax=g0
	  g(jmax+1)=g0
	go to 10
40	rombint=g0
	if (i.gt.MAXITER.and.abs(error).gt.tol)
     2    write(*,*) 'Rombint failed to converge; integral, error=',
     3    rombint,error
	return
	end
