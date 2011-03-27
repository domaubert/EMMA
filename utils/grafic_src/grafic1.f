cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	program grafic1
c  Generate initial conditions for cosmological N-body integration
c  as a Gaussian random field.
c  This version does not do constraints.  Instead, it produces output
c  compatible with grafic2 for multiscale initial conditions.
c  Disk use: lun=10 is output file, lun=11 is temp files.
c
	implicit none
	include 'grafic1.inc'
c
	real f(np1,np2,np3),slice(np1,np2)
	real dx,x1o,x2o,x3o,xoff,fm
	integer irand,iseed,itide,m1s,m2s,m3s,m1off,m2off,m3off
	integer m1t,m2t,m3t,m1offt,m2offt,m3offt
	integer i1,i2,i3,icomp,idim,nrefine,l1,l2,l3,m1,m2,m3
	double precision twopi,sigma,dsigma,rombin
	parameter (twopi=6.283185307179586d0)
	logical hanning
	character*80 filename
	real astart,omegam,omegav,h0,dladt,fomega,vfact,pbar,pcdm
	real asig,dxr,dpls,dplus,adp
	common /cosmoparms/ omegam,omegav,h0
	common /bigarray/ f
	common /dsig/ asig,dxr
	external pbar,pcdm,dsigma,rombin,dplus,adp
c
	print*
	print*,'Will generate initial conditions on grid of size ',
     &     np1,np2,np3
	print*
c
c  Initialize power spectrum.
	call pini
c  Initialize grid spacing.
	print*,'Enter dx (initial particle spacing in Mpc, not Mpc/h)'
	print*,'   or enter -boxlength in Mpc/h'
	print*,'   i.e. dx in Mpc if > 0, or -boxlength in Mpc/h if < 0'
	read(*,*) dx
	if (dx.lt.0.0) then
	  dx=-dx*100./(h0*np1)
	  print*,'  dx=',dx,' Mpc'
	end if
c  Set astart based on smallest grid spacing.
	print*,'Enter 1 if these initial conditions will not be ',
     &    'refined, otherwise'
	print*,'  enter the ultimate refinement factor for the ',
     &    'smallest grid spacing'
	print*,'  (This is used only to set astart.)'
	read(*,*) nrefine
	if (nrefine.lt.1) then
	  print*,'Error! refinement factor must be >= 1'
	  stop
	end if
	dxr=dx/nrefine
	asig=1.0d0
	sigma=2.0*twopi*rombin(dsigma,0.0d0,0.5d0*twopi/dxr,1.0d-7)
c  This is sigma at a=1.
	sigma=sqrt(sigma)
c  Normalize so that rms density flutuation=sigstart at starting
c  redshift scaling back the fluctuations as if they grew like cdm.
	dpls=sigstart/sigma*dplus(1.0,omegam,omegav)
	astart=adp(dpls,omegam,omegav)
	asig=astart
	sigma=2.0*twopi*rombin(dsigma,0.0d0,0.5d0*twopi/dxr,1.0d-7)
	sigma=sqrt(sigma)
	print*,'Scaling initial conditions to starting a=',astart
	print*,'  when sigma at ultimate refinement scale=',real(sigma)
	print*
c  velocity (proper km/s) =  Displacement (comoving Mpc at astart) * vfact.
c  vfact = dln(D+)/dtau where tau=conformal time.
	vfact=fomega(astart,omegam,omegav)
     &       *h0*dladt(astart,omegam,omegav)/astart
c
c  Now set output parameters.  There are two cases:
c  hanning=T: no further refinement.
c  hanning=F: prepare for further refinement.
c
	print*,'Enter 0 for final output or 1 for further refinement'
	read(*,*) irand
	if (irand.eq.1) then
	  hanning=.false.
	  print*,'Enter size (n1c,n2c,n3c) of the level-1 subgrid ',
     &      'to be extracted.'
	  print*,'  Sizes must be even numbers and be no larger than ',
     &      int(0.5*np1),int(0.5*np2),int(0.5*np3)
	  read(*,*) m1s,m2s,m3s
	  if (mod(m1s,2).ne.0.or.mod(m2s,2).ne.0.or.mod(m3s,2).ne.0) then
	    print*,'Error! Sizes must be even numbers!'
	    stop
	  end if
	  if (2*m1s.gt.np1.or.2*m2s.gt.np2.or.2*m3s.gt.np3) then
	    print*,'Error!  Subgrid is too large'
	    stop
	  end if
	  print*,'Enter offset of level-1 subgrid (m1off,m2off,m3off).'
	  print*,'  Offsets are relative to level-0 grid corner.'
	  print*,'  Offsets may be positive or negative, with absolute ',
     & ' values no larger than',int(0.5*np1),int(0.5*np2),int(0.5*np3)
	  read(*,*) m1off,m2off,m3off
c  Now get coordinates for tidal volume.
	  print*,'Enter size (m1t,m2t,m3t) of the final subvolume, ',
     &      'in units of top grid spacing'
	  read(*,*) m1t,m2t,m3t
	  if (m1t.gt.m1s.or.m2t.gt.m2s.or.m3t.gt.m3s) then
	    print*,'Error! Final subvolume cannot be larger than ',
     &        'level-1 subgrid'
	    stop
	  end if
	  print*,'Enter offset of final subvol. (m1offt,m2off2,m3offt).'
	  print*,'  Offsets are relative to level-0 grid corner.'
	  print*,'  Final subvolume must lie within level-1 subgrid'
	  read(*,*) m1offt,m2offt,m3offt
	  if (m1offt.lt.m1off.or.m1offt+m1t.gt.m1off+m1s.or.
     &        m2offt.lt.m2off.or.m2offt+m2t.gt.m2off+m2s.or.
     &        m3offt.lt.m3off.or.m3offt+m3t.gt.m3off+m3s) then
	    print*,'Error! Final subvolume isn''t contained within ',
     &             'level-1 subgrid'
	    stop
	  end if
c  Coordinates of the subgrid corner.
	  x1o=m1off*dx
	  x2o=m2off*dx
	  x3o=m3off*dx
	  print*,'Grafic1 will output level-1 subgrid as follows:'
	  print*,'  xmin, xmax = ',x1o,x1o+m1s*dx
	  print*,'  ymin, ymax = ',x2o,x2o+m2s*dx
	  print*,'  zmin, zmax = ',x3o,x3o+m3s*dx
	  print*,'Grafic1 will compute tides assuming final subvolume:'
	  print*,'  xmin, xmax = ',m1offt*dx,(m1offt+m1t)*dx
	  print*,'  ymin, ymax = ',m2offt*dx,(m2offt+m2t)*dx
	  print*,'  zmin, zmax = ',m3offt*dx,(m3offt+m3t)*dx
c  Open output file and write header.
	  open(10,file='grafic2.top',form='unformatted',status='unknown')
	  rewind 10
	  write(10)  2*m1s,2*m2s,2*m3s,dx,x1o,x2o,x3o,m1t,m2t,m3t,
     &      m1offt-m1off,m2offt-m2off,m3offt-m3off,hanning,astart,
     &      omegam,omegav,h0
	  print*,'Output file is grafic2.top use as input to subsequent',
     &           ' run of grafic2'
	else
	  hanning=.true.
	  m1s=np1
	  m2s=np2
	  m3s=np3
	  x1o=0.0
	  x2o=0.0
	  x3o=0.0
	  print*,'Setting output grid to ',m1s,m2s,m3s
	  print*,'Enter <RETURN> to skip output grid size'
	  read(*,*)
	  print*,'Enter <RETURN> to skip output grid offset'
	  read(*,*)
	  print*,'Enter <RETURN> to skip final grid size'
	  read(*,*)
	  print*,'Enter <RETURN> to skip final grid offset'
	  read(*,*)
	  print*,'Will produce output files ic_deltab, ic_velb[xyz],',
     &           ' ic_velc[xyz]'
	end if
c  Set parameters for subgrid noise.
	print*
	print*,'Subgrid white noise:'
	print*,'  Choose irand (1 or 2) from the following list:'
c	print*,'    irand=0 to generate new noise and don''t save it'
	print*,'    irand=1 to generate new noise and save to file'
	print*,'    irand=2 to read noise from existing file'
	print*,'Enter irand'
	read(*,*) irand
	if (irand.lt.0.or.irand.gt.2) then
	  print*,'Illegal value of irand'
	  stop
	end if
	print*,'  Enter random number seed (9-digit integer, ignored ',
     &    'if irand=2)'
	read(*,*) iseed
	print*,'  Enter filename of white noise file (or <RETURN> ',
     &    'if irand=0)'
	read(*,'(a80)') filename
c
c  Outer loop over components to be refined.
c
      do 10 icomp=0,12
c  0: baryon density.
c  1,2: inner,outer baryon x-velocity.
c  3,4: inner,outer baryon y-velocity.
c  5,6: inner,outer baryon z-velocity.
c  7,8: inner,outer CDM x-velocity.
c  9,10: inner,outer CDM y-velocity.
c  11,12: inner,outer CDM z-velocity.
c
c  If not generating grafic2.dat for further refinement, and
c    if icomp=odd, then skip.
	if (hanning.and.mod(icomp,2).eq.1) go to 10
c
c  Optional half-grid offset to velocity fields.
	if (hanning.and.icomp.gt.0.and.icomp.le.6) then
	   xoff=0.5*dx*offvelb
	else if (hanning.and.icomp.gt.6) then
	   xoff=0.5*dx*offvelc
	else
	  xoff=0.0
	end if
c  Don't need to recompute random numbers.
	if (icomp.gt.0.and.irand.eq.1) irand=2
c
c  Generate random fields.  idim=0,1,2,3 for density, displacement. 
	itide=0
	if (mod(icomp,2).eq.0) then
	  idim=icomp/2
c  Outer case.
	  if (.not.hanning.and.idim.gt.0) itide=1
	else
c  Inner case.
	  idim=(icomp+1)/2
	  if (.not.hanning.and.idim.gt.0) itide=-1
	end if
	if (idim.le.3) then
	  call ic4(idim,irand,iseed,itide,m1t,m2t,m3t,m1offt,m2offt,
     &      m3offt,hanning,filename,astart,pbar,dx,xoff,f,f,fm)
	else
	  idim=idim-3
	  call ic4(idim,irand,iseed,itide,m1t,m2t,m3t,m1offt,m2offt,
     &      m3offt,hanning,filename,astart,pcdm,dx,xoff,f,f,fm)

	end if
c
c  Prepare data for output.
c
	if (hanning) then
c  hanning=T, prepare final output files for simulation codes.
	  if (idim.gt.0) then
	    do i3=1,np3
	    do i2=1,np2
	    do i1=1,np1
c  Convert displacement to proper peculiar velocity in km/s.
	      f(i1,i2,i3)=f(i1,i2,i3)*vfact
	    end do
	    end do
	    end do
	  end if
c  Output files.
	  if (icomp.eq.0) then
	    open(11,file='ic_deltab',form='unformatted')
	  else if (icomp.eq.2) then
	    open(11,file='ic_velbx',form='unformatted')
	  else if (icomp.eq.4) then
	    open(11,file='ic_velby',form='unformatted')
	  else if (icomp.eq.6) then
	    open(11,file='ic_velbz',form='unformatted')
	  else if (icomp.eq.8) then
	    open(11,file='ic_velcx',form='unformatted')
	  else if (icomp.eq.10) then
	    open(11,file='ic_velcy',form='unformatted')
	  else if (icomp.eq.12) then
	    open(11,file='ic_velcz',form='unformatted')
	  end if
	  rewind 11
	  write(11) np1,np2,np3,dx,x1o+xoff,x2o+xoff,x3o+xoff,
     &      astart,omegam,omegav,h0
	  do i3=1,np3
	    write(11) ((f(i1,i2,i3),i1=1,np1),i2=1,np2)
	  end do
	  close(11)
	else
c  hanning=F, extract next-level subgrid and append to grafic2.dat.
c  First surround subvolume with 1/2-size buffer and wrap periodically.
	  sigma=0.0
	  do m3=1,2*m3s
	    i3=m3+m3off
c  Periodic boundary conditions on top grid.
	    if (i3.lt.1) i3=i3+np3
	    if (i3.gt.np3) i3=i3-np3
	    l3=i3
	    if (m3.gt.1.5*m3s) l3=l3-2*m3s
	    if (l3.lt.1) l3=l3+np3
	  do m2=1,2*m2s
	    i2=m2+m2off
	    if (i2.lt.1) i2=i2+np2
	    if (i2.gt.np2) i2=i2-np2
	    l2=i2
	    if (m2.gt.1.5*m2s) l2=l2-2*m2s
	    if (l2.lt.1) l2=l2+np2
	  do m1=1,2*m1s
	    i1=m1+m1off
	    if (i1.lt.1) i1=i1+np1
	    if (i1.gt.np1) i1=i1-np1
	    l1=i1
	    if (m1.gt.1.5*m1s) l1=l1-2*m1s
	    if (l1.lt.1) l1=l1+np1
	    slice(m1,m2)=f(l1,l2,l3)
	    sigma=sigma+f(l1,l2,l3)**2
	  end do
	  end do
	    write(10) ((slice(m1,m2),m1=1,2*m1s),m2=1,2*m2s)
	  end do
	  sigma=sqrt(sigma/(8*m1s*m2s*m3s))
	  print*,'After extraction, component ',icomp,
     &      ' has subvolume RMS=',real(sigma)
	end if
c  End loop over icomp.
10	continue
c
	if (.not.hanning) close(10)
	stop
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function dsigma(ak)
c  This function calculates the variance of density with a Hanning
c  filter at the grid Nyquist frequency.
c
	implicit none
	double precision dsigma,ak
	real asig,dx,p
	common /dsig/ asig,dx
	external p
c
	if (ak.le.0.0d0.or.ak.gt.3.1415926535898d0/dx) then
	  dsigma=0.0d0
	  return
	end if
c  Hanning filter.
	dsigma=ak*ak*p(real(ak),asig)*cos(0.5*ak*dx)
c
	return
	end
