cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	program grafic2
c  Generate multiscale gaussian initial conditions by applying mesh
c  refinement to an upper-level grid.
c  Before it can be run, the user must have a set of upper-grid data
c  files from either grafic1 or from a previous run of grafic2.
c  In the latter case, be sure to rename grafic2.dat.
c  grafic1 is used to produce data for the topmost periodic grid.
c  grafic2 is used to produce subgrid data.  It may be run on its
c  own output to produce nested refinements.
c
c  User must set parameters for the subgrid in include file grafic2.inc.
c
c  Memory requirement: a little more than 2*n1c*n2c*n3c*nref2**3.
c  If this is too much, then can split off wnsub and extrsub to be
c  stand-alone programs and reduce noisec to complex(n1c*nref,n2c*nref).
c  However, this will greatly increase the scratch disk requirement
c  because extrsub will have to write out all tidal fields onto disk
c  instead of returning them one at a time in noisec.
c  Disk use: lun=10 is input file, lun=12 is output file, lun=11 is temp files.
c
	implicit none
	include 'grafic2.inc'
c
	real transf(n1c*nref2,n2c*nref2,n3c*nref2)
	real noise(n1c*nref2,n2c*nref2,n3c*nref2)
	real slice(n1c*nref2,n2c*nref2)
	complex transfc((n1c*nref+1)*n2c*n3c*nref2*nref2)
	complex noisec(n1c*nref,n2c*nref2,n3c*nref2)
	complex noiseny(n2c*nref2,n3c*nref2)
	common /bigarrays/ transfc,noise,noiseny
	equivalence (transf,transfc),(noise,noisec,slice)
c
	integer m1s,m2s,m3s,irand,iseed,icomp,idim,m1off,m2off,m3off
	integer m1t,m2t,m3t,m1offt,m2offt,m3offt
	integer i1,i2,i3,k1,k2,k3,l1,l2,l3,m1,m2,m3,n1,n2,n3,k,k23
	logical hanning
	real dx,x1o,x2o,x3o,dk1,dk2,dk3,xoff,theta
	complex z
	double precision twopi,sigma
	parameter (twopi=6.283185307179586d0)
	character*80 filename
	real astart,omegam,omegav,h0,dladt,fomega,vfact,pbar,pcdm
	common /cosmoparms/ omegam,omegav,h0
	external pbar,pcdm
c
c  Structure of grafic2:
c   Loop over 13 fields: deltab, velinb(3), veloutb(3), velinc(3), veloutc(3).
c    N.B. b and c refer to baryons and CDM; in and out refer to tidal fields.
c    For each field, read the upper-level grid data, then convolve with
c     the appropriate anti-aliasing filter.
c     Then call wnsub to generate subgrid noise and convolve noise with the
c     appropriate transfer function.
c     Add long+short results and output to file.
c
c  Grafic2 has two different output modes:
c
c  The first mode produces output density and velocity fields after
c    refinement.  The output data files are ic_deltab, ic_velbn, ic_velcn.
c    The velocity fields are complete (not tidal fields).
c  The second mode produces intermediate output that will be further
c    refined before ic_* files are written.  In this case the output
c    is in one file, grafic2.dat, which provides the input for a
c    subsequent run of grafic2.  grafic2.dat must be renamed so that
c    it is not overwritten in a later run.
c  These two modes allow grafic2 to perform interative refinement
c  to multiple levels and grids.  The output mode is selected at runtime
c  when the user answers the query,
c    'Enter 0 for final output or 1 for further refinement'
c
	print*,'Will refine upper level in subgrid of size ',
     &    '(n1c,n2c,n3c)=',n1c,n2c,n3c
	print*,'Refinement factor=',nref
	print*
c
c  Open top grid file (lun=10).
c
	print*,'Enter input filename containing upper level grid data'
	read(*,'(a80)') filename
	if (filename.eq."grafic2.dat") then
	  print*,'You must rename your upper-level file!'
	  print*,'grafic2.dat is reserved for the output filename'
	  stop
	end if
	open(10,file=filename,form='unformatted',status='old')
	rewind 10
	read(10)  m1s,m2s,m3s,dx,x1o,x2o,x3o,m1t,m2t,m3t,
     &    m1offt,m2offt,m3offt,hanning,astart,omegam,omegav,h0
	if (m1s/2.ne.n1c.or.m2s/2.ne.n2c.or.m3s/2.ne.n3c) then
	  print*,'Error! Subgrid has wrong size!'
	  print*,filename,' has ',m1s,m2s,m3s
	  print*,'  Should be',2*n1c,2*n2c,2*n3c
	  print*,'  Change subgrid size in grafic2.inc and remake'
	  stop
	end if
	if (hanning) then
	  print*,'File: ',filename
	  print*,'  The coarse grid sample has been convolved with ',
     2      'a Hanning filter!'
	  print*,'  Error! Mesh refinement requires hanning=F!'
	  stop
	end if
c
	print*,'Reading upper-level grid data from file ',filename
	print*,' astart, Omegam, Omegav, H0=',astart,omegam,omegav,h0
	print*,' dx (Mpc), offsets/dx=',dx,x1o,x2o,x3o
	print*,' Tidal volume size/dx=',m1t,m2t,m3t
	print*,' Tidal volume offsets/dx=',m1offt,m2offt,m3offt
	print*
c  Refine grid spacing.
	dx=dx/nref
	m1t=m1t*nref
	m2t=m2t*nref
	m3t=m3t*nref
	m1offt=m1offt*nref
	m2offt=m2offt*nref
	m3offt=m3offt*nref
c
c  Initialize power spectrum and other parameters.
c
	call pini
c  velocity (proper km/s) =  Displacement (comoving Mpc at astart) * vfact.
c  vfact = dln(D+)/dtau where tau=conformal time.
	vfact=fomega(astart,omegam,omegav)
     &       *h0*dladt(astart,omegam,omegav)/astart
c
c  Now set output parameters.  There are output modes:
c  hanning=T: no further refinement.
c  hanning=F: prepare for further refinement.
c
	print*,'Enter 0 for final output or 1 for further refinement'
	read(*,*) irand
	if (irand.eq.1) then
	  hanning=.false.
	  print*,'Enter size (m1s,m2s,m3s) of the next-level subgrid ',
     &      'to be prepared for extraction.'
	  print*,'  Sizes must be even numbers and be no larger than ',
     &      int(0.5*nref*n1c),int(0.5*nref*n2c),int(0.5*nref*n3c)
	  print*,'      and ',m1t,m2t,m3t
	  read(*,*) m1s,m2s,m3s
	  if (mod(m1s,2).ne.0.or.mod(m2s,2).ne.0.or.mod(m3s,2).ne.0) then
	    print*,'Error! Sizes must be even numbers!'
	    stop
	  end if
	  if (2*m1s.gt.n1c*nref.or.2*m2s.gt.n2c*nref.or.
     &      2*m3s.gt.n3c*nref) then
	    print*,'Error!  Subgrid is too large'
	    stop
	  end if
	  if (m1s.lt.m1t.or.m2s.lt.m2t.or.m3s.lt.m3t) then
	    print*,'Error! Subgrid cannot be less than tidal subvolume!'
	    stop
	  end if
	  print*,'Enter offset of next-level subgrid (m1off,m2off,m3off).'
	  print*,'  m1off must lie in the range',
     &     int(0.5*m1s),int(n1c*nref-1.5*m1s)
	  print*,'  m2off must lie in the range',
     &     int(0.5*m2s),int(n2c*nref-1.5*m2s)
	  print*,'  m3off must lie in the range',
     &     int(0.5*m3s),int(n3c*nref-1.5*m3s)
	  read(*,*) m1off,m2off,m3off
	  if ((m1off.lt.0.5*m1s.or.m1off.gt.n1c*nref-1.5*m1s).or.
     &      (m2off.lt.0.5*m2s.or.m2off.gt.n2c*nref-1.5*m2s).or.
     &      (m3off.lt.0.5*m3s.or.m3off.gt.n3c*nref-1.5*m3s)) then
	    print*,'Error! (m1off,m2off,m3off) out of bounds'
	    stop
	  end if
	  if (m1offt.lt.m1off.or.m1offt+m1t.gt.m1off+m1s.or.
     &        m2offt.lt.m2off.or.m2offt+m2t.gt.m2off+m2s.or.
     &        m3offt.lt.m3off.or.m3offt+m3t.gt.m3off+m3s) then
	    print*,'Error! Final subvolume isn''t contained within ',
     &             'next-level subgrid'
	    stop
	  end if
c  Coordinates of the subgrid corner.
	  x1o=x1o+m1off*dx
	  x2o=x2o+m2off*dx
	  x3o=x3o+m3off*dx
c  Open output file and write header.
	  open(12,file='grafic2.dat',form='unformatted',status='new')
	  rewind 12
	  write(12) 2*m1s,2*m2s,2*m3s,dx,x1o,x2o,x3o,m1t,m2t,m3t,
     &      m1offt-m1off,m2offt-m2off,m3offt-m3off,hanning,astart,
     &      omegam,omegav,h0
	  print*
	  print*,'***Output file is grafic2.dat; rename it for use as ',
     &           ' input to subsequent run of grafic2'
	else
	  hanning=.true.
	  m1s=n1c*nref
	  m2s=n2c*nref
	  m3s=n3c*nref
	  m1off=0
	  m2off=0
	  m3off=0
	  m1t=m1s
	  m2t=m2s
	  m3t=m3s
	  m1offt=0
	  m2offt=0
	  m3offt=0
	  print*,'Setting output grid to ',m1s,m2s,m3s
	  print*,'Enter <RETURN> to skip output grid size'
	  read(*,*)
	  print*,'Enter <RETURN> to skip output grid offset'
	  read(*,*)
	  print*,'Will produce output files ic_deltab, ic_velb[xyz],',
     &           ' ic_velc[xyz]'
	end if
c  Grid frequencies for phase shifts in case hanning=T for velb and velc.
	dk1=twopi/(n1c*nref2*dx)
	dk2=twopi/(n2c*nref2*dx)
	dk3=twopi/(n3c*nref2*dx)
c  Set parameters for subgrid noise.
	print*,'Subgrid white noise:'
	print*,'  Choose irand (1 or 2) from the following list:'
c	print*,'    irand=0 to generate new noise and don''t save it'
	print*,'    irand=1 to generate new noise and save to file'
	print*,'    irand=2 to read existing noise from file'
	print*,'Enter irand'
	read(*,*) irand
	if (irand.lt.0.or.irand.gt.2) then
	  print*,'Illegal value of irand'
	  stop
	end if
	print*,'  Enter random number seed (9-digit integer, ignored ',
     &    'if irand=2)'
	read(*,*) iseed
	print*,'  Enter filename of wnsub noise file (or <RETURN> ',
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
c  Optional half-grid offset to velocity fields.
	if (hanning.and.icomp.gt.0.and.icomp.le.6) then
	   xoff=0.5*dx*offvelb
	else if (hanning.and.icomp.gt.6) then
	   xoff=0.5*dx*offvelc
	else
	  xoff=0.0
	end if
c
c  First do long-wavelength part from upper-level grid.
c
c  Read data from upper-level grid.
	do m3=0,2*n3c-1
	  read(10) ((noise(1+m1*nref,1+m2*nref,1+m3*nref),m1=0,2*n1c-1),
     &      m2=0,2*n2c-1)
	end do
c  Spread to subgrid and normalize for FFT.
	sigma=0.0
	do m3=0,2*n3c-1
	do m2=0,2*n2c-1
	do m1=0,2*n1c-1
	  theta=noise(1+m1*nref,1+m2*nref,1+m3*nref)
	  sigma=sigma+theta**2
	  do n3=1,nref
	  do n2=1,nref
	  do n1=1,nref
c  Spread and normalize for FFT.
	    noise(n1+m1*nref,n2+m2*nref,n3+m3*nref)=
     & theta/(n1c*n2c*n3c*nref2**3)
	  end do
	  end do
	  end do
	end do
	end do
	end do
	sigma=sqrt(sigma/(8*n1c*n2c*n3c))
	print*
	print*,'Before refinement: component ',icomp,
     &    ' has subvolume RMS=',real(sigma)
c  Transform upper-level data to k-space.
c  N.B. Do Not zero the k=0 mode!
	call fft3rk(noise,noiseny,n1c*nref2,n2c*nref2,n3c*nref2)
c  Get anti-aliasing filter.
	if (icomp.eq.0) then
c  Minimal k-space anti-aliasing filter for density.
	  idim=0
	  call filta(astart,pbar,dx,hanning,idim,transfc)
	else if (mod(icomp,2).eq.1) then
c  Minimal k-space anti-aliasing filter for inner velocity.
	  idim=(icomp+1)/2
	  if (idim.le.3) then
	    call filta(astart,pbar,dx,hanning,idim,transfc)
	  else
	    idim=idim-3
	    call filta(astart,pcdm,dx,hanning,idim,transfc)
	  end if
	else
c  Sharp k-space anti-aliasing filter for outer velocity.
	  idim=icomp/2
	  if (idim.gt.3) idim=idim-3
	  call filtak(hanning,transfc)
	end if
c  Convolve upper-level sample with anti-aliasing filter.
	sigma=0.0
	do k3=1,n3c*nref2
	  l3=k3-1
	  if (k3.gt.n3c*nref) l3=l3-n3c*nref2
	  do k2=1,n2c*nref2
	    k23=k2+(k3-1)*n2c*nref2
	    l2=k2-1
	    if (k2.gt.n2c*nref) l2=l2-n2c*nref2
	    do k1=1,n1c*nref
	      k=k1+(k23-1)*n1c*nref
	      l1=k1-1
C  For testing with no noise, uncomment the next line.
CCC      noisec(k1,k2,k3)=1.0/sqrt(1.0*n1c*n2c*n3c*nref2**3)
c  Shift using offsets, with care at the Brillouin zone boundaries.
	      theta=l1*xoff*dk1
	      if (k2.ne.n2c*nref+1) theta=theta+l2*xoff*dk2
	      if (k3.ne.n3c*nref+1) theta=theta+l3*xoff*dk3
	      z=cmplx(cos(theta),sin(theta))
c  These factors correctly average shifts at the Nyquist planes.
	      if (k2.eq.n2c*nref+1) z=z*cos(l2*xoff*dk2)
	      if (k3.eq.n3c*nref+1) z=z*cos(l3*xoff*dk3)
c  Convolution.
	      transfc(k)=transfc(k)*noisec(k1,k2,k3)*z
	      if (k1.eq.1) then
	        sigma=sigma+transfc(k)*conjg(transfc(k))
 	      else
	        sigma=sigma+2*transfc(k)*conjg(transfc(k))
	      end if
	    end do
c  k1=kny.
	    k1=n1c*nref+1
	    l1=k1-1
 	    k=n1c*n2c*n3c*nref*nref2*nref2+k23
C  For testing with no noise, uncomment the next line.
CCC      noiseny(k2,k3)=1.0/sqrt(1.0*n1c*n2c*n3c*nref2**3)
c  Shift using offsets, with care at the Brillouin zone boundaries.
	    theta=0.0
	    if (k2.ne.n2c*nref+1) theta=theta+l2*xoff*dk2
	    if (k3.ne.n3c*nref+1) theta=theta+l3*xoff*dk3
	    z=cmplx(cos(theta),sin(theta))
c  These factors correctly average shifts at the Nyquist planes.
	    z=z*cos(l1*xoff*dk1)
	    if (k2.eq.n2c*nref+1) z=z*cos(l2*xoff*dk2)
	    if (k3.eq.n3c*nref+1) z=z*cos(l3*xoff*dk3)
	    transfc(k)=transfc(k)*noiseny(k2,k3)*z
	    sigma=sigma+transfc(k)*conjg(transfc(k))
	  end do
	end do
c	print*,'Long wavelength sample after anti-aliasing filter ',
c     &    'has rms=',real(sqrt(sigma))
c  Transform back to real space.
	call fft3rinv(transf,n1c*nref2,n2c*nref2,n3c*nref2)
c  Save in a temporary file.
	if (mod(icomp,2).eq.0) then
c  Beware: temp.out is used for full density (icomp=0) case as well as
c  outer velocity field.
	  open (11,file='temp.out',form='unformatted',status='unknown')
	else
	  open (11,file='temp.in',form='unformatted',status='unknown')
	end if
	rewind 11
	write(11) icomp,n1c,n2c,n3c,nref,dx
c  Write long-wavelength data.
	do i3=1,n3c*nref2
	  write(11) ((transf(i1,i2,i3),i1=1,n1c*nref2),i2=1,n2c*nref2)
	end do
	close(11)
c
c  If not preparing for further refinement, then don't need
c  short-wavelength noise twice.
	if (hanning.and.mod(icomp,2).eq.1) go to 10
c
c  Now compute short-wavelength noise part on subgrid.
c
c  Get subgrid noise sample.
	if (icomp.eq.0) then
	  if (.not.hanning) then
c  Generate inner+outer tidal noise fields for later use if needed.
	    call wnsub(irand,iseed,.true.,m1t,m2t,m3t,m1offt,
     &        m2offt,m3offt,filename,noisec,noiseny)
c  Don't need to recompute random numbers.
	    if (irand.eq.1) irand=2
	  end if
c  Generate full noise sample.
	  call wnsub(irand,iseed,.false.,m1t,m2t,m3t,m1offt,
     &      m2offt,m3offt,filename,noisec,noiseny)
	  if (irand.eq.1) irand=2
c  Read white noise sample from disk if needed
	else
	  if (hanning) then
c  Read full noise sample from temp file written by wnsub.
	    open(11,file='temp.wnsub_full',form='unformatted')
c  Read inner or outer  noise sample from temp file.
	  else if (mod(icomp,2).eq.0) then
c  Read outer noise sample if hanning=F and icomp even (further refinement).
	    open(11,file='temp.wnsub_out',form='unformatted')
	  else 
c  Read inner noise sample if hanning=F and icomp odd (further refinement).
	    open(11,file='temp.wnsub_in',form='unformatted')
	  end if
	  rewind 11
	  read(11) n1,n2,n3,m1t,m2t,m3t,m1offt,m2offt,m3offt
	  do k3=1,n3c*nref2
	    read(11) ((noisec(k1,k2,k3),k1=1,n1c*nref),k2=1,n2c*nref2)
	  end do
	  read(11) ((noiseny(k2,k3),k2=1,n2c*nref2),k3=1,n3c*nref2)
	  close(11)
	end if
c  Get transfer function.
	if (idim.eq.0) then
c  Density.
	  if (hanning) then
c  Spherical transfer function.
	    call transfers(astart,pbar,dx,hanning,idim,transf)
	    call fft3r(transfc,n1c*nref2,n2c*nref2,n3c*nref2)
	  else
c  Minimal k-space transfer function.
	    call transfd(astart,pbar,dx,hanning,transfc)
	  end if
	else
c  Displacement.
	  if (icomp.le.6) then
	    call transfers(astart,pbar,dx,hanning,idim,transf)
	  else
	    call transfers(astart,pcdm,dx,hanning,idim,transf)
	  end if
	  call fft3r(transfc,n1c*nref2,n2c*nref2,n3c*nref2)
	end if
c  Convolve subgrid noise with transfer function.
	sigma=0.0
	do k3=1,n3c*nref2
	  l3=k3-1
	  if (k3.gt.n3c*nref) l3=l3-n3c*nref2
	  do k2=1,n2c*nref2
	    k23=k2+(k3-1)*n2c*nref2
	    l2=k2-1
	    if (k2.gt.n2c*nref) l2=l2-n2c*nref2
	    do k1=1,n1c*nref
	      k=k1+(k23-1)*n1c*nref
	      l1=k1-1
C  For testing with no noise, uncomment the next line.
CCC      noisec(k1,k2,k3)=1.0/sqrt(1.0*n1c*n2c*n3c*nref2**3)
c  Shift using offsets, with care at the Brillouin zone boundaries.
	      theta=l1*xoff*dk1
	      if (k2.ne.n2c*nref+1) theta=theta+l2*xoff*dk2
	      if (k3.ne.n3c*nref+1) theta=theta+l3*xoff*dk3
	      z=cmplx(cos(theta),sin(theta))
c  These factors correctly average shifts at the Nyquist planes.
	      if (k2.eq.n2c*nref+1) z=z*cos(l2*xoff*dk2)
	      if (k3.eq.n3c*nref+1) z=z*cos(l3*xoff*dk3)
c  Convolution.
	      transfc(k)=transfc(k)*noisec(k1,k2,k3)*z
	      if (k1.eq.1) then
	        sigma=sigma+transfc(k)*conjg(transfc(k))
 	      else
	        sigma=sigma+2*transfc(k)*conjg(transfc(k))
	      end if
	    end do
c  k1=kny.
	    k1=n1c*nref+1
	    l1=k1-1
 	    k=n1c*n2c*n3c*nref*nref2*nref2+k23
C  For testing with no noise, uncomment the next line.
CCC      noiseny(k2,k3)=1.0/sqrt(1.0*n1c*n2c*n3c*nref2**3)
c  Shift using offsets, with care at the Brillouin zone boundaries.
	    theta=0.0
	    if (k2.ne.n2c*nref+1) theta=theta+l2*xoff*dk2
	    if (k3.ne.n3c*nref+1) theta=theta+l3*xoff*dk3
	    z=cmplx(cos(theta),sin(theta))
c  These factors correctly average shifts at the Nyquist planes.
	    z=z*cos(l1*xoff*dk1)
	    if (k2.eq.n2c*nref+1) z=z*cos(l2*xoff*dk2)
	    if (k3.eq.n3c*nref+1) z=z*cos(l3*xoff*dk3)
	    transfc(k)=transfc(k)*noiseny(k2,k3)*z
	    sigma=sigma+transfc(k)*conjg(transfc(k))
	  end do
	end do
c	print*,'Short wavelength sample after transfer function ',
c     &    'has rms=',real(sqrt(sigma))
c  Transform back to real space.
	call fft3rinv(transf,n1c*nref2,n2c*nref2,n3c*nref2)
c
c  Combine long and short, in and out data for output.
c
	if (hanning) then
c  hanning=T, prepare final output files for simulation codes.
c  Only get here if icomp=even.
c
c  Add temp.out and transf
	  open (11,file='temp.out',form='unformatted',status='old')
	  rewind 11
	  read(11) i1,n1,n2,n3,m1,dx
c  Read long-wavelength data and add to short in transf.
	  sigma=0.0
	  do i3=1,n3c*nref2
	    read(11) ((slice(i1,i2),i1=1,n1c*nref2),i2=1,n2c*nref2)
	    do i2=1,n2c*nref2
	    do i1=1,n1c*nref2
	      transf(i1,i2,i3)=transf(i1,i2,i3)+slice(i1,i2)
	if (i1.le.n1c*nref.and.i2.le.n2c*nref.and.i3.le.n3c*nref) then
	      sigma=sigma+transf(i1,i2,i3)**2
	end if
	    end do
	    end do
	  end do
	  close(11)
	  sigma=sqrt(sigma/(1.0*n1c*n2c*n3c*nref**3))
c  Add temp.in for displacement.
	  if (icomp.gt.0) then
	    open (11,file='temp.in',form='unformatted',status='old')
	    rewind 11
	    read(11) i1,n1,n2,n3,m1,dx
	    sigma=0.0
	    do i3=1,n3c*nref2
	      read(11) ((slice(i1,i2),i1=1,n1c*nref2),i2=1,n2c*nref2)
	      do i2=1,n2c*nref2
	      do i1=1,n1c*nref2
	        theta=transf(i1,i2,i3)+slice(i1,i2)
	if (i1.le.n1c*nref.and.i2.le.n2c*nref.and.i3.le.n3c*nref) then
	        sigma=sigma+theta**2
	end if
c  Multiply by vfact to convert displacement to proper peculiar velocity
c  in km/s at astart.
	        transf(i1,i2,i3)=theta*vfact
	      end do
	      end do
	    end do
	    close(11)
	    sigma=sqrt(sigma/(n1c*n2c*n3c*nref**3))
	  end if
	  print*,'After refinement, component ',icomp,
     &      ' has subvolume RMS=',real(sigma)
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
	  write(11) n1c*nref,n2c*nref,n3c*nref,dx,x1o+xoff,
     &      x2o+xoff,x3o+xoff,astart,omegam,omegav,h0
	  do i3=1,n3c*nref
	    write(11) ((transf(i1,i2,i3),i1=1,n1c*nref),i2=1,n2c*nref)
	  end do
	  close(11)
	else
c  hanning=F, extract next-level subgrid and append to grafic2.dat.
c
c  icomp=0: output temp.out plus transf
c  icomp=odd: output temp.in plus transf
c  icomp=even>1: output temp.out plus transf
c
c  Add temp.in/out and transf for delta_b and vel_out.
	  if (mod(icomp,2).eq.0) then
	    open (11,file='temp.out',form='unformatted',status='old')
	  else
	    open (11,file='temp.in',form='unformatted',status='old')
	  end if
	  rewind 11
	  read(11) i1,n1,n2,n3,m1,dx
c  Read long-wavelength data and add to short in transf.
	  do i3=1,n3c*nref2
	    read(11) ((slice(i1,i2),i1=1,n1c*nref2),i2=1,n2c*nref2)
	    do i2=1,n2c*nref2
	    do i1=1,n1c*nref2
	      transf(i1,i2,i3)=transf(i1,i2,i3)+slice(i1,i2)
	    end do
	    end do
	  end do
	  close(11)
c  Surround subvolume with 1/2-size buffer and wrap periodically.
	  sigma=0.0
	  do m3=1,2*m3s
	    i3=m3+m3off
	    l3=i3
	    if (m3.gt.1.5*m3s) l3=l3-2*m3s
	  do m2=1,2*m2s
	    i2=m2+m2off
	    l2=i2
	    if (m2.gt.1.5*m2s) l2=l2-2*m2s
	  do m1=1,2*m1s
	    i1=m1+m1off
	    l1=i1
	    if (m1.gt.1.5*m1s) l1=l1-2*m1s
	    slice(m1,m2)=transf(l1,l2,l3)
	    sigma=sigma+transf(l1,l2,l3)**2
	  end do
	  end do
	    write(12) ((slice(m1,m2),m1=1,2*m1s),m2=1,2*m2s)
	  end do
	  sigma=sqrt(sigma/(8*m1s*m2s*m3s))
	  print*,'After refinement, component ',icomp,
     &      ' has subvolume RMS=',real(sigma)
	end if
c  End loop over icomp.
10	continue
c
	close(10)
	if (.not.hanning) close(12)
	stop
	end
