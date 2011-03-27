ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine wnsub(irand,iseed,tideflag,m1s,m2s,m3s,m1off,
     &    m2off,m3off,filename,noise,noiseny)
c  Generate a high-frequency sample of subgrid white noise, transform it
c  to k-space, and save it to file(s) as well as in arrays noise,noiseny.
c  Input: irand,iseed,tideflag,m1s,m2s,m3s,m1off,m2off,m3off,filename.
c  irand=0: use randg to generate white noise, don't save in real space.
c  irand=1: use randg to generate white noise, then save in real space
c    in filename.
c  irand=2: read filename to get random numbers in real space.
c  iseed: 9-digit integer random number seed.  Beware that rand8 does not
c    give the same random numbers on 32-bit and 64-bit machines!
c  tideflag=F to use full subvolume and output temp.wnsub_full.
c  tideflag=T to split subvolume into inner/outer and write temp.wnsub_in/out.
c  m?s = size of next-level subvolume to split if tideflag=T.
c  m?off = offset of next-level subvolume to split if tideflag=T.
c  N.B. noise and noiseny are returned with the k-space noise which
c  is also written out to temp.* files in the local directory.
c  Need to pass noise,noiseny for array storage.
c  
	implicit none
	include 'grafic2.inc'
c
	integer irand,iseed,m1s,m2s,m3s,m1off,m2off,m3off
	character*(*) filename
	logical tideflag
	real noise(n1c*nref2,n2c*nref2,n3c*nref2)
	complex noiseny(n2c*nref2,n3c*nref2)
c
	integer i1,i2,i3,k1,k2,k3,m1,m2,m3,n1,n2,n3,ndof,itime
	real xr,avg,fact,anu
	logical inflag
	double precision chisq
c
	if (irand.lt.0.or.irand.gt.2) then
	  print*,'Error in wnsub! Need irand=0,1,2'
	  stop
	end if
c
c  Counter to keep track of outer/inner.
	itime=0
c  Zero noise array.
10	do i3=1,n3c*nref2
	do i2=1,n2c*nref2
	do i1=1,n1c*nref2
	  noise(i1,i2,i3)=0.0
	end do
	end do
	end do
	if (irand.lt.2) then
	  if (itime.eq.0) print*,'Warning: Generating new random ',
     &      'numbers in wnsub!'
	  call randini(iseed)
c  Fill noise array.
	  do i3=1,n3c*nref
	  do i2=1,n2c*nref
	  do i1=1,n1c*nref
	    call randg(xr)
	    noise(i1,i2,i3)=xr
	  end do
	  end do
	  end do
c  Output subgrid noise field in real space.
	  if (irand.eq.1.and.itime.eq.0) then
	    print*,'Writing random numbers used in wnsub to ',filename
	    open(11,file=filename,form='unformatted')
	    rewind 11
	    write(11) n1c*nref,n2c*nref,n3c*nref,iseed
	    do i3=1,n3c*nref
	      write(11) ((noise(i1,i2,i3),i1=1,n1c*nref),i2=1,n2c*nref)
	    end do
	    close(11)
	  end if
	else
c  irand=2.
c	  if (itime.eq.0) print*,'Reading random numbers used in ',
c     &        'wnsub from ',filename
	  open(11,file=filename,form='unformatted')
	  rewind 11
	  read(11) n1,n2,n3,iseed
	  if (n1.ne.n1c*nref.or.n2.ne.n2c*nref.or.n3.ne.n3c*nref) then
	    print*,'Error in wnsub!'
	    print*,' File ',filename,' has n1,n2,n3=',n1,n2,n3
	    print*,' Expected ',n1c*nref,n2c*nref,n3c*nref
	    stop
	  end if
	  print*,'  Random numbers generated with iseed=',iseed
	  do i3=1,n3
	    read(11) ((noise(i1,i2,i3),i1=1,n1),i2=1,n2)
	  end do
	  close(11)
	end if
c  Normalization factor.
C	fact=sqrt(1.0*np1*np2*np3*nref*nref*nref)
c  Compensate for this change in normalization by increasing the
c  normalization of the transfer functions by the above factor.
c  This has the advantage of eliminating np1 everywhere.
	fact=1.0
c  Loop over coarse grid cells to be refined.
	do m3=0,n3c-1
	do m2=0,n2c-1
	do m1=0,n1c-1
c  Loop over refinements.
	  avg=0.0
	  do n3=1,nref
	  do n2=1,nref
	  do n1=1,nref
	    noise(n1+m1*nref,n2+m2*nref,n3+m3*nref)=
     &      noise(n1+m1*nref,n2+m2*nref,n3+m3*nref)*fact
	    avg=avg+noise(n1+m1*nref,n2+m2*nref,n3+m3*nref)
	  end do
	  end do
	  end do
	  avg=avg/(nref*nref*nref)
c  Now subtract mean over coarse grid cells.
	  do n3=1,nref
	  do n2=1,nref
	  do n1=1,nref
	    noise(n1+m1*nref,n2+m2*nref,n3+m3*nref)=
     &      noise(n1+m1*nref,n2+m2*nref,n3+m3*nref)-avg
	  end do
	  end do
	  end do
c  End loops over coarse grid.
	end do
	end do
	end do
c
c  Normalize noise as in eq. (18), compute chisq, and zero for tides.
	chisq=0.0
	do i3=1,n3c*nref2
	do i2=1,n2c*nref2
	do i1=1,n1c*nref2
c  Zero a portion if tideflag=T.
	  if (tideflag) then
	    inflag=((i1.gt.m1off.and.i1.le.m1off+m1s).and.
     &         (i2.gt.m2off.and.i2.le.m2off+m2s).and.
     &         (i3.gt.m3off.and.i3.le.m3off+m3s))
c  Outer.
	  if (itime.eq.0.and.inflag) noise(i1,i2,i3)=0.0
c  Inner.
	  if (itime.eq.1.and.(.not.inflag)) noise(i1,i2,i3)=0.0
	  end if
	  chisq=chisq+noise(i1,i2,i3)**2
	  noise(i1,i2,i3)=noise(i1,i2,i3)/(n1c*n2c*n3c*nref2**3)
	end do
	end do
	end do
	chisq=chisq*fact**2
	ndof=n1c*n2c*n3c*(nref**3-1)
	anu=(nref**3-1.0)/(nref**3)
	if (tideflag.and.itime.eq.0) then
	  ndof=ndof-m1s*m2s*m3s*anu
	else if (tideflag.and.itime.eq.1) then
	  ndof=m1s*m2s*m3s*anu
	end if
	anu=(chisq-ndof)/sqrt(float(ndof))
	print*,'wnsub white noise: chisq, dof, nu=',real(chisq),ndof,anu
c  Transform noise to k-space.
c  N.B. Do Not zero the k=0 mode!
	call fft3rk(noise,noiseny,n1c*nref2,n2c*nref2,n3c*nref2)
c  Write to disk.  This will be used later by grafic2.
c  Better to write to disk than to return, because will have to
c  reuse many times.
	if (tideflag.and.itime.eq.0) then
	  open(11,file='temp.wnsub_out',form='unformatted')
	else if (tideflag.and.itime.eq.1) then
	  open(11,file='temp.wnsub_in',form='unformatted')
	else
	  open(11,file='temp.wnsub_full',form='unformatted')
	end if
	rewind 11
	write(11) n1c*nref,n2c*nref2,n3c*nref2,m1s,m2s,m3s,m1off,
     &    m2off,m3off
        do k3=1,n3c*nref2
c  Assume that f77 stores complex numbers as pairs of adjacent reals.
	  write(11) ((noise(k1,k2,k3),k1=1,n1c*nref2),k2=1,n2c*nref2)
	end do
	write(11) ((noiseny(k2,k3),k2=1,n2c*nref2),k3=1,n3c*nref2)
	close(11)
c  Redo with outer if needed.
	if (tideflag.and.itime.eq.0) then
	  itime=1
	  go to 10
	end if
c
	return
	end
