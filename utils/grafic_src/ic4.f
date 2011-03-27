cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ic4(idim,irand,iseed,itide,m1s,m2s,m3s,m1off,
     &    m2off,m3off,hanning,filename,astart,pk,dx,xoff,f,fc,fm)
c  Generate an unconstrained sample of (rho,psi1,psi2,psi3,phi) for
c  idim=0,1,2,3,4.
c  Input: idim, irand, iseed, itide, m?s, m?off, hanning, filename,
c         astart, pk, dx, xoff
c  irand=0: use randg to generate white noise, don't save.
c  irand=1: use randg to generate white noise, then save in real space
c    in filename.
c  irand=2: read filename to get random numbers.
c  iseed: 9-digit integer random number seed.  Beware that rand8 does not
c    give the same random numbers on 32-bit and 64-bit machines!
c  itide=0 to use full subvolume for computing f.
c  itide=1 to set xi=0 inside subvolume so as to get outer field.
c  itide=-1 to set xi=0 outside subvolume so as to get inner field.
c  m?s = size of next-level subvolume to split if itide.ne.0.
c  m?off = offset of next-level subvolume to split if itide.ne.0
c  hanning=T to apply hanning filter to f.
c  hanning=F to not apply hanning filter to f.
c  filename = file containing random numbers in real space.
c  astart = expansion factor
c  pk(ak,astart) = power spectrum function for wavenumber ak
c  dx = grid spacing.
c  xoff = offset to evaluate fields (e.g. use to shift baryon or cdm fields).
c  Output: f=fc (sampled field in real space), fm (maximum absolute value of f).
c  N.B. f and fc must point to the same place in memory - they are listed
c  separately in the subroutine call because f77 will not allow equivalencing
c  pointers.  The calling routine must pass the same pointer for each.
c
	implicit none
	include 'grafic1.inc'
c
	integer n12,n22,n32,n23,ngr2,npow
	parameter (n12=np1/2,n22=np2/2,n32=np3/2,n23=np2*np3)
	parameter (ngr2=n23*n12)
	parameter (npow=30720)
c
	integer idim,irand,iseed,itide,m1s,m2s,m3s,m1off,m2off,m3off
	logical hanning
	character*(*) filename
	real astart,pk,dx,xoff,fm,f(np1,np2,np3)
	complex fc(n12,np2,np3)
c
	complex fcny(np2,np3),z
	double precision twopi,avg,sigma,chisq
	parameter (twopi=6.283185307179586d0)
	integer j,i1,i2,i3,k1,k2,k3,k23,jp,ndof,modes
	real dk1,dk2,dk3,d3k,akmax,akmaxf,ak,ak23,akk,fact,xr
	real ak1,ak2,ak3,ak33,anu,dp,tf,tsav(0:npow),theta
	logical inflag
	external pk
c
	print*
	if (idim.lt.0.or.idim.gt.4) then
	  print*,'Error in ic4! Need idim=0,1,2,3,4'
	  stop
	end if
	if (irand.lt.0.or.irand.gt.2) then
	  print*,'Error in ic4! Need irand=0,1,2'
	  stop
	end if
	if (itide.lt.-1.or.itide.gt.1) then
	  print*,'Error in ic4! Need itide=-1,0,1'
	  stop
	end if
	if (itide.ne.0.and.(m1s.gt.0.5*np1.or.m2s.gt.0.5*np2
     &                  .or.m3s.gt.0.5*np3)) then
	  print*,'Error in ic4! Subvolume must be no larger than half ',
     &      'the size of the top grid'
	  print*,'Top grid size=',np1,np2,np3
	  print*,'Subvolume size=',m1s,m2s,m3s
	  stop
	end if
c
	dk1=twopi/(np1*dx)
	dk2=twopi/(np2*dx)
	dk3=twopi/(np3*dx)
	d3k=dk1*dk2*dk3
	akmax=twopi/dx
c  Precompute transfer function table for interpolation.
c  N.B. must go to at least sqrt(3)*akmax/2 unless use Hanning filter,
c  in which case go to akmax/2.
	if (hanning) then
	  akmaxf=akmax/2.0
	else
	  akmaxf=akmax
	end if
	do j=0,npow
	  ak=j*akmaxf/npow
	  tsav(j)=sqrt(pk(ak,astart)*d3k)
	  if (hanning) tsav(j)=tsav(j)*cos(0.5*ak*dx)
	  if (idim.gt.0.and.j.gt.0) then
	    tsav(j)=tsav(j)/(ak*ak)
	  end if
	end do
c  Get white noise sample.
	if (irand.lt.2) then
	  print*,'Warning: Generating new random numbers in ic4!'
	  call randini(iseed)
	  do i3=1,np3
	  do i2=1,np2
	  do i1=1,np1
	    call randg(xr)
	    f(i1,i2,i3)=xr
	  end do
	  end do
	  end do
c  Output white noise field.
	  if (irand.eq.1) then
	    print*,'Writing random numbers used in ic4 to ',filename
	    open(11,file=filename,form='unformatted')
	    rewind 11
	    write(11) np1,np2,np3,iseed
	    do i3=1,np3
	      write(11) ((f(i1,i2,i3),i1=1,np1),i2=1,np2)
	    end do
	    close(11)
	  end if
	else
c  irand=2.
c	  print*,'Reading random numbers used in ic4 from ',filename
	  open(11,file=filename,form='unformatted')
	  rewind 11
	  read(11) i1,i2,i3,iseed
	  if (i1.ne.np1.or.i2.ne.np2.or.i3.ne.np3) then
	    print*,'Error in ic4! file has np1,np2,np3=',i1,i2,i3
	    print*,' Expected ',np1,np2,np3
	    stop
	  end if
	  print*,'  Random numbers generated with iseed=',iseed
	    do i3=1,np3
	      read(11) ((f(i1,i2,i3),i1=1,np1),i2=1,np2)
	    end do
	    close(11)
	end if
c  Compute mean.
	avg=0.0
	do i3=1,np3
	do i2=1,np2
	do i1=1,np1
	  avg=avg+f(i1,i2,i3)
	end do
	end do
	end do
	avg=avg/(np1*np2*np3)
	fact=sqrt(1.0*np1*np2*np3)
c  Enforce zero mean for simulations with periodic boundary conditions.
c  Subtracting the mean of a white noise field is equivalent to applying
c  the mean constraint using the Hoffman-Ribak method.
c  Compute chisq for this sample, too.
	chisq=0.0
	do i3=1,np3
	  k3=i3-m3off
c  Apply periodic boundary conditions on the top grid for the
c  purpose of determining whether the point is in the subvolume.
	  if (k3.gt.np3) k3=k3-np3
	  if (k3.lt.1) k3=k3+np3
	do i2=1,np2
	  k2=i2-m2off
	  if (k2.gt.np2) k2=k2-np2
	  if (k2.lt.1) k2=k2+np2
	do i1=1,np1
	  k1=i1-m1off
	  if (k1.gt.np1) k1=k1-np1
	  if (k1.lt.1) k1=k1+np1
c  Subtract mean.
	  f(i1,i2,i3)=f(i1,i2,i3)-avg
c  Zero a portion if itide.ne.1, with periodic boundary conditions.
	    inflag=((k1.gt.0.and.k1.le.m1s).and.
     &         (k2.gt.0.and.k2.le.m2s).and.
     &         (k3.gt.0.and.k3.le.m3s))
	  if (itide.ne.0) then
c  Outer.
	    if (itide.eq.1.and.inflag) f(i1,i2,i3)=0.0
c  Inner.
	    if (itide.eq.-1.and.(.not.inflag)) f(i1,i2,i3)=0.0
	  end if
	  chisq=chisq+f(i1,i2,i3)**2
c  Standard deviation is fact, but divide by np1*np2*np3=fact**2 to
c  normalize f for FFT
	  f(i1,i2,i3)=f(i1,i2,i3)/fact
	end do
	end do
	end do
	ndof=np1*np2*np3-1
	if (itide.eq.1) then
	  ndof=ndof-m1s*m2s*m3s
	else if (itide.eq.-1) then
	  ndof=m1s*m2s*m3s
	end if
	anu=(chisq-ndof)/sqrt(float(ndof))
	print*,'ic4 white noise: chisq, dof, nu=',real(chisq),ndof,anu
c  Transform noise to k-space.
	call fft3rk(f,fcny,np1,np2,np3)
	chisq=0.0
	sigma=0.0
c  Generate unconstrained sample in Fourier transform space.
	  do 30 k3=1,np3
	  ak3=(k3-1)*dk3
	  if (k3.gt.n32) ak3=ak3-akmax
	  ak33=ak3*ak3
	    do 20 k2=1,np2
	    ak2=(k2-1)*dk2
	    if (k2.gt.n22) ak2=ak2-akmax
	    ak23=ak2*ak2+ak33
	    k23=k2-1+(k3-1)*np2
	      do 10 k1=1,n12
c  Do k1=n12+1 separately below.
	      ak1=(k1-1)*dk1
	      akk=ak1*ak1+ak23
	      ak=sqrt(akk)
c  Evaluate transfer function.
	      dp=npow*ak/akmaxf
	      if (dp.ge.npow) then
	        tf=0.0
	      else
	        jp=int(dp)
	        dp=dp-jp
	        tf=(1.0-dp)*tsav(jp)+dp*tsav(jp+1)
	      end if
c  Shift using offsets, with care at the Brillouin zone boundaries.
	      theta=ak1*xoff
	      if (k2.ne.n22+1) theta=theta+ak2*xoff
	      if (k3.ne.n32+1) theta=theta+ak3*xoff
	      z=cmplx(cos(theta),sin(theta))
c  These factors correctly average shifts at the Nyquist planes.
	      if (k2.eq.n22+1) z=z*cos(ak2*xoff)
	      if (k3.eq.n32+1) z=z*cos(ak3*xoff)
c  Convolve white noise with transfer function.
	      fc(k1,k2,k3)=fc(k1,k2,k3)*z*tf
	      if (idim.eq.1) then
	        fc(k1,k2,k3)=cmplx(0.0,1.0)*ak1*fc(k1,k2,k3)
	      else if (idim.eq.2) then
	        fc(k1,k2,k3)=cmplx(0.0,1.0)*ak2*fc(k1,k2,k3)
	        if (k2.eq.n22+1) fc(k1,k2,k3)=0.0
	      else if (idim.eq.3) then
	        fc(k1,k2,k3)=cmplx(0.0,1.0)*ak3*fc(k1,k2,k3)
	        if (k3.eq.n32+1) fc(k1,k2,k3)=0.0
	      end if
c  Double the contribution to account for modes with k1 > n12+1 (k1 < 0).
	      modes=2
	      if (k1.eq.1) modes=1
	      sigma=sigma+modes*fc(k1,k2,k3)*conjg(fc(k1,k2,k3))
	      if (idim.eq.0.or.idim.eq.4) then
	        chisq=chisq+modes*tf*tf
	      else
	        chisq=chisq+modes*tf*tf*akk/3.0
	      end if
10	    continue
c  Do k1=n12+1.
	    ak1=0.5*akmax
	    akk=ak1*ak1+ak23
	    ak=sqrt(akk)
c  Evaluate transfer function.
	    dp=npow*ak/akmaxf
	    if (dp.ge.npow) then
	      tf=0.0
	    else
	      jp=int(dp)
	      dp=dp-jp
	      tf=(1.0-dp)*tsav(jp)+dp*tsav(jp+1)
	    end if
c  Shift using offsets, with care at the Brillouin zone boundaries.
	    theta=0.0
	    if (k2.ne.n22+1) theta=theta+ak2*xoff
	    if (k3.ne.n32+1) theta=theta+ak3*xoff
	    z=cmplx(cos(theta),sin(theta))
c  These factors correctly average shifts at the Nyquist planes.
	    z=z*cos(ak1*xoff)
	    if (k2.eq.n22+1) z=z*cos(ak2*xoff)
	    if (k3.eq.n32+1) z=z*cos(ak3*xoff)
c  Convolve white noise with transfer function.
	    fcny(k2,k3)=fcny(k2,k3)*z*tf
	    if (idim.eq.1) then
	      fcny(k2,k3)=0.0
	    else if (idim.eq.2) then
	      fcny(k2,k3)=cmplx(0.0,1.0)*ak2*fcny(k2,k3)
	      if (k2.eq.n22+1) fcny(k2,k3)=0.0
	    else if (idim.eq.3) then
	      fcny(k2,k3)=cmplx(0.0,1.0)*ak3*fcny(k2,k3)
	      if (k3.eq.n32+1) fcny(k2,k3)=0.0
	    end if
	    modes=1
	    sigma=sigma+modes*fcny(k2,k3)*conjg(fcny(k2,k3))
	    if (idim.eq.0.or.idim.eq.4) then
	      chisq=chisq+modes*tf*tf
	    else
	      chisq=chisq+modes*tf*tf*akk/3.0
	    end if
20	  continue
30	continue
c  Enforce zero mean.
	fc(1,1,1)=cmplx(0.0,0.0)
	chisq=sqrt(chisq)
	sigma=sqrt(sigma)
c  Transform to position space.
	call fft3kr(fc,fcny,np1,np2,np3)
        fm=0.0
        do i3=1,np3
	do i2=1,np2
	do i1=1,np1
          fm=max(fm,abs(f(i1,i2,i3)))
	end do
	end do
	end do
	print*,'Statistics of ic4 for idim, itide=',idim,itide
	print*,'   Mean sigma, sampled sigma, maximum=',real(chisq),
     2    real(sigma),fm
c
	return
	end
