c  Compute anti-aliasing filter in k-space using the
c  minimal k-space method.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine filta(astart,pk,dxf,hanning,idim,transfc)
c  Generate the anti-aliasing filter in k-space to correct the
c  coarse-grid sample for aliasing when it is spread to the subgrid.
c  The filter is stored in transfc; the real-space octant of size
c  (n1c*nref,n2c*nref,n3c*nref) is doubled in each dimension to handle
c  non-periodic BCs.  The k-space dimensions are effectively
c  complex(n1c*nref/2,n2c*nref,n3c*nref)+(1,n2c*nref,n3c*nref)
c  although the second part (the 1-direction Nyquist plane) is appended
c  onto transfc.
c  Input: astart = expansion factor
c  pk(ak,astart) = power spectrum function for wavenumber ak
c	     dxf = fine grid spacing (Mpc)
c        hanning = true or false for use of hanning filter on scale dxf
c           idim = component (0,1,2,3) to be returned (idim=0 for density,
c  idim=1,2,3 for three cartesian components of displacement)
c  Output: transfc (transfer function, k-space)
c  
	implicit none
	include 'grafic2.inc'
c
	integer n12,n23,ngr2,npow
	parameter (n12=nref*n1c,n23=nref2*n2c*nref2*n3c)
	parameter (ngr2=n23*n12)
	parameter (npow=30720)
c
	real astart,pk,dxf
	logical hanning
	integer idim
	complex transfc(ngr2+n23)
c
	real dk1,dk2,dk3,akmax,volumek,ak,ak0,ak1,ak2,ak3,ak01,ak02,ak03
	real dp,tf,tf0,tsav0(0:npow),tsav(0:npow)
	double precision twopi
	parameter (twopi=6.283185307179586d0)
	integer i,j,jp,k,kc,k1,k2,k3,k01,k02,k03,b1,b2,b3,k23,k23c
	external pk
c
c  N.B. The k-space here is specific to the subgrid of size
c  (n1c,n2c,n3c)*nref2*dxf.
	dk1=twopi/(n1c*nref2*dxf)
	dk2=twopi/(n2c*nref2*dxf)
	dk3=twopi/(n3c*nref2*dxf)
	akmax=twopi/dxf
	volumek=akmax**3
c  Initialize transfer function to zero.
	do j=1,ngr2+n23
	  transfc(j)=0.0
	end do
c  Precompute transfer function tables for interpolation of T(k0)
c  and T(k).  N.B. must go to at least sqrt(3)*akmax/2
c  unless use Hanning filter, in which case must go to akmax/2.
	do j=0,npow
	  ak0=j*akmax/npow
	  tsav0(j)=sqrt(pk(ak0,astart)*volumek)
	  if (hanning) then
	    ak=ak0/2.0
	    tsav(j)=sqrt(pk(ak,astart)*volumek)*cos(0.5*ak*dxf)
	  else
	    tsav(j)=tsav0(j)
	  end if
	end do
c  Loop over k-space.
	do k3=1,nref2*n3c
	  ak3=(k3-1)*dk3
	  if (k3.gt.nref*n3c) ak3=ak3-akmax
c  Project wavevector into the fundamental Brillouin zone.
	  b3=(k3-1)/(2*n3c)
	  k03=k3-b3*2*n3c
	  ak03=(k03-1)*dk3
	  if (k03.gt.n3c) ak03=ak03-akmax/nref
	do k2=1,nref2*n2c
	  k23=(k2-1+(k3-1)*nref2*n2c)*n12
	  ak2=(k2-1)*dk2
	  if (k2.gt.nref*n2c) ak2=ak2-akmax
	  b2=(k2-1)/(2*n2c)
	  k02=k2-b2*2*n2c
	  ak02=(k02-1)*dk2
	  if (k02.gt.n2c) ak02=ak02-akmax/nref
	do k1=2,n12
c  Do k1=1 and k1=n12+1 separately below.
	  k=k1+k23
	  ak1=(k1-1)*dk1
	  b1=(k1-1)/(2*n1c)
	  k01=k1-b1*2*n1c
	  ak01=(k01-1)*dk1
	  if (k01.gt.n1c) ak01=ak01-akmax/nref
	  ak=sqrt(ak1*ak1+ak2*ak2+ak3*ak3)
	  ak0=sqrt(ak01*ak01+ak02*ak02+ak03*ak03)
	  dp=npow*ak0/akmax
	  if (dp.ge.npow) then
	    tf0=0.0
	  else
	    jp=int(dp)
	    dp=dp-jp
	    tf0=(1.0-dp)*tsav0(jp)+dp*tsav0(jp+1)
	  end if
	  dp=npow*ak/akmax
	  if (hanning) dp=dp*2.0
	  if (dp.ge.npow) then
	    tf=0.0
	  else
	    jp=int(dp)
	    dp=dp-jp
	    tf=(1.0-dp)*tsav(jp)+dp*tsav(jp+1)
	  end if
c  Compute component of displacement from density, being careful
c  to zero it at Brillouin zone boundary.
	  if (idim.eq.1) then
	    tf=tf*ak1/(ak*ak)
	    if (k01.eq.1.or.k01.eq.n1c+1) then
	      tf0=0.0
	    else
	      tf0=tf0*ak01/(ak0*ak0)
	    end if
	  else if (idim.eq.2) then
	    if (k2.eq.1.or.k2.eq.nref*n2c+1) then
	      tf=0.0
	    else
	      tf=tf*ak2/(ak*ak)
	    end if
	    if (k02.eq.1.or.k02.eq.n2c+1) then
	      tf0=0.0
	    else
	      tf0=tf0*ak02/(ak0*ak0)
	    end if
	  else if (idim.eq.3) then
	    if (k3.eq.1.or.k3.eq.nref*n3c+1) then
	      tf=0.0
	    else
	      tf=tf*ak3/(ak*ak)
	    end if
	    if (k03.eq.1.or.k03.eq.n3c+1) then
	      tf0=0.0
	    else
	      tf0=tf0*ak03/(ak0*ak0)
	    end if
	  end if
c  Anti-aliasing filter.
c N.B. If fft this, then normalize by /(2*ngr2).
	  if (tf0.eq.0.0) then
	    transfc(k)=0.0
	  else
	    transfc(k)=tf/tf0
	  end if
C	TRANSFc(k)=TRANSFc(k)/(2*ngr2)
c  Close loop over k1.
	end do
c  Now do k1=1 and k1=n12+1.
	k23=k2-1+(k3-1)*nref2*n2c
	k23c=nref2*n2c+1-k2+(nref2*n3c+1-k3)*nref2*n2c
	if (k2.eq.1) k23c=k23c-nref2*n2c
	if (k3.eq.1) k23c=k23c-n23
	if (k23.le.k23c) then
	  do i=1,2
c  Do k1=1 and k1=n12+1.
	    if (i.eq.1) then
	      k1=1
	      k01=1
	      k=1+n12*k23
	      kc=1+n12*k23c
	      ak1=0.0
	      ak01=0.0
	    else
	      k1=n12+1
	      k=ngr2+1+k23
	      kc=ngr2+1+k23c
	      ak1=-0.5*akmax
	      b1=nref/2
	      k01=k1-b1*2*n1c
	      ak01=(k01-1)*dk1
	      if (k01.gt.n1c) ak01=ak01-akmax/nref
	    end if
	    ak=sqrt(ak1*ak1+ak2*ak2+ak3*ak3)
	    ak0=sqrt(ak01*ak01+ak02*ak02+ak03*ak03)
	    dp=npow*ak0/akmax
	    if (dp.ge.npow) then
	      tf0=0.0
	    else
	      jp=int(dp)
	      dp=dp-jp
	      tf0=(1.0-dp)*tsav0(jp)+dp*tsav0(jp+1)
	    end if
	    dp=npow*ak/akmax
	    if (hanning) dp=dp*2.0
	    if (dp.ge.npow) then
	      tf=0.0
	    else
	      jp=int(dp)
	      dp=dp-jp
	      tf=(1.0-dp)*tsav(jp)+dp*tsav(jp+1)
	    end if
c  Displacement=0 at k1=1, n12+1.
	    if (idim.eq.1) then
	      tf=0.0
	      if (k01.eq.1.or.k1.eq.n1c+1) then
	        tf0=0.0
	      else
	        tf0=tf0*ak01/(ak0*ak0)
	      end if
	    else if (idim.eq.2) then
	      if (k2.eq.1.or.k2.eq.nref*n2c+1) then
	        tf=0.0
	      else
	        tf=tf*ak2/(ak*ak)
	      end if
	      if (k02.eq.1.or.k02.eq.n2c+1) then
	        tf0=0.0
	      else
	        tf0=tf0*ak02/(ak0*ak0)
	      end if
	    else if (idim.eq.3) then
	      if (k3.eq.1.or.k3.eq.nref*n3c+1) then
	        tf=0.0
	      else
	        tf=tf*ak3/(ak*ak)
	      end if
	      if (k03.eq.1.or.k03.eq.n3c+1) then
	        tf0=0.0
	      else
	        tf0=tf0*ak03/(ak0*ak0)
	      end if
	    end if
c  Anti-aliasing filter.
c N.B. If fft this, then normalize by /(2*ngr2).
	    if (ak.eq.0.0) then
	      transfc(k)=1.0
	    else if (tf0.eq.0.0) then
	      transfc(k)=0.0
	    else
	      transfc(k)=tf/tf0
	    end if
C	TRANSFc(k)=TRANSFc(k)/(2*ngr2)
	    transfc(kc)=conjg(transfc(k))
c  Close loop over k1=1, nref*n1c+1.
	  end do
	end if
c  Close loops over k2, k3.
	end do
	end do

c  Fourier transform trans1c to spatial domain (testing only).
c  If do this, be sure to normalize as above.
C	call fft3rinv(transfc,nref2*n1c,nref2*n2c,nref2*n3c)

	return
	end
