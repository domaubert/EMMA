c  Compute density transfer function in k-space using the
c  minimal k-space method.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine transfd(astart,pk,dxf,hanning,transfc)
c  Generate the transfer function in k-space.
c  The transfer function is stored in transfc; the real-space octant of size
c  (n1c*nref,n2c*nref,n3c*nref) is doubled in each dimension to handle
c  non-periodic BCs.  The k-space dimensions are effectively
c  complex(n1c*nref/2,n2c*nref,n3c*nref)+(1,n2c*nref,n3c*nref)
c  although the second part (the 1-direction Nyquist plane) is appended
c  onto transfc.
c  Input: astart = expansion factor
c  pk(ak,astart) = power spectrum function for wavenumber ak
c	     dxf = fine grid spacing (Mpc)
c        hanning = true or false for use of hanning filter on scale dxf
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
	complex transfc(ngr2+n23)
c
	real dk1,dk2,dk3,akmax,volumek,akmaxf,ak,ak1,ak2,ak3,dp,tf
	real tsav(0:npow)
	integer j,jp,k,k1,k2,k3,k23
	double precision twopi
	parameter (twopi=6.283185307179586d0)
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
c  Precompute transfer function table for interpolation.
c  N.B. must go to at least sqrt{3}*akmax/2
c  unless use Hanning filter, in which case go to akmax/2.
	if (hanning) then
	  akmaxf=akmax/2.0
	else
	  akmaxf=akmax
	end if
	do j=0,npow
	  ak=j*akmaxf/npow
	  tsav(j)=sqrt(pk(ak,astart)*volumek)
	  if (hanning) tsav(j)=tsav(j)*cos(0.5*ak*dxf)
	end do
c  Loop over k-space.
	do k3=1,nref2*n3c
	  ak3=(k3-1)*dk3
	  if (k3.gt.nref*n3c) ak3=ak3-akmax
	do k2=1,nref2*n2c
	  k23=k2-1+(k3-1)*nref2*n2c
	  ak2=(k2-1)*dk2
	  if (k2.gt.nref*n2c) ak2=ak2-akmax
	do k1=1,n12
c  Do k1=n12+1 separately below.
	  k=k1+k23*n12
	  ak1=(k1-1)*dk1
	  ak=sqrt(ak1*ak1+ak2*ak2+ak3*ak3)
	  dp=npow*ak/akmaxf
	  if (dp.ge.npow) then
	    tf=0.0
	  else
	    jp=int(dp)
	    dp=dp-jp
	    tf=(1.0-dp)*tsav(jp)+dp*tsav(jp+1)
	  end if
c  Transfer function.
c N.B. If fft this, then normalize by /(2*ngr2).
	  transfc(k)=tf
C	TRANSFc(k)=TRANSFc(k)/(2*ngr2)
c  Close loop over k1.
	end do
c  Now do k1=n12+1.
	k1=n12+1
	k=ngr2+1+k23
	ak1=-0.5*akmax
	ak=sqrt(ak1*ak1+ak2*ak2+ak3*ak3)
	dp=npow*ak/akmaxf
	if (dp.ge.npow) then
	  tf=0.0
	else
	  jp=int(dp)
	  dp=dp-jp
	  tf=(1.0-dp)*tsav(jp)+dp*tsav(jp+1)
	end if
c  Transfer function.
c N.B. If fft this, then normalize by /(2*ngr2).
	transfc(k)=tf
C	TRANSFc(k)=TRANSFc(k)/(2*ngr2)
c  Close loops over k2, k3.
	end do
	end do

c  Fourier transform trans1c to spatial domain (testing only).
c  If do this, be sure to normalize as above.
C	call fft3rinv(transfc,nref2*n1c,nref2*n2c,nref2*n3c)

	return
	end
