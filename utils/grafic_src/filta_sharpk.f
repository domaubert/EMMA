cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine filtak(hanning,transfc)
c  Sharp k-space filter, with or without Hanning filter at the grid scale.
c  The filter is zero outside the primary Brillouin zone.
c  
	implicit none
	include 'grafic2.inc'
c
	integer n12,n23,ngr2
	parameter (n12=nref*n1c,n23=nref2*n2c*nref2*n3c)
	parameter (ngr2=n23*n12)
c
	logical hanning
	complex transfc(ngr2+n23)
c
	integer k,k2,k3,l1,l2,l3,k23
	real ak
	double precision twopi
	parameter (twopi=6.283185307179586d0)
c
c  Initialize transfer function to zero.
	do k=1,ngr2+n23
	  transfc(k)=0.0
	end do
c  Loop over primary Brillouin zone.
	do l3=-n3c,n3c-1
	  k3=l3+1
	  if (k3.lt.1) k3=k3+nref2*n3c
	  do l2=-n2c,n2c-1
	    k2=l2+1
	    if (k2.lt.1) k2=k2+nref2*n2c
	    k23=k2-1+(k3-1)*nref2*n2c
	    do l1=0,n1c-1
	      k=1+l1+k23*n12
	      ak=sqrt((1.*l1/n1c)**2+(1.*l2/n2c)**2+(1.*l3/n3c)**2)
	      if (hanning) then
	        if (ak.lt.nref) then
	          transfc(k)=cos(0.25*twopi*ak/nref)**2
	        else
	          transfc(k)=0.0
	        end if
	      else
	        transfc(k)=1.0
	      end if
C	TRANSFc(k)=TRANSFc(k)/(2*ngr2)
	    end do
c  Now do k1=n12+1.
	    if (nref.eq.1.and.(.not.hanning)) transfc(ngr2+1+k23)=1.0
C	TRANSFc(k)=TRANSFc(k)/(2*ngr2)
	  end do
	end do

c  Fourier transform trans1c to spatial domain (testing only).
C	call fft3rinv(transfc,nref2*n1c,nref2*n2c,nref2*n3c)

	return
	end
