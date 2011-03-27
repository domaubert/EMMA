c  Compute transfer functions for density and velocity in real space
c  assuming spherical symmetry.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine transfers(astart,pk,dxf,hanning,idim,transf)
c  Generate the density transfer function by one-dimensional FFT
c  quadrature assuming spherical symmetry.
c  The transfer function is stored in transf; the octant of size
c  (n1c*nref,n2c*nref,n3c*nref) is doubled in each dimension to handle
c  non-periodic BCs.
c  Input: astart = expansion factor
c  pk(ak,astart) = power spectrum function for wavenumber ak
c	     dxf = fine grid spacing (Mpc)
c        hanning = true or false for use of hanning filter on scale dxf
c           idim = component (0,1,2,3) to be returned (idim=0 for density,
c  idim=1,2,3 for three cartesian components of displacement)
c  Output: transf (transfer function, real space)
c  
	implicit none
	include 'grafic2.inc'
c
	integer nsav
	parameter (nsav=32768)
c
	real astart,pk,dxf,transf(n1c*nref2,n2c*nref2,n3c*nref2)
	logical hanning
	integer idim
c
	real drsav,fact,fr,akmax,volumek,akmaxh,dksav,d,trf
	real ak,aak,x1,x2,x3,r
	integer j,j1,j2,j3,j1o,j2o,j3o,ja,jb
	double precision tdelta(0:nsav/2),dtdelta(0:nsav/2)
	double precision tvr(0:nsav/2),dtvr(0:nsav/2),trf0
	common /transfar/ tdelta,dtdelta,tvr,dtvr
	save /transfar/
	complex trfc(0:nsav-1)
	double precision twopi
	parameter (twopi=6.283185307179586d0)
c  Increase the radius of the hanning filter in k-space by a factor fact
c  to account for the extra power with no hanning filter.
c  The results are insensitive to this for factors between 2 and 5.
c  fr must be at least equal to fact.
	parameter (fact=4.0,fr=4)
	external pk
c
	drsav=dxf/fr
	if (tdelta(0).eq.0) then
c  Initialize density.
	  akmax=twopi/dxf
	  volumek=akmax**3
	  akmaxh=0.5*akmax
	  if (.not.hanning) akmaxh=fact*akmaxh
	  dksav=akmax*fr/nsav
	  trf0=0.0
c  Precompute spherical transfer function tables for interpolation.
	  do j=0,nsav-1
	    ak=j*dksav
	    if (j.ge.nsav/2) ak=ak-akmax*fr
	    aak=ak
	    if (ak.lt.0.0) aak=-ak
	    if (aak.lt.akmaxh) then
	      trf=dksav*ak*sqrt(pk(aak,astart)*volumek)*cos(0.25*twopi*ak/akmaxh)
	    else
	      trf=0.0
	    end if
	    trfc(j)=cmplx(trf,0.0)
	    trf0=trf0+ak*trf
	  end do
	  call FFT1(trfc,nsav,1)
	  tdelta(0)=trf0*twopi/akmax**3
	  tvr(0)=0.0
	  do j=1,nsav/2
	    r=j*drsav
	    tdelta(j)=aimag(trfc(j))*twopi/(r*akmax**3)
	    tvr(j)=-r*r*drsav*tdelta(j)
	  end do
	  call splini
	  call splder(tdelta,dtdelta,nsav/2+1)
	  call splder(tvr,dtvr,nsav/2+1)
	  call spintn(tvr,dtvr,tvr,nsav/2+1)
	  tvr(0)=0.0
	  do j=1,nsav/2
	    r=j*drsav
	    tvr(j)=tvr(j)/(r*r)
	  end do
	end if
c  Initialize transfer function to zero.
	do j3=1,n3c*nref2
	do j2=1,n2c*nref2
	do j1=1,n1c*nref2
	  transf(j1,j2,j3)=0.0
	end do
	end do
	end do
c  Fill transfer function.
	do j3=1,n3c*nref+1
	  x3=(j3-1)*dxf
	  j3o=n3c*nref2-j3+2
	do j2=1,n2c*nref+1
	  x2=(j2-1)*dxf
	  j2o=n2c*nref2-j2+2
	do j1=1,n1c*nref+1
	  x1=(j1-1)*dxf
	  j1o=n1c*nref2-j1+2
	    r=sqrt(x1*x1+x2*x2+x3*x3)
c  Spline interpolation.
	    ja=r/drsav
	    jb=ja+1
	    if (jb.ge.nsav/2) then
	      print*,'nsav too small in transfd! nsav, jb=',nsav,jb
	      stop
	    end if
	    d=r/drsav-ja
	    if (idim.eq.0) then
C	      trf=(1.0-d)*tdelta(ja)+d*tdelta(jb)
	      trf=tdelta(ja)+d*(dtdelta(ja)+d*(3.0d0*(tdelta(jb)-
     2 tdelta(ja))-2.0d0*dtdelta(ja)-dtdelta(jb)+d*(dtdelta(ja)+
     3 dtdelta(jb)+2.0d0*(tdelta(ja)-tdelta(jb)))))
	    else
C	      trf=(1.0-d)*tvr(ja)+d*tvr(jb)
	      trf=tvr(ja)+d*(dtvr(ja)+d*(3.0d0*(tvr(jb)-
     2 tvr(ja))-2.0d0*dtvr(ja)-dtvr(jb)+d*(dtvr(ja)+
     3 dtvr(jb)+2.0d0*(tvr(ja)-tvr(jb)))))
c  Multiply the radial component of displacement by
c  x(idim)/r to get appropriate component.
	      if (r.ne.0.0) then
	        if (idim.eq.1) then
	          trf=trf*x1/r
	        else if (idim.eq.2) then
	          trf=trf*x2/r
	        else if (idim.eq.3) then
	          trf=trf*x3/r
	        end if
	      else
	        trf=0.0
	      end if
	    end if
c  (j1,j2,j3) label position of one octant in transf.
c  (j1o,j2o,j3o) label other octants; must exclude j?o=n?c*nref2 (j?=0).
c  Fill octants of transf by reflection.
	    transf(j1,j2,j3)=trf
	    if (j1.ne.1) transf(j1o,j2,j3)=trf
	    if (j2.ne.1) transf(j1,j2o,j3)=trf
	    if (j1.ne.1.and.j2.ne.1) transf(j1o,j2o,j3)=trf
	    if (j3.ne.1) then
	      transf(j1,j2,j3o)=trf
	      if (j1.ne.1) transf(j1o,j2,j3o)=trf
	      if (j2.ne.1) transf(j1,j2o,j3o)=trf
	      if (j1.ne.1.and.j2.ne.1) transf(j1o,j2o,j3o)=trf
	    end if
c  N.B. displacement window function is odd along dimension idim!
	    if (idim.eq.1.and.j1.ne.1) then
	      transf(j1o,j2,j3)=-transf(j1o,j2,j3)
	      if (j2.ne.1) transf(j1o,j2o,j3)=-transf(j1o,j2o,j3)
	      if (j3.ne.1) transf(j1o,j2,j3o)=-transf(j1o,j2,j3o)
	    else if (idim.eq.2.and.j2.ne.1) then
	      transf(j1,j2o,j3)=-transf(j1,j2o,j3)
	      if (j1.ne.1) transf(j1o,j2o,j3)=-transf(j1o,j2o,j3)
	      if (j3.ne.1) transf(j1,j2o,j3o)=-transf(j1,j2o,j3o)
	    else if (idim.eq.3.and.j3.ne.1) then
	      transf(j1,j2,j3o)=-transf(j1,j2,j3o)
	      if (j1.ne.1) transf(j1o,j2,j3o)=-transf(j1o,j2,j3o)
	      if (j2.ne.1) transf(j1,j2o,j3o)=-transf(j1,j2o,j3o)
	    end if
	    if (idim.ne.0.and.j1.ne.1.and.j2.ne.1.and.j3.ne.1)
     &        transf(j1o,j2o,j3o)=-transf(j1o,j2o,j3o)
	end do
	end do
	end do

	return
	end
