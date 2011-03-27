c  Make the 4-panel velx plot showing errors due to radial
c  transfer function and minimal filter.
c
	real f1(256,256,128),f2(128,128,128)
	logical hanning
	double precision sig1,sig2,sigma
	character*80 filename
	common /big/ f1,f2
c
	print*,'Enter *1024.256 filename'
	read(*,'(a80)') filename
	open(10,file=filename,form="unformatted")
	rewind 10
	read(10) np1,np2,np3,dx,astart,omegam,omegav,h0
	print*, np1,np2,np3,dx,astart,omegam,omegav,h0
	read(10) (((f1(i,j,k),i=1,256),j=1,256),k=1,128)
	close(10)
	print*,'Enter ic_* file to compare'
	read(*,'(a80)') filename
	open(10,file=filename,form="unformatted")
	rewind 10
	read(10) np1,np2,np3,dx,x1o,x2o,x3o,m1t,m2t,m3t,
     &    m1offt,m2offt,m3offt,hanning,astart,omegam,omegav,h0
	print*, np1,np2,np3,dx,x1o,x2o,x3o,m1t,m2t,m3t,
     &    m1offt,m2offt,m3offt,hanning,astart,omegam,omegav,h0
	do k=1,np3
	  read(10) ((f2(i,j,k),i=1,np1),j=1,np2)
	end do
	close(10)
c  Compare.
	sig1=0.0
	sig2=0.0
	sigma=0.0
	do k=1,128
	do j=1,128
	do i=1,128
	  sig1=sig1+f1(i,j,k)**2
	  sig2=sig2+f2(i,j,k)**2
	  f1(i,j,k)=f1(i,j,k)-f2(i,j,k)
	  sigma=sigma+f1(i,j,k)**2
	end do
	end do
	end do
	sig1=sqrt(sig1/(128*128*128))
	sig2=sqrt(sig2/(128*128*128))
	sigma=sqrt(sigma/(128*128*128))
	print*,'File 1 rms=',sig1
	print*,'File 2 rms=',sig2
	print*,'RMS difference=',sigma
c  Output
	open(10,file="diff2.dat",form="unformatted")
	rewind 10
	write(10) np1,np2,np3,n1c,n2c,n3c,nref,dx,
     &           astart,omegam,omegav,h0,hanning,xoff,yoff,zoff
	write(10) (((f1(i,j,k),i=1,128),j=1,128),k=1,128)
	close(10)
	stop
	end
