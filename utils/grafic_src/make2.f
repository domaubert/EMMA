c  n can be between 64 and 128.
	parameter (n=64)
	real f(256,256,256),f1(128,128,128),g(2*n,2*n,2*n)
	common /biga/ f,f1,g
c  Read top grid data.
	open(10,file='wntop.0',form='unformatted')
	rewind 10
	read(10) n1,n2,n3,iseed
	print*,'n1,n2,n3,iseed=',n1,n2,n3,iseed
	do i3=1,n3
	  read(10) ((f(i1,i2,i3),i1=1,n1),i2=1,n2)
	end do
	close(10)
c  Read subgrid data.
	open(10,file='wnsub.1',form='unformatted')
	rewind 10
	read(10) n1,n2,n3,iseed
	print*,'n1,n2,n3,iseed=',n1,n2,n3,iseed
	do i3=1,n3
	  read(10) ((f1(i1,i2,i3),i1=1,n1),i2=1,n2)
	end do
	close(10)
c  Now spread portion of top array, renormalizing for refinement.
	do m3=0,n-1
	  i3=1+m3-0.5*(n-32)
	  if (i3.lt.1) i3=i3+256
	do m2=0,n-1
	  i2=1+m2-0.5*(n-32)
	  if (i2.lt.1) i2=i2+256
	do m1=0,n-1
	  i1=1+m1-0.5*(n-32)
	  if (i1.lt.1) i1=i1+256
	  do n3=1,2
	  do n2=1,2
	  do n1=1,2
	    g(n1+2*m1,n2+2*m2,n3+2*m3)=f(i1,i2,i3)/sqrt(8.0)
	  end do
	  end do
	  end do
	end do
	end do
	end do
c  Average subgrid data over cells 2x2x2 and multiply by sqrt(8)
c  to adjust noise level.
	do m3=0,63
	do m2=0,63
	do m1=0,63
	  avg=0.0
	  do n3=1,2
	  do n2=1,2
	  do n1=1,2
	    avg=avg+f1(n1+2*m1,n2+2*m2,n3+2*m3)
	  end do
	  end do
	  end do
	  avg=avg/sqrt(8.0)
c  Insert data into g.
	  g(1+m1+n-32,1+m2+n-32,1+m3+n-32)=avg
	end do
	end do
	end do
c  Write data.
	open(10,file='wnsub.dat',form='unformatted')
	rewind 10
	write(10) 2*n,2*n,2*n,0
	do i3=1,2*n
	  write(10) ((g(i1,i2,i3),i1=1,2*n),i2=1,2*n)
	end do
	close(10)
	stop
	end
