c  Call randini(iseed) with a large integer seed (smaller than 2**32-5).
c  Then call randa(x) repeatedly to get uniform deviates x, 0.le.x.le.1.
c  Call randg(x) for standard normal deviates and randp(alam,n) for Poisson
c  deviates (with Poisson parameter alam).
c  Written by E. Bertschinger, 1992.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine randini(iseed)
c  randini is used to initialize the random number generators.  iseed is a
c  positive integer.  The basic random number generator is a
c  subtract-with-borrow lagged Fibonacci generator with base b=2**32-5 and
c  period b**r-b**s=10**414.
c  See Marsaglia and Zaman, Ann. Appl. Prob. 1, 462 (1991).
c  This generator is shuffled with an independent one to break up serial
c  correlations.
	integer iseed,r,s,c,nshuf,n,irptr,i,j,nmr,nms,inorm
	parameter (r=43,s=22,nshuf=157)
	real randar(nshuf),x,xnorm
	double precision xt(r),xn,b,binv
	common /randnos/ xt,randar,irptr,n,c
	common /randnor/ xnorm,inorm
	save /randnos/,/randnor/
	parameter (b=4294967291.0d0,binv=1.0d0/b)
c
	xn=iseed
	xn=dmod(xn,b)
	if (xn.lt.0) xn=xn+b
	if (xn.lt.0.or.xn.ge.b) then
	  write(*,*) 'Error with seed in randini'
	  stop
	end if
	xt(1)=xn
c  Initialize xt using multiplicative congruential generator.
	  do 10 i=2,r
	  xn=dmod(16807.0d0*xn,b)
	  xt(i)=xn
10	continue
	n=r
	c=0
c  Warm up generator.
	  do 30 j=1,5
	    do 20 i=1,nshuf
	    n=n+1
	    if (n.gt.r) n=1
	    nmr=n-r
	    if (nmr.lt.1) nmr=nmr+r
	    nms=n-s
	    if (nms.lt.1) nms=nms+r
	    xn=xt(nms)-xt(nmr)-c
	    if (xn.ge.0) then
	      c=0
	    else
	      c=1
	      xn=xn+b
	    end if
	    xt(n)=xn
c  Fill shuffling table.
	    randar(i)=xn*binv
20	  continue
30	continue
	irptr=1
c  Initialize shuffling generator.
	call randin1(iseed)
c  Run a few times to shuffle.
	  do 40 i=1,5*nshuf
	  call randa(x)
40	continue
	xnorm=0.0
	inorm=0
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine randa(x)
	integer r,s,c,nshuf,n,irptr,nmr,nms
	parameter (r=43,s=22,nshuf=157)
	real randar(nshuf),x,x1
	double precision xt(r),xn,b,binv
	common /randnos/ xt,randar,irptr,n,c
	save /randnos/
	parameter (b=4294967291.0d0,binv=1.0d0/b)
c
c  Extract random number from shuffling table.
10	x=randar(irptr)
c  Generate a new random number.
	n=n+1
	if (n.gt.r) n=1
	nmr=n-r
	if (nmr.lt.1) nmr=nmr+r
	nms=n-s
	if (nms.lt.1) nms=nms+r
	xn=xt(nms)-xt(nmr)-c
	if (xn.ge.0) then
	  c=0
	else
	  c=1
	  xn=xn+b
	end if
	xt(n)=xn
c  Refill shuffling table.
	randar(irptr)=xn*binv
c  Use independent generator to get pointer into shuffling table.
	call randa1(x1)
	irptr=int(nshuf*x1)+1
	irptr=min(irptr,nshuf)
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine randg(x)
c  randg generates standard normal deviates.  It uses randa, which must be
c  initialized by a call to randini.
c
	real x,xr,xnorm,twopi
	parameter (twopi=6.28318531)
	common /randnor/ xnorm,inorm
	save /randnor/
c
	if (inorm.eq.1) then
	  x=xnorm
	  inorm=0
	else
10	  call randa(x)
	  if (x.le.0.0) go to 10
	  xr=sqrt(-2.0*log(x))
	  call randa(x)
	  xnorm=xr*sin(twopi*x)
	  x=xr*cos(twopi*x)
	  inorm=1
	end if
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine randp(alam,np)
c  randp generates a sample np from a Poisson distribution with mean alam.
c  It uses randa, which must be initialized by a call to randini.
c
	real alam,x,t,tt
c
	if (alam.le.0.0) then
	  np=-1
	  return
	end if
c
	if (alam.gt.50.0) then
c  Use asymptotic Gaussian distribution for alam > 50.
	  call randg(x)
	  np=alam+sqrt(alam)*x
	  return
	end if
c
	np=0
	t=1.0
	tt=exp(-alam)
10	  call randa(x)
	  t=t*x
	  if (t.lt.tt) go to 20
	  np=np+1
	  go to 10
20	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine randin1(iseed)
c  randini is used to initialize the random number generators.  iseed is a
c  positive integer less than 1.0e9.  The basic random number generator is
c  taken from Press et al., Numerical Recipes, p. 199 and is based on
c  Knuth's suggestion for a portable random number generator.  mseed is
c  any large number less than m=1.0e9.
c
	parameter (m=1000000000,mseed=161803398,rm=1.0/m,nrand=157)
c  The dimension 55 is special and should not be modified; see Knuth.
	integer ma(55)
	real randar(nrand)
	common /randno1/ randar,irptr,ma,inext,inextp
	save /randno1/
c
	isee=mod(iseed,m)
c  Initialize ma(55).
	mj=mseed-isee
	mj=mod(mj,m)
	if (mj.lt.0) mj=mj+m
	ma(55)=mj
	mk=1
c  Now initialize the rest of the table, in a slightly random order,
c  with numbers that are not especially random.
	  do 10 i=1,54
	  ii=mod(21*i,55)
	  ma(ii)=mk
	  mk=mj-mk
	  if (mk.lt.0) mk=mk+m
	  mj=ma(ii)
10	continue
c  Randomize them by "warming up the generator."
	  do 30 k=1,4
	    do 20 i=1,55
	    ma(i)=ma(i)-ma(1+mod(i+30,55))
	    if (ma(i).lt.0) ma(i)=ma(i)+m
20	  continue
30	continue
	inext=0
	inextp=31
c  Exercise generator before storing in shuffling table.
	  do 40 i=1,nrand
	  inext=inext+1
	  if (inext.eq.56) inext=1
	  inextp=inextp+1
	  if (inextp.eq.56) inextp=1
	  mj=ma(inext)-ma(inextp)
	  if (mj.lt.0) mj=mj+m
	  ma(inext)=mj
40	continue
c  Now fill shuffling table.
	  do 50 i=1,nrand
	  inext=inext+1
	  if (inext.eq.56) inext=1
	  inextp=inextp+1
	  if (inextp.eq.56) inextp=1
	  mj=ma(inext)-ma(inextp)
	  if (mj.lt.0) mj=mj+m
	  ma(inext)=mj
	  randar(i)=mj*rm
50	continue
	irptr=1
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine randa1(x)
c  randa generates uniform random numbers in the interval [0,1).
c  It must be initialized with a call to randin1.
c
	parameter (m=1000000000,mseed=161803398,rm=1.0/m,nrand=157)
c  The dimension 55 is special and should not be modified; see Knuth.
	integer ma(55)
	real randar(157),x
	common /randno1/ randar,irptr,ma,inext,inextp
	save /randno1/
c
c  Extract random number from shuffling table.
	x=randar(irptr)
c  Generate a new random number.
	inext=inext+1
	if (inext.eq.56) inext=1
	inextp=inextp+1
	if (inextp.eq.56) inextp=1
	mj=ma(inext)-ma(inextp)
	if (mj.lt.0) mj=mj+m
	ma(inext)=mj
c  Shuffle.
	randar(irptr)=mj*rm
	irptr=int(nrand*x)+1
	irptr=min(irptr,nrand)
	return
	end
