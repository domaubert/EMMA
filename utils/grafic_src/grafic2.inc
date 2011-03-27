c  grafic2.inc
c  User must set parameters for each level of multiscale refinement,
c  then make grafic2
c
c  Parameters for offsetting baryon and CDM velocity fields by
c  0.5 grid spacing.
c  The origin of the baryon density field is (0,0,0).  The origin
c  of the baryon velocity field will be offvelb*(0.5,0.5,0.5)*dx
c  and similarly for offvelc.  offvelb/c should be +/-1 or 0.
	integer offvelb,offvelc
	parameter (offvelb=0,offvelc=+1)
c
c  Parameters for refined mesh.
c
	integer n1c,n2c,n3c,nref,nref2
	parameter (n1c=64,n2c=64,n3c=32,nref=2, nref2=2*nref)
c
c n[1-3]c:	size of subgrid to sample at higher density.
c		  must be even numbers, and must be no larger than
c                 half of the upper grid size.  They need not be equal.
c                 They must equal the values of (m1s,m2s,m3s) of the
c                 subgrid extracted in the previous level of refinement.
c nref:		refinement factor for subgrid.  Must be > 1.
c  The memory required by grafic2 is 2*n1c*n2c*n3c*nref2**3 words plus a few.
