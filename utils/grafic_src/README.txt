README.txt for grafic2 written by E. Bertschinger, March 18, 2001
  Please send all questions, comments, bugs, etc. to edbert@mit.edu

*** This file contains important information for users of grafic2.
*** Please read and retain for reference as you run grafic2.

This package produces multiscale gaussian random fields on nested
Cartesian lattices.  There are two programs: grafic1 and grafic2.
grafic1 produces the top-level periodic grid.  grafic2 is used for
refinement.  Refinement may be recursive to give subgrids within
subgrids like Russian matrioshka dolls.

The periodic top grid is called the level-0 grid.  The first refinement
is called the level-1 grid.  The level-n grid may be refined by grafic2
to produce the level-(n+1) grid.

Unavoidable numerical errors are introduced at each stage of refinement.
The errors are greatest for the velocity field and they arise from aliasing
of the tidal fields.  For typical grids and a CDM-like spectrum the rms
errors are less than 5 percent at the level-2 grid.

For an example of grafic1 and grafic2, type 'make' and then 'make test'
Compare test_results/test.out and test_results/*.gif with the files in
tests_sample.


Procedure for a single refinement level.
----------------------------------------

1. Edit grafic1.inc to set the number of grid points in each dimension
  for the level-0 grid (np1,np2,np3).  Also set offsets for baryon and
  cdm velocities (offvelb, offvelc) and the parameter sigstart that fixes
  the initial redshift by setting the standard deviation of the density
  fluctuation at the highest level of refinement.

2. Edit the Makefile to replace the f77 compiler name if necessary
  (e.g. for Linux, F77 = g77) and to use your favorite compiler
  flags for optimization, etc.  Then type 'make grafic1'

3. Run grafic1 interactively once or twice to see what parameters must
  be entered.  The memory requirement is a little more than (np1*np2*np3)
  words.  Look at grafic1.inp to see an example of the input parameters.

4. The first parameters entered specify the transfer function calculation
  for the baryon and CDM power spectra.  (These are allowed to be different,
  although they are set equal unless one uses LINGERS.)

  Transfer functions for baryons and CDM may be precomputed by running
  LINGERS (a modification of LINGER_SYN from the COSMICS package that
  outputs multiple times, allowing interpolation of the transfer function
  to an arbitrary redshift).  Alternatively, transfer functions for CDM
  models may be given by fitting formulae (BBKS or enter your own into
  power.f), or scale-free power-laws may be used.

***BEWARE*** the BBKS transfer function significantly underestimates
     the amount of large-scale power for currently favored models.
     It assumes omegab=0.  For realistic models one should run LINGERS,
     which is provided with this package.
  Scale-free power spectra are not set up in grafic to automatically
  normalize multiscale initial conditions.  For scale-free initial conditions
  the user will have to increase the normalization k^3P(k) (a parameter
  input to grafic1) at each refinement level by a factor nref**(3+n).

5. For realistic models, grafic1 can normalize the power spectrum either
  with the CMB quadrupole using Q_rms-PS or with sigma_8.  Whatever your
  choice, grafic1 will calculate the other choice and inform you.  For
  accurate transfer functions computed with LINGERS, I recommend CMB
  normalization with Q_rms-PS = 18 micro-K for a n=1 spectrum, following
  Hinshaw et al 1997, ApJ 464, L17.  Careful simulators will also run
  CMBFAST to get the shape of the CMB power spectrum.  Beware that CMBFAST
  normalizes using the Bunn-White prescription; however, you can determine
  Q_rms-PS from C_2.

6. For multiscale initial conditions, the user will need to generate
  density and velocity fields for each of the grids and refinement levels
  that are used.  In the case of a single refinement level, there are
  two sets of data files: one for the top (level-0) grid and one for the
  level-1 refinement.  grafic1 must be run once for each --- once to produce
  output for the level-0 grid and once to produce input to grafic2 which
  will produce the level-1 output.  When running grafic1 (and similarly
  grafic2), the user must inform the program whether to produce output
  fields at the current refinement level or to produce input for further
  refinement.  grafic1 asks the user to specify which case is being run.
  Note that output fields are smoothed with a Hanning filter (spherical
  cosine window function that goes to zero at the one-dimensional Nyquist
  frequency).

7. If further refinement will be performed, grafic1 prompts the user
  to enter the size of the subgrid (in units of the top grid spacing dx)
  to be extracted for refinement, and grid offsets to give the position
  of its corner.  These are all integers.  The subgrid must be no larger
  than half the top grid size (np1,np2,np3).  The offset may be positive
  or negative, from -0.5*(np1,np2,np3) to +0.5*(np1,np2,np3).  Periodic
  boundary conditions apply on the top grid.

8. In order to reduce numerical aliasing errors with multiple refinement
  levels, grafic1 then prompts the user to enter the final subvolume size
  and offsets (integers, with units given by the top grid spacing dx).
  For refinement to level 1, these are the same as the "next grid"
  values (see grafic1.inp).

9. The last step for grafic1 is to enter information for random numbers.
  It is strongly recommended that the user save the white noise files
  created and used by grafic1 (and grafic2).  For multiple refinement,
  the same random numbers must always be used for all refinements across
  a given level of the hierarchy; saving the files and reusing them
  can ensure this.  Also beware that the random number generator
  will not necessarily produce the same results (for a given seed)
  with 32-bit and 64-bit processors.  If the user wishes, s/he may
  provide the random number routine (by replacing randa with a uniform
  (0,1) pseudorandom generator) or the white noise files.  Note that
  the random fields are sampled in real space, not Fourier space.

10. grafic1 should be run once to produce final output files at the
  top grid level.  These will be called ic_deltab, ic_velbx, etc.
  They should be put into a separate directory labeled "level_0" to
  avoid confusion with the files that will be produced for refined
  grids by grafic2.  These output files provide the information
  needed for comsological simulation codes.  They are written as follows:
	open(11,file='ic_deltab',form='unformatted')
	rewind 11
	write(11) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
	do i3=1,np3
	  write(11) ((deltab(i1,i2,i3),i1=1,np1),i2=1,np2)
	end do
	close(11)
  The velocity files are written the same way.  The parameters are
  np[1-3] = grid size, dx = grid spacing in comoving Mpc,
  x[1-3]o = grid offsets in comoving Mpc (=0 unless offvelb/c is set nonzero
  in grafic1.inc), astart is the starting expansion factor (with a=1 today),
  omegam,omegav,h0 are self-explanatory (H0 has units of km/s/Mpc).
  Note that deltab is the dimensionless density fluctuation in real space.
  The velocity (for ic_velbx, etc.) is the peculiar velocity in units of
  proper km/s at astart.  grafic1 and grafic2 no longer output particle
  positions like their predecessor grafic does.  See the end of this file
  for information on how to compute them from the velocities.

11. grafic1 should be run a second time to produce intermediate output
  for the level-1 grid.  This time the user must specify the subvolume
  that will be extracted (size and offset).  Both of these quantities are
  given in terms of the top grid spacing.  Thus, if the top grid is a
  120 Mpc cube np1=np2=np3=240, and one wishes a 40 Mpc cubical subvolume
  with 4 times higher resolution centered at (0,30,60) Mpc, then the
  subgrid size is n1c=n2c=n3c=80 and the offsets are m1off=-40, m2off=20,
  m3off=80.  The refinement factor nref=4 isn't needed by grafic1.
  Save it for grafic2.  Note that the same random numbers must be used
  for both the first and second runs of grafic1.  They determine the
  structure in the box.

12. The second run of grafic1 will produce an output file grafic2.top.
  This should be moved to a new directory called "level_1".  It is the
  input data file needed by grafic2.

13. Standard output of grafic1 contains useful information about the
  random numbers and the standard deviation of each field it computes.
  When run in the "intermediate output" mode to produce grafic2.top,
  the velocity fields are split into "inner" and "outer" parts due to
  fluctuations inside and outside of the subvolume.  The user may
  generally ignore the standard output (except when warnings and errors
  are written), although the sophisticated user will find useful
  statistics.

14. To generate the level-1 subgrid, edit grafic2.inc to include the
  subgrid size (n1c,n2c,n3c) (exactly as entered into grafic1) and
  the refinement factor nref, i.e. the factor by which the grid spacing
  dx will be decreased.  The memory requirement for grafic2 is a little
  more than 16*n1c*n2c*n3c*nref**3.  Then make grafic2

15. Run grafic2 in the level_1 directory.  Its input is similar to grafic1,
  except that one must specify the grafic2.dat input filename from grafic1.
  Normally this would be grafic2.top unless one changed the name.  grafic2
  doesn't need information about the final subvolume since these were already
  given to grafic1 and saved in the grafic2.top file.  To produce output
  density and velocity fields at level 1, choose the correct mode.
  grafic2 has 2 modes (final output, intermediate output) just like grafic1.
  Finally, choose a different random number seed for grafic2.  These random
  numbers are used for the level-1 subgrid.  grafic2 will produce output
  files ic_deltab, etc exactly like grafic1.  Keep these in the level_1
  directory so that they won't be confused with the level_0 files.
  (The headers contain the grid spacing dx for each level, allowing them
  to be distinguished in case of confusion.)

16. grafic2 takes a large amount of cpu time; it performs 25-39 FFTs of
  size 8*n1c*n2c*n3c*nref**3.  For a size of 512**3, on the Origin 2000
  the runtime is about 6 hours for the final output mode (25 FFTs).

17. grafic2 produces several temporary files (temp.*) that may be removed
  after each grafic2 run.  A large amount of scratch disk space is needed
  for grafic2.  Including the input and output files, the scratch disk
  requirement can be as large as 52*n1c*n2c*n3c*nref**3 4-byte words.


Procedure for multiple refinement levels.
-----------------------------------------

1. If one wishes to produce initial conditions for more than two levels
  of refinement, grafic2 may be run recursively.   Producing output for
  refinement level n requires running grafic1 once and grafic2 n times.
  Thus, producing output at all levels (0, 1, ..., n) of an n-level
  multiscale simulation requires n+1 runs of grafic1 and n(n+1)/2 runs of
  grafic2.  The user will have to keep track of the files at each level
  and provide the correct subgrid information to grafic1 and 2 for each run.
  It might be worthwhile writing a perl script to automate this.

2. For example, suppose that one has already done 2-level (level_0
  and level_1) as described above, and now wishes to further refine
  to produce a level_2 subgrid.  The first step is to set the size and
  offset of this subgrid relative to the level_0 grid.  Subgrids of all
  refinement levels must be aligned on level_0 cells; all subvolumes
  of all refinement levels must have physical sizes that are a whole
  number of top grid spacings across.  Once the user has set the size
  and location of the level_2 subgrid, the first step is to run grafic1.
  grafic1 needs to know both the size and location of the level_1 subgrid
  and the level_2 subgrid.  (To produce output for level_n, grafic1 needs
  to know the level_n subgrid as the "final" subgrid.)  Note that the
  same random numbers must be used with grafic1 as were used in the
  earlier grafic runs for earlier levels of the grid hierarchy.

3. After running grafic1, the user edits grafic2.inc and sets
  (n1c,n2c,n3c) equal to the level_1 subgrid size that was entered
  into grafic1.  The user also sets the refinement factor for level 1
  and makes grafic2.  grafic2 is then run with the output of grafic1
  in the "intermediate output" mode, which requires specifying the size
  and offsets of the level_2 subgrid that will be extracted.  Be careful
  here.  The size is in units of the level_1 grid, and is therefore
  larger by a factor of the level_1 refinement factor than the level_2
  grid size entered into grafic1.  Similarly, the offsets given to
  grafic2 are relative to the origin of the level_1 grid and are in units
  of the level_1 grid spacing.  The same random numbers must be used for
  this run as were used in the previous grafic2 run at level_1.  grafic2
  will produce grafic2.dat.

4. The user then edits grafic2.inc again to enter the level_2 subgrid size
  that was entered into grafic2 in step 3 as well as the refinement factor
  for level_2, i.e. the factor by which level_1 is refined to yield level_2.
  grafic2 is run in final output mode, and different random numbers should
  be used than were used at level_0 and level_1.  (Just name the random number
  files wnsub_level0, wnsub_level1, etc., and keep careful track of the
  refinement level, and you should have no problem.)  The result will be
  ic_* for level_2.

5. To generate output for higher levels, one repeats this pattern, extending
  the chain of runs of grafic2 by one for each additional refinement level.


Computing particle positions from velocities.
---------------------------------------------

Particle positions follow readily from velocities, as illustrated by
the following f77 code fragment:

c  (np1,np2,np3) are the size of the level_n subgrid.
	real velx(np1,np2,np3),x(np1,np2,np3)
c
	print*,'Enter top grid spacing in Mpc'
	read(*,*) dx0
c  Change the filename for the appropriate component;
c  ic_velbx will give the x-coordinate of baryon particles.
	open(11,file='ic_velbx',form='unformatted')
	rewind 11
	read(11) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
c  velocity (proper km/s) =  Displacement (comoving Mpc at astart) * vfact.
c  vfact = dln(D+)/dtau where tau=conformal time.
c  These functions (fomega, dladt) are in time.f, so link them in.
	vfact=fomega(astart,omegam,omegav)
     &       *h0*dladt(astart,omegam,omegav)/astart
c  Offset the positions so that the (dx0/dx)**3 particles corresponding
c  to one level_0 particle are centered on the level_0 particle.
	offset=0.5*(dx0-dx)
c  Read each slice of the velocity field.
	do i3=1,np3
	  read(11) ((velx(i1,i2,i3),i1=1,np1),i2=1,np2)
c  Unperturbed grid position.
	  z0=(i3-1)*dx+offset
	  do i2=1,np2
	    y0=(i2-1)*dx+offset
	    do i1=1,np1
	      x0=(i1-1)*dx+offset
c  For the y-component, replace x0 with y0.
	      x(i1,i2,i3)=x0+velx(i1,i2,i3)/vfact
	    end do
	  end do
	end do
	close(11)



Additional Notes
----------------

1. The initial density fluctuations will be smaller at the coarser grid
  levels.  The user's multiscale evolution code may have to compute
  accurate gravitational forces when the density is almost uniform over
  the top grid yet highly deformed at the highest level of refinement.
  The user may wish to apply the Zel'dovich approximation at the top
  refinement level while using a Poisson solver only for the refinements,
  until the density fluctuations become sufficiently large.  You should
  check how well your code reproduces the Zel'dovich approximation at
  early times, using the accurate velocity and displacement fields computed
  by grafic1 and grafic2.

2. Multiple subvolumes may be used at any given refinement level.
  For example, if the top grid spacing is 1 Mpc, then non-overlapping
  subgrids of spacing 0.25 Mpc may be put at several places in the box.
  However, the small-scale fields within one will be ignored when
  computing the tides felt by the other.  Thus, the subgrids should be
  widely separated.

3. The refinement factor may be any whole number greater than 1.  Different
  refinement factors may be used between different levels.

4. A utility ic2gif is provided to produce images of the ic_* files.
  This can be helpful in checking.  It is created by typing 'make ic2gif'
  Then run it as follows: 'ic2gif ic_filename gif_filename.gif'
  where ic_filename is the data file to be images and gif_filename.gif
  is the output gif file.  ic2gif will image one slice.
