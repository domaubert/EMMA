#! /bin/csh -f
# Script for testing grafic1 and grafic2.
# Warning: the grafic2 run requires about 250 MB of scratch disk storage.
#
cd test_results
cp ../linger.lcdm .
#
# Testing grafic1 using linger.lcdm precomputed by lingers.
# N.B. grafic1.inc is set for 64*64*32 grid.
# It is shorter in the third direction simply to cut down on the
# time required to complete this test.
# First generate level_0 output.
#
../grafic1 >& test.out <<!
1		# use linger.dat file for transfer functions
linger.lcdm
1		# spectral index n
-0.9		# normalization (sigma_8=-0.9)
0 0		# don't output power.dat
0.5		# Grid spacing (Mpc)
2		# Ultimate refinement factor (used only to set astart)
0		# further refinement (1 or 0 for yes or no)
64 64 32	# next grid size (used only if further refinement)
32 32 16	# next grid offsets (used only if further refinement)
64 64 32	# final subvolume size (used only if further refinement)
32 32 16	# final subvolume offsets (used only if further refinement)
1		# irand =1 to generate new nos, =2 use old ones
314159265	# iseed
wnsub.lev0
!
time >>& test.out
# Make an image of density slice through middle of box.
../ic2gif ic_deltab deltab_lev0.gif >>& test.out <<!
32
0.4
!
# Save output in level_0 subdirectory
mkdir level_0
\mv -f ic_* level_0
#
# Now run grafic1 again to produce intermediate level_1 output.
#
../grafic1 >>& test.out <<!
1		# use linger.dat file for transfer functions
linger.lcdm
1		# spectral index n
-0.9		# normalization (sigma_8=-0.9)
0 0		# don't output power.dat
0.5		# Grid spacing (Mpc)
2		# Ultimate refinement factor (used only to set astart)
1		# further refinement (1 or 0 for yes or no)
64 64 32	# next grid size (used only if further refinement)
32 32 16	# next grid offsets (used only if further refinement)
64 64 32	# final subvolume size (used only if further refinement)
32 32 16	# final subvolume offsets (used only if further refinement)
2		# irand =1 to generate new nos, =2 use old ones
314159265	# iseed
wnsub.lev0
!
time >>& test.out
\mv -f grafic2.top grafic2.lev0
#
# Now run grafic2 to produce level_1 output.
# grafic2.inc is set up with nref=2 and 64**3 subgrid
#   (i.e. 32**3 "next grid size" from input to grafic1 above)
#
../grafic2 >>& test.out <<!
grafic2.lev0
1		# use linger.dat file for transfer functions
linger.lcdm
1		# spectral index n
-0.9		# normalization (sigma_8=-0.9)
0 0		# don't output power.dat
0		# further refinement (1 or 0 for yes or no)
64 64 32	# next grid size (used only if further refinement)
32 32 16	# next grid offsets (used only if further refinement)
1		# irand =1 to generate new nos, =2 use old ones
134519265	# iseed
wnsub.lev1
!
time >>& test.out
# Make an image of density slice through middle of box.
# Note 2x magnification relative to deltab_lev0.gif.
../ic2gif ic_deltab deltab_lev1.gif >>& test.out <<!
32
0.4
!
# Save output in level_1 subdirectory
mkdir level_1
\mv -f ic_* level_1
