# PRECISION should be DOUB or SING
PRECISION = DOUB

# Add -DADD0US or -DADD2US to CFLAGS if fortran objects have 0 or 2 underscores added
# Default is one underscore.
CC = mpicc
CFLAGS =

FC = mpif90
FFLAGS = -cpp  -I/usr/local/include -D$(PRECISION) #-fpe0 #-O3 -ipo  -axP #-static

ifeq ($(PRECISION),DOUB)
	LIBS = -L/usr/local/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
	LIBS = -L/usr/local/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif


OBJS = grafic_types.o parallel_io.o grafic_io.o transform.o \
      paste.o random.o util.o time.o mpnorm.o mpic4.o mppower.o


all:	mpgrafic


mpgrafic: $(OBJS) mpgrafic.o grafic1.inc
	$(FC) $(FFLAGS) -o $@ mpgrafic.o $(OBJS) $(LIBS)

clean:
	rm -rf $(OBJS) *.mod mpgrafic.o mpgrafic

#.f90.o:
#	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.f90 grafic1.inc
	$(FC) $(FFLAGS) -c $< -o $@

