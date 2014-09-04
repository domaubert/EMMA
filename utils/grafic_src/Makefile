
CC = gcc
F77 = ifort

#DEBUG = -W -Wall -O2 -fbounds-check
#DEBUG = -DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON:verbose_runtime=ON -g
DEBUG = -O2 -cpp

.c.o:
	$(CC) $(DEBUG) -c $<

.f.o:
	$(F77) $(DEBUG) -c $<

G2OBJ = grafic2.o wnsub.o filta_minimal.o filta_sharpk.o transf_spherical.o \
transfd_minimal.o power.o time.o util.o random8.o fft3r.o

G1OBJ = grafic1.o ic4.o power.o time.o util.o random8.o fft3r.o

LOBJ = lingers.o dverk.o

IOBJ = ic2gif.o compress.o

all:	grafic1 grafic2 lingers ic2gif

grafic1: $(G1OBJ)
	$(F77) $(DEBUG) $(G1OBJ) -o grafic1

grafic2: $(G2OBJ)
	$(F77) $(DEBUG) $(G2OBJ) -o grafic2

lingers: $(LOBJ)
	$(F77) $(DEBUG) $(LOBJ) -o lingers

ic2gif: $(IOBJ)
	$(CC) $(DEBUG) $(IOBJ) -o ic2gif

test:
	@if [ ! -d ./test_results ]; then ( mkdir -p ./test_results ) ; fi;
	@csh test.csh
	@echo "test results are in ./test_results"                             

clean:
	rm *.o grafic1 grafic2 lingers ic2gif temp.*

grafic1.o:		grafic1.inc
ic4.o:			grafic1.inc
grafic2.o:		grafic2.inc
transf_spherical.o:	grafic2.inc
transfd_minimal.o:	grafic2.inc
filta_minimal.o:	grafic2.inc
filta_sharpk.o:		grafic2.inc
wnsub.o:		grafic2.inc
ic2gif.o:		cmap.h
