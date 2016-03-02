FC = mpif90  
MPIFC = mpif90  
OPTS = -O3 -g -warn all  -assume byterecl -fpp -DLONGINT
#OPTS = -O3 -g -warn all -convert big_endian -assume byterecl -fpp -DLONGINT -I/usr/local/include
#OPTS = -C -g -warn all -convert big_endian -assume byterecl -fpp -DLONGINT -I/usr/local/include
#OPTS = -O3 -g -warn all -convert big_endian -assume byterecl -fpp -I/usr/local/include

MODS = 			\
	modconst.mod	\
	modvariable.mod	\
	modtiming.mod	\
	modparam.mod	\
	modmpicom.mod	\
	modio.mod	\
	modsort.mod	\
	modfofpara.mod	

OBJS = 			\
	modconst.o	\
	modvariable.o	\
	modtiming.o	\
	modparam.o	\
	modmpicom.o	\
	modio.o 	\
	modsort.o	\
	modfofpara.o	\
	fof.o


.SUFFIXES: .f90

.f90.o: $(MODS)
	$(FC) $(OPTS) -c $<

.f90.mod:
	$(FC) $(OPTS) -c $<

all : fof

fof : $(MODS) $(OBJS)
	$(MPIFC) $(OBJS) -o fof 

clean :
	rm -f *.mod *.o *~ fof

