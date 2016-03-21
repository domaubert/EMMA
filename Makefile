##########################################
ARCH = CPU

C_LIBS = -O2 -Wimplicit  -g -lm  # -fopenmp #-O3 -ftree-vectorize -ffast-math -fno-cx-limited-range  #-fopenmp # -lstdc++ -g -std=c11
HDF5_LIBS = -I /usr/local/hdf5/include/ -L /usr/local/hdf5/lib/
C_FLAGS =
C_OBJS= emma.o \
				hilbert.o \
				io.o \
				cic.o \
				oct.o \
				particle.o \
				tools.o \
				amr.o \
				segment.o \
				communication.o \
				friedmann.o \
				advanceamr.o \
				hydro_utils.o \
				poisson_utils.o \
				rad_utils.o \
				chem_utils.o \
				src_utils.o \
				stars.o \
				zoom.o \
				supernovae.o \
				movie.o \
				convert.o \
				parameters.o \
				restart.o \
				ic.o\


DEFINES =
OBJDIR = obj
SRCDIR = src
OBJ=$(patsubst %,$(OBJDIR)/%,$(C_OBJS))

#=========================================== CODE PARAMETERS =====================
include param.mk

ifeq ($(ARCH),GPU)
DEFINESGLOB= $(DEFINES) -DGPUAXL
EXECUTABLE = emmagpu
CUDA_OBJS= interface.o poisson_utils_gpu.o hydro_utils_gpu.o rad_utils_gpu.o chem_utils_gpu.o # cic_gpu.o
CUDA_LIBS =  -I/workdir/observatoire/aubert/cudpp_src_2.0/include -L/workdir/observatoire/aubert/cudpp_src_2.0/lib -L/usr/local/cuda-5.0/lib64 -lcudart  -I/usr/local/cuda-5.0/include -I/usr/lib/openmpi/include -L/usr/lib/openmpi/lib/ -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -ldl -lcudpp #-lcuda
else
DEFINESGLOB= $(DEFINES)
EXECUTABLE = emmacpu
CUDA_OBJS=
CUDA_LIBS =
endif
CUDA_OBJ=$(patsubst %,$(OBJDIR)/%,$(CUDA_OBJS))


NVCC= nvcc -lstdc++ --ptxas-options=-v #-g -G #--device-emulation
CC = mpicc
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(HDF5_LIBS) $(CUDA_LIBS) $(C_FLAGS) -c $< -o $@ -lgsl -lgslcblas -lhdf5
ifeq ($(ARCH),GPU)
$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(DEFINESGLOB) $(C_LIBS) $(CUDA_LIBS) -I$(SRCDIR) -arch=sm_35 -c $< -o $@
endif

all:$(OBJ) $(CUDA_OBJ)
	$(CC) $(OBJ)  $(CUDA_OBJ) $(HDF5_LIBS) $(C_LIBS) $(CUDA_LIBS) -I$(SRCDIR) -o $(EXECUTABLE) -lgsl -lgslcblas -lhdf5

oct2grid:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) utils/oct2grid.c -o utils/oct2grid -lm
field2grid:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) utils/field2grid.c -o utils/field2grid -lm
alloct:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/alloct utils/alloct.c -lm
part2cic:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) utils/part2cic.c -o utils/part2cic -lm
cube2silo:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/cube2silo utils/cube2silo.c utils/libsilo.a -lm
part2silo:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/part2silo utils/part2silo.c utils/libsilo.a -lm


clean:
	rm -f $(OBJDIR)/*.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) *~
