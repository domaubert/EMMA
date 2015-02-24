##########################################
ARCH = CPU
C_LIBS =  -O2 -lm -Wimplicit #-O3 -ftree-vectorize -ffast-math -fno-cx-limited-range  #-fopenmp # -lstdc++ -g
C_FLAGS =
C_OBJS= quartz.o hilbert.o io.o cic.o oct.o particle.o tools.o amr.o segment.o communication.o hydro_utils.o friedmann.o advanceamr.o poisson_utils.o rad_utils.o chem_utils.o src_utils.o stars.o zoom.o supernovae.o movie.o
DEFINES  =

#=========================================== CODE PARAMETERS =====================
include param.mk

ifeq ($(ARCH),GPU)
DEFINESGLOB= $(DEFINES) -DGPUAXL
EXECUTABLE = quartzgpu
CUDA_OBJS= interface.o poisson_utils_gpu.o hydro_utils_gpu.o rad_utils_gpu.o chem_utils_gpu.o # cic_gpu.o
CUDA_LIBS =  -I/workdir/observatoire/aubert/cudpp_src_2.0/include -L/workdir/observatoire/aubert/cudpp_src_2.0/lib -L/usr/local/cuda-5.0/lib64 -lcudart  -I/usr/local/cuda-5.0/include -I/usr/lib/openmpi/include -L/usr/lib/openmpi/lib/ -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -ldl -lcudpp #-lcuda
else
DEFINESGLOB= $(DEFINES)
EXECUTABLE = quartzcpu
CUDA_OBJS=
CUDA_LIBS =
endif

NVCC= nvcc -lstdc++ --ptxas-options=-v #-g -G #--device-emulation
CC = mpicc
OBJECTS = $(C_OBJS) $(CUDA_OBJS)
.c.o:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(CUDA_LIBS) $(C_FLAGS) -c $<
ifeq ($(ARCH),GPU)
%.o:%.cu
	$(NVCC) $(DEFINESGLOB) $(C_LIBS) $(CUDA_LIBS) -arch=sm_35 -c $*.cu
endif

all:$(C_OBJS) $(CUDA_OBJS)
	$(CC)  $(C_OBJS)  $(CUDA_OBJS) $(C_LIBS) $(CUDA_LIBS) -o $(EXECUTABLE)

oct2grid:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) utils/oct2grid.c -o utils/oct2grid -lm
alloct:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/alloct utils/alloct.c -lm
cube2silo:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/cube2silo utils/cube2silo.c utils/libsilo.a -lm

clean:
	rm -f *.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) *~
