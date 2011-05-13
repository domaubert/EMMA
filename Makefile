##########################################
ARCH = GPU
C_LIBS = -lm 
C_FLAGS = #-O1
C_OBJS= quartz.o hilbert.o vector.o io.o cic.o oct.o particle.o tools.o amr.o segment.o communication.o
DEFINES  = -DNEWASSIGN -DPART2 -DALLCELL  -DNEWCIC -DNEWFORCE  -DNEWJACK -DNEWJACK2  -DWMPI #-DNEWPART #-DWMPI # -DEGYCSV O2 #-DNEWCIC #-O2 -DWMPI
ifeq ($(ARCH),GPU)
DEFINESGLOB= $(DEFINES) -DWGPU
EXECUTABLE = quartz
CUDA_OBJS= interface.o vector_gpu.o # cic_gpu.o 
CUDA_LIBS = -lcuda -L/usr/local/cuda/lib64 -lcudart -lcudpp_x86_64 -I/usr/lib/openmpi/include -L/usr/lib/openmpi/lib/ -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -lm -ldl
else
DEFINESGLOB= $(DEFINES) 
EXECUTABLE = quartzcpu
CUDA_OBJS= 
CUDA_LIBS =
endif
NVCC= nvcc #--device-emulation
CC = mpicc
OBJECTS = $(C_OBJS) $(CUDA_OBJS)
.c.o:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -c $<
ifeq ($(ARCH),GPU)
%.o:%.cu
	$(NVCC) $(DEFINESGLOB) $(C_LIBS) $(CUDA_LIBS) -c $*.cu	
endif

all:$(C_OBJS) $(CUDA_OBJS)
	$(CC)  $(CUDA_OBJS) $(C_OBJS) $(C_LIBS) $(CUDA_LIBS) -o $(EXECUTABLE)

clean:
	rm -f *.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) *~
