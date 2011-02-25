##########################################
C_LIBS = -lm 
CUDA_LIBS =
C_OBJS= quartz.o hilbert.o vector.o io.o cic.o oct.o particle.o tools.o amr.o segment.o communication.o
CUDA_OBJS= 
DEFINES  = -g -DNBUFF=4096 -DNEWASSIGN -DTESTPLUM -DNDUMP=1 -DNSTEP=10 -DLCOARSE=5 -DLMAX=6 -DCIC2 -DNITER=4096 -DALLCELL -DSTRIDE=512 -DDT=1e-4  -DWMPI
NVCC= nvcc #--device-emulation
CC = mpicc
OBJECTS = $(C_OBJS) $(CUDA_OBJS)
EXECUTABLE = quartz
.c.o:
	$(CC) $(DEFINES) -c $<
%.o:%.cu
	$(NVCC) $(DEFINES) $(C_LIBS) $(CUDA_LIBS) -c $*.cu	

all:$(C_OBJS) $(CUDA_OBJS)
	$(CC)  $(CUDA_OBJS) $(C_OBJS) $(C_LIBS) $(CUDA_LIBS) -o $(EXECUTABLE)

clean:
	rm -f *.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) *~
