##########################################
C_LIBS = -lm 
CUDA_LIBS =
C_OBJS= quartz.o hilbert.o vector.o io.o cic.o oct.o particle.o tools.o amr.o segment.o communication.o
CUDA_OBJS= 
DEFINES  =-g  -DNEWASSIGN -DTESTCOSMO -DALLCELL  -DNEWCIC -DNEWFORCE  -DNEWJACK #-DWMPI # -DEGYCSV O2 #-DNEWCIC #-O2 -DWMPI
NVCC= nvcc #--device-emulation
CC = gcc
OBJECTS = $(C_OBJS) $(CUDA_OBJS)
EXECUTABLE = quartz
.c.o:
	$(CC) $(DEFINES) $(C_LIBS) -c $<
%.o:%.cu
	$(NVCC) $(DEFINES) $(C_LIBS) $(CUDA_LIBS) -c $*.cu	

all:$(C_OBJS) $(CUDA_OBJS)
	$(CC)  $(CUDA_OBJS) $(C_OBJS) $(C_LIBS) $(CUDA_LIBS) -o $(EXECUTABLE)

clean:
	rm -f *.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) *~
