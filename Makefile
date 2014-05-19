########################################## 
ARCH = CPU
C_LIBS =   -O2 # -lstdc++ -g
C_FLAGS =
C_OBJS= quartz.o hilbert.o io.o cic.o oct.o particle.o tools.o amr.o segment.o communication.o hydro_utils.o friedmann.o advanceamr.o poisson_utils.o rad_utils.o chem_utils.o src_utils.o stars.o
DEFINES  =  
#=========================================== CODE PARAMETERS =====================

#------------ MAIN OPTIONS --------------------
DEFINES  +=  -DPIC
DEFINES  +=  -DWHYDRO2
DEFINES  +=  -DWGRAV 
DEFINES  +=  -DWRAD
DEFINES  +=  -DWRADHYD
DEFINES  +=  -DTESTCOSMO
#DEFINES  +=  -DWDBG
DEFINES  +=  -DSTARS

#------------ MPI OPTIONS ---------------------

DEFINES  +=  -DWMPI
#DEFINES  +=  -DFLOORDT

#------------ CUDA OPTIONS ---------------------
#DEFINES  +=  -DWCUDA_ERR
#DEFINES  +=  -DNOCOMP

#------------ ICs OPTIONS ---------------------

DEFINES  +=  -DGRAFIC
#DEFINES  +=  -DZELDOVICH
#DEFINES  +=  -DEVRARD
#DEFINES  +=  -DEDBERT
#DEFINES  +=  -DTUBE
#DEFINES  +=  -DPARTN
#DEFINES  +=  -DPART2
#DEFINES  +=  -DWRADTEST  
#DEFINES  +=  -DTESTCLUMP # RADTEST MUST BE SET

#------------ PIC OPTIONS ----------------------

#DEFINES += -DPART_EGY
#DEFINES += -DPERFECT

#------------ GRAV OPTIONS ----------------------
#DEFINES  +=  -DFASTGRAV 
DEFINES  += -DONFLYRED

# ----------- HYDRODYNAMICS OPTIONS ------------
DEFINES  +=  -DRIEMANN_HLLC
#DEFINES  +=  -DRIEMANN_EXACT
DEFINES  +=  -DPRIMITIVE
DEFINES  +=  -DDUAL_E


# ----------- RADIATION OPTIONS ------------
DEFINES  += -DWCHEM 
DEFINES  += -DS_100000
DEFINES  += -DCOOLING
#DEFINES  += -DUVBKG
DEFINES  += -DSEMI_IMPLICIT 
DEFINES  += -DACCEL_RAD_STAR

# ---- BOUNDARY CONDITIONS (PERIODIC BY DEFAULT)--
#DEFINES  +=  -DTRANSZM
#DEFINES  +=  -DTRANSZP
#DEFINES  +=  -DTRANSYM
#DEFINES  +=  -DTRANSYP
#DEFINES  +=  -DTRANSXM
#DEFINES  +=  -DTRANSXP
#DEFINES  +=  -DREFXM # TRANS must be turned on too
#DEFINES  +=  -DREFYM # TRANS must be turned on too
#DEFINES  +=  -DREFZM # TRANS must be turned on too




#=================================================================================

ifeq ($(ARCH),GPU)
DEFINESGLOB= $(DEFINES) -DGPUAXL
EXECUTABLE = quartzgpu
CUDA_OBJS= interface.o poisson_utils_gpu.o hydro_utils_gpu.o rad_utils_gpu.o chem_utils_gpu.o # cic_gpu.o 
CUDA_LIBS =  -I/workdir/observatoire/aubert/cudpp_src_2.0/include -L/workdir/observatoire/aubert/cudpp_src_2.0/lib -L/usr/local/cuda5/lib64 -lcudart  -I/usr/local/cuda5/include -I/usr/local/openmpi-1.6.i13/include -L/usr/local/openmpi-1.6.i13/lib/ -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -ldl -lcudpp #-lcuda
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
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/oct2grid utils/oct2grid.c utils/silo/lib/libsilo.a	
oct2cell:
	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/oct2cell utils/oct2cell.c utils/silo/lib/libsilo.a	

clean:
	rm -f *.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) *~
