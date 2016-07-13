CC= mpicc -g3
C_LIBS= -lm

EXECUTABLE= test


#------------------------------------------------------------------------------
# UNIT_C_OBJS= test_grid.o

#------------------------------------------------------------------------------
UNIT_C_OBJS= test_mpi.o
DEFINES = -DWMPI

#------------------------------------------------------------------------------
# UNIT_C_OBJS= test_grafic_ic.o
# DEFINES = -DPIC -DGRAFIC -DTESTCOSMO -DWGRAV

#------------------------------------------------------------------------------
# UNIT_C_OBJS= test_uniform_grid.o
# DEFINES = -DPIC -DWGRAV

#------------------------------------------------------------------------------
# UNIT_C_OBJS= test_star_formation.o
# DEFINES = -DPIC -DWGRAV


OBJDIR= obj
UNIT_SRCDIR= src
I_DIR= -I../src

$(OBJDIR)/%.o: $(UNIT_SRCDIR)/%.c
	$(CC) $(C_FLAGS) $(DEFINES) $(I_DIR) -c $< -o $@

UNIT_OBJ=$(patsubst %,$(OBJDIR)/%,$(UNIT_C_OBJS))


C_OBJS= \
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

SRCDIR = ../src
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(C_FLAGS) $(DEFINES) $(I_DIR) -c $< -o $@

MAIN_OBJ=$(patsubst %,$(OBJDIR)/%,$(C_OBJS))

OBJ= $(MAIN_OBJ) $(UNIT_OBJ)

all:$(OBJ)
	$(CC) $(OBJ) $(LD_DIR) $(C_LIBS) -o $(EXECUTABLE) $(LD_FLAGS)

clean:
	rm -f $(OBJDIR)/*.o $(EXECUTABLE) *~
