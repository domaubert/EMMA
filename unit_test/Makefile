CC= mpicc
C_LIBS= -lm
C_FLAGS= -g3

EXECUTABLE= test
OBJDIR= obj
UNIT_SRCDIR= src
SRCDIR = ../src

# include test.mk
# include mk/test_grid.mk
# include mk/test_grafic_ic.mk
include mk/test_grid_mpi.mk
# include mk/test_mpi.mk
# include mk/test_star_formation

$(OBJDIR)/%.o: $(UNIT_SRCDIR)/%.c
	$(CC) $(C_FLAGS) $(DEFINES) -I$(SRCDIR) -c $< -o $@
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

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(C_FLAGS) $(DEFINES) -I$(SRCDIR) -c $< -o $@
MAIN_OBJ=$(patsubst %,$(OBJDIR)/%,$(C_OBJS))

all: $(MAIN_OBJ) $(UNIT_OBJ)
	$(CC) $(MAIN_OBJ) $(UNIT_OBJ) $(LD_DIR) $(C_LIBS) -o $(EXECUTABLE) $(LD_FLAGS)


clean:
	rm -f $(OBJDIR)/*.o $(EXECUTABLE) *~
