# Compiler
CC   = mpif90
OPTS = -g

# Project name
PROJECT = adios2_miniapp

# Libraries
LIBS = `adios2-config --fortran-libs` -L/opt/view/lib
INCS = `adios2-config --fortran-flags`

OBJS = settings.o heat_transfer.o main.o

$(PROJECT): $(OBJS)
	$(CC) $(OPTS) $^ $(INCS) -o $(PROJECT) $(LIBS) -lhdf5

%.o: %.F90
	$(CC) $(OPTS) $(INCS) -c $< $(LIBS) -o $@

clean:
	rm -f $(PROJECT) *.o *.mod
