OBJECTS = LAMMPS_EON.o
MPI_OBJECTS = mpi.o
#CXX =		g++
CXXFLAGS += -I../..
#LINK =		g++

all: liblammps.a libfakempi.a

LIB := liblammps.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@

MPI_LIB := libfakempi.a
$(MPI_LIB) : $(MPI_OBJECTS)
	ar cru $@ $^
	ranlib $@

clean :
	rm -f $(OBJECTS) $(LIB) $(MPI_OBJECTS) $(MPI_LIB)

.PHONY : clean
