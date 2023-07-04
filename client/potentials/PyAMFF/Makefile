#SRCDIR = /home/ryanciufo/pyamff/pyamff/fingerprints
LIBS += libAMFF.a

OBJECTS = PyAMFF.o  #$(SRCDIR)/fptypeF.o $(SRCDIR)/neighborlistMT.o $(SRCDIR)/behlerParrinelloMT.o $(SRCDIR)/io.o $(SRCDIR)/fingerprintsF.o $(SRCDIR)/pyamffCalcF.o
#CXXFLAGS += -I../..
FC=gfortran
# Note that the optimization flag is overruled  with 0
# Seems to be necessary for the interface between fortran and gcc
#%.o: %.f90
#	$(FC) $(FFLAGS) -c $< -o $@
all: libPyAMFF.a
LIB := libPyAMFF.a

$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
#	$(AR) $@ $^
#	$(RANLIB) $@


clean :
	rm -f $(OBJECTS) $(LIB)
	rm -f libAMFF.a
