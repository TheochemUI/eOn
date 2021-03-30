SRCDIR = /home/ryanciufo/pyamff/pyamff/fingerprints
OBJECTS = $(SRCDIR)/fptypeF.o $(SRCDIR)/neighborlistMT.o $(SRCDIR)/behlerParrinelloMT.o $(SRCDIR)/fingerprintsF.o $(SRCDIR)/pyamffCalcF.o
FC=gfortran
AR=ar cru
RANLIB=ranlib
# Note that the optimization flag is overruled  with 0
# Seems to be necessary for the interface between fortran and gcc
%.o: %.f90
	$(FC) -c $< -o $@

LIB := libAMFF.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
