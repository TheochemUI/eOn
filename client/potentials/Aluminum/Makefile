
OBJECTS = dfrhoDblexp.o  embedenergy.o  fofrhoDblexp.o  forces.o  gagafeDblexp.o  potinit.o  sumembforce.o Aluminum.o

FC ?= gfortran
FFLAGS ?= -O2
AR ?= ar cru

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libAL.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)
