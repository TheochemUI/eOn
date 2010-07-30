
OBJECTS = Lenosky.o lenoskyFortran.o

#FC ?= gfortran
#FFLAGS ?= -O2
#AR ?= ar cru

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libLenosky.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)
