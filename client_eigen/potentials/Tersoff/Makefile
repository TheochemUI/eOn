
OBJECTS = Tersoff.o tersoffFortran.o

#FC ?= gfortran
#FFLAGS ?= -O2
#AR ?= ar cru

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libTersoff.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)
