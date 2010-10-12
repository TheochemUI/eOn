
OBJECTS = SW.o SWFortran.o

#FC ?= gfortran
#FFLAGS ?= -O2
#AR ?= ar cru

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libSW.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)
