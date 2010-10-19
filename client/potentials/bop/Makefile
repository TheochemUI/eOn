
OBJECTS = bop.o 

#FC ?= gfortran
#FFLAGS ?= -O2
AR = ar cru

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libbop.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)
