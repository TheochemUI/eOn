
OBJECTS = Tip4p_H.o potH_H2O.o

# Note that the optimization flag is overruled  with 0
# Seems to be necessary for the interface between fortran and gcc
%.o: %.f90
	$(FC) $(FFLAGS) -O0 -c $< -o $@

LIB := libtip4p_h.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
