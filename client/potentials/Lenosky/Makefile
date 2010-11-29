
OBJECTS = Lenosky.o lenoskyFortran.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libLenosky.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
