OBJECTS = Tersoff.o tersoffFortran.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libTersoff.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
