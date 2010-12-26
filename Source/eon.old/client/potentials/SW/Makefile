
OBJECTS = SW.o SWFortran.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libSW.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
