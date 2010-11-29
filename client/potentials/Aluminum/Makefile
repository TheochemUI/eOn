OBJECTS = dfrhoDblexp.o  embedenergy.o  fofrhoDblexp.o  forces.o  gagafeDblexp.o  potinit.o  sumembforce.o Aluminum.o

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libAL.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
