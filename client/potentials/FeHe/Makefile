OBJECTS = feforce.o FeHe.o

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libFeHe.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
