
OBJECTS = bop.o 

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libbop.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
