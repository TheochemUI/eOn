OBJECTS = NewPotential.o
CXXFLAGS += -I../..
LIB := libnewpotential.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
