OBJECTS = GPRPotential.o
CXXFLAGS += -I../..
LIB := libGPRPot.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
