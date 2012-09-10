OBJECTS = TerminalPotential.o
CXXFLAGS += -I../..
LIB := libterminalpotential.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
