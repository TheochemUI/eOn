OBJECTS = Morse.o
CXXFLAGS += -I../..
LIB := libMorse.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
