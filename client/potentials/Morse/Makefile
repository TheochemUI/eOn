OBJECTS = Morse.o
CXXFLAGS += -I../..
LIB := libMorse.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
