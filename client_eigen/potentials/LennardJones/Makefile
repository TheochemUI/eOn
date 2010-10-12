OBJECTS = LJ.o
CXXFLAGS += -I../..
LIB := libLJ.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
