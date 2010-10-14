OBJECTS = VASP.o
CXXFLAGS += -I../.. -Wall -pedantic
LIB := libVASP.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
