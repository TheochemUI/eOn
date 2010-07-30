OBJECTS = VASP.o
CXXFLAGS += -I../..
LIB := libVASP.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
