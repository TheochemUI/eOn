OBJECTS = VASP.o
CXXFLAGS += -I../.. -Wall -pedantic
LIB := libVASP.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
