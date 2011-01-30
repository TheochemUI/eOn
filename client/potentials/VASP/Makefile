OBJECTS = VASP.o
CXXFLAGS += -I../.. -Wall
LIB := libVASP.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
