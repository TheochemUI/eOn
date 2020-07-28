OBJECTS = PyAMFF.o
CXXFLAGS += -I../.. -Wall
LIB := libPyAMFF.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
