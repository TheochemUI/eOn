CXX=mpic++
OBJECTS = GPAW.o
CXXFLAGS += -I../.. -Wall
LIB := libGPAW.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
