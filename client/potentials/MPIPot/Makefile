OBJECTS = MPIPot.o
CXXFLAGS += -I../.. -Wall -O2
LIB := libMPIPot.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
