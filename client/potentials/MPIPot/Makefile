OBJECTS = MPIPot.o
CXXFLAGS += -I../.. -Wall -O0 -g
LIB := libMPIPot.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
