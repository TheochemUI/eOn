OBJECTS = LJ.o
CXXFLAGS += -I../..
LIB := libLJ.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
