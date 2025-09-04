OBJECTS = NewPot.o
CXXFLAGS += -I../..
LIB := libnewpot.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
