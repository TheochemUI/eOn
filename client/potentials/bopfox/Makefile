OBJECTS = bopfox.o
CXXFLAGS += -I../..
LIB := libbopfox.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
