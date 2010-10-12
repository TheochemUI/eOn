OBJECTS = bopfox.o
CXXFLAGS += -I../..
LIB := libbopfox.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
