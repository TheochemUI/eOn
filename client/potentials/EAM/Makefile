OBJECTS = EAM.o
CXXFLAGS += -I../.. -Wall
LIB := libEAM.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
