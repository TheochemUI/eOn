OBJECTS = EAM.o
CXXFLAGS += -I../..
LIB := libEAM.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
