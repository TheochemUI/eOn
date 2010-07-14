OBJECTS = EAM.o Parameters.o
CXXFLAGS += -I../.. -Wall -pedantic -Wextra
LIB := libEAM.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
