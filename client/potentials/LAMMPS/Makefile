
OBJECTS = LAMMPS_EON.o
#CXX =		g++
CXXFLAGS += -I../.. -Isrc #-Isrc/STUBS
#LINK =		g++

LIB := liblammps.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
