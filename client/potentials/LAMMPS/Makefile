OBJECTS = LAMMPS.o
#CXX =		g++
CXXFLAGS += -I../..
#LINK =		g++

all: liblammps.a

LIB := liblammps.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@

clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
