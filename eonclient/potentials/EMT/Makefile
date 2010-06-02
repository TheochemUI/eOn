OBJECTS = EffectiveMediumTheory.o Asap/Atoms.o Asap/NeighborList.o Asap/CellList.o Asap/EMT.o Asap/SuperCell.o Asap/EMTDefaultParameterProvider.o Asap/Vec.o Asap/Exception.o Asap/GhostAtoms.o 

CXXFLAGS += -O3 -I../..
CFLAGS += -O3 -I../..
LIB := libEMT.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
