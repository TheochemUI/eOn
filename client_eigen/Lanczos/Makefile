OBJECTS := lanczos.o lanczos_for_eon.o tools.o sort.o tridiagv.o
#CXXFLAGS += -I..
LIB := libLanczos.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)
.PHONY : clean
