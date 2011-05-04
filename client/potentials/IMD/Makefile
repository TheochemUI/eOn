OBJECTS = IMD.o
CXXFLAGS += -I../.. -Wall
LIB := libIMD.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
