OBJECTS = AMS.o
CXXFLAGS += -I../.. -Wall
LIB := libAMS.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
