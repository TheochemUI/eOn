OBJECTS = AMS_IO.o
CXXFLAGS += -I../.. -Wall
LIB := libAMS_IO.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
