OBJECTS = lj_wales.o
CXXFLAGS += -I../..
LIB := liblj_wales.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
