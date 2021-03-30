OBJECTS = ExtPot.o
CXXFLAGS += -I../..
LIB := libextpot.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
