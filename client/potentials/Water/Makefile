OBJECTS := ccl.o potential_base.o spce_ccl.o tip4p_ccl.o Water.o
CXXFLAGS += -I..
LIB := libwater.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
.PHONY : clean
