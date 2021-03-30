OBJECTS = LJCluster.o
CXXFLAGS += -I../..
LIB := libLJCluster.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
