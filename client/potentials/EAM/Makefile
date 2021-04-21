OBJECTS = EAM.o
CXXFLAGS += -I../.. -Wall -pedantic -Wextra -O0 -g
LIB := libEAM.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
test : CXXFLAGS += -DEAM_STANDALONE
test : $(OBJECTS) test.o
	g++ -o test $^

.PHONY : clean
