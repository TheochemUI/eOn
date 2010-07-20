OBJECTS = QSC.o
CXXFLAGS += -I../.. -Wall -pedantic -Wextra -g -O0
LIB := libQSC.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)
test : CXXFLAGS += -DQSC_STANDALONE
test : $(OBJECTS) test.o
	g++ -llua -o test $^

.PHONY : clean
