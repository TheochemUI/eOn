OBJECTS = QSC.o
CXXFLAGS += -I../.. -Wall -pedantic -O2
LIB := libQSC.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
test : CXXFLAGS += -DQSC_STANDALONE
test : $(OBJECTS) test.o
	g++ -llua -o test $^

.PHONY : clean
