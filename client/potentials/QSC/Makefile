OBJECTS = QSC.o
CXXFLAGS += -I../.. -Wall -O2
LIB := libQSC.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
QSC.o: Parameters.h
clean :
	rm -f $(OBJECTS) $(LIB)
test : CXXFLAGS += -DQSC_STANDALONE
test : $(OBJECTS) test.o
	g++ -llua -o test $^

.PHONY : clean
