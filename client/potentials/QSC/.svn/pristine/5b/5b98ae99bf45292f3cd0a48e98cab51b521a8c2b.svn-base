OBJECTS = QSC.o
CXXFLAGS += -I../.. -Wall -O2
LIB := libQSC.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
QSC.o: Parameters.h
clean :
	rm -f $(OBJECTS) $(LIB)

.PHONY : clean
