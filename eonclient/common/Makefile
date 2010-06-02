OBJECTS = ConjugateGradients.o Matter.o Quickmin.o Constants.o Parameters.o HelperFunctions.o \
    Potentials.o PotentialsInterface.o debug.o
CXXFLAGS += -I..
ifdef NO_POT_EMT
   CXXFLAGS += -DNO_POT_EMT
endif
ifdef NO_POT_ZPICE
   CXXFLAGS += -DNO_POT_ZPICE
endif
libCommon.a : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean:
	rm -f $(OBJECTS) libCommon.a

.PHONY : clean