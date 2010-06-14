
OBJECTS = dfrhoDblexp.o  embedenergy.o  fofrhoDblexp.o  forces.o  gagafeDblexp.o  potinit.o  sumembforce.o Aluminum.o

FC = gfortran

%.o: %.f
	${FC} -c $< -o $@

CXXFLAGS += -O3 -I../..
CFLAGS += -O3 -I../..

LIB := libAL.a

$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@
clean :
	rm -f $(OBJECTS) $(LIB)




