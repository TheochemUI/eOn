CXX=g++
CXXFLAGS=-Wall -Wfatal-errors -O3
CPPFLAGS=-I../../client -I.
OBJS=mastereqn.o mpreal.o
LDFLAGS=-lgmp -lmpfr
mastereqn: $(OBJS)
	$(CXX)  -o mastereqn $(OBJS) $(LDFLAGS)

clean:
	rm mastereqn *.o
