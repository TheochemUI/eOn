CXX=g++
CXXFLAGS=-Wall -Wfatal-errors -O3
CPPFLAGS=-I../../client -I.
OBJS=mastereqn.o mpreal.o
LDFLAGS=-L/usr/local/lib -lgmp -lmpfr
mastereqn: $(OBJS)
	$(CXX) $(LDFLAGS) -o mastereqn $(OBJS)

clean:
	rm mastereqn *.o
