CXXFLAGS=-O3 -fPIC
#can be compiled with threading
#CXXFLAGS=-O3 -fPIC -fopenmp
CPPFLAGS=-Ilibqd -I../../client
LDFLAGS=libqd/libqd.a

libmcamc.so: libmcamc.cpp libqd/libqd.a
	$(CXX) -shared -o $@ libmcamc.cpp $(CXXFLAGS) $(LDFLAGS) $(CPPFLAGS)

libqd/libqd.a:
	$(MAKE) -C libqd

.PHONY: clean

clean:
	-rm -f libmcamc.so libqd/{*.a,*.o}
