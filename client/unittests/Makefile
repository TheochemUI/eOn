#utFiles
UTEXEC := Test.out
UTSRC  := $(wildcard ./*.cpp)
UTOBJ  := $(patsubst %.cpp, %.o, $(UTSRC))
DOTD   := $(patsubst %.cpp, %.d, $(UTSRC))

#options
INCLUDE := -I./..

#rules
$(UTEXEC): $(UTOBJ)
	@echo "$(UTSRC)"
	$(CXX) $(INCLUDE) -o $@ $^ $(DEPOBJECTS) $(DEPLIBS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(INCLUDE) $(DEPFLAGS) $(CXXFLAGS) -c $<

$(UTOBJ): $(UTSRC)

#targets
.PHONY: clean clobber
clobber:
	rm -r */ 
clean:
	rm -f $(UTEXEC) $(UTOBJ) $(DOTD)  *~ .*~
