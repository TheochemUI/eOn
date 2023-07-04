CXX=g++
CXXFLAGS=-I. -O3 -fPIC
OBJS= bits.o c_qd.o dd_real.o qd_const.o util.o c_dd.o dd_const.o \
fpu.o qd_real.o

libqd.a: $(OBJS)
	ar rcs libqd.a $(OBJS)

.PHONY: clean

clean:
	-rm -f *.o libqd.a
