import os
import numpy
import ctypes

solver = ctypes.CDLL(os.path.abspath('solver.so'))

Q = numpy.array([ [0.0, 0.5, 0.0],
                  [0.5, 0.0, 0.5],
                  [0.0, 0.5, 0.0] ])

R = numpy.array([ [0.5, 0.0],
                  [0.0, 0.0],
                  [0.0, 0.5] ])

c = numpy.array([1.0, 2.0, 0.5])

Qflat = list(Q.ravel())
Qflat = (ctypes.c_double * len(Qflat))(*Qflat)
Rflat = list(R.ravel())
Rflat = (ctypes.c_double * len(Rflat))(*Rflat)
cflat = list(c)
cflat = (ctypes.c_double * len(cflat))(*cflat)

Bflat = (ctypes.c_double * len(Rflat))()
tflat = (ctypes.c_double * len(cflat))()

# void solve(int Qsize, double *Qflat, int Rcols, double *Rflat, double *c_in, double *B, double *t)

solver.solve(Q.shape[0], Qflat, R.shape[1], Rflat, cflat, Bflat, tflat)

print "B:\n", numpy.array(list(Bflat)).reshape(R.shape)
print "t:\n", list(tflat)





