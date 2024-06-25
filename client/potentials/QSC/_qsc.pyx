import numpy as np
cimport numpy as np

cdef extern from "QSC.h": 
    cdef cppclass QSC:
        QSC()

        long vlist_updates
        void set_verlet_skin(double dr)
        void set_cutoff(double c)
        double get_cutoff()
        void set_qsc_parameter(int Z, double n, double m, double epsilon,
                               double c, double a)

        void force(long N, const double *R, const int *atomicNrs,
                   double *F, double *U, double *variance, const double *box) except +

class QSCMissingParameter(Exception):
    pass

cdef class PyQSC:
    cdef QSC *thisptr

    def __cinit__(self):
        self.thisptr = new QSC()

    def __dealloc__(self):
        del self.thisptr

    def set_parameter(self, Z, n, m, epsilon, c, a):
        self.thisptr.set_qsc_parameter(Z,n,m,epsilon,c,a)

    def get_vlist_updates(self):
        return self.thisptr.vlist_updates

    def set_verlet_skin(self, dr):
        self.thisptr.set_verlet_skin(dr) 

    def set_cutoff(self, c):
        self.thisptr.set_cutoff(c)

    def get_cutoff(self):
        return self.thisptr.get_cutoff()
        
    def force(self,
              N, 
              np.ndarray[np.double_t,ndim=1] R,
              np.ndarray[np.int,ndim=1] atomicNrs,
              np.ndarray[np.double_t,ndim=1] box):
        #XXX: UGLY TYPE HACK
        cdef np.ndarray[np.int32_t]  Z=atomicNrs.astype(np.int32)
        cdef np.ndarray[np.double_t] U=np.zeros((1,),dtype=np.double) 
        cdef np.ndarray[np.double_t] F=np.zeros((3*N,),dtype=np.double)
        cdef np.ndarray[np.double_t] variance=np.zeros((3*N,),dtype=np.double)
        try:
            self.thisptr.force(N,<double*>R.data,<int*>Z.data,<double*>F.data,
                               <double*>U.data,<double*>variance.data,<double*>box.data)
        except RuntimeError:
            raise QSCMissingParameter("Missing Parameter")
        return U[0],F

