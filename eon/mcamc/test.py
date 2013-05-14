#!/usr/bin/env python
import os
import numpy
import ctypes
from time import time
from mcamc import mcamc


def random_chain(t,r,p):
    Q = numpy.random.random((t,t))
    R = p*numpy.random.random((t,r))
    c = numpy.ones(t)
    return Q,R,c

#def auto_precision(Q,R):
#    N = len(Q)
#    for i in range(Q
#    Q = numpy.divide(Q,numpy.sum(Q,1)+numpy.sum(R,1))
#    print Q
#
#    cond_number = 0.0
#    for r in range(N):
#        idx = numpy.argsort(Q[r])
#        k = 1.0
#        for i in range(len(Q[r])):
#            k -= Q[r,i]
#        #cond_number = max(cond_number, -numpy.log10(abs(k)))
#        print k
#
#    if cond_number < 6:
#        prec = 'f'
#    elif cond_number < 12:
#        prec = 'd'
#    elif cond_number < 24:
#        prec = 'dd'
#    elif cond_number < 48:
#        prec = 'qd'
#    else:
#        raise RuntimeWarning("Problems is too ill-conditioned")
#    print 'chose precision: %s' % prec

def main():

    Ntrans = 200
    Nabs = 50
    print 'PRECISION TESTING'
    print '-----------------'
    print 'Solving MCAMC problem with %i transient states and %i absorbing states' % (Ntrans, Nabs)

    print '%5s %5s %5s %5s %5s' % ('p', 'f', 'd', 'dd', 'qd')
    for p in [1e-5,1e-10, 1e-15, 1e-25, 1e-30, 1e-40, 1e-60, 1e-80]:
        errors = []
        for prec in ['f','d','dd','qd']:
            Q,R,c = random_chain(Ntrans,Nabs,p)
            t, B = mcamc(Q,R,c,prec)
            errors.append(abs(1.0-max(abs(numpy.sum(B,1)))))
        print '%5.0e %5.0e %5.0e %5.0e %5.0e' % (p,errors[0],errors[1],errors[2],errors[3])

    print '\nPERFORMANCE TESTING'
    print '-------------------'

    print '%8s %8s %8s %8s %8s' % ('N', 'f', 'd', 'dd', 'qd')
    for Ntrans in [10,100,200,500,800,1000,1200,1500]:
        times = []
        for prec in ['f','d','dd','qd']:
            p = 1e-5
            Nabs = 50
            Q,R,c = random_chain(Ntrans,Nabs,p)
            t1 = time()
            t, B = mcamc(Q,R,c,prec)
            t2 = time()
            times.append(t2-t1)
        print '%8i %8.3f %8.3f %8.3f %8.3f' % (Ntrans,times[0],times[1],times[2],times[3])

main()
