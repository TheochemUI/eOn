#!/usr/bin/env python
import numpy
from os.path import isdir, join
from os import listdir
from os.path import basename
from sys import argv, exit

if len(argv) < 2:
    print '%s: state [state ...]' % basename(argv[0])
    print 'must be run in top level of the simulation directory'
    exit(2)

numpy.set_printoptions(precision=5, suppress=True, linewidth=200, threshold=10000)

def read_processtable(path):
    f = open(path)
    f.readline()
    products = []
    rates = []
    for line in f:
        line = line.strip()
        fields = line.split()
        product = int(fields[3])
        if product == -1: continue
        rate = float(fields[7])
        products.append(product)
        rates.append(rate)
    f.close()
    return products, rates


states = [ int(d) for d in listdir('states') if isdir(join('states', d)) ]
states = sorted(states)

N = len(states)
rate_matrix = numpy.zeros((N,N))

for i in states:
    products, rates = read_processtable(join('states',str(i),'processtable'))
    for j in range(len(products)):
        if i == j: continue
        rate_matrix[i,i] += rates[j]
        rate_matrix[products[j],i] -= rates[j]



numpy.savetxt('rate_matrix', rate_matrix)

ew, ev = numpy.linalg.eig(rate_matrix)
idx = numpy.argsort(ew)
ew = ew[idx]
ew[0] = 0.0
ev = ev[:,idx]

for i in range(N):
    ev[:,i] /= sum(ev[:,i])

p0 = numpy.zeros(N)
p0[0] = 1.0
c0 = numpy.linalg.solve(ev, p0)

states = [ int(a) for a in argv[1:] ]
print 'set logscale x'
plots = []
for j in states:
    plots.append('"-" w l t "state %i"' % j)
print 'plot ' + ','.join(plots)

for j in states:
    for t in numpy.logspace(-15, 0, 100):
        print t, abs(sum([ c0[i]*numpy.exp(-ew[i]*t)*ev[:,i] for i in xrange(N) ])[j])
    if j != states[-1]:
        print 'e'
