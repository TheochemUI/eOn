#!/usr/bin/env python
import numpy
from os.path import isdir, join
from os import listdir
from sys import argv

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


fh = open(argv[1], 'w')
fh.write(str(N)+'\n')
numpy.savetxt(fh, rate_matrix)
fh.close()
