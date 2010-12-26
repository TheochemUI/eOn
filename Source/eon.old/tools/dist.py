#!/usr/bin/env python

import sys
import numpy

import pathfix
import io, atoms

p1 = io.loadcon(sys.argv[1])
p2 = io.loadcon(sys.argv[2])

print     "\n%5s    %10s    %10s    %10s    %10s" % ("atom", "x", "y", "z", "norm")
print "-------------------------------------------------------------"

for i in range(len(p1)):
    diff = atoms.pbc(p2.r[i] - p1.r[i], p1.box)
    print "%5d    %10f    %10f    %10f    %10f" % (i, diff[0], diff[1], diff[2], numpy.linalg.norm(diff))

print "\n total norm:", numpy.linalg.norm(atoms.pbc(p2.r - p1.r, p1.box)), "\n"
