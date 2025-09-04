#!/usr/bin/env python

import sys
import numpy

import pathfix
import fileio as io
import atoms

p1 = io.loadcon(sys.argv[1])
p2 = io.loadcon(sys.argv[2])

distances = []

for i in range(len(p1)):
    diff = atoms.pbc(p2.r[i] - p1.r[i], p1.box)
    distances.append((i, diff[0], diff[1], diff[2], numpy.linalg.norm(diff)))

distances.sort(key=lambda d: d[4])

print("\n%5s    %10s    %10s    %10s    %10s" % ("atom", "x", "y", "z", "norm"))
print("-------------------------------------------------------------")

for d in distances:
    print("%5d    %10f    %10f    %10f    %10f" % (d[0], d[1], d[2], d[3], d[4]))

print("\n total norm:", numpy.linalg.norm(atoms.pbc(p2.r - p1.r, p1.box)), "\n")
