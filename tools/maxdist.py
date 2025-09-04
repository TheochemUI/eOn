#!/usr/bin/env python

import sys

import pathfix
import atoms
import fileio as io

p1 = io.loadcon(sys.argv[1])
for file2 in sys.argv[2:]:
    p2 = io.loadcon(file2)
    distances = atoms.per_atom_norm(p1.r - p2.r, p1.box)

    max_i=0
    max_d=0.0
    for i in range(len(distances)):
        if distances[i]>max_d:
            max_d = distances[i]
            max_i = i

    print("%s: max distance: %f atom index: %i" % (file2, max_d, max_i))
