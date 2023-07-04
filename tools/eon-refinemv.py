#!/usr/bin/env python

import sys

from eon import fileio as io
from eon import atoms
import numpy

def abort():
    print("\n   usage: eon-refinemv.pl [input movie file] [ouput movie file] <intermediate frames>")
    sys.exit()

try:
    p = io.loadposcars(sys.argv[1])
    if len(p) == 0:
        p = io.loadcons(sys.argv[1])
    if len(p) == 0:
        abort()
except:
    abort()

try:
    outfile = sys.argv[2]
except:
    abort()

try:
    N = int(sys.argv[3])
except:
    N = 1

f = open(outfile, 'w')
f.close()

for i in range(len(p) - 1):
    v = atoms.pbc(p[i+1].r - p[i].r, p[i].box)
    d = numpy.linalg.norm(v)
    v /= d
    io.savecon(outfile, p[i], 'a')
    for j in range(N):
        temp = p[i].copy()
        temp.r += j * (v * (d/(N-1)))
        io.savecon(outfile, temp, 'a')
io.savecon(outfile, p[len(p) - 1], 'a')
