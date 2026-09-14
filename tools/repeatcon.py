#!/usr/bin/env python

import sys
import numpy

import pathfix
import io, atoms

confile = sys.argv[1]
xscale = int(sys.argv[2])
yscale = int(sys.argv[3])
zscale = int(sys.argv[4])
try:
    outcon = sys.argv[5]
except:
    outcon = "repeat-" + confile

p0 = io.loadcon(sys.argv[1])

t0 = atoms.Atoms(0)
for x in range(xscale):
    for a in range(len(p0)):
        t0.append(p0.r[a] + p0.box[0] * x, p0.free[a], p0.names[a], p0.mass[a])
t1 = atoms.Atoms(0)
for y in range(yscale):
    for a in range(len(t0)):
        t1.append(t0.r[a] + p0.box[1] * y, t0.free[a], t0.names[a], t0.mass[a])
t2 = atoms.Atoms(0)
for z in range(zscale):
    for a in range(len(t1)):
        t2.append(t1.r[a] + p0.box[2] * z, t1.free[a], t1.names[a], t1.mass[a])

t2.box[0] = p0.box[0] * xscale
t2.box[1] = p0.box[1] * yscale
t2.box[2] = p0.box[2] * zscale

io.savecon(outcon, t2)
