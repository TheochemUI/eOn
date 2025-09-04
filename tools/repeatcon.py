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
        t0.r = numpy.append(t0.r, [p0.r[a] + p0.box[0] * x], 0)
        t0.free = numpy.append(t0.free, p0.free[a])
        t0.names.append(p0.names[a])
        t0.mass = numpy.append(t0.mass, p0.mass[a])
t1 = atoms.Atoms(0)
for y in range(yscale):
    for a in range(len(t0)):
        t1.r = numpy.append(t1.r, [t0.r[a] + p0.box[1] * y], 0)
        t1.free = numpy.append(t1.free, t0.free[a])
        t1.names.append(t0.names[a])
        t1.mass = numpy.append(t1.mass, t0.mass[a])
t2 = atoms.Atoms(0)
for z in range(zscale):
    for a in range(len(t1)):
        t2.r = numpy.append(t2.r, [t1.r[a] + p0.box[2] * z], 0)
        t2.free = numpy.append(t2.free, t1.free[a])
        t2.names.append(t1.names[a])
        t2.mass = numpy.append(t2.mass, t1.mass[a])

t2.box[0] = p0.box[0] * xscale
t2.box[1] = p0.box[1] * yscale
t2.box[2] = p0.box[2] * zscale

io.savecon(outcon, t2)
