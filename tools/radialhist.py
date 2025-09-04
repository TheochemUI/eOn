#!/usr/bin/env python

#this script outputs a gnuplot file to plot the
#radial histogram

import sys
import os
import numpy

import pathfix
import fileio as io, atoms

def usage():
    print("usage: %s poscarfile" % (os.path.basename(sys.argv[0])))
    sys.exit(1)

if len(sys.argv) < 2:
    usage()

try:
    p = io.loadcon(sys.argv[1])
except IOError as xxx_todo_changeme:
    (errno, strerrno) = xxx_todo_changeme.args
    print("%s: %s" % (sys.argv[1], strerrno))
    usage()

print("set boxwidth 0.10")
print("set xlabel 'Distance (Angstrom)'")
print("set ylabel 'Number of Atoms'")
print("set key off")
print("plot [0:8] '-' with boxes")

hist = {}
for i in range(0, len(p)):
    for j in range(i+1, len(p)):
        d = numpy.linalg.norm(p.r[i]-p.r[j])
        d = round(d, 1)
        if d not in hist:
            hist[d] = 1
        else:
            hist[d] += 1

for k,v in hist.items():
    print(k,v)
