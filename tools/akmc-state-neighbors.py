#!/usr/bin/env python

import os
import sys
import numpy

import pathfix
import fileio as io
import atoms

stpath = os.path.join('states', sys.argv[1])

lines = open(os.path.join(stpath, 'info'), 'r').readlines()

reacenergy = 0.0

for line in lines:
    if "reactant energy" in line:
        reacenergy = float(line.split()[-1])

reac = io.loadcon(os.path.join(stpath, 'reactant.con'))

print("%12s %12s %12s" % ('pid', 'max dist', 'e diff'))

lines = open(os.path.join(stpath, 'processtable'), 'r').readlines()
for line in lines[1:]:
    split = line.strip().split()
    pid = split[0]
    pu = float(split[4])
    pr = io.loadcon(os.path.join(stpath, 'procdata', 'product_%s.con' % pid))
    maxd = 0.0
    for i in range(len(reac)):
        d = numpy.linalg.norm(atoms.pbc(reac.r[i] - pr.r[i], reac.box))
        maxd = max(maxd, d)
    print("%12s %12.6f %12.6f" % (pid, maxd, pu - reacenergy))
