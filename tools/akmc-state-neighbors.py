#!/usr/bin/env python

import sys
from pathlib import Path

import numpy

import pathfix
import fileio as io
import atoms

stpath = Path("states") / sys.argv[1]

lines = (stpath / "info").read_text().splitlines()

reacenergy = 0.0

for line in lines:
    if "reactant energy" in line:
        reacenergy = float(line.split()[-1])

reac = io.loadcon(str(stpath / "reactant.con"))

print("%12s %12s %12s" % ('pid', 'max dist', 'e diff'))

lines = (stpath / "processtable").read_text().splitlines()
for line in lines[1:]:
    split = line.strip().split()
    pid = split[0]
    pu = float(split[4])
    pr = io.loadcon(str(stpath / "procdata" / ("product_%s.con" % pid)))
    maxd = 0.0
    for i in range(len(reac)):
        d = numpy.linalg.norm(atoms.pbc(reac.r[i] - pr.r[i], reac.box))
        maxd = max(maxd, d)
    print("%12s %12.6f %12.6f" % (pid, maxd, pu - reacenergy))
