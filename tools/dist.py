#!/usr/bin/env python

import sys

import pathfix
import io, atoms

p1 = io.loadcon(sys.argv[1])
p2 = io.loadcon(sys.argv[2])

diff = atoms.per_atom_norm(p2.r-p1.r, p1.box)

if len(sys.argv) > 3:
    cutoff = float(sys.argv[3])
    for r in diff:
        if r >= cutoff:
            print r
else:
    print diff
