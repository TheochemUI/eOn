#!/usr/bin/env python

import sys

import pathfix
import io, atoms

p1 = io.loadcon(sys.argv[1])
p2 = io.loadcon(sys.argv[2])
print max(atoms.per_atom_norm(p1.r - p2.r, p1.box))

