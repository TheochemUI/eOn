#!/usr/bin/env python

import sys

import pathfix
import atoms
import fileio as io

p = io.loadcon(sys.argv[1])

if len(sys.argv) == 3:
    cut = float(sys.argv[2])
else:
    cut = 3.3

print(atoms.coordination_numbers(p, cut))
