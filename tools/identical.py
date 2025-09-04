#!/usr/bin/env python
import sys

import pathfix
import fileio as io
import atoms
import config
config.comp_eps_r=float(sys.argv[3])

p1 = io.loadcon(sys.argv[1])
p2 = io.loadcon(sys.argv[2])
print(atoms.identical(p1,p2))
