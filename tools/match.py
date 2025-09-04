#!/usr/bin/env python
import sys

import pathfix
import atoms
import fileio as io

p1 = io.loadcon(sys.argv[1])
p2 = io.loadcon(sys.argv[2])
ret = atoms.get_mappings(p1, p2, float(sys.argv[3]), float(sys.argv[4]))
if ret:
    print(sys.argv[1], sys.argv[2])
