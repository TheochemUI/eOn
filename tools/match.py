#!/usr/bin/env python
import sys

import pathfix
import atoms
import fileio as io

p1 = io.loadcon(sys.argv[1])
p2 = io.loadcon(sys.argv[2])
print atoms.match(p1,p2,float(sys.argv[3]),3.3)
