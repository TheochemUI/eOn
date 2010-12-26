#!/usr/bin/env python
import sys

import pathfix
import io, atoms

p1 = io.loadcon(sys.argv[1])
p2 = io.loadcon(sys.argv[2])
print atoms.identical(p1,p2,float(sys.argv[3]))
