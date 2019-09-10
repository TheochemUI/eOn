#!/usr/bin/env python

import sys

import pathfix
import io, atoms

for filename in sys.argv[1:]:
    p = io.loadcon(filename)

    cn=atoms.coordination_numbers(p, 3.3)

    min=99999
    for i in range(len(cn)):
        if cn[i]<min:
            min=cn[i]
            min_index=i

    print("%s: min coordination %i atom index %i" % (filename, min, min_index))
