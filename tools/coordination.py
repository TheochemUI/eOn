#!/usr/bin/env python

import sys

import pathfix
import io, atoms

p = io.loadcon(sys.argv[1])

print atoms.coordination_numbers(p, 3.3)
