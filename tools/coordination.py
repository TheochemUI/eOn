#!/usr/bin/env python

import sys

import pathfix
import atoms
import fileio as io

p = io.loadcon(sys.argv[1])

print atoms.coordination_numbers(p, 3.3)
