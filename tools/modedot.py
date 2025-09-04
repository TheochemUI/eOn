#!/usr/bin/env python

import sys
import numpy

import pathfix
import fileio as io

m1 = io.load_mode(sys.argv[1])
m2 = io.load_mode(sys.argv[2])

print((m1 * m2).sum())
