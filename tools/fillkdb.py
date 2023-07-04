#!/usr/bin/env python

# Assumes kdbinsert.py is in the path, that the user is inside an eon akmc
# directory, and that the states directory is named 'states'.

import os
import sys

args = ' '.join(sys.argv[1:])

n = 0
while os.path.exists(os.path.join('states', '%d' % n)):
	lines = open(os.path.join('states', '%d' % n, 'processtable'), 'r').readlines()[1:]
	procs = [line.split()[0] for line in lines]
	for proc in procs:
		path = os.path.join('states', '%d' % n, 'procdata')
		reactant_path = os.path.join('states', '%d' % n, 'procdata', 'reactant_%s.con' % proc)
		saddle_path = os.path.join('states', '%d' % n, 'procdata', 'saddle_%s.con' % proc)
		product_path = os.path.join('states', '%d' % n, 'procdata', 'product_%s.con' % proc)
		mode_path = os.path.join('states', '%d' % n, 'procdata', 'mode_%s.dat' % proc)
		os.system('kdbinsert.py %s %s %s -o %s %s' % (reactant_path, saddle_path, product_path, mode_path, args))
	n += 1
