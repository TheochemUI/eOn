#!/usr/bin/env python

# Assumes kdbinsert.py is in the path, that the user is inside an eon akmc
# directory, and that the states directory is named 'states'.

import os
import sys
from pathlib import Path

args = ' '.join(sys.argv[1:])

n = 0
while (Path("states") / ("%d" % n)).exists():
	state = Path("states") / ("%d" % n)
	lines = (state / "processtable").read_text().splitlines()[1:]
	procs = [line.split()[0] for line in lines]
	for proc in procs:
		procdata = state / "procdata"
		reactant_path = procdata / ("reactant_%s.con" % proc)
		saddle_path = procdata / ("saddle_%s.con" % proc)
		product_path = procdata / ("product_%s.con" % proc)
		mode_path = procdata / ("mode_%s.dat" % proc)
		os.system('kdbinsert.py %s %s %s -o %s %s' % (reactant_path, saddle_path, product_path, mode_path, args))
	n += 1
