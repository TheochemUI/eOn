#!/usr/bin/env python

import os
import sys

if not os.path.isdir('states'):
	print('Abort: this script should be executed in an akmc directory with a states subdirectory present')
	sys.exit()


print('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s' % ('state', 'good', 'NC', 'energy', 'iterations', 'nonlocal', 'barrier', 'nonneg', 'prefactor', 'unknown'))
i = 0
while os.path.isdir(os.path.join('states', str(i))):
	results = open(os.path.join('states', str(i), 'search_results.txt'), 'r').readlines()[2:]
	good = 0
	nc = 0
	energy = 0
	iterations = 0
        nonlocalt = 0
        barrier = 0
	unknown = 0
        nonneg = 0
        prefac = 0
	for r in results:
		if 'good' in r:
			good += 1
		elif 'repeat' in r:
			good += 1
		elif 'Not Connected' in r:
			nc += 1
		elif 'High Energy' in r:
			energy += 1
		elif 'Total Iterations' in r:
			iterations += 1
		elif 'reverse' in r:
			pass
                elif 'Nonlocal' in r:
                        nonlocalt += 1
                elif 'barrier >' in r:
                        barrier += 1
                elif 'Nonneg' in r:
                        nonneg += 1
                elif 'Failed Prefactor Calculation' in r:
                        prefac +=1
		else:
			unknown += 1
	print('%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d' % (i, good, nc, energy, iterations, nonlocalt, barrier, nonneg, prefac, unknown))
	i += 1
