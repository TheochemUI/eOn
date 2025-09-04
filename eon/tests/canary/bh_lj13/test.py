#!/usr/bin/env python
import os
import sys
sys.path.insert(0, "../../")
from fileio import parse_results

test_path = os.path.split(os.path.realpath(__file__))[0]
test_name = os.path.basename(test_path)

retval = os.system("../../client/client > stdout.dat")
if retval:
    print("%s: problem running eon" % test_name)
    sys.exit(1)

run_gmin = parse_results('results.dat')['minimum_energy']

global_min = -44.326801

error = abs(global_min-run_gmin)

if error == 0.0:
    print("%s: passed error of %.6e"%(test_name, error))
    sys.exit(0)
else:
    print("%s: failed error of %.6e"%(test_name,error))
    sys.exit(1)
