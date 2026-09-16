#!/usr/bin/env python
import os
import sys
from pathlib import Path
sys.path.insert(0, "../../")
from fileio import parse_results

test_path = str(Path(__file__).resolve().parent)
test_name = Path(test_path).name

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
