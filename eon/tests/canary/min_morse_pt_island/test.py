#!/usr/bin/env python

import os
import sys
from pathlib import Path

test_path = str(Path(__file__).resolve().parent)
test_name = Path(test_path).name

# this first command will echo the output to stdout as well
#os.system("../../client/client | tee minimization.txt")
os.system("../../../../client/client > minimization.txt")

for line in open("minimization.txt"):
    if "Final Energy" in line:
        energydata = line.strip().split()
        u = float(energydata[2])

for line in open("minimization.test"):
    if "Final Energy" in line:
        energydata = line.strip().split()
        r = float(energydata[2])

#print "Test energy: ",u
#print "Ref energy : ",r

rel_err = abs(u-r)/abs(u)

if rel_err > 0.01:
    print("%s: failed relative error of %.3f exeeds tolerence" % (test_name,rel_err))
    sys.exit(1)
else:
    #passed, so delete files
    print("%s: passed relative error %.3e" % (test_name, rel_err))
    os.system("rm -f minimization.txt min_0.xyz reactant.con")
