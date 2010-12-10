#!/usr/bin/env python

import os
import sys
import commands

failstr = "\n\nXXXXX XXXXX FAILED TEST 2 XXXXX XXXXX\n\n"
passstr = "\n\n+++++ +++++ PASSED TEST 2 +++++ +++++\n\n"

# this first command will echo the output to stdout as well
#os.system("../../client/client | tee minimization.txt")
print "\nRunning min_morse_pt_island test\n";
os.system("../../client/client > minimization.txt")

for line in open("minimization.test"):
    if "Final Energy" in line:
        energydata = line.strip().split()
        u = float(energydata[2])

for line in open("minimization.test"):
    if "Final Energy" in line:
        energydata = line.strip().split()
        r = float(energydata[2])

print "Test energy: ",u
print "Ref energy : ",r

if abs(u-r)/u > 0.01:
    print failstr
    sys.exit()

#passed, so delete files"
os.system("rm -f minimization.txt min_0.xyz")
print passstr

