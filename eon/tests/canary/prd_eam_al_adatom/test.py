#!/usr/bin/env python

import os
import sys

failstr = "\n\nXXXXX XXXXX FAILED TEST 1 XXXXX XXXXX\n\n"
passstr = "\n\n+++++ +++++ PASSED TEST 1 +++++ +++++\n\n"
print("\nRunning prd_qsc_al_adatom test\n");

os.system("echo 'y' | ../../../parallelreplica.py --reset")
for i in range(5):
    retval = os.system("../../../parallelreplica.py")
    if retval:
        print(failstr)
        sys.exit(1)

unit = open("dynamics.test", 'r').readlines()[2:]
result = open("dynamics.txt", 'r').readlines()[2:]

if len(unit) != len(result):
    print(failstr)
    sys.exit(1)

for i in range(len(unit)):
    unitLine = unit[i].strip().split()
    resultLine = result[i].strip().split()
    if unitLine[1] != resultLine[1]:
        print(failstr)
        sys.exit()
    if unitLine[2] != resultLine[2]:
        print(failstr)
        sys.exit()
    if unitLine[3] != resultLine[3]:
        print(failstr)
        sys.exit()
    u = float(unitLine[4])
    r = float(resultLine[4])
    if abs(u-r)/u > 0.01:
        print(failstr)
        sys.exit()
    u = float(unitLine[5])
    r = float(resultLine[5])
    if abs(u-r)/u > 0.01:
        print(failstr)
        sys.exit()

print(passstr)
