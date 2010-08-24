#!/usr/bin/env python

import os
import sys

os.system("yes | ../../akmc.py --reset")
for i in range(103):
    os.system("../../akmc.py")

unit = open("dynamics.unit", 'r').readlines()
result = open("dynamics.txt", 'r').readlines()

if len(unit) != len(result):
    print "unit test 1 failed"
    sys.exit()

for i in range(len(unit)):
    unitLine = unit[i].strip().split()
    resultLine = result[i].strip().split()
    if unitLine[0] != resultLine[0]:
        print "unit test 2 failed"
        sys.exit()
    if unitLine[1] != resultLine[1]:
        print "unit test 2 failed"
        sys.exit()
    u = float(unitLine[2])
    r = float(resultLine[2])
    if abs(u-r)/u > 0.01:
        print "unit test 2 failed"
        sys.exit()
        
print "\n***** ***** TEST PASSED ***** *****\n"

