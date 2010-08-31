#!/usr/bin/env python

import os
import sys

failstr = "\n\nXXXXX XXXXX FAILED TEST 1 XXXXX XXXXX\n\n"
passstr = "\n\n+++++ +++++ PASSED TEST 1 +++++ +++++\n\n"

os.system("yes | ../../akmc.py --reset")
for i in range(143):
    os.system("../../akmc.py")

unit = open("dynamics.test", 'r').readlines()
result = open("dynamics.txt", 'r').readlines()

if len(unit) != len(result):
    print failstr
    sys.exit()

for i in range(len(unit)):
    unitLine = unit[i].strip().split()
    resultLine = result[i].strip().split()
    if unitLine[0] != resultLine[0]:
        print failstr
        sys.exit()
    if unitLine[1] != resultLine[1]:
        print failstr
        sys.exit()
    u = float(unitLine[2])
    r = float(resultLine[2])
    if abs(u-r)/u > 0.01:
        print failstr
        sys.exit()
        
print passstr

