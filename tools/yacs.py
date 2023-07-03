#!/usr/bin/env python

import pathfix
import atoms
import fileio as io
import sys

try:
    p = io.loadcon(sys.argv[1])
except:
    print("\nusage: yacs input.con")
    sys.exit()

codes = atoms.cnar(p, 3.2, False)

newcodes = [[] for i in range(len(p))]

for i in range(len(p)):
    for code in codes[i]:
        newcodes[i].append(code.split(',') + [codes[i][code]])

for i in range(len(p)):
    newcodes[i].sort(key= lambda x: x[3], reverse=True)
    print('%6d' % i, end=' ')
    for code in newcodes[i]:
        print('%3d%10s' % (code[3], '(%s,%s,%s)' % (code[0], code[1], code[2])), end=' ')
    print()
