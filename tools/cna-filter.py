#!/usr/bin/env python

import sys

if 'help' in sys.argv or len(sys.argv) < 4:
    print('\nusage: cna-filter.py input.con output.con [bcc,a15]\n')
    sys.exit()

import pathfix
import fileio as io
import atoms

traj = io.loadcons(sys.argv[1])
filtered = []

for p in traj:

    rejects = []

    if 'bcc' in sys.argv:
        ids = atoms.cnat(p, 3.5)
        for i in range(len(ids)):
            if ids[i] == 5:
                rejects.append(i)

    if 'a15' in sys.argv:
        ids = atoms.cnat(p, 3.5)
        for i in range(len(ids)):
            if ids[i] in [1,2]:
                rejects.append(i)

    rejects = list(set(rejects))

    newp = atoms.Atoms(len(p) - len(rejects))
    newp.box = p.box
    index = 0
    for i in range(len(p)):
        if i not in rejects:
            newp.r[index] = p.r[i]
            newp.free[index] = p.free[i]
            newp.names[index] = p.names[i]
            newp.mass[index] = p.mass[i]
            index += 1
    filtered.append(newp)

io.savecon(sys.argv[2], filtered[0], 'w')

for p in filtered[1:]:
    io.savecon(sys.argv[2], p, 'a')
