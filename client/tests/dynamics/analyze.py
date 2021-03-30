#!/usr/bin/env python
import ase
import tsase

print 'loading dynamics.con'
traj = tsase.io.read_con('dynamics.con')

print 'averaging distances'
total_distance = 0.0
N = 0

start = 100
skip = 10
for atoms in traj[start::skip]:
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            r = atoms.get_distance(i,j,True)
            if r < 3.3:
                total_distance += 1
                N += 1
print total_distance/N
