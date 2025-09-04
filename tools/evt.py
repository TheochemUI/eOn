#!/usr/bin/env python


import sys

import pathfix
import io, atoms
import config
config.init(sys.argv[1])

import akmc

states = akmc.get_statelist(config.akmc_temperature/11604.5)

e = []
for i in range(states.get_num_states()):
    e.append(states.get_state(i).get_energy())

f = open('evt.plt', 'w')
print("set title \"Energy of Minima versus Time\"", f)
print("set xlabel \"time(s)\"", f)
print("set ylabel \"Energy (eV)\"", f)
print("set nokey", f)
print("plot \"-\" with lines", f)

dyn = io.Dynamics('dynamics_fixed.txt')
for i in dyn.get():
    print(i['totaltime'], e[i['reactant']], i['reactant'], f)
