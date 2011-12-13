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
print >> f, "set title \"Energy of Minima versus Time\""
print >> f, "set xlabel \"time(s)\""
print >> f, "set ylabel \"Energy (eV)\""
print >> f, "set nokey"
print >> f, "plot \"-\" with lines" 

dyn = io.Dynamics('dynamics_fixed.txt')
for i in dyn.get():
    print >> f, i['totaltime'], e[i['reactant']], i['reactant']

