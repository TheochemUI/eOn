#!/usr/bin/env python

from state import State
import os
import time
from kmc import kmc
import sys


s = State.load(sys.argv[1])


t = 0
time_i = 0
time_f = 0
nsteps = 0
for i in xrange(100000):
    time_i = time.time()
    s, dt = kmc(s)
    time_f = time.time()

    t+=dt
    nsteps+=1
    os.system('clear')
    print s
    print "energy:",s.energy
    print "time:",t
    print "dt:", dt
    print "nsteps:", nsteps
    print "nstates:", len(State.states)
    if time_f - time_i > 0:
        print "step rate:", 1/(time_f - time_i)
