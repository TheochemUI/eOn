import random
from math import log
from state import State
import superbasinscheme

from config import *

if use_sb:
    superbasining = superbasinscheme.TransitionCounting(50)

def kmc(s):
    ratesum = 0.0

    dosb = use_sb and superbasining.get_containing_superbasin(s)
    if dosb:
        rate_table, dt, exit_state = superbasining.get_containing_superbasin(s).step(s)
    else:
        rate_table = s.get_rate_table()
        for proc in rate_table:
            ratesum += proc['rate']
        dt = -log(random.random())/ratesum
    u = random.random()
    p = 0.0
    for proc in rate_table:
        p += proc['rate']/ratesum
        if p>u:
            newst = State.get_state(proc['product'].grid)
            if use_sb:
                if dosb:
                    superbasining.register_transition(exit_state, newst)
                else:
                    superbasining.register_transition(s, newst)
            return newst, dt
    else:
        print "Failed to choose process"
        return None
