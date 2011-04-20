import random
from math import log
from state import State

def kmc(s):
    ratesum = 0.0
    for proc in s.get_rate_table():
        ratesum += proc['rate']

    u = random.random()
    p = 0.0
    for proc in s.get_rate_table():
        p += proc['rate']/ratesum
        if p>u:
            dt = -log(random.random())/ratesum
            return State.get_state(proc['product'].grid), dt
            break
    else:
        print "Failed to choose process"
        return None


