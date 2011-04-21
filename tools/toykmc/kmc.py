import random
from math import log
from state import State

def kmc(s):
    ratesum = 0.0
    
    if s in superbasin:
        rate_table, dt = #superbasin that s is in /step() 
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
            return State.get_state(proc['product'].grid), dt
            break
    else:
        print "Failed to choose process"
        return None


