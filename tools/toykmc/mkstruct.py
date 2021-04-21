import random
from state import State
import sys

g = []
w = 20
h = 20

density = .05
for i in range(h):
    u = []
    for j in range(w):
        if random.random()<density:
            u.append(True)
        else:
            u.append(False)
    g.append(u)
State(g).save(sys.argv[1])
