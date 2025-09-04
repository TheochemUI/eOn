from __future__ import division
import copy
from math import exp
import odict

from config import *

class State:
    states = odict.OrderedDict()

    @staticmethod
    def get_state(grid):
        gs = State.gridhash(grid)
        if gs in State.states:
            return State.states[gs]
        else:
            s = State(grid)
            State.states[gs] = s
            return s

    @staticmethod
    def saddle_energy(s1, s2):
        return max(s1.energy, s2.energy) + dEsaddle/(max(abs(s1.energy-s2.energy),1))

    #XXX: Ugly
    @staticmethod
    def gridhash(grid):
        #XXX: Hyper-lazy hash function. Still results in speed-up once things get all superbasiney.
        gs = ""
        for i in grid:
            for j in i:
                gs += "T" if j else "F"
        return gs

    def __init__(self, grid, energy = None):
        self.grid = grid
        self.w = len(self.grid[0])
        self.h = len(self.grid)

        if not energy:
            self.energy = self.calc_energy()
        else:
            self.energy = energy
        self.rate_table = None

    def calc_energy(self):
        e = 0
        for i in range(self.h):
            for j in range(self.w):
                if self.grid[i][j]:
                    e += self.calc_energy_at(i,j)
        return e

    def calc_energy_at(self, i, j):
        e = 0
        for z in energy_neighbors:
            e -= self.grid[(i+z[0])%self.h][(j+z[1])%self.w]
        for z in energy_neighbors_2:
            e -= self.grid[(i+z[0])%self.h][(j+z[1])%self.w]*energy_2
        return e

    def __eq__(self, other):
        return self.grid == other.grid

    def get_rate_table(self):
        if self.rate_table:
            return self.rate_table

        self.rate_table = []
        for i in range(self.h):
            for j in range(self.w):
                if self.grid[i][j]:
                    for z in move_neighbors:
                        m,n = (i+z[0])%self.h, (j+z[1])%self.w
                        if not self.grid[m][n]:
                            newgrid = copy.deepcopy(self.grid)
                            newgrid[i][j] = False
                            newgrid[m][n] = True
                            dE = self.calc_energy_at(m,n) - self.calc_energy_at(i,j)

                            proc = {}
                            proc['product'] = State(newgrid, self.energy+2*dE)
                            proc['barrier'] = State.saddle_energy(self, proc['product']) - self.energy
                            proc['rate'] = exp(-proc['barrier']/.01)
                            self.rate_table.append(proc)
        return self.rate_table

    def save(self, filename):
        f = open(filename, 'w')
        print >> f, self
        f.close()

    def __hash__(self):
        return hash(State.gridhash(self.grid))

    @staticmethod
    def load(filename):
        f = open(filename, 'r')
        grid = []
        for i in f:
            if i[0]=='+':
                continue
            gl = []
            for j in i[1:-2]:
                gl.append(False if j==' ' else True)
            grid.append(gl)
        f.close()

        return State(grid)

    def __str__(self):
        out = ""
        out+="+"
        for i in range(self.w):
            out+="-"
        out+="+\n"
        for i in range(self.h):
            out+="|"
            for j in range(self.w):
                if self.grid[i][j]:
                    out+="O"
                else:
                    out+=" "
            out+="|\n"
        out+="+"
        for i in range(self.w):
            out+="-"
        out+="+\n"

        return out

if __name__ == '__main__':
    import random
    g = []
    for i in range(20):
        u = []
        for j in range(70):
            u.append(random.choice([True, False]))
        g.append(u)

    s = State(g)
    print s
    print s.energy
