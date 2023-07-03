#!/usr/bin/env python

import os
import sys

import numpy

import pygmin.storage
from pygmin.landscape import Graph
from pygmin.utils import disconnectivity_graph

import matplotlib.pyplot as plt

nlevels = 20
if len(sys.argv) > 1:
    try:
        nlevels = int(sys.argv[1])
    except:
        print()
        print("usage: run dg.py in an eon akmc directory with the number of energy bins as an optional argument.")
        print()
        sys.exit()

print('Creating a disconnectivity graph with %d energy levels.' % nlevels)


minima_energies = [float(line.strip().split()[1]) for line in open("states/state_table", 'r').readlines()]
num_minima = len(minima_energies)

build_database = True
if os.path.exists('tree.db'):
    print("Using existing tree.db file. If you want to see more recent information, delete tree.db first.")
    build_database = False

db = pygmin.storage.Database('tree.db', accuracy=0)

if build_database:
    print('Loading states...')
    minima = []
    for i in range(num_minima):
        minima.append(db.addMinimum(minima_energies[i], numpy.array([i])))

    print('Loading transitions...')
    transitions = []
    for i in range(num_minima):
        lines = open('states/%d/processtable' % i, 'r').readlines()[1:]
        reac = minima[i]
        for j in range(len(lines)):
            line = lines[j].split()
            prod = int(line[3])
            if prod == -1:
                continue
            prod = minima[prod]
            e = float(line[1])
            proc_id = int(line[0])
            coords = int('%d%d' % (i+1,j))
            if [reac,prod] not in transitions and [prod,reac] not in transitions:
                db.addTransitionState(e, coords, reac, prod)
                transitions.append([reac,prod])

graphwrapper = Graph(db)
graph = graphwrapper.graph

dg = disconnectivity_graph.DisconnectivityGraph(graph,subgraph_size=0,nlevels=nlevels,center_gmin=True)
dg.calculate()

print('Plotting (this can take a while)...')
dg.plot()
print('Saving tree.pdf...')
plt.savefig("tree.pdf")
print('Displaying plot...')
plt.show()
print('Done.')
