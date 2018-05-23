"""
This code is used to construct trajectory file with states
for help:
states_traj.py --h
"""
#!/usr/bin/env python

import sys
import os
import re
import glob
import ast
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from pele.storage import Database
from pele.landscape import database2graph
from pele.utils.disconnectivity_graph import DisconnectivityGraph
import numpy
from ase.neighborlist import neighbor_list as nl
from ase.io import read
from ase.io.trajectory import Trajectory

current = os.getcwd()
state_main_dir = current+"/states/"
#os.chdir(state_main_dir)
#read in states with energy as dictionary
#with open('state_table') as f:
#   states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])


atoms = None
parser = argparse.ArgumentParser()
parser.add_argument('--states', type=int, nargs='+', metavar='StateNumbers', 
        default=None,
        help='States that used to constructe trajectories')

parser.add_argument('--start', type=int, metavar='StateNumbers', 
        default=None,
        help='the start point of the trajectory')

parser.add_argument('--end', type=int, metavar='StateNumbers', 
        default=None,
        help='the end point of the trajectory')

parser.add_argument('--akmc_step', type=int, metavar='AKMCstepnumber', 
        default=None,
        help='generate trajectories based on dynamics.txt by assigning a start akmc_step and an end point(product ID)')

args = parser.parse_args()
if args.start is not None and args.end is not None:
   log_structures = Trajectory(str(args.start)+'_'+str(args.end)+'.traj',
                            'w', atoms)
   for dir in range(args.start, args.end):
       os.chdir(state_main_dir+str(dir))
       atoms = read('reactant.con',index=0)
       log_structures.write(atoms)

if args.states:
   log_structures = Trajectory(str(args.states[0])+'_'+str(args.states[-1])+'.traj',
                            'w', atoms)
   for dir in args.states:
       os.chdir(state_main_dir+str(dir))
       atoms = read('reactant.con',index=0)
       log_structures.write(atoms)

if args.akmc_step is not None and args.end is not None:
   rs = []
   log_structures = Trajectory(str(args.akmc_step)+'.traj',
                            'w', atoms)
   dynamics = pd.read_table('dynamics.txt', delimiter = r'\s+', skiprows = [0,1], names=['step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy'])
   selected_dynamics=dynamics[dynamics['step-number']>=args.akmc_step]
   for i in range(len(selected_dynamics['step-number'])):
      rs.append(selected_dynamics['reactant-id'].iloc[i])
      print "======="
      print i, rs
      try:
        index = rs.index(selected_dynamics['product-id'].iloc[i])
        del rs[index:]
        print 'del', index
      except:
        pass
      if selected_dynamics['product-id'].iloc[i]==args.end:
         print "ending"
         rs.append(selected_dynamics['product-id'].iloc[i])
         break
   print "final:", rs
   for dir in rs:
       os.chdir(state_main_dir+str(dir))
       atoms = read('reactant.con',index=0)
       log_structures.write(atoms)
