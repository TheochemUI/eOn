#!/usr/bin/env python
"""
This code is used to constructure disconnectivity graph using data (/states) from eon_akmc runs
Use 'pele' package from Wales group
Inputs example:
accuracy = 1e-10
nlevels  = 350
max_state_n = 100
min_state_n = 0
emphasize_state_start = None
emphasize_state_end = None
emax = -249.5
n_Au =  19
symbol = o
minima_to_draw = 1 2 3 4 5 6
color = green brown purple orange pink red
symbol = ^  v  8   s  p  *   h
max_color = blue
check_structure = False
check_e = -249.8
fig_name = ./run9_40.png
fig_width = 6.0
fig_height = 7.0
ytick_fontsize = 14.0
marker_size = 50.0
"""

import sys
import os
import re
import glob
import ast
import pandas as pd
import matplotlib.pyplot as plt
from pele.storage import Database
from pele.landscape import database2graph
from pele.utils.disconnectivity_graph import DisconnectivityGraph
import numpy
from ase.neighborlist import neighbor_list as nl
from ase.io import read

default_paras = dict(
    accuracy = 1e-5,
    nlevels  = 350,
    work_dir = None,
    fig_name = None,
    fig_width = 6.0,
    fig_height = 7.0,
    ytick_fontsize = 14.0,
    marker_size = 50.0,
    #number of Au in the system
    n_Au =  19,
    #highest energy
    emax = None,
    #state range to draw
    max_state_n = 4000,
    min_state_n = 0,
    #define a range of state that will be emphasized
    emphasize_state_start = None,
    emphasize_state_end = None,
    #used to locate a specific state
    specified_state = None,
    #enable to draw states based on number of atoms on surface
    draw_symbol = False,
    minima_to_draw = [1, 2, 3, 4, 5, 6],
    color = ['green', 'brown', 'purple', 'orange', 'pink', 'red'],
    symbol = ['^',  'v',  '8',   's',  'p',  '*',   'h'],
    #color used for state with max number of atoms
    max_color = 'blue',
    #used for checking specifice states
    check_structure = False,
    check_e = -249.8,
    )

type_define=dict(
    accuracy = 1e-5,
    nlevels  = 350,
    work_dir = '/dir',
    fig_name = '../disccon_graph.png',
    fig_width = 6.0,
    fig_height = 7.0,
    ytick_fontsize = 14.0,
    marker_size = 50.0,
    #number of Au in the system
    n_Au =  19,
    #highest energy
    emax = -249.5,
    #state range to draw
    max_state_n = 700,
    min_state_n = 600,
    #define a range of state that will be emphasized
    emphasize_state_start = 662,
    emphasize_state_end = 667,
    #used to locate a specific state
    specified_state = 665,
    #enable to draw states based on number of atoms on surface
    draw_symbol = False,
    minima_to_draw = [1, 2, 3, 4, 5, 6],
    color = ['green', 'brown', 'purple', 'orange', 'pink', 'red'],
    symbol = ['^',  'v',  '8',   's',  'p',  '*',   'h'],
    #color used for state with max number of atoms
    max_color = 'blue',
    #used for checking specifice states
    check_structure = False,
    check_e = -249.8
    )

def readinputs(filename):
    f=open(filename, 'r')
    parameters = {}
    lines=f.readlines()
    for line in lines:
      if line.startswith('#'):
         continue
      fields = line.split('=')
      parameters[fields[0].strip()]=fields[1].replace("\n","").strip()
    return parameters

def type_convertion(para, value):
    if type(para) is int:
       return int(value)
    if type(para) is float:
       return float(value)
    if type(para) is str:
       return value
    if type(para) is list:
       return value.split()
    if type(para) is bool:
       return ast.literal_eval(value)

def compareStru():
#    if min1.cords == min2.cords:
#       return True
    return False

def get_coord(atoms):
    surf_Au = 0
    i = nl('i', atoms,
                     {('Au','Au'):3.3,
                     ('Au','Pd'):3.3,
                     ('Pd','Pd'):3.3
                     })
    coord = numpy.bincount(i)
    index_Au = [ atom.index for atom in atoms if atom.symbol=='Au']
    for i in range(len(coord)):
       if i in index_Au:
          if coord[i] < 10:
            surf_Au += 1
    return surf_Au

paras = readinputs('inputs.ini')
for para in default_paras:
   try:
     paras[para] = type_convertion(type_define[para], paras[para])
   except:
     paras[para] = default_paras[para]
Au_seg = {}
for i in range (paras['n_Au']+1):
    Au_seg[i] = []

if paras['work_dir'] is not None:
   current = paras['work_dir']
   output_dir = os.getcwd()
else:
   current = os.getcwd()
   output_dir = current

state_listdir=[]

#if(len(states) == 0):
state_main_dir = current+"/states/"
os.chdir(state_main_dir)
db = Database(accuracy=paras['accuracy'],compareMinima=compareStru())
#read in states with energy as dictionary
with open('state_table') as f:
   states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])

for f in glob.glob('*'):
   try:
      state_listdir.append(int(f))
   except:
      continue
state_listdir.sort()

process_names=['proc', 'saddle energy', 'prefactor', 'productID', 'product energy', 'product prefactor', 'barrier', 'rate', 'repeats']
emphasize = []
specified_min = None
for dir in state_listdir:
   if dir < paras['min_state_n']:
      continue
   if dir == paras['max_state_n']:
      break
   print("state:", dir)
   try:
     rs = db.addMinimum(states_e[dir], [dir])
   except:
     continue
   if dir==paras['min_state_n']:
      start_state = rs
   if paras['emphasize_state_start']:
      if dir >= paras['emphasize_state_start'] and dir <= paras['emphasize_state_end']:
         emphasize.append(rs)
   if dir == paras['specified_state']:
      specified_min = rs
   os.chdir(state_main_dir+str(dir))

   atoms = read('reactant.con',index=0)
   Au_seg[get_coord(atoms)].append(rs)

   process_table = pd.read_table('processtable', delimiter = r'\s+', skiprows=[0], names=process_names)
   #only forward process stored to avoid duplicated process
   if paras['min_state_n'] > 0:
       selected_procs = process_table[process_table['productID']>=dir]
   else:
       selected_procs = process_table[process_table['productID']>=0]
   for i in range(len(selected_procs['productID'])):
       ps_ID = selected_procs['productID'].iloc[i]
       #print ps_ID
       if ps_ID >= paras['max_state_n']:
          continue
       try:
         ps = db.addMinimum(states_e[ps_ID], [ps_ID], max_n_minima=-1)
         ts = db.addTransitionState(selected_procs['saddle energy'].iloc[i], [dir, ps_ID], rs, ps)
       except:
         continue
print("#of states:",len(db.minima()))
#for mini in db.minima():
#   print mini.coords
Emax = paras['emax']
if Emax is None:
   Emax = -1e20
   for ts in db.transition_states():
      if ts.energy > Emax:
         Emax = ts.energy
   print('max ts:', Emax)
#check the structures with energy larger than check_e
#print paras['check_structure']
if paras['check_structure']:
   for ts in db.transition_states():
       if ts.energy > paras['check_e']:
          print(ts.coords)
   for rs in db.minima():
       if ts.energy > paras['check_e']:
          print(rs.coords)
graph = database2graph(db)
dg = DisconnectivityGraph(graph,nlevels=paras['nlevels'],Emax=Emax+0.05,node_offset=0)
dg.calculate()

fig = plt.figure(figsize=(paras['fig_width'],paras['fig_height']))
#fig = plt.figure(figsize=(7, 8))
fig.set_facecolor('white')
ax = fig.add_subplot(111, adjustable='box')
#ax = fig.add_subplot(111)

dg.plot(linewidth=1.0, axes=ax)

plt.yticks(fontsize=paras['ytick_fontsize'], weight=20)
plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.1)

if paras['min_state_n'] == 0:
   dg.draw_minima([start_state],c='tab:gray', s=50)

if len(emphasize) > 0:
   dg.draw_minima(emphasize,marker='o',c='tab:red')

if specified_min is not None:
   dg.draw_minima([specified_min], marker='8',c='tab:green')
#find max number of atoms on surface
max_key = 0
for key in Au_seg:
   if len(Au_seg[key]) >0 and key > max_key:
      max_key = key
print('Max # of surface Au:',max_key)

if paras['draw_symbol']:
   for i in range(len(paras['minima_to_draw'])):
      print(paras['minima_to_draw'][i],paras['symbol'][i],paras['color'][i])
#      try:
      dg.draw_minima(Au_seg[int(paras['minima_to_draw'][i])],marker=paras['symbol'][i],c='tab:'+paras['color'][i], s=paras['marker_size'])
#      except:
#         continue

   if str(max_key) not in paras['minima_to_draw']:
      dg.draw_minima(Au_seg[max_key],marker='<',c='tab:'+paras['max_color'], s=paras['marker_size'])
if paras['fig_name']:
   plt.savefig(output_dir+'/'+paras['fig_name'])
   plt.show()
else:
   plt.show()
