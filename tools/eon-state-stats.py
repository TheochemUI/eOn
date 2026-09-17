#!/usr/bin/env python

import glob
import os
import re
import sys
from pathlib import Path

states = sys.argv[1:]
num_states = len(sys.argv)-1

state_listdir=[]
current_cwd=os.getcwd()
if(len(states) == 0):
    state_main_dir = "states/"
    os.chdir(state_main_dir)
    for f in glob.glob('*'):
       try:
          state_listdir.append(int(f))
       except:
          continue
    state_listdir.sort()
    for dir in state_listdir:
        state_dir = str(dir)
        if Path(state_dir).is_dir():
            states.append(str(dir))

print("#{:>11s}:{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>20s}{:>20s}{:>20s}".format(
                                  "state", "num_search", "num_good", "num_unique", "frac_good", "fcalls_sad",
                                  "fcalls_min","fcalls_dyn", "fcalls_sad_total", "fcalls_sad_bad", "fcalls_sad_abort"))
t_num_search    = 0
t_num_good      = 0
t_num_unique    = 0
t_frac_good     = 0
t_fcalls_sad    = 0
t_fcalls_min    = 0
t_fcalls_dyn    = 0
t_fc_sad_tot    = 0
t_fc_sad_abort  = 0
t_fc_sad_bad    = 0

for state in states:

    state_dir = str(state)

    # read the search data for the state
    if Path(state_dir).is_dir():
        state_file = open(state_dir+"/search_results.txt",'r')
        state_data = state_file.read()
        state_data_lines = state_data.splitlines()[2:]

    # initialize
    num_search = 0
    num_good = 0
    num_unique = 0
    fc_sad = 0
    fc_sad_tot = 0
    fc_sad_abort = 0
    fc_sad_bad = 0
    fc_min = 0
    fc_dyn = 0
    e_sad = 0

    # parse the search data
    for line in state_data_lines:
        num_search += 1
        line_item = line.split()
        search_status = line_item[7]
        fc_sad_tot += int(line_item[4])
        if "Nonnegative Displacement Abort" in line:
            fc_sad_abort += int(line_item[4])
        if (re.search("good",search_status) or re.search("repeat",search_status)):
            num_good += 1
            if (re.search("good",search_status)):
                num_unique += 1
            e_sad += float(line_item[2])
            fc_sad += float(line_item[4])
            fc_min += float(line_item[5])
            fc_dyn += float(line_item[6])
        else:
            fc_sad_bad += int(line_item[4])

    if(num_search>0):
        frac_good = float(num_good)/float(num_search)
    if(num_good>0):
        fcalls_sad = float(fc_sad)/float(num_good)
        fcalls_min = fc_min/num_good
        fcalls_dyn = fc_dyn/num_good
    print("{:>11s}:{:>12d}{:>12d}{:>12d}{:>12.2f}{:>12.2f}{:>12.2f}{:>12.2f}{:>20d}{:>20d}{:>20d}".format(
          state, num_search, num_good, num_unique, frac_good, fcalls_sad, fcalls_min, fcalls_dyn, fc_sad_tot, fc_sad_bad, fc_sad_abort))
    t_num_search   += num_search
    t_num_good     += num_good
    t_num_unique   += num_unique
    t_frac_good    += frac_good
    t_fcalls_sad   += fcalls_sad
    t_fcalls_min   += fcalls_min
    t_fcalls_dyn   += fcalls_dyn
    t_fc_sad_tot   += fc_sad_tot
    t_fc_sad_abort += fc_sad_abort
    t_fc_sad_bad   += fc_sad_bad

n_s = len(states)
os.chdir(current_cwd)
print("{:>11s}:{:>12d}{:>12d}{:>12d}{:>12.2f}{:>12.2f}{:>12.2f}{:>12.2f}{:>20d}{:>20d}{:>20d}".format(
      "#average", int(t_num_search/n_s), int(t_num_good/n_s), int(t_num_unique/n_s), t_frac_good/n_s, t_fcalls_sad/n_s,
      t_fcalls_min/n_s, t_fcalls_dyn/n_s, int(t_fc_sad_tot/n_s), int(t_fc_sad_bad/n_s), int(t_fc_sad_abort/n_s)))
print("")
