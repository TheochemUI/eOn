#!/usr/bin/env python

import sys
import os
import re

states = sys.argv[1:]
num_states = len(sys.argv)-1
#print "num_states: ",num_states

if(len(states) == 0):
    state_main_dir = "states/"
    state_listdir = os.listdir(state_main_dir)
    for dir in state_listdir:
        state_dir = state_main_dir+dir
        if(os.path.isdir(state_dir)):
            states.append(dir)

for state in states:

    print "\nstate:",state
    state_dir = state_main_dir+str(state)

    # read the search data for the state
    if(os.path.isdir(state_dir)):
        state_file = open(state_dir+"/search_results.txt",'r')
        state_data = state_file.read()
        state_data_lines = state_data.splitlines()[2:]
#        print "\n".join(state_data_lines)

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
#            print " ".join(line_item)
            num_good += 1
            if (re.search("good",search_status)):
                num_unique += 1
            e_sad += float(line_item[2])
            fc_sad += float(line_item[4])
            fc_min += float(line_item[5])
            fc_dyn += float(line_item[6])
        else:
            fc_sad_bad += int(line_item[4])

    # print stats
    print " num search :",num_search
    print " num good   :",num_good
    print " num unique :",num_unique
    if(num_search>0):   
        print " frac good  : {0:.2f}".format(float(num_good)/float(num_search))
    if(num_good>0):
        print " fcalls sad : {0:.2f}".format(float(fc_sad)/float(num_good))
        print " fcalls min : {0:.2f}".format(fc_min/num_good)
        print " fcalls dyn : {0:.2f}".format(fc_dyn/num_good)
        print " fcalls sad total : %d" % fc_sad_tot
        print " fcalls sad abort : %d" % fc_sad_abort
        print " fcalls sad bad : %d" % fc_sad_bad

print ""
