""" The statelist module. """

import logging
logger = logging.getLogger('statelist')
import math
import os
import shutil
import sys

from ConfigParser import SafeConfigParser 

import atoms
import config
import akmcstate
import statelist


class AKMCStateList(statelist.StateList):
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, state_path, kT, thermal_window, max_thermal_window, 
                 epsilon_e, epsilon_r, use_identical, initial_state = None, 
                 list_search_results = False, filter_hole = False):
        # aKMC data.
        self.kT = kT
        self.thermal_window = thermal_window
        self.max_thermal_window = max_thermal_window
        self.list_search_results = list_search_results
        self.filter_hole = filter_hole
        statelist.StateList.__init__(self, state_path, epsilon_e, epsilon_r,
                                     use_identical, akmcstate.AKMCState, initial_state)
                 
        ''' Check to see if state_path exists and that state zero exists.
            Initializes state zero when passed a initial_state only if state
            zero doesn't already exist. '''

