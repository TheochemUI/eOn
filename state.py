##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##
##-----------------------------------------------------------------------------------
""" The state module. """

import os
import shutil
import math
import logging
logger = logging.getLogger('state')
from ConfigParser import SafeConfigParser 
import StringIO

import numpy

import io
import atoms


class State:
    """ The State super class. """
    def __init__(self, statepath, statenumber, statelist, previous_state_num = -1, 
                 reactant_path = None):                 
        """ Creates a new State, with lazily loaded data. """

        # The parent statelist.
        self.statelist = statelist

        # The path to and number of this state.
        self.path = statepath
        self.number = statenumber

        # Lazily loaded data. Should use the get/set methods for these.
        self.info = None
        self.procs = None
        self.proc_repeat_count = None

        self.info_path = os.path.join(self.path, "info")
        self.procdata_path = os.path.join(self.path, "procdata")
        self.reactant_path = os.path.join(self.path, "reactant.con")
        self.proctable_path = os.path.join(self.path, "processtable")
        self.search_result_path = os.path.join(self.path, "search_results.txt")

        # If this state does not exist on disk, create it.
        if not os.path.isdir(self.path):
            if reactant_path == None:
                raise IOError("State needs a reactant_path when it is being instantiated to disk.")
            os.mkdir(self.path)
            os.mkdir(self.procdata_path)
            shutil.copy(reactant_path, self.reactant_path)
            config = SafeConfigParser()
            config.add_section("MetaData")
            config.set("MetaData", "previous state", str(previous_state_num))            
            config.write(open(self.info_path, 'w'))
            f = open(self.proctable_path, 'w')
            f.write(self.processtable_header)
            f.close()
            if self.statelist.list_search_results:
                f = open(self.search_result_path, 'w')
                f.write(self.search_result_header)
                f.close()

    def __eq__(self, other):
        return self.number == other.number

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "State #%i" % self.number

    def get_energy(self):
        """ Loads the info file if it is not already loaded and returns the energy, or None
            if it is not there. """
        self.load_info()
        try:        
            return self.info.getfloat("MetaData", "reactant energy")
        except:
            return None       

    def set_energy(self, e):
        """ Loads the info file if it has not been loaded and sets the reactant energy 
            variable. """
        self.load_info()
        self.info.set("MetaData", "reactant energy", "%f" % e)
        self.save_info()        

    def get_reactant(self):
        """ Loads the reactant.con into a point and returns it. """
        return io.loadcon(self.reactant_path)

    def load_info(self):
        """ Loads the info file if it has not been loaded. """
        if self.info == None:
            self.info = SafeConfigParser()
            self.info.read(self.info_path)

    def save_info(self):
        """ Saves the info object if it exists. """
        if self.info == None:
            return
        self.info.write(open(self.info_path, 'w'))

    # Utility functions for loading process .con and mode files.
    def get_process_reactant(self, id):
        return io.loadcon(self.proc_reactant_path(id))
    def get_process_product(self, id):
        return io.loadcon(self.proc_product_path(id))

    # Utility functions for compiling procdata paths, whether the files exist or not.    
    def proc_reactant_path(self, id):
        return os.path.join(self.procdata_path, "reactant_%d.con" % id)
    def proc_product_path(self, id):
        return os.path.join(self.procdata_path, "product_%d.con" % id)
    def proc_results_path(self, id):
        return os.path.join(self.procdata_path, "results_%d.dat" % id)

    def get_process(self, id):
        self.load_process_table()
        return self.procs[id]        

    def get_process_ids(self):
        """ Returns the list of ids in the rate table. """
        return [b[0] for b in self.get_ratetable()]

    def get_previous_state(self):
        self.load_info()
        try:
            return self.info.getint("MetaData", "previous state")
        except:
            return -1

    def get_num_procs(self):
        """ Loads the process table if it is not already loaded and returns the length of it """
        self.load_process_table()
        return len(self.procs)

    def load_info(self):
        """ Loads the info file if it has not been loaded. """
        if self.info == None:
            self.info = SafeConfigParser()
            self.info.read(self.info_path)

    def get_process_table(self):
        self.load_process_table()
        return self.procs
