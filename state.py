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
    """ The State class. """


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
