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
import prstate
import statelist


class PRStateList(statelist.StateList):
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, state_path, epsilon_e, epsilon_r, use_identical, 
                 list_search_results=None, initial_state = None):
        statelist.StateList.__init__(self, state_path, epsilon_e, epsilon_r,
                                     use_identical, prstate.PRState, list_search_results,
                                     initial_state)

    def register_process(self, reactant_number, product_number, process_id):
        # Get the reactant and product state objects.
        reactant = self.get_state(reactant_number)
        product = self.get_state(product_number)
        reactant.load_process_table()
        product.load_process_table()

        # Make the reactant process point to the product state number.
        reactant.procs[process_id]["product"] = product_number
        reactant.save_process_table()

        product.set_energy(reactant.procs[process_id]["product_energy"])
        # Update the reactant state to point at the new state id.
        reactant.procs[process_id]['product'] = product_number
        reactant.save_process_table()
