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
