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
        # Create the reverse process in the new state.
        barrier = reactant.procs[process_id]['saddle_energy'] - reactant.procs[process_id]['product_energy']
        shutil.copy(reactant.proc_saddle_path(process_id), product.proc_saddle_path(0))
        shutil.copy(reactant.proc_reactant_path(process_id), product.proc_product_path(0))
        shutil.copy(reactant.proc_product_path(process_id), product.proc_reactant_path(0))
        shutil.copy(reactant.proc_results_path(process_id), product.proc_results_path(0))
        shutil.copy(reactant.proc_mode_path(process_id), product.proc_mode_path(0))
        product.append_process_table(id = 0, 
                                     saddle_energy = reactant.procs[process_id]['saddle_energy'], 
                                     prefactor = reactant.procs[process_id]['product_prefactor'], 
                                     product = reactant_number, 
                                     product_energy = reactant.get_energy(),
                                     product_prefactor = reactant.procs[process_id]['prefactor'],
                                     barrier = barrier, 
                                     rate = reactant.procs[process_id]['product_prefactor'] * math.exp(-barrier / self.kT), 
                                     repeats = 0)
        product.save_process_table()

