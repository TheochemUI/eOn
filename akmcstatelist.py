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
import shutil


import atoms
import akmcstate
import statelist


class AKMCStateList(statelist.StateList):
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, kT, thermal_window, max_thermal_window, 
                  initial_state = None, filter_hole = False):
        statelist.StateList.__init__(self, akmcstate.AKMCState, initial_state)
        # aKMC data.
        self.kT = kT
        self.thermal_window = thermal_window
        self.max_thermal_window = max_thermal_window
        self.filter_hole = filter_hole

    def register_process(self, reactant_number, product_number, process_id):

        # Get the reactant and product state objects.
        reactant = self.get_state(reactant_number)
        product = self.get_state(product_number)
        reactant.load_process_table()

        # Forward process (reac->prod)
        # Make the reactant process point to the product state number.
        reactant.procs[process_id]['product'] = product_number
        reactant.save_process_table()

        # Reverse process (prod->reac)
        # Have any processes been determined for the product state
        if product.get_num_procs() != 0:
            product.load_process_table()
            reverse_procs = product.get_process_table()

            saddle_energy = reactant.procs[process_id]['saddle_energy']
            candidates = []

            # An alike reverse process might already exist
            for id in reverse_procs.keys():
                proc = reverse_procs[id]
                if abs(proc['saddle_energy'] - saddle_energy) < self.epsilon_e:
                    # Returns if an alike process is already identified:
                    if (proc['product'] == reactant_number): 
                        return

                    # It might be unidentified
                    elif (proc['product'] == -1):
                        candidates.append(id)

            # Not an exact match, check the candidates
            if len(candidates):
                reactant_conf = reactant.get_reactant()
                for id in candidates:
                    conf = product.get_process_product(id)
                    dist = max(atoms.per_atom_norm(reactant_conf.r - conf.r, reactant_conf.box))
                    if dist < self.epsilon_r:
                        # Remember we are now looking at the reverse processes 
                        reverse_procs[id]['product'] = reactant_number
                        product.save_process_table()
                        return

            # The process is not in the list! 
            reverse_process_id = product.get_num_procs()
            # BUG return SHOULD BE REMOVED, HOWEVER, IF THE PROCESS IS ADDED THE META FILE IT NOT UP TO DATE
            # HOW COULD THIS BE SOLVED 
            return
        else:
            # This must be a new state
            reverse_process_id = 0
            product.set_energy(reactant.procs[process_id]['product_energy'])

        # Add the reverse process in the product state
        barrier = reactant.procs[process_id]['saddle_energy'] - reactant.procs[process_id]['product_energy']
        shutil.copy(reactant.proc_saddle_path(process_id), product.proc_saddle_path(0))
        shutil.copy(reactant.proc_reactant_path(process_id), product.proc_product_path(0))
        shutil.copy(reactant.proc_product_path(process_id), product.proc_reactant_path(0))
        shutil.copy(reactant.proc_results_path(process_id), product.proc_results_path(0))
        shutil.copy(reactant.proc_mode_path(process_id), product.proc_mode_path(0))
        product.append_process_table(id = reverse_process_id, 
                                     saddle_energy = reactant.procs[process_id]['saddle_energy'], 
                                     prefactor = reactant.procs[process_id]['product_prefactor'], 
                                     product = reactant_number, 
                                     product_energy = reactant.get_energy(),
                                     product_prefactor = reactant.procs[process_id]['prefactor'],
                                     barrier = barrier, 
                                     rate = reactant.procs[process_id]['product_prefactor'] * math.exp(-barrier / self.kT), 
                                     repeats = 0)
        product.save_process_table()
            
        return


    def connect_states(self, states):
        for i in states:
            proc_tab = i.get_process_table()
            for j in proc_tab:
                if proc_tab[j]['product'] != -1:
                    continue
                enew = proc_tab[j]['product_energy']
                energetically_close = []
                for state in states:
                    if abs(state.get_energy() - enew) < self.epsilon_e:
                        energetically_close.append(state)

                # Perform distance checks on the energetically close configurations.
                if len(energetically_close) > 0:
                    pnew = i.get_process_product(j)
                    for state in energetically_close:
                        p = state.get_reactant()
                        if atoms.match(p, pnew, True):
                            # Update the reactant state to point at the new state id.
                            self.register_process(i.number, state.number, j)                            
