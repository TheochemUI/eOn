""" The statelist module. """


import os
import shutil
import sys
import math
from ConfigParser import SafeConfigParser 

import logging
logger = logging.getLogger('statelist')

import state
import atoms


class StateList:
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """



    def __init__(self, state_path, kT, thermal_window, max_thermal_window, epsilon_e, epsilon_r, use_identical, initial_state = None):
        ''' Check to see if state_path exists and that state zero exists.
            Initializes state zero when passed a initial_state only if state
            zero doesn't already exist. '''

        # aKMC data.
        self.path = state_path
        self.kT = kT
        self.thermal_window = thermal_window
        self.max_thermal_window = max_thermal_window
        self.epsilon_e = epsilon_e
        self.epsilon_r = epsilon_r
        self.use_identical = use_identical

        # Paths
        self.state_table_path = os.path.join(self.path, "state_table")

        # Create the statelist directory if it does not exist.
        if not os.path.isdir(self.path):
            logger.warning("state list path does not exist, creating %s" % self.path)
            os.makedirs(self.path)
            open(self.state_table_path, 'w').close()

        # Create the zero state directory if it does not exist.
        if not os.path.isdir(os.path.join(self.path, "0")):
            if initial_state == None:
                raise IOError("Missing zeroth state directory and no reactant provided.")
            state.State(os.path.join(self.path, "0"), 0, kT, thermal_window, max_thermal_window, epsilon_e, epsilon_r, initial_state)

        # Other class variables.
        self.states = {}
        self.state_table = None



    def get_num_states(self):
        """ Returns the number of lines in the state_table file. """
        self.load_state_table()
        return len(self.state_table)


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
                        if self.use_identical:
                            if atoms.identical(p, pnew, self.epsilon_r):
                                # Update the reactant state to point at the new state id.
                                self.register_process(i.number, state.number, j)                            
                        else:
                            dist = max(atoms.per_atom_norm(p.r - pnew.r, p.box))
                            if dist < self.epsilon_r:
                                self.register_process(i.number, state.number, j)                            


    def get_product_state(self, state_number, process_id):
        ''' Returns a State object referenced by state_number and process_id. '''
        #TODO: Compare configuration of product with existing states.
        # If the number of states in state_table is zero, we need to add the zero state and energy to the state table.
        if self.get_num_states() == 0:
            zst = self.get_state(0)
            self.append_state_table(zst.get_energy())
        # Load the state object containing the process we want the product for.
        st = self.get_state(state_number)
        st.load_process_table()
        # Get the state number for the product.
        newstnr = st.procs[process_id]['product']
        # If the product id is not initialized, make sure it is not a copy of an existing state.
        # Otherwise, create it, connect it to st, and return it.
        if newstnr == -1:
            # Make a list of states for which we need to compare configurations.
            enew = st.procs[process_id]['product_energy']
            energetically_close = []
            for id in range(self.get_num_states()):
                if abs(self.get_state(id).get_energy() - enew) < self.epsilon_e:
                    energetically_close.append(id)
            # Perform distance checks on the energetically close configurations.
            if len(energetically_close) > 0:
                pnew = st.get_process_product(process_id)
                for id in energetically_close:
                    p = self.get_state(id).get_reactant()
                    if self.use_identical:
                        if atoms.identical(p, pnew, self.epsilon_r):
                            # Update the reactant state to point at the new state id.
                            self.register_process(st.number, id, process_id)                            
                            return self.get_state(id)
                    else:
                        dist = max(atoms.per_atom_norm(p.r - pnew.r, p.box))
                        if dist < self.epsilon_r:
                            # Update the reactant state to point at the new state id.
                            self.register_process(st.number, id, process_id)                            
                            return self.get_state(id)
            # The id for the new state is the number of states.
            newstnr = self.get_num_states()
            # Create the new state object.
            newst = state.State(self.state_path(newstnr), newstnr, self.kT, 
                    self.thermal_window, self.max_thermal_window, self.epsilon_e, 
                    self.epsilon_r, st.proc_product_path(process_id))
            self.register_process(st.number, newstnr, process_id)
            # Append the new state to the state table.
            self.append_state_table(st.procs[process_id]['product_energy'])
        # The product state already exists, so get it.
        else:
            newst = self.get_state(newstnr)
        # Return the product state.
        return newst


    def register_process(self, reactant_number, product_number, process_id):
        # Get the reactant and product state objects.
        reactant = self.get_state(reactant_number)
        product = self.get_state(product_number)
        reactant.load_process_table()
        product.load_process_table()
        reverse_procs = product.get_process_table()
        # Make the reactant process point to the product state number.
        reactant.procs[process_id]["product"] = product_number
        reactant.save_process_table()
        if product.get_num_procs() != 0:
            # Check to see if the reverse process is already identified (not -1). If so, return.
            for id in reverse_procs.keys():
                proc = reverse_procs[id]
                if proc["product"] == reactant_number:
                    return
            # The reverse process might already exist, but is unidentified (-1). See if this is the case and fix it.
            candidates = []
            for id in reverse_procs.keys():
                if reverse_procs[id]["product"] == -1:
                    candidates.append(id)
            esaddle = proc["saddle_energy"]
            energetically_close = []
            for id in candidates:
                proc = reverse_procs[id]
                if abs(proc["saddle_energy"] - esaddle) < self.epsilon_e:
                    energetically_close.append(id)
            if len(energetically_close) > 0:
                saddle_config = reactant.get_reactant()
                for id in energetically_close:
                    temp_config = product.get_process_saddle(id)
                    dist = max(atoms.per_atom_norm(saddle_config.r - temp_config.r, saddle_config.box))
                    if dist < self.epsilon_r:
                        reverse_procs[id][product] = product_number
                        product.save_process_table()
                        return
        else:
            # There are no processes, so this must be a new state. Add the reverse process.
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
                                         product_energy = reactant.procs[process_id]['product_energy'],
                                         product_prefactor = reactant.procs[process_id]['prefactor'],
                                         barrier = barrier, 
                                         rate = reactant.procs[process_id]['product_prefactor'] * math.exp(-barrier / self.kT), 
                                         repeats = 0)
            product.save_process_table()


    def get_state(self, state_number):
        ''' Returns a state object. '''
        if state_number in self.states:
            return self.states[state_number]
        st = state.State(os.path.join(self.path, str(state_number)), state_number, self.kT, self.thermal_window, self.max_thermal_window, self.epsilon_e, self.epsilon_r)
        self.states[state_number] = st
        return st


    def load_state_table(self):
        if self.state_table == None:
            f = open(self.state_table_path, 'r')
            lines = f.readlines()
            f.close()
            self.state_table = []
            for l in lines:
                self.state_table.append(float(l.strip().split()[1]))

    def save_state_table(self):
        if self.state_table != None:
            f = open(self.state_table_path, 'w')
            for i in range(len(self.state_table)):
                f.write("% 7d %16.5f\n" % (i, self.state_table[i]))
            f.close()

    def append_state_table(self, energy):
        number = self.get_num_states()
        f = open(self.state_table_path, 'a')
        f.write("% 7d %16.5f\n" % (number, energy))
        f.close()
        if self.state_table != None:    
            self.state_table.append(energy)
        
                                            
        

    

    def state_path(self, state_number):
        """ Utility function to return the compiled path of a state, whether it exists or not. """
        return os.path.join(self.path, str(state_number))



if __name__ == "__main__":
    pass

    


