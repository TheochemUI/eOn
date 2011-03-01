##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

""" The statelist module. """

import logging
logger = logging.getLogger('statelist')
import os
import math
import shutil


import atoms
import config
import state
import tables
import parsers

class StateList:
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, StateClass, initial_state = None):
        ''' Check to see if state_path exists and that state zero exists.
            Initializes state zero when passed a initial_state only if state
            zero doesn't already exist. '''

        self._StateClass = StateClass

        # Paths
        self._path = config.path_states
        self._state_table_path = os.path.join(self._path, config.file_state_table)

        # Create the statelist directory if it does not exist.
        if not os.path.isdir(self._path):
            logger.warning("state list path does not exist, creating %s" % self._path)
            os.makedirs(self._path)

        # Create the zero state directory if it does not exist.
        if not os.path.isdir(os.path.join(self._path, "0")):
            if initial_state == None:
                raise IOError("Missing zeroth state directory and no reactant provided.")
            self._StateClass(statepath = os.path.join(self._path, "0"), 
                        statenumber = 0, 
                        statelist = self,
                        reactant_path = initial_state)

        # Other class variables.
        self._states = {}
        self._state_table = None

    def __len__(self):
        return self.get_num_states()

    def get_state(self, state_nr):
        ''' Returns a state object. '''
        state = None
        # It is already in the table
        if state_nr in self._states:
            state = self._states[state_nr]
        # It is added to the table
        else:
            state = self._StateClass(statepath = os.path.join(self._path, str(state_nr)), 
                                statenumber = state_nr, 
                                statelist = self)
            self._states[state_nr] = state
            
        return state

    def state_path(self, state_nr):
        """ Utility function to return the compiled path of a state, whether it exists or not. """
        return os.path.join(self._path, str(state_nr))

    def initialize_state_table(self):
        # Loads table if nessesary
        if self._state_table == None:
            self._state_table = tables.StateTable(self._state_table_path) 
        return
        
    def get_num_states(self):
        """ Returns the number of lines in the state_table file. """
        self.initialize_state_table()
        return len(self._state_table)

    def get_energy_state_nr(self, nr):
        """ Returns the state nr at key in state_table file. """
        self.initialize_state_table()        
        return self._state_table.get_energy_state_nr[nr]

    def get_min_energy_state(self):
        self.initialize_state_table()        
        min_energy = float('inf')
        min_nr = None
        for nr in self.get_state_nrs():
            energy = self.get_energy_state_nr(nr)
            if energy < min_energy:
                min_energy = energy
                min_nr = nr
        return min_energy, min_nr

    def get_state_nrs(self):
        """ Returns the indices for the states stored. """
        self.initialize_state_table()
        return self._state_table.keys()

    def add_state(self, nr, energy):
        self.initialize_state_table()
        self._state_table.add_state(nr, energy)
        
    def if_known_return_state_nr(self, state_energy, state_structure):
        match_nr = -1
        # Make a list of states for which we need to compare configurations.
        candidates = []
        for trial_nr in range(self.get_num_states()):
            trial_energy = self.get_state(trial_nr).get_energy()
            if abs(trial_energy - state_energy) < config.comp_eps_e:
                candidates.append(trial_nr)

        # Perform distance checks on the energetically close configurations.
        if len(candidates) > 0:
            for trial_nr in candidates:
                trial_structure = self.get_state(trial_nr).get_reactant()

                if atoms.match(trial_structure, state_structure, True):
                    match_nr = trial_nr
                    break
        return match_nr
        

        

class AKMCStateList(StateList):
    """ The AKMCStateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, initial_state = None, filter_hole = False):
        StateList.__init__(self, state.AKMCState, initial_state)
        
        # AKMC data.
        self.filter_hole = filter_hole

    def register_choosen_process(self, state_nr, process_nr):
        ''' Returns a State object referenced by state_nr and process_nr. '''
        # If the number of states in state_table is zero
        # we need to add the zero state and energy to the state table.
        self.initialize_state_table()
        if self.get_num_states() == 0:
            zst = self.get_state(0)
            self.add_state(0, zst.get_energy())

        # Load the state object containing the process we want the product for.
        st = self.get_state(state_nr)

        product_nr = st.get_process_product_nr(process_nr)

        # If the product nr is not initialized, make sure it is not a copy of an existing state.
        if product_nr == -1:
            product_energy = st.get_process_product_energy(process_nr)
            product_structure = st.get_process_product(process_nr)
            product_nr = self.if_known_return_state_nr(product_energy, product_structure)
            
        # It is not known, create it.
        if product_nr == -1:
            # The id for the new state is the number of states.
            product_nr = self.get_num_states()
            # Create the new state object.
            product = self._StateClass(statepath = self.state_path(product_nr), 
                                    statenumber = product_nr, 
                                    statelist = self,
                                    reactant_path = st.proc_product_path(process_nr),
                                    energy = product_energy) 
            # Append the new state to the state table.
            self.add_state(product_nr, st.get_process_product_energy(process_nr))

        # The product state is already known, so get it.
        else:
            product = self.get_state(product_nr)

        self._register_in_process_tables(state_nr, product_nr, process_nr)

        # Return the product state.
        return product

    def _register_in_process_tables(self, reactant_number, product_number, for_nr):
        # Forward process (reac->prod)
        # Make the reactant process point to the product state number.
        reactant = self.get_state(reactant_number)
        reactant.set_process_product_nr(for_nr, product_number)

        # Reverse process (prod->reac)
        # Have any processes been determined for the product state
        rev_known = False
        product = self.get_state(product_number)

        if product.get_num_procs() != 0:
            for_saddle_energy = reactant.get_process_saddle_energy(for_nr)

            # An alike reverse process might already exist
            candidates = []
            for rev_nr in xrange(product.get_num_procs()):
                rev_saddle_energy = product.get_process_saddle_energy(rev_nr)
                if abs(rev_saddle_energy-for_saddle_energy) < config.comp_eps_e:
                    if product.get_process_product_nr(rev_nr) == reactant_number:
                        rev_known = True
                    else:
                        candidates.append(rev_nr)
            
            if len(candidates):
                reactant_conf = reactant.get_reactant()
                for rev_nr in candidates:
                    conf = product.get_process_product(rev_nr)
                    dist = max(atoms.per_atom_norm(reactant_conf.r - conf.r, reactant_conf.box))
                    if dist < config.comp_eps_r:
                        # Remember we are now looking at the reverse processes 
                        product.set_process_product_nr(rev_nr, reactant_number)
                        rev_known = True

        if not rev_known:
            process = reactant.get_process(for_nr)
            product.register_reverse_process(reactant, process, for_nr)
        return

    def connect_states(self, states):
        for i in states:
            proc_tab = i.get_process_table()
            for j in proc_tab:
                if proc_tab[j]['product'] != -1:
                    continue
                enew = proc_tab[j]['product_energy']
                candidates = []
                for state in states:
                    if abs(state.get_energy() - enew) < config.comp_eps_e:
                        candidates.append(state)

                # Perform distance checks on the energetically close configurations.
                if len(candidates) > 0:
                    pnew = i.get_process_product(j)
                    for state in candidates:
                        p = state.get_reactant()
                        if atoms.match(p, pnew, True):
                            # Update the reactant state to point at the new state nr.
                            self.register_process(i.get_nr(), state.get_nr(), j)
        return

class CoarseStateList(StateList):
    """ The CoarseStateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, nr, state_nrs):
        StateList.__init__(self, state.AKMCState)
        
        self._nr = nr
        self._state_nrs = state_nrs
        self._info = None
    
#        # Creates new info file for the coarse state list
#        if not os.path.isfile(self._path_coarse, self._nr):
#            self._info = parsers.CoarseStateInfo(self._path_coarse, self._nr)
#            self._info.new()
#        return
        
#    def add_state_nr(self, nr, energy):
#        StateList.add_state(self, nr, energy)
#        self._initialize_info()
#        self._info.add_state(nr)
#
#        state_nrs = self._info.get_state_nrs()
#        
#        lowest_barrier = float("inf")
#        lowest_barrier_process_nr = None
#        lowest_barrier_state_nr = None
#
#        for state_nr in state_nrs:
#            st = self.get_state(state_nr)
#            barrier, process_nr = st.get_lowest_barrier_excluding_product_nrs(state_nrs)
#            if barrier < lowest_barrier:
#                lowest_barrier = barrier
#                lowest_barrier_process_nr = process_nr
#                lowest_barrier_state_nr = state_nr                
#
#        self._info.set_lowest_barrier(lowest_barrier, lowest_barrier_process_nr, lowest_barrier_state_nr)
        
    #------------------
    # Info file related
    def _initialize_info(self):
        if self._info == None:
            self._info = parsers.CoarseStateListInfo(config.path_coarse_states, self._nr)
            state_nrs = self._info.get_state_nrs()
            # State has been added
            # Info file is outdated
            if len(state_nrs) != len(self._state_nrs):
                self._info.new()
                self._info.set_state_nrs(self._state_nrs)
                self.determine_lowest_escape_barrier()
                self.determine_lowest_minima()
                self.determine_escape_rate()
                self._info.save()
#                print "needed update"
#                print state_nrs
#                print self._state_nrs
#                raw_input()
#        print "get_state_nrs():",self._info.get_state_nrs()
#        print self._info.get_state_nrs()
        return

    def get_lowest_minima(self):
        """Return the lowest minima."""
        self._initialize_info()
        return self._info.get_lowest_minima()

    def get_lowest_escape_barrier(self):
        """Return the lowest exit barrier."""
        self._initialize_info()
        return self._info.get_lowest_barrier()

    def get_escape_rate(self):
        """Return the lowest exit barrier."""
        self._initialize_info()
        return self._info.get_escape_rate()


    def determine_lowest_minima(self):
        self._initialize_info()
        if self._info.get_lowest_minima() is None:
            state_nrs = self._info.get_state_nrs()
            lowest_energy = float("inf")
            for state_nr in state_nrs:
                st = self.get_state(state_nr)
                energy = st.get_energy()
                if energy < lowest_energy:
                    lowest_energy = energy
            self._info.set_lowest_minima(lowest_energy)
            self._info.save()
        return

    def determine_lowest_escape_barrier(self):
        self._initialize_info()
        if self._info.get_lowest_barrier() is None:
            state_nrs = self._info.get_state_nrs()
            lowest_barrier = float("inf")
            for state_nr in state_nrs:
                st = self.get_state(state_nr)
                barrier, process_nr, product = st.get_lowest_barrier_excluding_product_nrs(state_nrs)
                if barrier < lowest_barrier:
                    lowest_barrier = barrier
                    lowest_barrier_process_nr = process_nr
                    lowest_barrier_state_nr = state_nr                
            self._info.set_lowest_barrier(lowest_barrier, lowest_barrier_process_nr, lowest_barrier_state_nr)
            self._info.save()
        return

    def determine_escape_rate(self):
        self._initialize_info()
        if self._info.get_escape_rate() is None:
            state_nrs = self._info.get_state_nrs()
            # Assuming a constant prefactor of 1 ps
            prefactor = 1e12
            sum_saddles = 0.
            sum_minima = 0.
            kBT = config.main_temperature * 8.617e-5
            reference_energy = self.get_lowest_minima()
            
            for state_nr in state_nrs:
                st = self.get_state(state_nr)
                sum_minima += math.exp(-(st.get_energy()-reference_energy)/kBT) 
                saddle_energies = st.get_saddle_energies_excluding_product_nrs(state_nrs)
                for saddle_energy in saddle_energies:
                    sum_saddles += math.exp(-(saddle_energy-reference_energy)/kBT) 
            ratio = sum_saddles / sum_minima
            self._info.set_escape_rate(prefactor * ratio)
            self._info.save()
        return 

        
     
        
        
