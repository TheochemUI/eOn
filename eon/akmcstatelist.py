""" The statelist module. """

from eon.config import config
import logging
logger = logging.getLogger('statelist')
import math
import shutil

from eon import atoms
from eon import akmcstate
from eon import statelist

from eon.config import config as EON_CONFIG
from eon.config import ConfigClass # Typing


class AKMCStateList(statelist.StateList):
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(
        self,
        kT,
        thermal_window,
        max_thermal_window,
        initial_state=None,
        filter_hole=False,
        config: ConfigClass = EON_CONFIG,
    ):
        self.config = config
        statelist.StateList.__init__(self, akmcstate.AKMCState, initial_state, self.config)
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

        # If the process is of type reac->reac no reverse process needs to be added and we can return.
        # This can be the case if for example a water molecule rotates to mirrored conf.
        if reactant_number == product_number:
            return

        # Loads 'reference' energy for the reactant state as defined in meta-data.
        reactant_energy = reactant.get_energy()
        saddle_energy = reactant.procs[process_id]['saddle_energy']

        # Reverse process (prod->reac).
        # Have any processes been determined for the product state.
        if product.get_num_procs() != 0:
            print("register_process: found process in product state")
            product.load_process_table()
            reverse_procs = product.get_process_table()
            candidates = []

            # An alike reverse process might already exist
            for id in list(reverse_procs.keys()):
                proc = reverse_procs[id]
                if ( abs(proc['saddle_energy'] - saddle_energy) < self.epsilon_e ) and (proc['product']==-1):
                    candidates.append(id)

            if len(candidates):
                print("register_process: some candidate reverse processes found")
                reactant_conf = reactant.get_reactant()
                for id in candidates:
                    conf = product.get_process_product(id)

                    # The process is known but has not been accepted yet.
                    if atoms.match(reactant_conf, conf, self.config.comp_eps_r, self.config.comp_neighbor_cutoff, False):

                        # Reverse process table should be updated to ensure that the two processes (reac->prod & proc->reac) are symmetric.
                        reactant.load_process_table()

                        # Set maximum rate, if defined
#                        cur_rate = reactant.procs[process_id]['product_prefactor'] * math.exp( - ( saddle_energy - product.get_energy() ) /self.kT)
#                        if self.config.akmc_max_rate > 0 and cur_rate > self.config.akmc_max_rate:
#                            # print "max rate exceeded: ", cur_rate
#                            cur_rate = self.config.akmc_max_rate

                        # Set equilibrium rate, if defined
                        print("register_process: into eq rate test")
                        reverse_rate = reactant.procs[process_id]['product_prefactor'] * math.exp(-(saddle_energy - product.get_energy()) / self.kT)
                        forward_barrier = saddle_energy - reactant.get_energy()
                        forward_rate = reactant.procs[process_id]['prefactor'] * math.exp(-forward_barrier / self.kT)

                        eq_rate_flag = False
                        if self.config.akmc_eq_rate > 0 and forward_rate > self.config.akmc_eq_rate and reverse_rate > self.config.akmc_eq_rate:
                            eq_rate_flag = True
                            print("eq_rate exceeded, forward:", forward_rate, " reverse: ", reverse_rate)
                            if forward_rate < reverse_rate:
                                forward_eq_rate = self.config.akmc_eq_rate
                                reverse_eq_rate = self.config.akmc_eq_rate * (reverse_rate / forward_rate)
                            else:
                                forward_eq_rate = self.config.akmc_eq_rate * (forward_rate / reverse_rate)
                                reverse_eq_rate = self.config.akmc_eq_rate
                            print("new eq forward rate:", forward_eq_rate, " reverse: ", reverse_eq_rate)

                        # Remember we are now looking at the reverse processes
                        reverse_procs[id]['product'] = reactant_number
                        reverse_procs[id]['saddle_energy'] = saddle_energy
                        reverse_procs[id]['prefactor'] = reactant.procs[process_id]['product_prefactor']
                        reverse_procs[id]['product_energy'] = reactant.get_energy()
                        reverse_procs[id]['product_prefactor'] = reactant.procs[process_id]['prefactor']
                        reverse_procs[id]['barrier'] = saddle_energy - product.get_energy()
                        reverse_procs[id]['rate'] = reverse_rate
                        product.save_process_table()

                        # If equilibrium rate, change the forward and reverse rate
                        if eq_rate_flag:
                            print("register_process: setting eq rates")
                            reactant.procs[process_id]['rate'] = forward_eq_rate
                            reverse_procs[id]['rate'] = reverse_eq_rate
                            reactant.save_process_table()
                            product.save_process_table()

                        # We are done.
                        return

            # The process is not in the list and must be added as a new process.
            reverse_process_id = product.get_num_procs()

        else:
            # This must be a new state.
            print("register_process: new product state")
            product.set_energy(reactant.procs[process_id]['product_energy'])
            reverse_process_id = 0

        # The product state does not know the reverse process yet.
        # Reverse id has been determined above. (0 if it is a new state, else the last element in the proc table + 1)
        shutil.copy(reactant.proc_saddle_path(process_id), product.proc_saddle_path(reverse_process_id))
        shutil.copy(reactant.proc_reactant_path(process_id), product.proc_product_path(reverse_process_id))
        shutil.copy(reactant.proc_product_path(process_id), product.proc_reactant_path(reverse_process_id))
        shutil.copy(reactant.proc_results_path(process_id), product.proc_results_path(reverse_process_id))
        shutil.copy(reactant.proc_mode_path(process_id), product.proc_mode_path(reverse_process_id))

        # Add the reverse process in the product state
        barrier = saddle_energy - product.get_energy()

        # Set maximum rate, if defined
        # print "max rate code"
#        cur_rate = reactant.procs[process_id]['product_prefactor'] * math.exp(-barrier / self.kT)
#        if self.config.akmc_max_rate > 0 and cur_rate > self.config.akmc_max_rate:
#            # print "max rate exceeded: ", cur_rate
#            cur_rate = self.config.akmc_max_rate

        product.append_process_table(id = reverse_process_id,
                                     saddle_energy = saddle_energy,
                                     prefactor = reactant.procs[process_id]['product_prefactor'],
                                     product = reactant_number,
                                     product_energy = reactant_energy,
                                     product_prefactor = reactant.procs[process_id]['prefactor'],
                                     barrier = barrier,
                                     rate = reactant.procs[process_id]['product_prefactor'] * math.exp(-barrier / self.kT),
#                                     rate = cur_rate,
                                     repeats = 0)

        # Set equilibrium rate, if defined
        print("register_process: into eq rate test")
        forward_barrier = saddle_energy - reactant.get_energy()
        forward_rate = reactant.procs[process_id]['prefactor'] * math.exp(-forward_barrier / self.kT)
        reverse_rate = reactant.procs[process_id]['product_prefactor'] * math.exp(-(saddle_energy - product.get_energy()) / self.kT)

        eq_rate_flag = False
        if self.config.akmc_eq_rate > 0 and forward_rate > self.config.akmc_eq_rate and reverse_rate > self.config.akmc_eq_rate:
            eq_rate_flag = True
            print("eq_rate exceeded, forward:", forward_rate, " reverse: ", reverse_rate)
            if forward_rate < reverse_rate:
                forward_eq_rate = self.config.akmc_eq_rate
                reverse_eq_rate = self.config.akmc_eq_rate * (reverse_rate / forward_rate)
            else:
                forward_eq_rate = self.config.akmc_eq_rate * (forward_rate / reverse_rate)
                reverse_eq_rate = self.config.akmc_eq_rate
            print("new eq forward rate:", forward_eq_rate, " reverse: ", reverse_eq_rate)
            reactant.procs[process_id]['rate'] = forward_eq_rate
            product.procs[reverse_process_id]['rate'] = reverse_eq_rate

        #GH: added this first line
        reactant.save_process_table()
        product.save_process_table()

        # Update the metadata.
        # If this the forward process was a random proc, increase the repeat count.
        if(process_id in reactant.get_proc_random_count() ):
            product.inc_proc_random_count(reverse_process_id)
        product.set_unique_saddle_count( product.get_unique_saddle_count() + 1 )
        product.update_lowest_barrier( barrier )

        # Register the process in the search result file.
        result_fake = {'barrier_reactant_to_product' : barrier,
                       'displacement_saddle_distance' : 0.0,
                       'force_calls_saddle' : 0,
                       'force_calls_minimization' : 0,
                       'force_calls_prefactors' : 0}
        if self.config.akmc_server_side_process_search:
            first_column = "search_id"
        else:
            first_column = "wuid"
        result = { first_column : 0,
                   'type' : 'reverse',
                   'results' : result_fake}
        product.append_search_result(result, 'reverse from '+str(reactant_number), None)
        return

    def connect_states(self, states):
        '''
        This function goes through the process tables of all states in the argument and checks if any of the
        unregistered processes connect these states. It thus tries to connect update the processtables of the states.
        '''
        # print "connect_states, states: ",states
        for i in states:
            proc_tab = i.get_process_table()
            for j in proc_tab:
                # print "checking state, process: ",i," ",j
                if proc_tab[j]['product'] != -1:
                    continue
                enew = proc_tab[j]['product_energy']
                energetically_close = []
                for state in states:
                    if abs(state.get_energy() - enew) < self.epsilon_e:
                        energetically_close.append(state)

                # Perform distance checks on the energetically close configurations.
                if len(energetically_close) > 0:
                    # print "energetically close: ",energetically_close
                    pnew = i.get_process_product(j)
                    for state in energetically_close:
                        p = state.get_reactant()
                        # print "atoms.match between state, process, state: ",i," ",j," ",state
                        if atoms.match(p, pnew, self.config.comp_eps_r, self.config.comp_neighbor_cutoff, True):
                            # Update the reactant state to point at the new state id.
                            # print "structures match"
                            self.register_process(i.number, state.number, j)

    def connect_state_sets(self, states1, states2):
        '''
        This function goes through the process tables of all states in states1 checks if any of the unregistered
        processes connect to a state in state2. It thus tries to connect update the processtables of the states.
        '''
        # print "connect_state_sets: ",states1," ",states2
        for i in states1:
            proc_tab = i.get_process_table()
            for j in proc_tab:
                # print "checking state, process: ",i," ",j
                if proc_tab[j]['product'] != -1:
                    continue
                enew = proc_tab[j]['product_energy']
                energetically_close = []
                for state in states2:
                    if abs(state.get_energy() - enew) < self.epsilon_e:
                        energetically_close.append(state)

                # Perform distance checks on the energetically close configurations.
                if len(energetically_close) > 0:
                    # print "energetically close: ",energetically_close
                    pnew = i.get_process_product(j)
                    for state in energetically_close:
                        p = state.get_reactant()
                        # print "atoms.match between state, process, state: ",i," ",j," ",state
                        if atoms.match(p, pnew, self.config.comp_eps_r, self.config.comp_neighbor_cutoff, True):
                            # Update the reactant state to point at the new state id.
                            # print "structures match"
                            self.register_process(i.number, state.number, j)

        for i in states2:
            proc_tab = i.get_process_table()
            for j in proc_tab:
                # print "checking state, process: ",i," ",j
                if proc_tab[j]['product'] != -1:
                    continue
                enew = proc_tab[j]['product_energy']
                energetically_close = []
                for state in states1:
                    if abs(state.get_energy() - enew) < self.epsilon_e:
                        energetically_close.append(state)

                # Perform distance checks on the energetically close self.configurations.
                if len(energetically_close) > 0:
                    # print "energetically close: ",energetically_close
                    pnew = i.get_process_product(j)
                    for state in energetically_close:
                        p = state.get_reactant()
                        # print "atoms.match between state, process, state: ",i," ",j," ",state
                        if atoms.match(p, pnew, self.config.comp_eps_r, self.config.comp_neighbor_cutoff, True):
                            # Update the reactant state to point at the new state id.
                            # print "structures match"
                            self.register_process(i.number, state.number, j)
