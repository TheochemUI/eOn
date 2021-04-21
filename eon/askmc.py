import math
import os
import copy

import numpy
import logging
logger = logging.getLogger('askmc')

ID, ENERGY, PREFACTOR, PRODUCT, PRODUCT_ENERGY, PRODUCT_PREFACTOR, BARRIER, RATE, REPEATS = list(range(9))
processtable_head_fmt = "%7s %16s %11s %9s %16s %17s %8s %12s %7s\n"
mod_processtable_head_fmt = "%7s %16s %11s %9s %16s %17s %8s %12s %12s\n"
processtable_header = processtable_head_fmt % ("proc #", "saddle energy", "prefactor", "product", "product energy", "product prefactor", "barrier", "rate", "repeats")
mod_processtable_header = processtable_head_fmt % ("proc #", "saddle energy", "prefactor", "product", "product energy", "product prefactor", "barrier", "rate", "view count")
processtable_line = "%7d %16.5f %11.5e %9d %16.5f %17.5e %8.5f %12.5e %7d\n"


class ASKMC:
    """ This is a class to keep track of things associated with
        performing the Chatterjee & Voter accelerated superbasin KMC method. """

    def __init__(self, kT, states, confidence, alpha, gamma, barrier_test_on, connection_test_on, sb_recycling_on, path_root, thermal_window, recycle_path):
        self.kT = kT
        self.states = states
        self.Beta = 1/(kT)
        self.sb_recycling_on = sb_recycling_on
        self.path = path_root
        self.thermal_window = thermal_window
        self.recycle_path = recycle_path

        # --- Defining the "acceleration" parameters ---

        # "delta" is more or less a measure of maximum error encountered as a result of using this method.
        # Smaller values of delta correspond to higher assurance that the superbasin exiting will be correct
        #  with respect to exit time and direction (and smaller computational speed-up).
        self.delta = (1 - confidence)
        # alpha is the factor by which the rate constants are lowered when the superbasin criterion is passed.
        # In general Chatterjee & Voter found 2 to be an acceptable value, but if there are "overlapping timescales",
        # it should be less than the square root of gamma.
        # Default in the config:  1.5.
        self.alpha = alpha
        # gamma determines what is a "high" barrier in the superbasin criterion.
        # Chatterjee and Voter found 2 to be an acceptable value.
        # Default in the config:  2.
        self.gamma = gamma
        # Nf is the number of times to have seen a process before performing the superbasin criterion.
        # This comes from eq. 10 from the Chatterjee & Voter paper.
        self.Nf = math.ceil((self.alpha - 1)/self.delta*math.log(1/self.delta))

        # --- Extra check(s) on the superbasin criterion? ---

        # "checksb_barriers" will check to see if there are any low-barrier, unvisited processes
        # originating from states in the superbasin.
        # Default in the config:  On
        if barrier_test_on:
            self.checksb_barriers = True
        else:
            self.checksb_barriers = False
        # "checksb_connections" will implement a method to determine whether there are any missed
        # connections between states identified as "in the superbasin".
        # Default in the config:  Off
        if connection_test_on:
            self.checksb_connections = True
        else:
            self.checksb_connections = False

        if sb_recycling_on:
            self.sb_recycling_on = True
        else:
            self.sb_recycling_on = False

    def compile_process_table(self, current_state):
        """ Load the normal process table, and replace all the rates edited by AS-KMC. """
        current_state_procs_compiled = self.get_real_process_table(current_state)
        current_state_mod_procs = self.get_modified_process_table(current_state)
        for process_id in list(current_state_mod_procs.keys()):
            try:
                current_state_procs_compiled[process_id] = current_state_mod_procs[process_id]
            except KeyError:
                # An exception is raised:
                # If the process is somehow in the modified process list, but it's not
                # in the normal list, there's a problem.
                raise Exception("Somehow the modified process table has a process that\
                    the 'comprehensive' process list does not.  That's a problem...")
        return current_state_procs_compiled

    def get_ratetable(self, current_state):
        """ Generate a rate table according to the kT and thermal_window. """
        current_state_procs = self.compile_process_table(current_state)
        # Find the lowest barrier process for the state.
        lowest = 1e300
        for process_id in list(current_state_procs.keys()):
            barrier = current_state_procs[process_id]["barrier"]
            if barrier < lowest:
                lowest = barrier
        # Make the rate table.
        table = []
        for process_id in list(current_state_procs.keys()):
            proc = current_state_procs[process_id]
            if proc["barrier"] < lowest + (self.kT * self.thermal_window):
                table.append((process_id, proc["rate"]))
        return table

    def get_askmc_metadata(self):
        """ Check to see if the AS-KMC data has been saved.
            If it hasn't, simply use the correct starting values.
            Otherwise, read the current values. """
        # "sb_check_count" is a method to reduce the wasteful calls of working through the superbasin criterion.
        # Rather than checking (calling raiseup) each time a process is seen more than Nf times (see register_transition),
        # it is checked when a process is seen more than (2^sb_check_count)*Nf times with sb_check_count starting at 0
        # and increasing by one each time the superbasin criterion fails.
        # In order to keep track of it, it will be written to disk, along with "num_rate_changes".
        # "num_rate_changes" keeps up with how many times barriers were raised in the entire simulation.
        if not os.path.isfile(os.path.join(self.path,"askmc_data.txt")):
            sb_check_count = 0
            num_rate_changes = 0
        else:
            fi = open(os.path.join(self.path,"askmc_data.txt"),"r")
            lines = fi.readlines()
            sb_check_count = int(lines[1].strip().split()[-1])
            num_rate_changes = int(lines[2].strip().split()[-1])
        return sb_check_count, num_rate_changes

    def save_askmc_metadata(self, sb_check_count, num_rate_changes):
        """ Save the current values of the AS-KMC data. """
        fo = open(os.path.join(self.path,"askmc_data.txt"),"w")
        fo.write("[Info for the Chatterjee & Voter AS-KMC method]\n")
        fo.write("sb_check_count = %d\n" % (sb_check_count))
        fo.write("num_rate_changes = %d\n" % (num_rate_changes))

    def get_real_process_table(self, current_state):
        """ Return the real process table. """
        current_state.load_process_table()
        return copy.deepcopy(current_state.procs)

    def get_modified_process_table(self, current_state):
        """ Return the table of modified processes. If it doesn't exist yet, it's hopefully because the system is in the first state.
            These will be "substituted" in place of their corresponding rates in the normal rate table. """
        mod_proctable_path = os.path.join(current_state.path,"askmc_processtable")
        if not os.path.isfile(mod_proctable_path):
            return {}
        fi = open(mod_proctable_path)
        lines = fi.readlines()
        fi.close()
        procs = {}
        for l in lines[1:]:
            l = l.strip().split()
            procs[int(l[ID])] = {"saddle_energy":     float(l[ENERGY]),
                                 "prefactor":         float(l[PREFACTOR]),
                                 "product":           int  (l[PRODUCT]),
                                 "product_energy":    float(l[PRODUCT_ENERGY]),
                                 "product_prefactor": float(l[PRODUCT_PREFACTOR]),
                                 "barrier":           float(l[BARRIER]),
                                 "rate":              float(l[RATE]),
                                 "view_count":        int  (l[REPEATS])}
        return procs

    def save_modified_process_table(self, current_state, current_state_mod_procs):
        """ Write the modified process table for the current state to disk. """
        mod_proctable_path = os.path.join(current_state.path,"askmc_processtable")
        fo = open(mod_proctable_path,"w")
        fo.write(mod_processtable_header)
        for process_id in list(current_state_mod_procs.keys()):
            proc = current_state_mod_procs[process_id]
            fo.write(processtable_line % (process_id,
                                         proc['saddle_energy'],
                                         proc['prefactor'],
                                         proc['product'],
                                         proc['product_energy'],
                                         proc['product_prefactor'],
                                         proc['barrier'],
                                         proc['rate'],
                                         proc['view_count']))
        fo.close()

    def append_modified_process_table(self, current_state, process_id, saddle_energy, prefactor, product, product_energy, product_prefactor, barrier, rate, view_count):
        """ Append a single line to the modified process table on disk. """
        mod_proctable_path = os.path.join(current_state.path,"askmc_processtable")
        # If the file doesn't exist yet, save a ready copy with the header.
        if not os.path.isfile(mod_proctable_path):
            self.save_modified_process_table(current_state, {})
        fo = open(mod_proctable_path, 'a')
        fo.write(processtable_line % (process_id, saddle_energy, prefactor, product, product_energy, product_prefactor, barrier, rate, view_count))
        fo.close()

    def get_process_id(self, current_state_procs, next_state_num, flag):
        """ Return the process id of the process going from the current state
            to the next state (number). """
        # If in "try" mode, return "None" if there is no process from current_state to next_state_num.
        if flag == "try":
            next_state_process_id = None
            for process_id in list(current_state_procs.keys()):
                if current_state_procs[process_id]["product"] == next_state_num:
                    next_state_process_id = process_id
                    break
        # If in "find" mode, an exception will be raised if the next_state_num is never found.
        elif flag == "find":
            for process_id in list(current_state_procs.keys()):
                if current_state_procs[process_id]["product"] == next_state_num:
                    next_state_process_id = process_id
                    break
        return next_state_process_id

    def register_transition(self, current_state, next_state):
        """ Whenever there is a move (in KMC), update the view_count, and potentially
            lower rate constants and raise the corresponding transition state energies.
            If superbasin recycling is on, and a superbasin has just been exited from,
            Perform the superbasin recycling. """
        # If for some reason, there wasn't actually a move, don't do anything.
        if current_state == next_state:
            return
        sb_check_count, num_rate_changes = self.get_askmc_metadata()
        current_state_mod_procs = self.get_modified_process_table(current_state)
        # Determine if the process is new -- Try to find the process_id in the modified process table.
        next_state_process_id = self.get_process_id(current_state_mod_procs, next_state.number, "try")
        # If we've never been along this path before, add it to our modified list
        if next_state_process_id is None:
            # Find the find the process id from the real process table
            current_state_real_procs = self.get_real_process_table(current_state)
            next_state_process_id = self.get_process_id(current_state_real_procs, next_state.number, "find")
            # Get the info from the actual process table
            procs = current_state_real_procs[next_state_process_id]
            saddle_energy = procs["saddle_energy"]
            prefactor = procs["prefactor"]
            product = procs["product"]
            product_energy = procs["product_energy"]
            product_prefactor = procs["product_prefactor"]
            barrier = procs["barrier"]
            rate = procs["rate"]
            self.append_modified_process_table(current_state, next_state_process_id, saddle_energy, prefactor, product, product_energy, product_prefactor, barrier, rate, 1)
        # Otherwise it's been seen before, and one simply needs to be added to the view_count
        else:
            current_state_mod_procs[next_state_process_id]["view_count"] += 1
            self.save_modified_process_table(current_state, current_state_mod_procs)
            # And if it's been seen enough, check for raising the barriers
            # in the potential superbasin
            if current_state_mod_procs[next_state_process_id]["view_count"] >= (2**sb_check_count)*self.Nf:
                self.raiseup(current_state, next_state, sb_check_count, num_rate_changes)

    def in_array(self, test, array):
        """ Determine whether two arrays are equivalent.
            Because of the use of numpy, the natural python 'in' test does not work properly.  So this is a workaround.
            "test" should be a 1Xn (1X2) list, and "array" should be a mXn (mX2) array. """
        test = numpy.array(test)
        testvar = 0
        for i in array:
            if numpy.array_equal(test, i):
                testvar = 1
                break
        return testvar

    def is_equal(self, a, b):
        """ Determine whether two floats are 'close'. """
        if abs(a - b) < 1e-5:
            return True
        else:
            return False

    def edgelist_to_statelist(self):
        """ Compile a list of state numbers from the current edgelist with no repeats. """
        state_nums = []
        for edge in self.edgelist:
            ref_state_a_num = edge[0]
            ref_state_b_num = edge[1]
            if ref_state_a_num not in state_nums:
                state_nums.append(ref_state_a_num)
            if ref_state_b_num not in state_nums:
                state_nums.append(ref_state_b_num)
        return state_nums

    def missed_connections(self, origEtrans):
        """ Determine whether or not there are low barrier processes which
            connect states in the proposed superbasin which have been missed. """
        # testvar indicates whether or not a problem has arisen.
        # A value of zero indicates "all clear"
        testvar = 0
        # Compile a list of the state *numbers* with no repeats.
        state_nums_in_basin = self.edgelist_to_statelist()
        # Make a list of state objects.
        states_in_basin = [self.states.get_state(state_num) for state_num in state_nums_in_basin]
        # Convert all the "-1's" in the process tables to the correct
        # product state number (if it's in the superbasin)
        self.states.connect_states(states_in_basin)
        # See if any of the states have "low-barrier" processes
        # to other states in the superbasin which haven't been caught
        for state in states_in_basin:
            state_real_procs = self.get_real_process_table(state)
            for process_id in list(state_real_procs.keys()):
                product_state_num = state_real_procs[process_id]["product"]
                # If (1) the process goes somewhere in the basin, (2) it's not in the basin edgelist, and (3) it has a low barrier, that's bad.
                # Note that the "real" process table may be used in this case because if it were in the "modified" process table,
                # then this process *really* should have been caught in the "locsearch" function
                if product_state_num in state_nums_in_basin \
                    and not (self.in_array([product_state_num, state.number], self.edgelist) or self.in_array([state.number, product_state_num], self.edgelist)) \
                    and state_real_procs[process_id]["saddle_energy"] < origEtrans + math.log(self.gamma)/self.Beta:
                        testvar = 1
                        return testvar
        return testvar

    def missed_low_barriers(self, origEtrans):
        """ Determine whether or not there are low barrier processes originating
            from any of the states in identified superbasin that *were missed*. """
        # testvar indicates whether or not a problem has arisen.
        # A value of zero indicates "all clear"
        testvar = 0
        # Compile a list of the state *numbers* with no repeats.
        state_nums_in_basin = self.edgelist_to_statelist()
        # Make a list of state objects from the basin.
        states_in_basin = [self.states.get_state(state_num) for state_num in state_nums_in_basin]
        # See if any of the states have "low-barrier" processes from them which have not been seen.
        # Only processes not in the modified list (unvisited processes) will be considered.
        for state in states_in_basin:
            state_real_procs = self.get_real_process_table(state)
            state_mod_procs = self.get_modified_process_table(state)
            mod_proc_ids = [process_id for process_id in list(state_mod_procs.keys())]
            for process_id in list(state_real_procs.keys()):
                if process_id not in mod_proc_ids and state_real_procs[process_id]["saddle_energy"] < origEtrans + math.log(self.gamma)/self.Beta:
                    testvar = 1
                    return testvar
        return testvar

    def locsearch(self, current_state, origEtrans):
        """ Act recursively to find if the superbasic criterion is met.
            However, ---Note--- that in this implementation, not all neighbors are necessarily checked --
            just those which have been seen so far (which is *probably* all low barrier processes in a superbasin).
            This can be remedied by using the checksb_barriers option.
            origEtrans should be transition state energy (eV) of the process initially used for the call.
            Immediately before each call, self.welltest should be set to 1 and self.edgelist should be any 2X1 numpy array. """
        # "self.welltest" will serve as a flag to stop the search if it has failed.
        # A value of 1 will indicate it has not yet failed.
        # "self.edgelist" is the list of edges that have been examined.
        if self.welltest == 1:
            current_state_mod_procs = self.get_modified_process_table(current_state)
            # scan all the neighbor states
            for next_state_id in list(current_state_mod_procs.keys()):
                next_state_num = current_state_mod_procs[next_state_id]["product"]
                if self.welltest == 1:
                    # If this process has not already been traversed (either direction),
                    # and the barrier is considered "low", its number of sightings (view_count),
                    # in both the forward and backward direction will be checked.
                    if not (self.in_array([current_state.number, next_state_num], self.edgelist) or self.in_array([next_state_num, current_state.number], self.edgelist)) \
                        and current_state_mod_procs[next_state_id]["saddle_energy"] < origEtrans + math.log(self.gamma)/self.Beta:
                            next_state = self.states.get_state(next_state_num)
                            next_state_mod_procs = self.get_modified_process_table(next_state)
                            # Try to find the "process_id" of the current state relative to the next state
                            current_state_id = self.get_process_id(next_state_mod_procs, current_state.number, "try")
                            # If the number of sightings is below the cut-off, the superbasin criterion fails.
                            # Both forward and backward reaction counts are checked.
                            if current_state_mod_procs[next_state_id]["view_count"] < self.Nf or current_state_id is None or next_state_mod_procs[current_state_id]["view_count"] < self.Nf:
                                self.welltest = 0
                            # Otherwise, add this process to "edgelist" to note it's been viewed, and continue the search
                            else:
                                self.edgelist = numpy.vstack((self.edgelist,[current_state.number, next_state_num]))
                                self.locsearch(next_state, origEtrans)

    def raiseup(self, current_state, next_state, sb_check_count, num_rate_changes):
        """ Determine (with the aid of locsearch) if a superbasin is present.
            If a superbasin is found to be present, raise the rate constants,
            lower the saddle energies and the barriers, and save the states in the basin. """
        # Preparing for the "locsearch" function.
        self.welltest = 1
        self.edgelist = numpy.array([0,0])
        sb_check_count += 1
        current_state_mod_procs = self.get_modified_process_table(current_state)
        # Find the "process_id" of the next state
        next_state_id = self.get_process_id(current_state_mod_procs, next_state.number, "find")
        origEtrans = current_state_mod_procs[next_state_id]["saddle_energy"]
        # Perform "locsearch" to determine if the current state is in a superbasin,
        # and if so, the "bottom" of edgelist will be the contained states.
        self.locsearch(current_state, origEtrans)
        if self.welltest:
            self.edgelist = self.edgelist[1:,:]
            # "continue_flag" will indicate if it is safe to proceed (in case one of the following tests is performed).
            continue_flag = 1
            # If the extra check to find missed connections within the discovered superbasin is on,
            # and the test indicates there was a missed process (returns 1), that's problematic.
            if self.checksb_connections and self.missed_connections(origEtrans):
                # NOTE - Matt/Rye, do you think this should be just something to log, or should it raise an exception?
                continue_flag = 0
                logger.info("Based on 'connections-check', default superbasin criterion apparently missed a barrier - Not raising the barriers")
            # If the extra check to find missed low barriers from superbasin states is on,
            # and the test indicates there was a missed low barrier (returns 1), that's problematic.
            if self.checksb_barriers and self.missed_low_barriers(origEtrans):
                continue_flag = 0
                logger.info("Based on 'Barrier-check', default superbasin criterion apparently missed a barrier - Not raising the barriers")
            if continue_flag:
                # Reset the check count, because a barrier raise is being performed
                sb_check_count = 0
                # Performing adjustments
                for state_pair in self.edgelist:
                    # The state numbers
                    state_a_num = state_pair[0]
                    state_b_num = state_pair[1]
                    logger.info('Raising barrier between states %d and %d' % (state_a_num, state_b_num))
                    # The state objects
                    state_a = self.states.get_state(state_a_num)
                    state_b = self.states.get_state(state_b_num)
                    # The state process tables
                    state_a_mod_procs = self.get_modified_process_table(state_a)
                    state_b_mod_procs = self.get_modified_process_table(state_b)
                    # The state energies
                    state_a_energy = state_a.get_energy()
                    state_b_energy = state_b.get_energy()
                    # And the state process id's from the other's perspective
                    state_b_id = self.get_process_id(state_a_mod_procs, state_b_num, "find")
                    state_a_id = self.get_process_id(state_b_mod_procs, state_a_num, "find")
                    # The new rates
                    a_b_new_rate = state_a_mod_procs[state_b_id]["rate"]/self.alpha
                    b_a_new_rate = state_b_mod_procs[state_a_id]["rate"]/self.alpha
                    # The new saddle energies. Both calculated then compared for assurance
                    a_b_new_saddle = math.log(a_b_new_rate/state_a_mod_procs[state_b_id]["prefactor"])/(-self.Beta) + state_a_energy
                    b_a_new_saddle = math.log(b_a_new_rate/state_b_mod_procs[state_a_id]["prefactor"])/(-self.Beta) + state_b_energy
                    if not self.is_equal(a_b_new_saddle, b_a_new_saddle):
                        # NOTE - this exception is primarily for ensuring proper functioning of the method.
                        # If it's never raised, then only one of the above saddles needs to be calculated,
                        # and this "if" clause may be removed. (Note that below, *_*_new_saddle must then be replaced)
                        logger.warning("When recalculating saddle energies from changed rate constants,\
                                        different saddle energies were obtained from the two direction")
                        raise ValueError("When recalculating saddle energies from changed rate constants,\
                                         different saddle energies were obtained from the two directions!")
                    # The new "barrier" values.
                    a_b_new_barrier = a_b_new_saddle - state_a_energy
                    b_a_new_barrier = b_a_new_saddle - state_b_energy
                    # Saving the new data
                    state_a_mod_procs[state_b_id]["rate"] = a_b_new_rate
                    state_b_mod_procs[state_a_id]["rate"] = b_a_new_rate
                    state_a_mod_procs[state_b_id]["saddle_energy"] = a_b_new_saddle
                    state_b_mod_procs[state_a_id]["saddle_energy"] = b_a_new_saddle
                    state_a_mod_procs[state_b_id]["barrier"] = a_b_new_barrier
                    state_b_mod_procs[state_a_id]["barrier"] = b_a_new_barrier
                    state_a_mod_procs[state_b_id]["view_count"] = 0
                    state_b_mod_procs[state_a_id]["view_count"] = 0
                    self.save_modified_process_table(state_a, state_a_mod_procs)
                    self.save_modified_process_table(state_b, state_b_mod_procs)
                # If using super-basin recycling, write out the list of states that are present in this 'superbasin'.
                if self.sb_recycling_on:
                    if not os.path.isdir(self.recycle_path):
                        os.mkdir(self.recycle_path)
                    sb_data_path = os.path.join(self.recycle_path, "current_sb_states")
                    fo = open(sb_data_path, "w")
                    fo.write("Most recent 'superbasin' states:\n")
                    state_list = self.edgelist_to_statelist()
                    fo.write(repr(state_list))
                    fo.close()
                # Update the number of times rate constants have been adjusted.
                num_rate_changes += 1
            # Finally, save the 'tallies' being kept track of
        self.save_askmc_metadata(sb_check_count, num_rate_changes)
