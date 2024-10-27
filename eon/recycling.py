
import os
import numpy
from eon import atoms
from eon import fileio as io
from eon.config import config as EON_CONFIG
from eon.config import ConfigClass # Typing

class SB_Recycling:
    """ Constructs a super-basin recycling object.
        The basic methods this undergoes are as follows:

        *First, determine whether this is a continuation of a previous recycling effort.
        *If it is a continuation, simply read in the metadata, and be ready to make suggestions
            for saddles from where we left off previously.
        *If it is a new recycling process, first determine whether or not we just left a superbasin.
        *Then, if we did just leave a superbasin, obtain the list of states we will be recycling from,
            (the list of states in the superbasin which was just exited)
        *Next, make educated guesses about what the corresponding states in the new superbasin should look like.
        *Then, using the states in the exited superbasin, find the process for each which leads to the corresponding state
            in the new superbasin.
        *Use the StateList object to make an actual state object of the corresponding states.
        *Finally, it is ready to begin making saddle-search suggestions, using the corresponding pairs
            of states, one from the old superbasin and its corresponding state in the new superbasin.
        *In either case, this will use "Recycling" to actually generate the suggestions for each state pair.

        states - the StateList object
        previous_state - the state object just moved from
        current_state - the state object just moved to
        move_didstance - distance an atom must move to be in the "hole"
        path - the path to the superbasin recycling directory.
        sb_scheme - either "mcamc" or "askmc"
        superbasining - if sb_scheme is "mcamc", then this is the superbasining object; otherwise it is "None"
    """


    def __init__(self, states, previous_state, current_state, move_distance, recycle_save, path, sb_scheme, superbasining):
        """ Initialize the data for the super-basin recycling object. """
        # Establish the base data
        self.states = states
        self.previous_state = previous_state
        self.current_state = current_state
        self.path = path
        self.move_distance = move_distance
        self.recycle_save = recycle_save
        self.sb_scheme = sb_scheme
        self.superbasining = superbasining

        # Establish the working directory.
        if not os.path.isdir(self.path):
            os.mkdir(self.path)

        # Read in the metadata from before.  This sets the following:
        # self.sb_state_nums, self.on_state_num, and self.in_progress
        self.read_metadata()
        # If the passed current_state is not in the current sb_state_nums from metadata,
        # check to see if this "previous_state" was in a superbasin, and if it was,
        # this is a new superbasin recycling process.
        if (self.current_state.number not in [pair[1] for pair in self.sb_state_nums]):
            self.in_progress = True
            # Get sb_state_nums and the sb_states
            if self.sb_scheme is None:
                self.in_progress = False
            elif self.sb_scheme == "mcacm":
                # Getting the superbasin that the previous state was in (if it was in one)
                containing_sb = self.superbasining.get_containing_superbasin(self.previous_state)
                # If there was a return value other than None,
                if containing_sb:
                    # Get the states in the superbasin
                    self.sb_states = containing_sb.states
                    # Find the "real" previous state -- the state in the basin which has
                    # the fastest process leading to the current state.
                    # NOTE -- this may be problematic if multiple states in the basin can lead
                    # to the current state or if multpile processes can go the same state from
                    # a given superbasin state... Not really sure of the best way to deal with that.
                    fastest = 1e-300
                    for state in self.sb_states:
                        state.load_process_table()
                        pid = self.get_process_id(state.procs, self.current_state.number)
                        if pid and state.procs[pid]["rate"] > fastest:
                            self.previous_state = state
                            prev_state_index = self.sb_states.index(state)
                    # Reorganize the states so that the 'corresponding' state is not defined.
                    # This is rather ugly, but will be convenient for keeping track of the corresponding states.
                    self.sb_states = [[state, None] for state in self.sb_states]
                    # Get the state numbers as well.
                    self.sb_state_nums = [[pair[0].number, None] for pair in self.sb_states]
                    # Place the current state into its 'corresponding' position.
                    self.sb_state_nums[prev_state_index][1] = self.current_state.number
                    # Make the states that we hope to find in the next superbasin.
                    self.generate_corresponding_states()
                # Otherwise, the previous state was not in a superbasin -- we're done for now.
                else:
                    self.in_progress = False
            elif self.sb_scheme == "askmc":
                if os.path.isfile(os.path.join(self.path, "current_sb_states")):
                    fi = open(os.path.join(self.path, "current_sb_states"), "r")
                    fi.readline() # The header
                    self.sb_state_nums = eval(fi.readline())
                    fi.close()
                    # If the previous state was not in the last superbasin
                    if self.previous_state.number not in self.sb_state_nums:
                        self.sb_state_nums = [[number, None] for number in self.sb_state_nums]
                        self.in_progress = False
                    # The previous state *was* in the last superbasin
                    else:
                        # Find the index of the previous state.
                        for i in range(len(self.sb_state_nums)):
                            if self.sb_state_nums[i] == self.previous_state.number:
                                prev_state_index = i
                                break
                        # Reorganize the state numbers so that the 'corresponding' state is not defined.
                        # This is rather ugly, but will be convenient for keeping track of the corresponding states.
                        self.sb_state_nums = [[number, None] for number in self.sb_state_nums]
                        # Place the current state into its 'corresponding' position.
                        self.sb_state_nums[prev_state_index][1] = self.current_state.number
                        # Make the state objects of these states.
                        self.sb_states = [[self.states.get_state(pair[0]), None] for pair in self.sb_state_nums]
                        # Make the states that we hope to find in the next superbasin.
                        self.generate_corresponding_states()
                # If askmc hasn't even made a basin yet, we won't be doing recycling.
                else:
                    self.in_progress = False
        # If we're still in the middle of a previous superbasin recycling process,
        # use the established data and remake the list of state objects.
        elif self.in_progress:
            self.sb_states = [[self.states.get_state(pair[0]), self.states.get_state(pair[1])]
                               for pair in self.sb_state_nums]

    def get_process_id(self, current_state_procs, next_state_num):
        """ Return the process id of the process going from the current state
            to the next state (number). """
        next_state_process_id = None
        for process_id in list(current_state_procs.keys()):
            if current_state_procs[process_id]["product"] == next_state_num:
                next_state_process_id = process_id
                break
        return next_state_process_id

    def make_suggestion(self):
        """ This will use the recycling class to make saddle search suggestions. """
        if not self.in_progress:
            self.write_metadata()
            return None, None
        ref_state_index = [pair[1] for pair in self.sb_state_nums].index(self.current_state.number)
        ref_state = self.sb_states[ref_state_index][0]
        recycler = Recycling(self.states, ref_state, self.current_state, self.move_distance, self.recycle_save, from_sb = True)
        sugg_saddle, sugg_mode = recycler.make_suggestion()
        # Write the data before we send the search recommendation, because akmc.py *may* be about to terminate.
        self.write_metadata()
        return sugg_saddle, sugg_mode

    def write_metadata(self):
        """ Write to disk the superbasin state list, and the new possibilities to consider.
            Also write how far along in this superbasin recycling process we are. """
        data_path = os.path.join(self.path, "recycling_data.txt")
        fo = open(data_path, "w")
        fo.write("Superbasin Recycling Metadata\n")
        fo.write("'Prev_Superbasin_state_list' = %s\n"          % repr(self.sb_state_nums))
        fo.write("'in_progress' = %s\n"                         % self.in_progress)
        fo.close()

    def read_metadata(self):
        """ Read the metadata and return the superbasin state list, the "current/previous"
            states, and how far along in the superbasin recycling process we are. """
        data_path = os.path.join(self.path, "recycling_data.txt")
        # If the metadata is not there yet,
        if not os.path.isfile(data_path):
            self.sb_state_nums = []
            return None
        else:
            fi = open(data_path, "r")
            fi.readline() # The header
            self.sb_state_nums  = eval(fi.readline().split("=")[1])
            self.in_progress    = eval(fi.readline().strip().split()[2])
            fi.close()

    def generate_corresponding_states(self):
        """ Generate the list of reactants expected as part of the new superbasin.
            Then, use the StateList object to create the actual states. """
        # First, create the expected states.
        # Only start/continue this process if there hasn't already been a stop flag.
        if not self.in_progress:
            return
        # Only generate states from states in the superbasin which are *not* the previous state.
        indices_to_gen_from = [i for i in range(len(self.sb_states))
                               if self.sb_states[i][0].number != self.previous_state.number]
        current_reactant = self.current_state.get_reactant()
        previous_reactant = self.previous_state.get_reactant()
        # First generate a list of atoms which have moved enough to be considered "in the hole"
        num_atoms = len(current_reactant)
        diff = atoms.per_atom_norm(current_reactant.r - previous_reactant.r,
                                   current_reactant.box)
        moved = []
        for i in range(num_atoms):
            if diff[i] > self.move_distance:
                moved.append(i)
        # Now, generate what we hope will look like the new states.
        # There should be as many as are in the sb state list less the original, premade one.
        state_possibilities = []
        for k in indices_to_gen_from:
            sb_state = self.sb_states[k][0]
            sb_state_reactant = sb_state.get_reactant()
            # Take the previous reactant as a base.
            new_state_reactant = sb_state_reactant.copy()
            # Try to take all the atoms that moved in the trigger process and move them in each of the sb_states.
            # If any of them cannot be moved as they did in the trigger process (i.e. because they're not in the same place as in "previous_state"),
            # then one of the superbasin atoms must have moved during the triggering process.
            # Thus, scrap this whole super-basin recycling effort!
            # (If any of the followng do not get changed to a 1, toss the sb recycling.)
            all_clear = [0]*len(moved)
            # Move the atoms that moved in the process that initiated this recycling.
            for i in moved:
                # If the any of the atoms has the same basic location and the same name,
                # it will be considered the same atom, and will be moved as it did in the triggering process.
                for j in range(num_atoms):
                    if (numpy.linalg.norm(sb_state_reactant.r[j] - previous_reactant.r[i]) < self.move_distance
                      and sb_state_reactant.names[j] == previous_reactant.names[i]):
                        # Move it to where it moved in the trigger process
                        new_state_reactant.r[j] = current_reactant.r[i]
                        all_clear[moved.index(i)] = 1
                        break
            # If we were able to move all of them successfully
            if 0 not in all_clear:
                state_possibilities.append(new_state_reactant)
            # Otherwise, this should no longer be "in progress",
            # and there's no use going on (sniffle).
            else:
                self.in_progress = False
                return

        # Next, having created what we expect the corresponding states in
        # the superbasin to look like, create the actual states.
        # First, find the process that got us from the previous state to the current state.
        self.previous_state.load_process_table()
        ref_pid = self.get_process_id(self.previous_state.procs, self.current_state.number)
        ref_rate = self.previous_state.procs[ref_pid]["rate"]
        ref_barrier = self.previous_state.procs[ref_pid]["barrier"]
        # Then, if we find a similar process from the other sb_states, which leads to a
        # state that is similar to one of our generated possibilities, it's a keeper.
        for i in indices_to_gen_from:
            sb_state = self.sb_states[i][0]
            sb_state.load_process_table()
            state_path = sb_state.path
            product_con = None
            for process_id in list(sb_state.procs.keys()):
                # If the process "looks" similar -- it has a rate less than an order of magnitude different, and a barrier less than 0.2 eV different.
                if (max(sb_state.procs[process_id]["rate"] / ref_rate, ref_rate / sb_state.procs[process_id]["rate"]) < 10
                  and abs(sb_state.procs[process_id]["barrier"] - ref_barrier) < 0.2):
                    # Manually load the product.con for the process id and see if it's similar to the "state_possibility"
                    product_path = os.path.join(state_path, "procdata", "product_%d.con" % process_id)
                    fi = open(product_path, "r")
                    product_con = io.loadcon(fi)
                    fi.close()
                    if atoms.identical(product_con, state_possibilities[indices_to_gen_from.index(i)], self.move_distance):
                        break
                    else:
                        product_con = None
            if product_con is None:
                self.in_progress = False
                return
            # Now that we know which process from the reference state
            # goes to the desired state, make that state.
            self.sb_states[i][1] = self.states.get_product_state(sb_state.number, process_id)
            self.sb_state_nums[i][1] = self.sb_states[i][1].number


class Recycling:
    """ Constructs a saddle point recycling object.

        states - the StateList object
        suggested_ref_state - state to make suggestions from
        curr_state - state to make suggestions for
        move_distance - distance an atom must move to be in the "hole"
    """


    def __init__(self, states, suggested_ref_state, new_state, move_distance, save=False, from_sb=False, config: ConfigClass = EON_CONFIG):
        """ Initialize the data for the recycling object.
            If there is a file containing the data for the current state,
            use this file.  Otherwise, the previous state will be the reference. """
        self.config = config
        self.states = states
        self.ref_state = suggested_ref_state
        self.current_state = new_state
        self.metadata_path = os.path.join(self.current_state.path, "recycling_info")
        self.from_sb = from_sb
        self.save = save
        # If this state has already used the recycling process, we know
        # most of what we need about the state.
        if os.path.isfile(self.metadata_path):
            # Establish (1) self.process_number, (2) self.num_procs,
            # (3) self.in_hole, (4) and self.not_in_hole.
            # Overwrite self.ref_state if not called by SB_recycling
            self.read_recycling_metadata()
            # Re-load the current and reference reactants.
            self.curr_reactant = self.current_state.get_reactant()
            self.ref_reactant = self.ref_state.get_reactant()
        # Otherwise, this is a new recycling process;
        # start with the first process, find the total number of processes,
        # and determine what is and isn't in the hole
        else:
            self.process_number = 0
            # Load the process table to determine the number
            # of processes to recycle.
            # Can't use the rate-table because that skips some processes.
            self.ref_state.load_process_table()
            self.num_procs = len(self.ref_state.procs)

            # Load the reference and current reactants.
            self.curr_reactant = self.current_state.get_reactant()
            self.ref_reactant = self.ref_state.get_reactant()

            # GH: using the active region for all searches
            self.process_atoms = atoms.get_process_atoms(self.curr_reactant, self.ref_reactant,self.config.comp_eps_r,self.config.recycling_active_region)
            #if self.config.saddle_method == 'dynamics':
            #    self.process_atoms = atoms.get_process_atoms(self.curr_reactant, self.ref_reactant,self.config.comp_eps_r,self.config.recycling_active_region)
            #else:
            #    self.process_atoms = atoms.get_process_atoms(self.curr_reactant, self.ref_reactant,self.config.comp_eps_r)

            # Make a vector of distances between previous
            # current positions for each atom in the state.
            diff = atoms.per_atom_norm(self.curr_reactant.r - self.ref_reactant.r,
                                       self.curr_reactant.box)

            # The saddle will be taken as the reactant and will be modified
            # (along with the mode) for each suggested process
            # and then recommended as a search.
            self.moved = []
            self.unmoved = []
            for i in range(len(self.curr_reactant)):
                if diff[i] < move_distance:
                    self.unmoved.append(i)
                else:
                    self.moved.append(i)
        # Set up the mode to be modified for process suggestions.
        self.mode = numpy.zeros((len(self.curr_reactant), 3))

    def make_suggestion(self):
        """ Makes a saddle suggestion and returns True.
            If no more saddle suggestions are possible
            (end of process table has been reached), returns False. """
        # If we've gone through all the possible processes
        if self.process_number >= self.num_procs:
            self.write_recycling_metadata()
            return None, None

        # Make a fresh copy of the "saddle" we're going to send,
        # based on the current reactant.
        saddle = self.curr_reactant.copy()
        # Determine what happens in the reference state saddle
        # for the current process number.
        process_saddle = self.ref_state.get_process_saddle(self.process_number)
        process_mode = self.ref_state.get_process_mode(self.process_number)

        # Now, for all the things that did *not* move getting to this state,
        # suggest this particular process's position to them.

        for i in self.unmoved:
            saddle.r[i] = process_saddle.r[i]
            self.mode[i] = process_mode[i]

        # And, for all the things that *did* move getting to this state,
        # suggest the *motion* that they had available before to the
        # current position they're in.

        for i in self.moved:
            movement = process_saddle.r[i] - self.ref_reactant.r[i]
            saddle.r[i] += movement
            self.mode[i] = process_mode[i]

        # Save suggestions?
        if self.save:
            save_path = os.path.join(self.current_state.path, "saddle_suggestions")
            if not os.path.isdir(save_path):
                os.mkdir(save_path)
            fo = open(os.path.join(save_path, "proc_%d" %self.process_number), "w")
            io.savecon(fo, saddle)
            fo.close()

        # Make a note of the fact that we've tried to recycle another saddle.
        self.process_number += 1
        self.write_recycling_metadata()
        # Note: Uncomment the final return values to also return the list of indices of atoms
        # which are not in the hole and in the hole.  No change should need to be made to the
        # line in akmc.py which calls this function (unless akmc.py wants them as well).
        return saddle.copy(), self.mode.copy()

    def read_recycling_metadata(self):
        """ Open the recycling metadata file located in the current state's directory.
            Return the state from which suggestions are being made.
            Return the process number current up for recycling consideration. """
        fi = open(self.metadata_path, "r")
        fi.readline() # The header
        lines = fi.readlines()
        # If called from superbasin recycling, use its
        # recommendation for reference state.
        if not self.from_sb:
            ref_state_num = int(lines[0].strip().split()[2])
            self.ref_state = self.states.get_state(ref_state_num)
        self.process_number = int(lines[1].split('=')[1].strip())
        self.num_procs = int(lines[2].split('=')[1].strip())
        self.moved = eval(lines[3].split('=')[1].strip())
        self.unmoved = eval(lines[4].split('=')[1].strip())
        self.process_atoms = eval(lines[5].split('=')[1].strip())

    def write_recycling_metadata(self):
        """ Write the recycling metadata file located in the current state's directory. """
        fo = open(self.metadata_path, "w")
        fo.write("Recycling Metadata\n")
        fo.write("Reference State:  %d\n" %(self.ref_state.number))
        fo.write("Current Process Number = %d\n" %(self.process_number))
        fo.write("Number of processes = %d\n" %(self.num_procs))
        fo.write("Indices of 'moved' atoms = %s\n" %(repr(self.moved)))
        fo.write("Indices of 'unmoved' atoms = %s\n" %(repr(self.unmoved)))
        fo.write("Indices of 'process' atoms = %s\n" %(repr(self.process_atoms)))

    def get_moved_indices(self):
        """ Return the indices of atoms that moved in the process getting
            to the current state from the reference state. """
        return self.moved
