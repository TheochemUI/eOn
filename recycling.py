import atoms
import io
import os.path
import copy
import atexit

import numpy

#class SB_Recycling:
#    """ Constructs a super-basin recycling object.
#        
#        states - the statelist object
#        sb_state_list - the list of states in the superbasin exited from
#        current_state - the state just moved to
#        move_didstance - distance an atom must move to be in the "hole"
#    """
#
#    def __init__(self, states, sb_state_list, previous_state, current_state, move_distance):
#        self.sb_state_list = sb_state_list
##        self.previous_state = previous_state
#        self.current_state = current_state
#        self.move_distance = move_distance
##        self.epsilon_r = epsilon_r
##        recycling = Recycling(states, previous_state, self.current_state, move_distance)
#        self.state_possibilities = self.generate_possibilities(self.sb_state_list, self.current_state)
#
#    def generate_possibilities(self, sb_state_list, current_state):
#        for state_number in sb_state_list:
            

class Recycling:
    """ Constructs a saddle point recycling object.
    
        states - the statelist object
        suggested_ref_state - state to make suggestions from
        curr_state - state to make suggestions for
        move_distance - distance an atom must move to be in the "hole"
    """

    def __init__(self, states, suggested_ref_state, curr_state, move_distance):
        """ Initialize the data for the recycling object.
            If there is a file containing the data for the current state,
            use this file.  Otherwise, the previous state will be the reference. """

        self.current_state = curr_state
        self.metadata_path = os.path.join(self.current_state.path,"recycling_info")
        # If this state has already used the recycling process,
        # we know what the reference state and current process number are.
        if os.path.isfile(self.metadata_path):
            self.ref_state_num, self.process_number = self.read_recycling_metadata()
            self.ref_state = states.get_state(self.ref_state_num)
        # Otherwise, this is a new recycling process.
        # Start with the first process and assume that
        # the reference state is the 'suggested'.
        # This is the previous state if called from akmc.
        else:
            self.process_number = 0
            self.ref_state = suggested_ref_state
        self.id_list = self.ref_state.get_process_ids()
        self.ref_state_number = self.ref_state.number

        # Only save data upon exit, to minimize disk io
        atexit.register(self.write_recycling_metadata)

        curr_reactant = self.current_state.get_reactant()
        ref_reactant = self.ref_state.get_reactant()
        # Make a vector of distances between previous
        # current positions for each atom in the state.
        diff = atoms.per_atom_norm(curr_reactant.r - ref_reactant.r, 
                curr_reactant.box)
        
        # The saddle and mode are taken as the reactant
        # and will be modified for each suggested process
        # and then recommended as a search.
        self.saddle = copy.deepcopy(curr_reactant)
        self.mode = numpy.zeros((len(curr_reactant), 3))
        self.in_hole = []
        self.not_in_hole = []
        for i in range(len(self.saddle)):
            if diff[i] < move_distance:
                self.not_in_hole.append(i)
            else:
                self.in_hole.append(i)

    def make_suggestion(self):
        """ Makes a saddle suggestion and returns True.
            If no more saddle suggestions are possible
            (end of rate table has been reached), returns False. """
        # If we've gone through all the possible processes
        if self.process_number >= len(self.id_list):
            return None, None
        #XXX: If only things in the hole move, we will get a bad saddle.
        if len(self.not_in_hole) == 0:
            return None, None
        
        process_saddle = self.ref_state.get_process_saddle(self.process_number)
        process_mode = self.ref_state.get_process_mode(self.process_number)
        
        for i in self.not_in_hole:
            self.saddle.r[i] = process_saddle.r[i]
            self.mode[i] = process_mode[i]
        
        #XXX: These should be returned, but akmc.py depends on the true/false return

        self.process_number += 1
        # Note: Uncomment the final return values to also return the list of indices of atoms
        # which are not in the hole and in the hole.  No change should need to be made to the
        # line in akmc.py which calls this function.
        return copy.deepcopy(self.saddle), copy.deepcopy(self.mode)#, self.not_in_hole, self.in_hole
    
    def read_recycling_metadata(self):
        """ Open the recycling metadata file located in the current state's directory.
            Return the state from which suggestions are being made.
            Return the process number current up for recycling consideration. """
        fi = open(self.metadata_path,"r")
        fi.readline() # The header
        line = fi.readlines()
        ref_state = int(line[0].strip().split()[2])
        process_number = int(line[1].strip().split()[4])
        return ref_state, process_number

    def write_recycling_metadata(self):
        """ Write the recycling metadata file located in the current state's directory. """
        fi = open(self.metadata_path,"w")
        fi.write("Recycling metadata\n")
        fi.write("Reference State:  %d\n" %(self.ref_state.number))
        fi.write("Current Process Number = %d\n" %(self.process_number))
