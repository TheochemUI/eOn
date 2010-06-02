import atoms
import io
import os.path
import copy

import numpy

class Recycling:
    def __init__(self, prev_state, curr_state, move_distance, start_process_number = 0):
        '''
        Constructs a saddle point recycling object.
        
        prev_state - state to make suggestions from
        curr_state - state to make suggestions for
        move_distance - distance an atom must to be in the "hole"
        start_process_number - number of the process to start from in prev_state
        '''
        
        self.prev_state = prev_state
        
        self.process_num = start_process_number
        self.id_list = prev_state.get_process_ids()
        
        prev_reactant = prev_state.get_reactant()
        curr_reactant = curr_state.get_reactant()
        diff = atoms.per_atom_norm(curr_reactant.r - prev_reactant.r, 
                curr_reactant.box)
        
        self.saddle = copy.deepcopy(curr_reactant)
        self.mode = numpy.zeros((len(curr_reactant), 3))
        self.not_in_hole = []
        for i in range(len(self.saddle)):
            if diff[i] < move_distance:
                self.not_in_hole.append(i)

    def make_suggestion(self, path):
        '''Makes a saddle suggestion and returns True. If no more saddle suggestions are possible (end of rate table has been reached), returns False'''
        if self.process_num >= len(self.id_list):
            return False
        
        process_id = self.id_list[self.process_num]
        process_saddle = self.prev_state.get_process_saddle(process_id)
        process_mode   = self.prev_state.get_process_mode(process_id)
        
        #XXX: If only things in the hole move, we will get a bad saddle. 
        for i in self.not_in_hole:
            self.saddle.r[i] = process_saddle.r[i]
            self.mode[i] = process_mode[i]

        io.savecon(os.path.join(path, "displacement_passed.con"), self.saddle) 
        io.save_mode(os.path.join(path, "mode_passed.dat"), self.mode, self.saddle)

        self.process_num += 1
        return True
