##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

from ConfigParser import SafeConfigParser 
import os

import config

class AKMCInfo():
    def __init__(self, path_results):
        self._results_path = path_results
        self._akmc_info_path = os.path.join(path_results, config.file_akmc)
        self._info = None

        self.initialize()
        self._changed = False

    def __del__(self):
        """ Saves to disk if it has been changed """
        if self._changed:
            self.save()
        return    

    def initialize(self):
        if self._info == None:
            self._info = SafeConfigParser()

            # This is a new simulation
            if not os.path.isfile(self._akmc_info_path):
                # Make directory for results if it does not exist
                if not os.path.isdir(self._results_path):
                    os.makedirs(self._results_path)
                self.new()
            else:
                self._info.read(self._akmc_info_path)
        return
        
    def new(self):
        self._info.add_section("Simulation Information")
        self._info.set("Simulation Information", "time", str(0))
        self._info.set("Simulation Information", "start_state_num", str(0))
        self._info.set("Simulation Information", "previous_state_num", str(-1))
        self._info.set("Simulation Information", "current_state_num", str(0))

        self._info.add_section("Metadata")
        self._info.set("Metadata", "wuid", str(0))
        self._info.set("Metadata", "first_run", "True")
#        self._info.set("Metadata", "searchdata", repr([]))
        return

    def save(self):
        self._info.write(open(self._akmc_info_path, 'w'))
        return

    def set_time(self, time):
        self._info.set("Simulation Information", "time", str(time))
        self._changed = True
        return
    def get_time(self):
        return self._info.getfloat("Simulation Information", "time")

    def set_current_state_nr(self, current_state_num):
        self._info.set("Simulation Information", "current_state_num", str(current_state_num))
        self._changed = True
        return
    def get_current_state_nr(self):
        return self._info.getint("Simulation Information", "current_state_num")
        
    def set_previous_state_nr(self, previous_state_num):
        self._info.set("Simulation Information", "previous_state_num", str(previous_state_num))
        self._changed = True
        return
    def get_previous_state_nr(self):
        return self._info.getint("Simulation Information", "previous_state_num")        
    
    def set_first_run(self, first_run):
        self._info.set("Metadata", "first_run", "False")
        self._changed = True
        return
    def get_first_run(self):
        first_run = False
        try:
            first_run = self._info.getboolean("Metadata", "first_run")
        except:
            first_run = False
        return first_run
        
    def increase_wuid(self):
        num = self.get_wuid() + 1
        self._info.set("Metadata", "wuid", str(num))
        self._changed = True
        return
    def get_wuid(self):
        return self._info.getint("Metadata", "wuid")
    def set_wuid(self, wuid):
        self._info.set("Metadata", "wuid", str(wuid))
        self._changed = True
        return
            
    def set_searchdata(self, searchdata):
        self._info.set("Metadata", "searchdata", repr(searchdata))    
        self._changed = True    
    def get_searchdata(self):
        searchdata = []
        try:
            searchdata = eval(self._info.get("Metadata", "searchdata"))
        except:
            pass
        return searchdata
        
    def get_random_seed(self):
        seed = None
        try:
            seed = eval(self._info.get("Metadata", "random seed"))
        except:
            pass
        return seed
    def set_random_seed(self, seed):
        self._info.set("Metadata", "random seed", repr(seed))
        self._changed = True
        return



        

class StateInfo():

    def __init__(self, state_path, state_id):
        self._state_info_path = os.path.join(state_path, config.file_state)
        self._state_id = state_id
        self._info = None
        self.initialize()
        self._changed = False

    def __del__(self):
        """ Saves to disk if it has been changed """
        if self._changed:
            self.save()
        return
        
    def initialize(self):
        if self._info == None:
            self._info = SafeConfigParser()
            self._info.read(self._state_info_path)
        return

    def new(self):
        # Make an empty file
        f = open(self._state_info_path, "w")
        f.close()
        self._info = SafeConfigParser()
        self._info.add_section("Metadata")
        self.save()
        
    def save(self):
        self._info.write(open(self._state_info_path, 'w'))
        return
        
    def register_process(self, process_nr):
        new_process = self.increase_proc_repeat_count(process_nr)
        self.increase_good_saddle_count()

        if new_process:
            self.increase_unique_saddle_count()

    # Single values
#    def set_previous_state(self, id):
#        self._info.set("Metadata", "previous state", str(id))
#        self._changed = True
#        return

    def get_previous_state(self):
        previous = -1
        try:
            previous = self._info.getint("MetaData", "previous state")
        except:
            pass
        return previous

    def set_reactant_energy(self, energy):
        self._info.set("Metadata", "reactant energy", str(energy))
        self._changed = True
        return
        
    def get_reactant_energy(self):
        energy = None
        try:        
            energy = self._info.getfloat("Metadata", "reactant energy")
        except:
            pass
        return energy      
        
    def set_lowest_barrier(self, barrier):
        self._info.set("Metadata", "lowest barrier", str(barrier))
        self._changed = True
        return

    def get_lowest_barrier(self):
        barrier = float("inf")
        try:        
            barrier = self._info.getfloat("Metadata", "lowest barrier")
        except:
            pass
        return barrier      

    def increase_unique_saddle_count(self):
        num = self.get_unique_saddle_count() + 1
        self._info.set("Metadata", "unique saddles", str(num))
        self._changed = True
        return
    def get_unique_saddle_count(self):
        saddles = 0
        try:
            saddles = self._info.getint("Metadata", "unique saddles")
        except:
            pass
        return saddles
    def _set_unique_saddle_count(self, num):
        self._info.set("Metadata", "unique saddles", str(num))
        self._changed = True
        return

    def increase_good_saddle_count(self):
        num = self.get_good_saddle_count() + 1
        self._info.set("Metadata", "good saddles", str(num))
        self._changed = True
        return
    def get_good_saddle_count(self):
        saddles = 0
        try:
            saddles = self._info.getint("Metadata", "good saddles")
        except:
            pass
        return saddles
    def _set_good_saddle_count(self, num):
        self._info.set("Metadata", "good saddles", str(num))
        self._changed = True
        return

    def increase_bad_saddle_count(self):
        num = self.get_bad_saddle_count() + 1
        self._info.set("Metadata", "bad saddles", str(num))
        self._changed = True
        return
    def get_bad_saddle_count(self):
        saddles = 0
        try:
            saddles = self._info.getint("Metadata", "bad saddles")
        except:
            pass
        return saddles
    def _set_bad_saddle_count(self, num):
        self._info.set("Metadata", "bad saddles", str(num))
        self._changed = True
        return
        
    def get_total_saddle_count(self):
        return self.get_good_saddle_count() + self.get_bad_saddle_count()

    # Lists of values
    def set_procs_not_in_hole(self, procs):
        self._info.set("Metadata", "procs not in hole", repr(procs))
        self._changed = True

    def get_procs_not_in_hole(self):
        procs = []
        try:
            procs = eval(self._info.get("Metadata", "procs not in hole"))
        except:
            pass
        return procs
        
    def set_procs_in_hole(self, procs):
        self._info.set("Metadata", "procs in hole", repr(procs))    
        self._changed = True    

    def get_procs_in_hole(self):
        procs = []
        try:
            procs = eval(self._info.get("Metadata", "procs in hole"))
        except:
            pass
        return procs
    
    # Dynamic list of values
    def increase_proc_repeat_count(self, proc_nr):
        is_process_new = False
        proc_repeat_count = self.get_proc_repeat_count()

        if proc_nr == len(proc_repeat_count):
            proc_repeat_count.append(1)
            is_process_new = True
        elif (proc_nr < len(proc_repeat_count)) and proc_nr >= 0:
            proc_repeat_count[proc_nr] += 1
            is_process_new = False

        self._info.set("Metadata", "proc repeat count", repr(proc_repeat_count))
        self._changed = True
        return is_process_new

    def get_proc_repeat_count(self):
        repeat_count = []
        try:
            repeat_count = eval(self._info.get("Metadata", "proc repeat count"))
        except:
            pass
        return repeat_count



class CoarseStateListInfo():
    def __init__(self, coarse_path, coarse_nr):
        self._coarse_state_info_path = os.path.join(coarse_path, str(coarse_nr))
        self._coarse_nr = coarse_nr
        self._info = None
        self.initialize()
        self._changed = False

#    This parser should not save itself when destructed
#    Problem araises when merging coarse states into one
#    def __del__(self):
#        """ Saves to disk if it has been changed """
#        if self._changed:
#            self.save()
#        return
        
    def initialize(self):
        if self._info == None:
            self._info = SafeConfigParser()
            self._info.read(self._coarse_state_info_path)
        return

    def new(self):
        # Make an empty file
        f = open(self._coarse_state_info_path, "w")
        f.close()
        self._info = SafeConfigParser()
        self._info.add_section("Metadata")
        self.save()

    def save(self):
        self._info.write(open(self._coarse_state_info_path, 'w'))
        return
        
    def set_lowest_barrier(self, lowest_barrier, 
            lowest_barrier_process_nr, lowest_barrier_state_nr):
        self._info.set("Metadata", "lowest barrier", str(lowest_barrier))
        self._info.set("Metadata", "lowest barrier process nr", str(lowest_barrier_process_nr))
        self._info.set("Metadata", "lowest barrier state nr", str(lowest_barrier_state_nr))
        self._changed = True
        return
    def get_lowest_barrier(self):
        lowest_barrier_data = {}
        try:
            lowest_barrier_data["barrier"] = self._info.getfloat("Metadata", "lowest barrier")
            lowest_barrier_data["process_nr"] = self._info.getint("Metadata", "lowest barrier process nr")
            lowest_barrier_data["state_nr"] = self._info.getint("Metadata", "lowest barrier state nr")
        except:
            lowest_barrier_data = None
        return lowest_barrier_data

    def set_lowest_minima(self, energy):
        self._info.set("Metadata", "lowest minima", str(energy))
        self._changed = True
        return
    def get_lowest_minima(self):
        try:
            energy  = self._info.getfloat("Metadata", "lowest minima")
        except:
            energy = None
        return energy

    def set_escape_rate(self, rate):
        self._info.set("Metadata", "escape rate", str(rate))
        self._changed = True
        return
    def get_escape_rate(self):
        try:
            rate  = self._info.getfloat("Metadata", "escape rate")
        except:
            rate = None
        return rate

    # Dynamic list of values
    def set_state_nrs(self, state_nrs):
        self._info.set("Metadata", "state nrs", repr(state_nrs))
        self._changed = True
        return
    def get_state_nrs(self):
        state_nrs = []
        try:
            state_nrs = eval(self._info.get("Metadata", "state nrs"))
        except:
            pass
        return state_nrs


class CoarseInfo():
    def __init__(self, coarse_path):
        self._coarse_path = coarse_path
        self._coarse_info_path = os.path.join(coarse_path, config.file_coarse)
        self._info = None
        self.initialize()
        self._changed = False

    def __del__(self):
        """ Saves to disk if it has been changed """
        if self._changed:
            self.save()
        return
        
    def initialize(self):
        if self._info == None:
            self._info = SafeConfigParser()
            self._info.read(self._coarse_info_path)
        return

    def new(self):
        # Makes an empty file
        f = open(self._coarse_info_path, "w")
        f.close()
        self._info = SafeConfigParser()
        self._info.add_section("Metadata")
        self.save()

    def save(self):
        self._info.write(open(self._coarse_info_path, 'w'))
        return
        
    def set_clusters(self, cluster_list):
        self._info.set("Metadata", "cluster list", repr(cluster_list))    
        self._changed = True    
    def get_clusters(self):
        cluster_list = []
        try:
            cluster_list = eval(self._info.get("Metadata", "cluster list"))
        except:
            pass
        return cluster_list

    def set_cg_states(self, cg_states_list):
        self._info.set("Metadata", "cg state list", repr(cg_states_list))    
        self._changed = True    
    def get_cg_states(self):
        cg_states_list = []
        try:
            cg_states_list = eval(self._info.get("Metadata", "cg state list"))
        except:
            pass
        return cg_states_list



#class SuperbasinStateInfo():
#    def __init__(self, state_path):
#        self._superbasin_info_path = os.path.join(state_path, config.file_superbasin)
#        self._info = None
#        self.initialize()
#        self._changed = False
#
#    def __del__(self):
#        """ Saves to disk if it has been changed """
#        if self._changed:
#            self.save()
#        return
#        
#    def initialize(self):
#        if self._info == None:
#            self._info = SafeConfigParser()
#            self._info.read(self._superbasin_info_path)
#        return
#
#    def new(self):
#        # Make an empty file
#        f = open(self._superbasin_info_path, "w")
#        f.close()
#        self._info = SafeConfigParser()
#        self._info.add_section("Metadata")
#        self.save()
#        
#    def save(self):
#        self._info.write(open(self._superbasin_info_path, 'w'))
#        return
#
#    def increase_transition_count(self):
#        num = self.get_transition_count() + 1
#        self._info.set("Metadata", "transition count", str(num))
#        self._changed = True
#        return
#    def get_transition_count(self):
#        transitions = 0
#        try:
#            transitions = self._info.getint("Metadata", "transition count")
#        except:
#            pass
#        return transitions
#
#    def set_reference_energy(self, energy):
#        self._info.set("Metadata", "reference energy", str(energy))
#        self._changed = True
#        return
#    def get_reference_energy(self):
#        energy = 0
#        try:
#            energy = self._info.getint("Metadata", "reference energy")
#        except:
#            pass
#        return energy


