""" The state module. """

import os
import shutil
import math
import logging
logger = logging.getLogger('state')
from ConfigParser import SafeConfigParser 

import numpy

import io
import atoms



ID, ENERGY, PREFACTOR, PRODUCT, PRODUCT_ENERGY, PRODUCT_PREFACTOR, BARRIER, RATE, REPEATS = range(9)
processtable_head_fmt = "%7s %16s %11s %9s %16s %17s %8s %12s %7s\n"
processtable_header = processtable_head_fmt % ("proc #", "saddle energy", "prefactor", "product", "product energy", "product prefactor", "barrier", "rate", "repeats")
processtable_line = "%7d %16.5f %11.5e %9d %16.5f %17.5e %8.5f %12.5e %7d\n"
#GH search_result_header = "%8s %16s %8s %16s %16s %16s %16s    %s\n" % ("wuid", "type", "barrier", "max-atom-dist", "saddle-fcs", "minimize-fcs", "prefactor-fcs", "result")
search_result_header = "%8s %10s %10s %10s %10s %10s %10s    %s\n" % ("wuid", "type", "barrier", "max-dist", "sad-fcs", "mins-fcs", "pref-fcs", "result")
search_result_header += "-" * len(search_result_header) + '\n'



class State:
    """ The State class. """


    def __init__(self, statepath, statenumber, kT, thermal_window, max_thermal_window, epsilon_e, epsilon_r, reactant_path = None):
        """ Creates a new State, with lazily loaded data. """

        # The path to and number of this state.
        self.path = statepath                               
        self.number = statenumber

        self.info_path = os.path.join(self.path, "info")
        self.procdata_path = os.path.join(self.path, "procdata")
        self.bad_procdata_path = os.path.join(self.path, "badprocdata")
        self.reactant_path = os.path.join(self.path, "reactant.con")
        self.proctable_path = os.path.join(self.path, "processtable")
        self.search_result_path = os.path.join(self.path, "search_results.txt")

        # Variables for comparing energies and distances.
        self.epsilon_e = epsilon_e
        self.epsilon_r = epsilon_r

        # Rate table variables.
        self.kT = kT
        self.thermal_window = thermal_window
        self.max_thermal_window = max_thermal_window

        # Lazily loaded data. Should use the get/set methods for these.
        self.info = None
        self.procs = None

        # If this state does not exist on disk, create it.
        if not os.path.isdir(self.path):
            if reactant_path == None:
                raise IOError("State needs a reactant_path when it is being instantiated to disk.")
            os.mkdir(self.path)
            os.mkdir(self.procdata_path)
            shutil.copy(reactant_path, self.reactant_path)
            config = SafeConfigParser()
            config.add_section("MetaData")
            config.write(open(self.info_path, 'w'))
            f = open(self.proctable_path, 'w')
            f.write(processtable_header)
            f.close()
            f = open(self.search_result_path, 'w')
            f.write(search_result_header)
            f.close()

        # Statistics
        self.good_saddle_count = None
        self.bad_saddle_count = None
        self.unique_saddle_count = None

    def __repr__(self):
        return "State #%i" % self.number

    def append_search_result(self, result, comment):
        try:
            f = open(self.search_result_path, 'a')
#GH            f.write("%8d %16s %8.5f %16.5e %16d %16d %16d    %s\n" % (resultdata["search_id"], 
            resultdata = result['results']
            f.write("%8d %10s %10.5f %10.5f %10d %10d %10d    %s\n" % (result["wuid"], 
                     result["type"], 
                     resultdata["potential_energy_saddle"] - resultdata["potential_energy_reactant"],
                     resultdata["displacement_saddle_distance"],
                     resultdata["force_calls_saddle_point_concave"] + resultdata["force_calls_saddle_point_convex"],
                     resultdata["force_calls_prefactors"],
                     resultdata["force_calls"] - resultdata["force_calls_saddle_point_concave"] - resultdata["force_calls_saddle_point_convex"] - resultdata["force_calls_prefactors"],
                     comment))
            f.close()
        except:
            logger.warning("Failed to append search result %s", result['id'])
        


    def add_process(self, result):
        """ Adds a process to this State. """

        self.set_good_saddle_count(self.get_good_saddle_count() + 1) 
        
        resultdata = result["results"] #The information from the result.dat file

        # We may not already have the energy for this State.  If not, it should be in the result data.
        if self.get_energy() == None:
            self.set_energy(resultdata["potential_energy_reactant"])

        # Calculate the forward barrier for this process, and abort if the energy is too high.
        barrier = resultdata["potential_energy_saddle"] - resultdata["potential_energy_reactant"]
        lowest = self.update_lowest_barrier(barrier)
        ediff = (barrier - lowest) - (self.kT * self.max_thermal_window)
        if ediff > 0.0:
            self.append_search_result(result, "barrier > max_thermal_window")
            return None

        # Determine the number of processes in the process table that have a similar energy.
        self.load_process_table()
        energetically_close = []
        for id in self.procs.keys():
            if abs(self.procs[id]['barrier'] - barrier) < self.epsilon_e:
                energetically_close.append(id)

        # If the number of energetically similar saddles is > 0, we need to do distance checks on them.
        if len(energetically_close) > 0:
            p0 = result["saddle"]
            for id in energetically_close:
                p1 = io.loadcon(self.proc_saddle_path(id))
                for dist in atoms.per_atom_norm_gen(p1.r - p0.r, p1.box):
                    if dist > self.epsilon_r:
                        break
                else:
                    ediff = (barrier - lowest) - (self.kT * self.thermal_window)
                    if ediff < 0.0:
                        self.set_sequential_redundant(self.get_sequential_redundant() + 1)
                    self.procs[id]['repeats'] += 1
                    self.save_process_table()
                    self.append_search_result(result, "repeat-%d" % id)
                    return None

        # This appears to be a unique process.
        self.set_unique_saddle_count(self.get_unique_saddle_count() + 1)
        if barrier == lowest:
            logger.info("found new lowest barrier %f for state %i", lowest, self.number)

        
        # If this barrier is within the thermal window, reset sequential_redundant.
        ediff = (barrier - lowest) - (self.kT * self.thermal_window)        
        if ediff < 0.0:    
            self.append_search_result(result, "good-%d" % self.get_num_procs())
            self.set_sequential_redundant(0)
        else:
            self.append_search_result(result, "barrier > thermal_window")


        # The id of this process is the number of processes.
        id = self.get_num_procs()

        # Move the relevant files into the procdata directory.
        io.savecon(self.proc_saddle_path(id), result['saddle'])
        io.save_results_dat(self.proc_results_path(id), result['results']) 
        open(self.proc_reactant_path(id), 'w').writelines(result['reactant'])
        open(self.proc_mode_path(id), 'w').writelines(result['mode'])
        open(self.proc_product_path(id), 'w').writelines(result['product'])

        # Append this barrier to the process table (in memory and on disk).
        self.append_process_table(id = id, 
                                  saddle_energy =     resultdata["potential_energy_saddle"], 
                                  prefactor =         resultdata["prefactor_reactant_to_product"], 
                                  product =           -1, 
                                  product_energy =    resultdata["potential_energy_product"],
                                  product_prefactor = resultdata["prefactor_product_to_reactant"], 
                                  barrier =           barrier, 
                                  rate =              resultdata["prefactor_reactant_to_product"] * math.exp(-barrier / self.kT), 
                                  repeats =           0)

        # This was a unique process, so return True
        return id



    def load_process_table(self):
        """ Load the process table.  If the process table is not loaded, load it.  If it is loaded, do nothing. """
        if self.procs == None:
            f = open(self.proctable_path)
            lines = f.readlines()
            f.close()
            self.procs = {}
            for l in lines[1:]:
                l = l.strip().split()
                self.procs[int(l[ID])] = {"saddle_energy":     float(l[ENERGY]), 
                                          "prefactor":         float(l[PREFACTOR]), 
                                          "product":           int  (l[PRODUCT]), 
                                          "product_energy":    float(l[PRODUCT_ENERGY]), 
                                          "product_prefactor": float(l[PRODUCT_PREFACTOR]), 
                                          "barrier":           float(l[BARRIER]), 
                                          "rate":              float(l[RATE]), 
                                          "repeats":           int  (l[REPEATS])}

    
    def save_process_table(self):
        """ If the processtable is present in memory, writes it to disk. """
        if self.procs != None:
            f = open(self.proctable_path, 'w')
            f.write(processtable_header)
            for id in self.procs.keys():
                proc = self.procs[id]
                f.write(processtable_line % (id, 
                                             proc['saddle_energy'], 
                                             proc['prefactor'], 
                                             proc['product'], 
                                             proc['product_energy'],
                                             proc['product_prefactor'], 
                                             proc['barrier'],
                                             proc['rate'],
                                             proc['repeats']))
            f.close() 



    def append_process_table(self, id, saddle_energy, prefactor, product, product_energy, product_prefactor, barrier, rate, repeats):
        """ Append to the process table.  Append a single line to the process table file.  If we 
            have loaded the process table, also append it to the process table in memory. """
        f = open(self.proctable_path, 'a')
        f.write(processtable_line % (id, saddle_energy, prefactor, product, product_energy, product_prefactor, barrier, rate, repeats))
        f.close()
        if self.procs != None:
            self.procs[id] = {"saddle_energy":    saddle_energy, 
                              "prefactor":        prefactor, 
                              "product":          product, 
                              "product_energy":   product_energy, 
                              "product_prefactor":product_prefactor, 
                              "barrier":          barrier, 
                              "rate":             rate, 
                              "repeats":          repeats}
        


    def get_ratetable(self):
        """ Loads the process table if it has not been loaded and generates a rate table according to kT and thermal_window. """
        self.load_process_table()
        lowest = self.get_lowest_barrier()
        table = []
        for id in self.procs.keys():
            proc = self.procs[id]
            if proc['barrier'] < lowest + (self.kT * self.thermal_window):
                table.append((id, proc['rate']))
        return table



    def get_process_ids(self):
        """ Returns the list of ids in the rate table. """
        return [b[0] for b in self.get_ratetable()]

        

    def update_lowest_barrier(self, barrier):
        """ Compares the parameter barrier to the lowest barrier stored in info. Updates the lowest
            barrier stored in info if the barrier parameter is lower and returns the (possibly new)
            lowest barrier. """
        self.load_info()
        try:
            lowest = self.info.getfloat("MetaData", "lowest barrier")
        except:
            lowest = 1e300
        if barrier < lowest:
            lowest = barrier
            self.info.set("MetaData", "lowest barrier", "%f" % lowest)
            self.save_info()        
        return lowest

    def get_lowest_barrier(self):
        self.load_info()
        try:
            return self.info.getfloat("MetaData", "lowest barrier")
        except:
            return 1e300
    


    def get_num_procs(self):
        """ Loads the process table if it is not already loaded and returns the length of it """
        self.load_process_table()
        return len(self.procs)


    def get_unique_saddle_count(self):
        if self.unique_saddle_count is None:
            self.load_info()
            try:
                return self.info.getint("MetaData", "unique_saddles")
            except:
                return 0
        else:
            return self.unique_saddle_count


    def set_unique_saddle_count(self, num):
        self.unique_saddle_count = num
        self.load_info()
        self.info.set("MetaData", "unique_saddles", "%d" % num)
        self.save_info()        

    def get_good_saddle_count(self):
        if self.good_saddle_count is None:
            self.load_info()
            try:
                return self.info.getint("MetaData", "good_saddles")
            except:
                return 0
        else:
            return self.good_saddle_count

    def set_good_saddle_count(self, num):
        self.good_saddle_count = num
        self.load_info()
        self.info.set("MetaData", "good_saddles", "%d" % num)
        self.save_info()        

    def get_total_saddle_count(self):
        return self.get_good_saddle_count() + self.get_bad_saddle_count()
    
    def get_bad_saddle_count(self):
        if self.bad_saddle_count is None:
            self.load_info()
            try:
                return self.info.getint("MetaData", "bad_saddles")
            except:
                return 0
        else:
            return self.bad_saddle_count

    def set_bad_saddle_count(self, num):
        self.bad_saddle_count = num
        self.load_info()
        self.info.set("MetaData", "bad_saddles", "%d" % num)
        self.save_info()        
    
    def get_sequential_redundant(self):
        """ Returns the number of sequential_redundant in a row for the current state. """
        self.load_info()
        try:        
            return self.info.getint("MetaData", "sequential_redundant")
        except:
            return 0       

    def set_sequential_redundant(self, repeats):
        """ Loads the info file if it has not been loaded and sets the sequential_redundant variable. """
        self.load_info()
        self.info.set("MetaData", "sequential_redundant", "%d" % repeats)
        self.save_info()

    def get_confidence(self):
        return 1.0 - (1.0 / max(self.get_sequential_redundant(), 1.0))

    def get_energy(self):
        """ Loads the info file if it is not already loaded and returns the energy, or None
            if it is not there. """
        self.load_info()
        try:        
            return self.info.getfloat("MetaData", "reactant energy")
        except:
            return None       

    def set_energy(self, e):
        """ Loads the info file if it has not been loaded and sets the reactant energy 
            variable. """
        self.load_info()
        self.info.set("MetaData", "reactant energy", "%f" % e)
        self.save_info()        
                
    def get_reactant(self):
        """ Loads the reactant.con into a point and returns it. """
        return io.loadcon(self.reactant_path)

    def get_process_reactant(self, id):
        return io.loadcon(self.proc_reactant_path(id))

    def get_process_saddle(self, id):
        return io.loadcon(self.proc_saddle_path(id))

    def get_process_product(self, id):
        return io.loadcon(self.proc_product_path(id))

    def get_process_mode(self, id):
        return io.load_mode(self.proc_mode_path(id))

    def get_reactant(self):
        """ Loads the reactant.con into a point and returns it. """
        return io.loadcon(self.reactant_path)

    def register_bad_saddle(self, result, store = False):
        """ Registers a bad saddle. """
        result_state_code = ["Good",
                             "Init",
                             "Saddle Search No Convex Region",
                             "Saddle Search Terminated Barrier",
                             "Saddle Search Terminated Total Iterations",
                             "Saddle Search Terminated Concave Iterations",
                             "Not Connected",
                             "Bad Prefactor",
                             "Bad Barrier"]
        self.set_bad_saddle_count(self.get_bad_saddle_count() + 1)
        self.append_search_result(result, result_state_code[result["results"]["termination_reason"]])
        if store:
            if not os.path.isdir(self.bad_procdata_path):
                os.mkdir(self.bad_procdata_path)
            io.savecon(os.path.join(self.bad_procdata_path, "saddle_%d.con" % result['wuid']), result['saddle'])
            io.save_results_dat(os.path.join(self.bad_procdata_path, "results_%d.dat" % result['wuid']), result['results'])
            open(os.path.join(self.bad_procdata_path, "reactant_%d.con" % result['wuid']), 'w').writelines(result['reactant'])
            open(os.path.join(self.bad_procdata_path, "product_%d.con" % result['wuid']), 'w').writelines(result['product'])
            open(os.path.join(self.bad_procdata_path, "mode_%d.dat" % result['wuid']), 'w').writelines(result['mode'])

    def load_info(self):
        """ Loads the info file if it has not been loaded. """
        if self.info == None:
            self.info = SafeConfigParser()
            self.info.read(self.info_path)

    def save_info(self):
        """ Saves the info object if it exists. """
        if self.info == None:
            return
        self.info.write(open(self.info_path, 'w'))

    # Utility functions for compiling procdata paths, whether the files exist or not.    
    def proc_saddle_path(self, id):
        return os.path.join(self.procdata_path, "saddle_%d.con" % id)

    def proc_reactant_path(self, id):
        return os.path.join(self.procdata_path, "reactant_%d.con" % id)

    def proc_product_path(self, id):
        return os.path.join(self.procdata_path, "product_%d.con" % id)

    def proc_mode_path(self, id):
        return os.path.join(self.procdata_path, "mode_%d.dat" % id)

    def proc_results_path(self, id):
        return os.path.join(self.procdata_path, "results_%d.dat" % id)
        
    def get_process_table(self):
        self.load_process_table()
        return self.procs

if __name__ == "__main__":
    pass
