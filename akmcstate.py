""" The state module. """

import os
import shutil
import math
import logging
logger = logging.getLogger('state')
from ConfigParser import SafeConfigParser 
import StringIO

import numpy

import atoms
import io
import state



ID, ENERGY, PREFACTOR, PRODUCT, PRODUCT_ENERGY, PRODUCT_PREFACTOR, BARRIER, RATE, REPEATS = range(9)

processtable_head_fmt = "%7s %16s %11s %9s %16s %17s %8s %12s %7s\n"

processtable_header = processtable_head_fmt % ("proc #", "saddle energy", "prefactor", 
                                               "product", "product energy", "product prefactor",
                                               "barrier", "rate", "repeats")

processtable_line = "%7d %16.5f %11.5e %9d %16.5f %17.5e %8.5f %12.5e %7d\n"

search_result_header = "%8s %10s %10s %10s %10s %10s %10s    %s\n" % ("wuid", "type", "barrier",
                                                                      "max-dist", "sad-fcs", 
                                                                      "mins-fcs", "pref-fcs", 
                                                                      "result")
search_result_header += "-" * len(search_result_header) + '\n'



class AKMCState(state.State):
    def __init__(self, statepath, statenumber, statelist, previous_state_num = -1, 
                 reactant_path = None):                 
        """ Creates a new State, with lazily loaded data. """
        state.State.__init__(self,statepath, statenumber,statelist, previous_state_num,
                    reactant_path)

        self.info_path = os.path.join(self.path, "info")
        self.procdata_path = os.path.join(self.path, "procdata")
        self.bad_procdata_path = os.path.join(self.path, "badprocdata")
        self.reactant_path = os.path.join(self.path, "reactant.con")
        self.proctable_path = os.path.join(self.path, "processtable")
        self.search_result_path = os.path.join(self.path, "search_results.txt")

        # Lazily loaded data. Should use the get/set methods for these.
        self.info = None
        self.procs = None
        self.proc_repeat_count = None

        # If this state does not exist on disk, create it.
        if not os.path.isdir(self.path):
            if reactant_path == None:
                raise IOError("State needs a reactant_path when it is being instantiated to disk.")
            os.mkdir(self.path)
            os.mkdir(self.procdata_path)
            shutil.copy(reactant_path, self.reactant_path)
            config = SafeConfigParser()
            config.add_section("MetaData")
            config.set("MetaData", "previous state", str(previous_state_num))            
            config.write(open(self.info_path, 'w'))
            f = open(self.proctable_path, 'w')
            f.write(processtable_header)
            f.close()
            if self.statelist.list_search_results:
                f = open(self.search_result_path, 'w')
                f.write(search_result_header)
                f.close()

        # Statistics
        self.good_saddle_count = None
        self.bad_saddle_count = None
        self.unique_saddle_count = None

    def add_process(self, result):
        """ Adds a process to this State. """
        
        self.set_good_saddle_count(self.get_good_saddle_count() + 1) 
         
        resultdata = result["results"] #The information from the result.dat file
        
        # We may not already have the energy for this State.  If not, it should be in the result data.
        if self.get_energy() == None:
            self.set_energy(resultdata["potential_energy_reactant"])

        # Calculate the forward barrier for this process, and abort if the energy is too high.
        oldlowest = self.get_lowest_barrier()
        barrier = resultdata["potential_energy_saddle"] - resultdata["potential_energy_reactant"]
        lowest = self.update_lowest_barrier(barrier)
        ediff = (barrier - lowest) - (self.statelist.kT * self.statelist.max_thermal_window)
        if ediff > 0.0:
            self.append_search_result(result, "barrier > max_thermal_window")
            return None

        # Determine the number of processes in the process table that have a similar energy.
        self.load_process_table()
        energetically_close = []
        for id in self.procs.keys():
            if abs(self.procs[id]['barrier'] - barrier) < self.statelist.epsilon_e:
                energetically_close.append(id)

        # If the number of energetically similar saddles is > 0, we need to do distance checks on them.
        if len(energetically_close) > 0:
            #load the saddle
            result["saddle"] = io.loadcon(result["saddle.con"])
            p0 = result["saddle"]
            ibox = numpy.linalg.inv(p0.box)
            for id in energetically_close:
                p1 = io.loadcon(self.proc_saddle_path(id))
                
                for dist in atoms.per_atom_norm_gen(p1.free_r() - p0.free_r(), p1.box, ibox):
                    if dist > self.statelist.epsilon_r:
                        break
                else:
                    self.procs[id]['repeats'] += 1
                    self.save_process_table()
                    self.append_search_result(result, "repeat-%d" % id)
                    self.inc_proc_repeat_count(id)
                    return None

        # This appears to be a unique process.
        # Check if the mode, reactant, saddle, and product are legit
        try:
            if 'mode' not in result:
                io.load_mode(result['mode.dat'])
            if 'reactant' not in result:
                io.loadcon(result['reactant.con'])
            if 'saddle' not in result:
                io.loadcon(result['saddle.con'])
            if 'product' not in result:
                io.loadcon(result['product.con'])
        except:
            logger.exception("Mode, reactant, saddle, or product has incorrect format")
            return None
        
        self.set_unique_saddle_count(self.get_unique_saddle_count() + 1)
        if barrier == lowest and barrier < oldlowest - self.statelist.epsilon_e:
            logger.info("found new lowest barrier %f for state %i", lowest, self.number)
        
        # Update the search result table.
        self.append_search_result(result, "good-%d" % self.get_num_procs())

        # The id of this process is the number of processes.
        id = self.get_num_procs()

        # Keep track of the number of searches, Ns.
        self.inc_proc_repeat_count(id)

        # Move the relevant files into the procdata directory.
        open(self.proc_reactant_path(id), 'w').writelines(result['reactant.con'].getvalue())
        open(self.proc_mode_path(id), 'w').writelines(result['mode.dat'].getvalue())
        open(self.proc_product_path(id), 'w').writelines(result['product.con'].getvalue())
        open(self.proc_saddle_path(id), 'w').writelines(result['saddle.con'].getvalue())
        open(self.proc_results_path(id), 'w').writelines(result['results.dat'].getvalue())

        # Append this barrier to the process table (in memory and on disk).
        self.append_process_table(id =                id, 
                                  saddle_energy =     resultdata["potential_energy_saddle"], 
                                  prefactor =         resultdata["prefactor_reactant_to_product"], 
                                  product =           -1, 
                                  product_energy =    resultdata["potential_energy_product"],
                                  product_prefactor = resultdata["prefactor_product_to_reactant"], 
                                  barrier =           barrier, 
                                  rate =              resultdata["prefactor_reactant_to_product"] * math.exp(-barrier / self.statelist.kT), 
                                  repeats =           0)

        # This was a unique process, so return the id.
        return id



    def append_search_result(self, result, comment):
        if self.statelist.list_search_results:
            try:
                f = open(self.search_result_path, 'a')
                resultdata = result['results']
                f.write("%8d %10s %10.5f %10.5f %10d %10d %10d    %s\n" % (result["wuid"], 
                         result["type"], 
                         resultdata["potential_energy_saddle"] - resultdata["potential_energy_reactant"],
                         resultdata["displacement_saddle_distance"],
                         resultdata["force_calls_saddle"] ,
                         resultdata["force_calls_minimization"] ,
                         resultdata["force_calls_prefactors"],
                         comment))
                f.close()
            except:
                logger.warning("Failed to append search result.")
        


    def get_ratetable(self):
        """ Loads the process table if it has not been loaded and generates a rate table 
            according to kT and thermal_window. """
        self.load_process_table()
        lowest = self.get_lowest_barrier()
        table = []
        for id in self.procs.keys():
            proc = self.procs[id]
            if proc['barrier'] < lowest + (self.statelist.kT * self.statelist.thermal_window):
                table.append((id, proc['rate']))
        return table
        
        

    def get_confidence(self):
        """ The confidence is a function of the ratio Nf/Ns, where Nf is the number of unique 
            processes and Ns is the number of searches performed.  
            
            When Nf or Ns are zero, there is no confidence.  Zero is returned in this case. 
            
            As the ratio Nf/Ns decreases, the confidence increases from 0.0 to a limit of 1.0.
            This is roughly equivalent to the statement: I am more confident that I have found
            all relevant processes if I have found 10 processes after 1000 searches than if I 
            have found 900 processes after 1000 searches.
            
            Nf is calculated to be the number of processes in the rate table. Ns is calculated
            to be the number of searches that resulted in a process on the rate table.
            
            When using recycling or kdb, it is useful to ignore processes that occur outside
            the hole, or the region in which the last process took place.  Focusing on this
            region means you can do possibly far fewer searches to reach confidence. When using
            the hole to filter processes, Nf and Ns only take into account processes that 
            intersect the hole. """
        prc = self.get_proc_repeat_count()
        rt = self.get_ratetable()
        Nf = 0.0
        Ns = 0.0
        for r in rt:
            if self.statelist.filter_hole:
                if self.proc_in_hole(r[0]):
                    Ns += prc[r[0]]
                    Nf += 1
            else:
                Ns += prc[r[0]]
                Nf += 1
        if Nf < 1 or Ns < 1:
            return 0.0
        return 1.0 + (Nf/Ns) * lambertw(-math.exp(-1.0 / (Nf/Ns))/(Nf/Ns))



    def proc_in_hole(self, procid):
        """ Returns True if the given process intersects the previous state/current state hole,
            False if not. """
        if self.get_previous_state() == -1:
            return True
        if procid in self.get_procs_in_hole():
            return True
        if procid in self.get_procs_not_in_hole():
            return False
        ps = self.statelist.get_state(self.get_previous_state()).get_reactant()
        cs = self.get_reactant()
        pan = atoms.per_atom_norm(ps.r - cs.r, ps.box)
        hole_atoms = []
        for i in range(len(pan)):
            if pan[i] > 0.2:
                hole_atoms.append(i)
        pp = self.get_process_product(procid)
        pan = atoms.per_atom_norm(pp.r - cs.r, cs.box)
        for i in range(len(pan)):
            if pan[i] > 0.2:
                if i in hole_atoms:
                    pih = self.get_procs_in_hole()
                    pih.append(procid)
                    self.set_procs_in_hole(pih)
                    return True
        pnih = self.get_procs_not_in_hole()
        pnih.append(procid)
        self.set_procs_not_in_hole(pnih)
        return False
                
        

    def get_procs_in_hole(self):
        self.load_info()
        try:
            return eval(self.info.get("MetaData", "procs in hole"))
        except:
            return []



    def get_procs_not_in_hole(self):
        self.load_info()
        try:
            return eval(self.info.get("MetaData", "procs not in hole"))
        except:
            return []



    def set_procs_in_hole(self, procs):
        self.load_info()
        self.info.set("MetaData", "procs in hole", repr(procs))
        self.save_info()



    def set_procs_not_in_hole(self, procs):
        self.load_info()
        self.info.set("MetaData", "procs not in hole", repr(procs))
        self.save_info()



    def load_process_table(self):
        """ Load the process table.  If the process table is not loaded, load it.  If it is 
            loaded, do nothing. """
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
                f.write(processtable_line % (id, proc['saddle_energy'], proc['prefactor'], 
                                             proc['product'], proc['product_energy'],
                                             proc['product_prefactor'], proc['barrier'],
                                             proc['rate'], proc['repeats']))
            f.close() 



    def append_process_table(self, id, saddle_energy, prefactor, product, product_energy, 
                             product_prefactor, barrier, rate, repeats):
        """ Append to the process table.  Append a single line to the process table file.  If we 
            have loaded the process table, also append it to the process table in memory. """
        f = open(self.proctable_path, 'a')
        f.write(processtable_line % (id, saddle_energy, prefactor, product, product_energy, 
                                     product_prefactor, barrier, rate, repeats))
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
    
    def get_process(self, id):
        self.load_process_table()
        return self.procs[id]        


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
        
        

    def get_previous_state(self):
        self.load_info()
        try:
            return self.info.getint("MetaData", "previous state")
        except:
            return -1



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



    def get_proc_repeat_count(self):
        self.load_info()
        if self.proc_repeat_count is None:
            try:
                self.proc_repeat_count = eval(self.info.get("MetaData", "proc repeat count"))
            except:
                self.proc_repeat_count = [0]
        return self.proc_repeat_count



    def inc_proc_repeat_count(self, procid):
        self.get_proc_repeat_count()
        if procid == len(self.proc_repeat_count):
            self.proc_repeat_count.append(1)
        elif procid < len(self.proc_repeat_count) and procid >= 0:
            self.proc_repeat_count[procid] += 1                
        self.info.set("MetaData", "proc repeat count", repr(self.proc_repeat_count))
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
                             "Bad Barrier",
                             "Minimum Not Converged"]
        self.set_bad_saddle_count(self.get_bad_saddle_count() + 1)
        self.append_search_result(result, result_state_code[result["results"]["termination_reason"]])
        if store:
            if not os.path.isdir(self.bad_procdata_path):
                os.mkdir(self.bad_procdata_path)
            open(os.path.join(self.bad_procdata_path, "reactant_%d.con" % result['wuid']), 'w').writelines(result['reactant.con'].getvalue())
            open(os.path.join(self.bad_procdata_path, "product_%d.con" % result['wuid']), 'w').writelines(result['product.con'].getvalue())
            open(os.path.join(self.bad_procdata_path, "mode_%d.dat" % result['wuid']), 'w').writelines(result['mode.dat'].getvalue())
            open(os.path.join(self.bad_procdata_path, "results_%d.dat" % result['wuid']), 'w').writelines(result['results.dat'].getvalue())
            open(os.path.join(self.bad_procdata_path, "saddle_%d.con" % result['wuid']), 'w').writelines(result['saddle.con'].getvalue())

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

    # Utility functions for loading process .con and mode files.
    def get_process_saddle(self, id):
        return io.loadcon(self.proc_saddle_path(id))
    def get_process_mode(self, id):
        return io.load_mode(self.proc_mode_path(id))

    # Utility functions for compiling procdata paths, whether the files exist or not.    
    def proc_saddle_path(self, id):
        return os.path.join(self.procdata_path, "saddle_%d.con" % id)
    def proc_mode_path(self, id):
        return os.path.join(self.procdata_path, "mode_%d.dat" % id)
    def get_process_table(self):
        self.load_process_table()
        return self.procs

# Lambert-W function for calculating 
# Source: http://keithbriggs.info/software/LambertW.py
from math import log, sqrt, exp

def lambertw(z):
    """Lambert W function, principal branch"""
    eps = 1.0e-12
    em1 = 0.3678794411714423215955237701614608
    if z < -em1:
        logger.error("Tried to evaluate Lambert W function @ < -1/e")
        raise ValueError()
    if 0.0 == z: 
        return 0.0
    if z < -em1 + 1e-4:
        q = z + em1
        r = sqrt(q)
        q2 = q * q
        q3 = q2 * q
        return\
         -1.0\
         +2.331643981597124203363536062168 * r\
         -1.812187885639363490240191647568 * q\
         +1.936631114492359755363277457668 * r * q\
         -2.353551201881614516821543561516 * q2\
         +3.066858901050631912893148922704 * r * q2\
         -4.175335600258177138854984177460 * q3\
         +5.858023729874774148815053846119 * r * q3\
         -8.401032217523977370984161688514 * q3 * q
    if z < 1.0:
        p = sqrt(2.0 * (2.7182818284590452353602874713526625 * z + 1.0))
        w = -1.0 + p * (1.0 + p * (-0.333333333333333333333 + p * 0.152777777777777777777777))
    else:
        w = log(z)
    if z > 3.0: w-=log(w)
    for i in xrange(10):
        e = exp(w)
        t = w * e - z
        p = w + 1.0
        t /= e * p - 0.5 * (p + 1.0) * t / p
        w -= t
        if abs(t) < eps * (1.0 + abs(w)): 
            return w
    logger.error("Failed to converge Lambert W function.")
    raise ValueError()

if __name__ == "__main__":
    pass
