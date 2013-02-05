##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

""" The state module. """

import os
import math
import logging
logger = logging.getLogger('state')

import numpy

import atoms
import config
import fileio as io
import state

class AKMCState(state.State):
    ID, ENERGY, PREFACTOR, PRODUCT, PRODUCT_ENERGY, PRODUCT_PREFACTOR, BARRIER, RATE, REPEATS = range(9)
    processtable_head_fmt = "%7s %16s %11s %9s %16s %17s %8s %12s %7s\n"
    processtable_header = processtable_head_fmt % ("proc #", "saddle energy", "prefactor", 
                                                   "product", "product energy", "product prefactor",
                                                   "barrier", "rate", "repeats")
    processtable_line = "%7d %16.5f %11.5e %9d %16.5f %17.5e %8.5f %12.5e %7d\n"

    def __init__(self, statepath, statenumber, statelist, previous_state_num = -1,
                 reactant_path = None):
        """ Creates a new State, with lazily loaded data. """
        if config.akmc_server_side_process_search:
            self.search_result_header = "%8s %10s %10s %10s %10s %10s %10s    %s\n" % ("searchid", "type", "barrier",
                                                                                  "max-dist", "sad-fcs",
                                                                                  "mins-fcs", "pref-fcs",
                                                                                  "result")
        else:
            self.search_result_header = "%8s %10s %10s %10s %10s %10s %10s    %s\n" % ("wuid", "type", "barrier",
                                                                                  "max-dist", "sad-fcs",
                                                                                  "mins-fcs", "pref-fcs",
                                                                                  "result")
        self.search_result_header += "-" * len(self.search_result_header) + '\n'
        state.State.__init__(self,statepath, statenumber,statelist, previous_state_num,
                    reactant_path)

        self.bad_procdata_path = os.path.join(self.path, "badprocdata")

        self.con_cache = {}
        

    def find_repeat(self, saddle_file, barrier):
        self.load_process_table()
        energy_a = barrier
        p1 = io.loadcon(saddle_file)
        for id in self.procs.keys():
            energy_b = self.procs[id]['barrier']
            if abs(energy_a - energy_b) > config.comp_eps_e:
                continue

            if id in self.con_cache:
                p2 = self.con_cache[id]
            else:
                p2 = io.loadcon(self.proc_saddle_path(id)) 
                self.con_cache[id] = p2

            if atoms.match(p1, p2, config.comp_eps_r, config.comp_neighbor_cutoff, False):
                return id
        return None

    def add_process(self, result):
        """ Adds a process to this State. """
        state.State.add_process(self, result)

        self.set_good_saddle_count(self.get_good_saddle_count() + 1)

        resultdata = result["results"] #The information from the result.dat file

        if 'simulation_time' in resultdata:
            self.increment_time(resultdata['simulation_time'])

        # We may not already have the energy for this State.  If not, it should be placed in the result data.
        if self.get_energy() == None:
            # This energy now defines the reference energy for the state
            self.set_energy(resultdata["potential_energy_reactant"])
        reactant_energy = self.get_energy()

        # Calculate the forward barrier for this process, and abort if the energy is too high.
        oldlowest = self.get_lowest_barrier()
        barrier = resultdata["potential_energy_saddle"] - reactant_energy

        lowest = self.update_lowest_barrier(barrier)
        ediff = (barrier - lowest) - (self.statelist.kT *
                (self.statelist.thermal_window+self.statelist.max_thermal_window))
        if ediff > 0.0:
            self.append_search_result(result, "barrier > max_thermal_window")
            return None

        # Determine the number of processes in the process table that have a similar energy.
        id = self.find_repeat(result["saddle.con"], barrier)
        if id != None:
            self.append_search_result(result, "repeat-%d" % id)
            self.procs[id]['repeats'] += 1
            self.save_process_table()
            if result['type'] == "random" or result['type'] == "dynamics":
                self.inc_proc_random_count(id)
                if id in self.get_relevant_procids():
                    self.inc_repeats()
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

        # Reset the repeat count.
        self.reset_repeats()

        # Respond to finding a new lowest barrier.
        self.set_unique_saddle_count(self.get_unique_saddle_count() + 1)
        if barrier == lowest and barrier < oldlowest - self.statelist.epsilon_e:
            logger.info("found new lowest barrier %f for state %i", lowest, self.number)
        elif config.saddle_method == 'dynamics':
            logger.info("found new barrier %f for state %i", barrier, self.number)


        # Update the search result table.
        self.append_search_result(result, "good-%d" % self.get_num_procs())

        # The id of this process is the number of processes.
        id = self.get_num_procs()

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

        # If this is a random search type, add this proc to the random proc dict.
        if result['type'] == "random" or result['type'] == "dynamics":
            self.inc_proc_random_count(id)

        # This was a unique process, so return the id.
        return id

    def append_search_result(self, result, comment):
        #try:
        f = open(self.search_result_path, 'a')
        resultdata = result['results']

        if config.akmc_server_side_process_search:
            first_column = "search_id"
        else:
            first_column = "wuid"
        f.write("%8d %10s %10.5f %10.5f %10d %10d %10d    %s\n" % (result[first_column], 
                 result["type"], 
                 resultdata["barrier_reactant_to_product"],
                 resultdata["displacement_saddle_distance"],
                 resultdata["force_calls_saddle"],
                 resultdata["force_calls_minimization"],
                 resultdata["force_calls_prefactors"],
                 comment))
        f.close()
        #except:
        #    logger.warning("Failed to append search result.")

    def get_ratetable(self):
        """ Loads the process table if it has not been loaded and generates a rate table 
            according to kT and thermal_window. """
        self.load_process_table()
        lowest = self.get_lowest_barrier()
        table = []

        for id in self.procs.keys():
            proc = self.procs[id]
            if proc['barrier'] > lowest + (self.statelist.kT * self.statelist.thermal_window):
                continue
            table.append((id, proc['rate'], proc['prefactor']))
        return table
 
 
    def get_relevant_procids(self):
        rt = self.get_ratetable()
        rps = []
        for r in rt:
            rps.append(r[0])
        return rps


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
            the hole, the region in which the last process took place.  Focusing on this
            region means you can do possibly far fewer searches to reach confidence. When using
            the hole to filter processes, Nf and Ns only take into account processes that 
            intersect the hole. """
        alpha = 1.0
        if config.akmc_confidence_correction:
            rt = self.get_ratetable(False)
            prc = self.get_proc_random_count()
            mn = 1e300
            mx = 0
            for r in rt:
                if r[0] in prc:
                    mn = min(mn, prc[r[0]])
                    mx = max(mx, prc[r[0]])
            if mx < 1:
                alpha = 1.0
            else:
                alpha = float(mn)/mx
        if config.akmc_confidence_scheme == "new":
            rt = self.get_ratetable(False)
            prc = self.get_proc_random_count()
            Nf = 0.0
            Ns = 0.0
            for r in rt:
                if r[0] in prc:
                    Nf += 1
                    Ns += prc[r[0]]
            if Ns < 1:
                return 0.0
            if Nf < 1:
                Nf = 1.0
            return 1.0 + (Nf/(alpha*Ns)) * lambertw(-math.exp(-1.0 / (Nf/(alpha*Ns)))/(Nf/(alpha*Ns)))
        elif config.akmc_confidence_scheme == 'sampling':
            all_repeats = self.get_proc_random_count()
            rt = self.get_ratetable()
            repeats = {}
            for event in rt:
                id = event[0]
                repeats[id] = all_repeats[id]
            #number of events
            m = len(repeats)
            #number of searches
            n = sum(repeats.values())
            if n < 10: return 0.0

            #probabilities
            ps = numpy.array(repeats.values(),dtype=numpy.float)
            ps /= sum(ps)

            C = numpy.zeros(m)
            for i in xrange(n):
                C += (1.0-C)*ps

            return sum(C)/float(m)

        elif config.akmc_confidence_scheme == 'dynamics':
            if self.get_time() == 0.0: return 0.0
            rt = self.get_ratetable()
            T1 = config.main_temperature
            T2 = config.saddle_dynamics_temperature

            #rates are at T1
            rates = numpy.array([ p[1] for p in rt ])
            prefactors = numpy.array([ p[2] for p in rt ])
            if len(rates) == 0: return 0.0
            #extrapolate to T2
            rates_md = prefactors*(rates/prefactors)**(T1/T2)
            
            time = self.get_time()*1e-15
            C = 1.0-numpy.exp(-time*rates_md*1.0)
            total_rate = sum(rates)

            conf = sum(C*rates)/total_rate
            return conf

        else:
            Nr = self.get_repeats()
            if Nr < 1:
                return 0.0
            else:
                return max(0.0, 1.0 - 1.0/(alpha*Nr))

    def get_proc_random_count(self):
        return eval(self.info.get("MetaData", "proc repeat count", "{}"))

    def inc_proc_random_count(self, procid):
        prc = self.get_proc_random_count()
        if procid not in prc:
            prc[procid] = 1
        else:
            prc[procid] += 1
        self.info.set("MetaData", "proc repeat count", repr(prc))

    def reset_repeats(self):
        self.info.set("MetaData", "repeats", 0)

    def get_repeats(self):
        return self.info.get("MetaData", "repeats", 0)

    def inc_repeats(self):
        self.info.set("MetaData", "repeats", self.get_repeats() + 1)

    def load_process_table(self):
        """ Load the process table.  If the process table is not loaded, load it.  If it is
            loaded, do nothing. """
        if self.procs != None:
            return
        f = open(self.proctable_path)
        lines = f.readlines()
        f.close()
        self.procs = {}
        for l in lines[1:]:
            l = l.strip().split()
            self.procs[int(l[self.ID])] = {"saddle_energy":     float(l[self.ENERGY]),
                                           "prefactor":         float(l[self.PREFACTOR]),
                                           "product":           int  (l[self.PRODUCT]),
                                           "product_energy":    float(l[self.PRODUCT_ENERGY]),
                                           "product_prefactor": float(l[self.PRODUCT_PREFACTOR]),
                                           "barrier":           float(l[self.BARRIER]),
                                           "rate":              float(l[self.RATE]),
                                           "repeats":           int  (l[self.REPEATS])}
        try:
            kT = self.info.get('MetaData', 'kT')
        except NameError:
            self.info.set('MetaData', 'kT', self.statelist.kT)
            return
        if kT != self.statelist.kT:
            for id, proc in self.procs.items():
                proc['rate'] = proc['prefactor'] * math.exp(-proc['barrier'] / self.statelist.kT)
            self.save_process_table()            
                

    def save_process_table(self):
        """ If the processtable is present in memory, writes it to disk. """
        if self.procs != None:
            f = open(self.proctable_path, 'w')
            f.write(self.processtable_header)
            for id in self.procs.keys():
                proc = self.procs[id]
                f.write(self.processtable_line % (id, proc['saddle_energy'], proc['prefactor'],
                                                  proc['product'], proc['product_energy'],
                                                  proc['product_prefactor'], proc['barrier'],
                                                  proc['rate'], proc['repeats']))
            f.close() 


    def append_process_table(self, id, saddle_energy, prefactor, product, product_energy,
                             product_prefactor, barrier, rate, repeats):
        """ Append to the process table.  Append a single line to the process table file.  If we
            have loaded the process table, also append it to the process table in memory. """
        f = open(self.proctable_path, 'a')
        f.write(self.processtable_line % (id, saddle_energy, prefactor, product, product_energy,
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


    def update_lowest_barrier(self, barrier):
        """ Compares the parameter barrier to the lowest barrier stored in info. Updates the lowest
            barrier stored in info if the barrier parameter is lower and returns the (possibly new)
            lowest barrier. """
        lowest = self.info.get("MetaData", "lowest barrier", 1e300)
        if barrier < lowest:
            lowest = barrier
            self.info.set("MetaData", "lowest barrier", lowest)
        return lowest

    def get_lowest_barrier(self):
        return self.info.get("MetaData", "lowest barrier", 1e300)

    def get_unique_saddle_count(self):
        return self.info.get("MetaData", "unique_saddles", 0)

    def set_unique_saddle_count(self, num):
        self.info.set("MetaData", "unique_saddles", num)

    def get_good_saddle_count(self):
        return self.info.get("MetaData", "good_saddles", 0)

    def set_good_saddle_count(self, num):
        self.info.set("MetaData", "good_saddles", num)

    def increment_time(self, dt):
        self.info.set("MetaData", "time", self.get_time() + dt)

    def get_time(self):
        return self.info.get("MetaData", "time", 0.0)

    def get_total_saddle_count(self):
        return self.get_good_saddle_count() + self.get_bad_saddle_count()

    def get_bad_saddle_count(self):
        return self.info.get("MetaData", "bad_saddles", 0)

    def set_bad_saddle_count(self, num):
        self.info.set("MetaData", "bad_saddles", num)

    def register_bad_saddle(self, result, store = False):
        """ Registers a bad saddle. """
        result_state_code = ["Good",
                             "Init",
                             "Saddle Search No Convex Region",
                             "Saddle Search Terminated High Energy",
                             "Saddle Search Terminated Concave Iterations",
                             "Saddle Search Terminated Total Iterations",
                             "Not Connected",
                             "Bad Prefactor",
                             "Bad Barrier",
                             "Minimum Not Converged",
                             "Failed Prefactor Calculation",
                             "Potential Failed",
                             "Nonnegative Displacement Abort",
                             "Nonlocal abort",
                             "Negative barrier"]
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


    # Utility functions for loading process .con and mode files.
    def get_process_reactant(self, id):
        return io.loadcon(self.proc_reactant_path(id))
    def get_process_saddle(self, id):
        return io.loadcon(self.proc_saddle_path(id))
    def get_process_product(self, id):
        return io.loadcon(self.proc_product_path(id))
    def get_process_mode(self, id):
        return io.load_mode(self.proc_mode_path(id))


    # Utility functions for compiling procdata paths, whether the files exist or not.
    def proc_saddle_path(self, id):
        return os.path.join(self.procdata_path, "saddle_%d.con" % id)
    def proc_mode_path(self, id):
        return os.path.join(self.procdata_path, "mode_%d.dat" % id)


# Lambert-W function: http://keithbriggs.info/software/LambertW.py
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
        r = math.sqrt(q)
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
        p = math.sqrt(2.0 * (2.7182818284590452353602874713526625 * z + 1.0))
        w = -1.0 + p * (1.0 + p * (-0.333333333333333333333 + p * 0.152777777777777777777777))
    else:
        w = math.log(z)
    if z > 3.0: w-=math.log(w)
    for i in xrange(10):
        e = math.exp(w)
        t = w * e - z
        p = w + 1.0
        t /= e * p - 0.5 * (p + 1.0) * t / p
        w -= t
        if abs(t) < eps * (1.0 + abs(w)): 
            return w
    logger.error("Failed to converge Lambert W function.")
    raise ValueError()
