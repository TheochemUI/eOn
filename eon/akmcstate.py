
""" The state module. """

import os
import math
import configparser
import logging
logger = logging.getLogger('state')

import numpy

from eon.config import config as EON_CONFIG
from eon.config import ConfigClass # Typing
from eon import atoms
from eon import fileio as io
from eon import state

class AKMCState(state.State):
    ID, ENERGY, PREFACTOR, PRODUCT, PRODUCT_ENERGY, PRODUCT_PREFACTOR, BARRIER, RATE, REPEATS = list(range(9))
    processtable_head_fmt = "%7s %16s %11s %9s %16s %17s %8s %12s %7s\n"
    processtable_header = processtable_head_fmt % ("proc #", "saddle energy", "prefactor",
                                                   "product", "product energy", "product prefactor",
                                                   "barrier", "rate", "repeats")
    processtable_line = "%7d %16.5f %11.5e %9d %16.5f %17.5e %8.5f %12.5e %7d\n"

    def __init__(self, statepath, statenumber, statelist, previous_state_num = -1,
                 reactant_path = None, config: ConfigClass = EON_CONFIG):
        """ Creates a new State, with lazily loaded data. """
        self.config = config
        if self.config.akmc_server_side_process_search:
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
                             reactant_path, self.config)

        self.bad_procdata_path = os.path.join(self.path, "badprocdata")

        self.con_cache = {}


    def find_repeat(self, saddle_file, barrier):
        self.load_process_table()
        energy_a = barrier
        p1 = io.loadcon(saddle_file)
        for id in list(self.procs.keys()):
            energy_b = self.procs[id]['barrier']
            if abs(energy_a - energy_b) > self.config.comp_eps_e:
                continue

            if id in self.con_cache:
                p2 = self.con_cache[id]
            else:
                p2 = io.loadcon(self.proc_saddle_path(id))
                self.con_cache[id] = p2

            if atoms.match(p1, p2, self.config.comp_eps_r, self.config.comp_neighbor_cutoff, False):
                return id
        return None

    def add_process(self, result, superbasin=None):
        """ Adds a process to this State. """
        state.State.add_process(self, result)

        self.set_good_saddle_count(self.get_good_saddle_count() + 1)

        resultdata = result["results"] #The information from the result.dat file

        if 'simulation_time' in resultdata:
            self.increment_time(resultdata['simulation_time'], resultdata['md_temperature'])

        # We may not already have the energy for this State.  If not, it should be placed in the result data.
        if self.get_energy() is None:
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
            self.append_search_result(result, "barrier > max_thermal_window", superbasin)
            return None

        # Determine the number of processes in the process table that have a similar energy.
        id = self.find_repeat(result["saddle.con"], barrier)
        if id != None:
            self.append_search_result(result, "repeat-%d" % id, superbasin)
            self.procs[id]['repeats'] += 1
            self.save_process_table()
            if result['type'] == "random" or result['type'] == "dynamics":
                self.inc_proc_random_count(id)
                # Do not increase repeats if we are currently in a
                # superbasin and the process does not lead out of it;
                # or if the process barrier is outside the thermal window.
                if id in self.get_relevant_procids(superbasin):
                    self.inc_repeats()
            if 'simulation_time' in resultdata:
                current_time = self.get_time()
                logger.debug("event %3i found at time %f fs" % (id, current_time))
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
            logger.info("Found new lowest barrier %f for state %i (type: %s)", lowest, self.number, result['type'])
        logger.info("Found new barrier %f for state %i (type: %s)", barrier, self.number, result['type'])


        # Update the search result table.
        self.append_search_result(result, "good-%d" % self.get_num_procs(), superbasin)

        # The id of this process is the number of processes.
        id = self.get_num_procs()

        if 'simulation_time' in resultdata:
            current_time = self.get_time()
            logger.debug("new event %3i found at time %f fs" % (id, current_time))

        # Move the relevant files into the procdata directory.
        open(self.proc_reactant_path(id), 'w').writelines(result['reactant.con'].getvalue())
        open(self.proc_mode_path(id), 'w').writelines(result['mode.dat'].getvalue())
        open(self.proc_product_path(id), 'w').writelines(result['product.con'].getvalue())
        open(self.proc_saddle_path(id), 'w').writelines(result['saddle.con'].getvalue())
        open(self.proc_results_path(id), 'w').writelines(result['results.dat'].getvalue())

        # Set maximum rate, if defined
#        forward_rate = resultdata["prefactor_reactant_to_product"] * math.exp(-barrier / self.statelist.kT)
#        if self.config.akmc_max_rate > 0 and cur_rate > self.config.akmc_max_rate:
#            # print "max rate exceeded: ", cur_rate
#            forward_rate = self.config.akmc_max_rate

        # Set equilibrium rate, if defined
        forward_rate = resultdata["prefactor_reactant_to_product"] * math.exp(-barrier / self.statelist.kT)
        reverse_barrier = barrier - (resultdata["potential_energy_product"] - reactant_energy)
        reverse_rate = resultdata["prefactor_product_to_reactant"] * math.exp(-reverse_barrier / self.statelist.kT)

        eq_rate_flag = False
        if self.config.akmc_eq_rate > 0 and forward_rate > self.config.akmc_eq_rate and reverse_rate > self.config.akmc_eq_rate:
            eq_rate_flag = True
            #print "eq_rate exceeded, forward:", forward_rate, " reverse: ", reverse_rate
            if forward_rate < reverse_rate:
                forward_eq_rate = self.config.akmc_eq_rate
                reverse_eq_rate = self.config.akmc_eq_rate * (reverse_rate / forward_rate)
            else:
                forward_eq_rate = self.config.akmc_eq_rate * (forward_rate / reverse_rate)
                reverse_eq_rate = self.config.akmc_eq_rate
            #print "new eq forward rate:", forward_eq_rate, " reverse: ", reverse_eq_rate

        # Append this barrier to the process table (in memory and on disk).
        self.append_process_table(id =                id,
                                  saddle_energy =     resultdata["potential_energy_saddle"],
                                  prefactor =         resultdata["prefactor_reactant_to_product"],
                                  product =           -1,
                                  product_energy =    resultdata["potential_energy_product"],
                                  product_prefactor = resultdata["prefactor_product_to_reactant"],
                                  barrier =           barrier,
                                  rate =              resultdata["prefactor_reactant_to_product"] * math.exp(-barrier / self.statelist.kT),
                                  # rate =              forward_rate,
                                  repeats =           0)

        # If equilibrium rate, change the forward rate as well
        if eq_rate_flag:
            self.procs[id]['rate'] = forward_eq_rate
            # reverse_procs[id]['rate'] = reverse_eq_rate

        # If this is a random search type, add this proc to the random proc dict.
        if result['type'] == "random" or result['type'] == "dynamics":
            self.inc_proc_random_count(id)

        # This was a unique process, so return the id.
        return id

    def append_search_result(self, result, comment, superbasin):
        #try:
        f = open(self.search_result_path, 'a')
        resultdata = result['results']

        if self.config.akmc_server_side_process_search:
            first_column = "search_id"
        else:
            first_column = "wuid"

        if superbasin:
            comment += " [%i]" % superbasin.id

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

    def get_ratetable(self, superbasin=None):
        """ Loads the process table if it has not been loaded and generates a rate table
            according to kT and thermal_window. """
        self.load_process_table()
        lowest = self.get_lowest_barrier()
        # Maximum barrier according to thermal window.
        max_barrier = lowest + self.statelist.kT * self.statelist.thermal_window
        rt = [(id, proc['rate'], proc['prefactor'])
              for id, proc in list(self.procs.items())
              if proc['barrier'] <= max_barrier]
        if not superbasin:
            return rt
        else:
            # Filter out processes that lead to another state in the superbasin.
            return [entry
                    for entry in rt
                    if self.procs[entry[0]]["product"] not in superbasin.state_dict]

    def get_process_table(self):
        """Return dictionary of processes (id->proc) inside thermal window."""
        return dict([id, self.procs[id]]
                    for id in self.get_relevant_procids())

    def get_relevant_procids(self, superbasin=None):
        """Get process IDs inside thermal window."""
        return [entry[0] for entry in self.get_ratetable(superbasin)]

    def get_confidence(self, superbasin=None):
        """Confidence that all relevant processes for this state were found.

        The confidence is a function of the ratio Nf/Ns, where Nf is
        the number of unique processes and Ns is the number of
        searches performed.

        When Nf or Ns are zero, there is no confidence.  Zero is
        returned in this case.

        As the ratio Nf/Ns decreases, the confidence increases from
        0.0 to a limit of 1.0.  This is roughly equivalent to the
        statement: I am more confident that I have found all relevant
        processes if I have found 10 processes after 1000 searches
        than if I have found 900 processes after 1000 searches.

        Nf is calculated to be the number of processes in the rate
        table. Ns is calculated to be the number of searches that
        resulted in a process on the rate table.

        When using recycling or kdb, it is useful to ignore processes
        that occur outside the hole, the region in which the last
        process took place.  Focusing on this region means you can do
        possibly far fewer searches to reach confidence. When using
        the hole to filter processes, Nf and Ns only take into account
        processes that intersect the hole.

        If superbasin is passed, that means the state is in that
        superbasin. The confidence calculation is adjusted so that
        only processes that lead out of the superbasin are
        counted. This may be disabled by the
        self.config.sb_superbasin_confidence option.

        """
        # Possibly disable superbasin feature.
        #if not self.config.sb_superbasin_confidence:
        try:
            if not self.config.sb_superbasin_confidence:
                superbasin = None
        except AttributeError:
            superbasin = None

        # Checking to see if all recycling jobs are complete
        if self.config.recycling_on and self.config.disp_moved_only:
            job_table_path = os.path.join(self.config.path_root, "jobs.tbl")
            job_table = io.Table(job_table_path)
            if any([ t == 'recycling' for t in job_table.get_column('type') ]):
                return 0.0

        # Load the rate table. If we are in a superbasin, we filter
        # out all processes leading to another state in the same
        # superbasin.
        rt = self.get_ratetable(superbasin)
        # Load the repeat counts and filter if we are in a superbasin.
        prc = self.get_proc_random_count()
        if superbasin:
            prc = dict([proc, count]
                       for proc, count in list(prc.items())
                       if self.procs[proc]["product"] not in superbasin.state_dict)
        alpha = 1.0
        if self.config.akmc_confidence_correction:
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

        if self.config.akmc_confidence_scheme == 'new':
            Nf = 0.0
            Ns = 0.0
            for r in rt:
                if r[0] in prc:
                    Nf += 1.0
                    Ns += prc[r[0]]
            if Ns < 1:
                return 0.0
            if Nf < 1:
                Nf = 1.0
            return 1.0 + (Nf/(alpha*Ns)) * lambertw(-math.exp(-1.0 / (Nf/(alpha*Ns)))/(Nf/(alpha*Ns)))

        elif self.config.akmc_confidence_scheme == 'sampling':
            all_repeats = prc
            repeats = {}
            for event in rt:
                id = event[0]
                repeats[id] = all_repeats[id]
            # number of events
            m = len(repeats)
            # number of searches
            n = sum(repeats.values())
            if n < 10: return 0.0

            # probabilities
            ps = numpy.array(list(repeats.values()),dtype=float)
            ps /= sum(ps)

            C = numpy.zeros(m)
            for i in range(n):
                C += (1.0-C)*ps

            return sum(C)/float(m)

        elif self.config.akmc_confidence_scheme == 'dynamics':
            #print("into dynamics confidence")
            # filter out recycled saddles if displace_moved_only is true
            dyn_saddles = set()
            if self.config.disp_moved_only:
                f = open(self.search_result_path)
                f.readline()
                f.readline()

                for line in f:
                    status = line.split()[7]
                    if not (status.startswith('good') or status.startswith('repeat')):
                        continue
                    search_type = line.split()[1]
                    if search_type != 'dynamics':
                        continue
                    pid = int(status.split('-')[1])
                    dyn_saddles.add(pid)
                f.close()

                rt = [entry
                      for entry in rt
                      if entry[0] in dyn_saddles]

            conf = 0.0
            total_rate_estimator = 0.0
            total_rate_found = 0.0
            T1 = self.config.main_temperature
            for T2, T2_time in list(self.get_time_by_temp().items()):
                #print("out of get_time_by_temp")
                if T2_time == 0.0:
                    continue

                # rates are at T1
                rates = numpy.array([ p[1] for p in rt ])
                prefactors = numpy.array([ p[2] for p in rt ])
                if len(rates) == 0: return 0.0

                # extrapolate to T2
                rates_md = prefactors*(rates/prefactors)**(T1/T2)

                time = T2_time*1e-15
                C = 1.0-numpy.exp(-time*rates_md)
                total_rate_found = sum(rates)

                # Chill confidence
                #conf += sum(C*rates)/total_rate_found

                # Lelievre-Jourdain confidence
                total_rate_estimator += sum(rates/C)

            if(total_rate_estimator > 0):
                conf = total_rate_found / total_rate_estimator
            return conf

        else:
            # No superbasin-specific code needed here, the number of
            # repeats is adjusted for that in the add_process()
            # method!
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
        if abs(kT - self.statelist.kT) > 1e-8:
            for id, proc in list(self.procs.items()):
                proc['rate'] = proc['prefactor'] * math.exp(-proc['barrier'] / self.statelist.kT)
            self.save_process_table()


    def save_process_table(self):
        """ If the processtable is present in memory, writes it to disk. """
        if self.procs != None:
            f = open(self.proctable_path, 'w')
            f.write(self.processtable_header)
            for id in list(self.procs.keys()):
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

    def increment_time(self, dt, T_search):
        """Increment MD search time by dt at temperature T_search."""
        self.info.set("MetaData", "time", self.get_time() + dt)
        temp_str = "%.0f" % T_search # rounded to nearest Kelvin
        if self.info.has_section("SearchTime"):
            # Increase time spent at current search temperature.
            new_time_at_current_temp = \
                self.info.get("SearchTime", temp_str, 0.0) + dt
        else:
            # The info file must come from an older version of EON,
            # which didn't have this section, yet. In this case we
            # assume all prior searches were done at the current
            # temperature and add the SearchTime section. This makes
            # this codes backwards compatible and the correctness is
            # not worse than before.
            new_time_at_current_temp = self.get_time() # time is already updated, so no +dt!
        self.info.set("SearchTime", temp_str, new_time_at_current_temp)

    def get_time(self):
        return self.info.get("MetaData", "time", 0.0)

    def get_time_by_temp(self):
        #print("into get_time_by_temp")
        try:
            return dict([int(temp), float(time)] for temp, time in self.info.items("SearchTime"))
        except configparser.NoSectionError:
            # The "info" file seems to have been produced by an old
            # version of EON which didn't have the SearchTime
            # section. We simply upgrade and try again (no recursion
            # to avoid endless recursion if the exception is raised
            # again).
            self.increment_time(0.0, self.config.saddle_dynamics_temperature)
            #try:
            #    print(self.info.items('SearchTime', raw=True))
            #except Exception as e:
            #    print("exception: " + str(e))
            return dict([int(temp), float(time)] for temp, time in self.info.items("SearchTime", raw=True))

    def get_number_of_searches(self):
        # TODO: this is inefficient!
        f = open(self.search_result_path)
        try:
            n = 0
            f.readline()
            f.readline()
            for line in f:
                n += 1
        finally:
            f.close()
        return n

    def get_total_saddle_count(self):
        return self.get_good_saddle_count() + self.get_bad_saddle_count()

    def get_bad_saddle_count(self):
        return self.info.get("MetaData", "bad_saddles", 0)

    def set_bad_saddle_count(self, num):
        self.info.set("MetaData", "bad_saddles", num)

    def register_bad_saddle(self, result, store=False, superbasin=None):
        """ Registers a bad saddle. """
        #print ("bad saddle ",result["results"]["termination_reason"])
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
                             "Nonlocal Abort",
                             "Negative Barrier",
                             "MD Trajectory Too Short",
                             "No Negative Mode at Saddle",
                             "No Forward Barrier in Minimized Band",
                             "MinMode Zero Mode Abort",
                             "Optimizer Error"
                             ]
        self.set_bad_saddle_count(self.get_bad_saddle_count() + 1)
        self.append_search_result(result, result_state_code[result["results"]["termination_reason"]], superbasin)

        # If a MD saddle search is too short add it to the total clock time.
        if 'simulation_time' in result['results']:
            if result['results']['termination_reason'] == 15: #too short
                self.increment_time(result['results']['simulation_time'], result['results']['md_temperature'])

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
    for i in range(10):
        e = math.exp(w)
        t = w * e - z
        p = w + 1.0
        t /= e * p - 0.5 * (p + 1.0) * t / p
        w -= t
        if abs(t) < eps * (1.0 + abs(w)):
            return w
    logger.error("Failed to converge Lambert W function")
    raise ValueError()
