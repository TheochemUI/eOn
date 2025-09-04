
""" The PRState module. """

import logging
logger = logging.getLogger('state')

from eon import fileio as io
from eon import state

class PRState(state.State):
    ID, PRODUCT, PRODUCT_ENERGY, TIME = list(range(4))
    processtable_head_fmt = "%7s %9s %16s %12s\n"
    processtable_header = processtable_head_fmt % ("proc #", "product", "product energy",
                                                   "time")
    processtable_line = "%7d %9d %16.5f %12.5e\n"
    search_result_header = "%8s %10s\n" % ("wuid", "result")
    search_result_header += "-" * len(search_result_header) + '\n'
    def __init__(self, statepath, statenumber, statelist, previous_state_num = -1,
                 reactant_path = None):
        """ Creates a new state, with lazily loaded data. """
        state.State.__init__(self,statepath, statenumber,statelist, previous_state_num,
                    reactant_path)

    def add_process(self, result):
        """ Adds a process to this state. """
        state.State.add_process(self, result)

        resultdata = result["results"] # The information from the result.dat file

        # We may not already have the energy for this State.  If not, it should be in the result data.
        if self.get_energy() is None:
            self.set_energy(resultdata["potential_energy_reactant"])

        # Check if the reactant, and product are legit
        try:
            if 'reactant' not in result:
                io.loadcon(result['reactant.con'])
            if 'product' not in result:
                io.loadcon(result['product.con'])
        except:
            logger.exception("Reactant or product has incorrect format")
            return None

        # Update the search result table.
        #self.append_search_result(result, "good-%d" % self.get_num_procs())

        # The id of this process is the number of processes.
        id = self.get_num_procs()

        # Keep track of the number of searches, Ns.
        #self.inc_proc_repeat_count(id)

        # Move the relevant files into the procdata directory.
        open(self.proc_reactant_path(id), 'w').writelines(result['reactant.con'].getvalue())
        open(self.proc_product_path(id), 'w').writelines(result['product.con'].getvalue())
        open(self.proc_results_path(id), 'w').writelines(result['results.dat'].getvalue())

        # Append this barrier to the process table (in memory and on disk).
        self.append_process_table(id =                id,
                                  product =           -1,
                                  product_energy =    resultdata["potential_energy_product"],
                                  time =              resultdata["transition_time_s"])

        # This was a unique process, so return the id.
        return id

    def load_process_table(self):
        """ Load the process table.  If the process table is not loaded, load it.  If it is
            loaded, do nothing. """
        if self.procs is None:
            f = open(self.proctable_path)
            lines = f.readlines()
            f.close()
            self.procs = {}
            for l in lines[1:]:
                l = l.strip().split()
                self.procs[int(l[self.ID])] = {
                                          "product":           int  (l[self.PRODUCT]),
                                          "product_energy":    float(l[self.PRODUCT_ENERGY]),
                                          "time":              float(l[self.TIME]),
                                         }

    def save_process_table(self):
        """ If the processtable is present in memory, writes it to disk. """
        if self.procs != None:
            f = open(self.proctable_path, 'w')
            f.write(self.processtable_header)
            for id in list(self.procs.keys()):
                proc = self.procs[id]
                f.write(self.processtable_line % (id, proc['product'], proc['product_energy'],
                                                  proc['time']))
            f.close()

    def append_process_table(self, id, product, product_energy, time):
        """ Append to the process table.  Append a single line to the process table file.  If we
            have loaded the process table, also append it to the process table in memory. """
        f = open(self.proctable_path, 'a')
        f.write(self.processtable_line % (id, product, product_energy, time))
        f.close()
        if self.procs != None:
            self.procs[id] = {
                              "product":           product,
                              "product_energy":    product_energy,
                              "time":              time
                             }

    def get_time(self):
        return self.info.get("MetaData", "accumulated_time", 0.0)

    def inc_time(self, timeinc):
        time = self.get_time()
        self.info.set("MetaData", "accumulated_time", time + timeinc)

    def zero_time(self):
        self.info.set("MetaData", "accumulated_time", "0.0")
