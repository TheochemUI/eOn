
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
                 reactant_path = None, config=None):
        """ Creates a new state, with lazily loaded data. """
        # Match StateList / AKMCState: accept config= (PR previously dropped it).
        kwargs = {}
        if config is not None:
            kwargs["config"] = config
        state.State.__init__(
            self,
            statepath,
            statenumber,
            statelist,
            previous_state_num,
            reactant_path,
            **kwargs,
        )

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

        # Content-addressed id from product geometry (not table length).
        product_bytes = result['product.con'].getvalue()
        id = self.allocate_process_id(b"pr-product", product_bytes)

        # Keep track of the number of searches, Ns.
        #self.inc_proc_repeat_count(id)

        # Move the relevant files into the procdata directory.
        for path, content in ((self.proc_reactant_path(id), result['reactant.con'].getvalue()),
                              (self.proc_product_path(id), product_bytes),
                              (self.proc_results_path(id), result['results.dat'].getvalue())):
            with io.atomic_write(path) as f:
                f.writelines(content)

        # Append this barrier to the process table (in memory and on disk).
        self.append_process_table(id =                id,
                                  product =           -1,
                                  product_energy =    resultdata["potential_energy_product"],
                                  time =              resultdata["transition_time_s"])

        # This was a unique process, so return the id.
        return id

    def load_process_table(self, force=False):
        """Load the process table from disk (see AKMCState.load_process_table)."""
        if self.procs is not None and not force:
            return
        f = open(self.proctable_path)
        lines = f.readlines()
        f.close()
        self.procs = {}
        n_rows = 0
        for l in lines[1:]:
            parts = l.strip().split()
            if not parts:
                continue
            n_rows += 1
            pid = int(parts[self.ID])
            self.procs[pid] = {
                "product": int(parts[self.PRODUCT]),
                "product_energy": float(parts[self.PRODUCT_ENERGY]),
                "time": float(parts[self.TIME]),
            }
        if n_rows > len(self.procs):
            logger.warning(
                "State %s processtable has %d rows but only %d distinct ids "
                "(duplicates collapsed; last row per id kept).",
                self.number,
                n_rows,
                len(self.procs),
            )

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
        self.load_process_table(force=True)
        if id in self.procs:
            raise RuntimeError(
                "refusing to clobber process id %d in state %s (already registered); "
                "use allocate_process_id() for a content-addressed free id"
                % (id, self.number)
            )
        f = open(self.proctable_path, 'a')
        f.write(self.processtable_line % (id, product, product_energy, time))
        f.close()
        self.procs[id] = {
            "product": product,
            "product_energy": product_energy,
            "time": time,
        }

    def get_time(self):
        return self.info.get("MetaData", "accumulated_time", 0.0)

    def inc_time(self, timeinc):
        time = self.get_time()
        self.info.set("MetaData", "accumulated_time", time + timeinc)

    def zero_time(self):
        self.info.set("MetaData", "accumulated_time", "0.0")
