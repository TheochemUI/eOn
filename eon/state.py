
""" The state module. """

import os
import shutil
import logging
logger = logging.getLogger('state')
from configparser import ConfigParser
import tarfile

from eon import fileio as io

from eon.config import config as EON_CONFIG
from eon.config import ConfigClass # Typing

class State:
    """ The state super class. """
    def __init__(self, statepath, statenumber, statelist, previous_state_num = -1,
                 reactant_path = None, config: ConfigClass = EON_CONFIG):
        """ Creates a new state, with lazily loaded data. """

        self.config = config
        # The parent statelist.
        self.statelist = statelist

        # The path to and number of this state.
        self.path = statepath
        self.number = statenumber

        # Lazily loaded data. Should use the get/set methods for these.
        self.procs = None
        self.proc_repeat_count = None

        self.procdata_path = os.path.join(self.path, "procdata")
        self.reactant_path = os.path.join(self.path, "reactant.con")
        self.proctable_path = os.path.join(self.path, "processtable")
        self.search_result_path = os.path.join(self.path, "search_results.txt")
        self.tar_path = os.path.join(self.path, "procdata.tar")
        self.info = io.ini(os.path.join(self.path, "info"))

        # If this state does not exist on disk, create it.
        if not os.path.isdir(self.path):
            if reactant_path is None:
                raise IOError("State needs a reactant_path when it is being instantiated to disk")
            os.mkdir(self.path)
            os.mkdir(self.procdata_path)
            shutil.copy(reactant_path, self.reactant_path)
            self.info.set("MetaData", "previous state", previous_state_num)
            f = open(self.proctable_path, 'w')
            f.write(self.processtable_header)
            f.close()
            f = open(self.search_result_path, 'w')
            f.write(self.search_result_header)
            f.close()

    def __repr__(self):
        return "State #%i" % self.number

    def __eq__(self, other):
        return isinstance(other, State) and self.number == other.number

    def __hash__(self):
        return hash(self.number)

    def add_process(self, result):
        if 'stdout.dat' in result:
            id = self.get_num_procs()
            open(self.proc_stdout_path(id), 'w').writelines(result['stdout.dat'].getvalue())

    def get_energy(self):
        return self.info.get("MetaData", "reactant energy", None)

    def set_energy(self, e):
        self.info.set("MetaData", "reactant energy", e)

    def get_reactant(self):
        """ Loads the reactant.con into a point and returns it. """
        return io.loadcon(self.reactant_path)

    # Utility functions for loading process .con and mode files.
    def get_process_reactant(self, id):
        return io.loadcon(self.proc_reactant_path(id))
    def get_process_product(self, id):
        return io.loadcon(self.proc_product_path(id))

    # Utility functions for compiling procdata paths, whether the files exist or not.
    def proc_reactant_path(self, id):
        return os.path.join(self.procdata_path, "reactant_%d.con" % id)
    def proc_product_path(self, id):
        return os.path.join(self.procdata_path, "product_%d.con" % id)
    def proc_results_path(self, id):
        return os.path.join(self.procdata_path, "results_%d.dat" % id)
    def proc_stdout_path(self, id):
        return os.path.join(self.procdata_path, "stdout_%d.dat" % id)

    def tar_procdata(self):
        if not self.procdata_tarred:
            tar = tarfile.TarFile(self.tar_path, 'w')
            for i in os.listdir(self.procdata_path):
                tar.add(i)
                os.unlink(i)
            tar.close()
        else:
            logger.warning("Attempted to tar an already tarred procdata")
        self.procdata_tarred = True

    def get_process(self, id):
        self.load_process_table()
        return self.procs[id]

    def get_process_ids(self):
        """ Returns the list of ids in the rate table. """
        return [b[0] for b in self.get_ratetable()]

    def get_previous_state(self):
        return self.info.get("MetaData", "previous state", -1)

    def get_num_procs(self):
        """ Loads the process table if it is not already loaded and returns the length of it """
        self.load_process_table()
        return len(self.procs)

    def get_process_table(self):
        self.load_process_table()
        return self.procs
