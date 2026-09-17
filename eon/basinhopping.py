#!/usr/bin/env python

from io import StringIO
import logging
import logging.handlers
logger = logging.getLogger('basinhopping')
import numpy
numpy.seterr(divide="raise", over="raise", under="print", invalid="raise")
import optparse
import os
import shutil
import sys
import random
from pathlib import Path

from eon import version
from eon import atoms
from eon import communicator
from eon.config import ConfigClass
from eon import fileio as io
from eon import locking

class BHStates:
    def __init__(self, config: ConfigClass):
        self.config = config
        states = Path(self.config.path_states)
        states.mkdir(exist_ok=True)
        state_table_path = str(states / "state_table")
        self.energy_table = io.Table(state_table_path, ['state', 'energy', 'repeats'])

    def get_random_minimum(self):
        if len(self.energy_table.rows) == 0: return None
        N = self.config.bh_initial_state_pool_size
        self.energy_table.rows.sort(key=lambda r:-r['energy'])
        lowest_N = self.energy_table.rows[:N]
        i = random.choice(lowest_N)['state']
        minimum = Path(self.config.path_states) / str(i) / "minimum.con"
        return StringIO(minimum.read_text())

    def add_state(self, result_files, result_info):
        energy = result_info['minimum_energy']

        energetically_close = []
        added = True

        if len(self.energy_table) != 0:
            for row in self.energy_table:
                if abs(energy-row['energy']) < self.config.comp_eps_e:
                    energetically_close.append(row['state'])
            if len(energetically_close) != 0:
                a1 = io.loadcon(result_files['min.con'])
                for state_number in energetically_close:
                    state_con_path = str(
                        Path(self.config.path_states)
                        / str(state_number)
                        / "minimum.con"
                    )
                    a2 = io.loadcon(state_con_path)
                    if atoms.match(a1, a2, self.config.comp_eps_r, self.config.comp_neighbor_cutoff, True, check_rotation=self.config.comp_check_rotation, use_identical=self.config.comp_use_identical):
                        logger.info("Found a repeat of state %i", state_number)
                        added = False
                        for row in self.energy_table.rows:
                            if row['state'] == state_number:
                                row['repeats'] += 1
                                self.energy_table.write()
                                break

        if added:
            state_number = len(self.energy_table)

            row = { 'state':state_number, 'energy':energy, 'repeats':0 }
            self.energy_table.add_row(row)
            self.energy_table.rows.sort(key=lambda r:-r['energy'])
            self.energy_table.write()

            state_path = Path(self.config.path_states) / str(state_number)
            state_path.mkdir()

            result_files['minimum.con'] = result_files['min.con']
            del result_files['min.con']

            for fn, fh in list(result_files.items()):
                if hasattr(fh, 'getvalue') == False:
                    continue
                (state_path / fn).write_text(fh.getvalue())

        return added

def basinhopping(config: ConfigClass = None):
    if config is None:
        raise TypeError("basinhopping requires a ConfigClass instance")
    logger.info('Eon version: %s', version)
    # First of all, does the root directory even exist?
    if not Path(config.path_root).is_dir():
        logger.critical("Root directory does not exist")
        sys.exit(1)

    # load metadata
    bhstates = BHStates(config)

    wuid_path = Path("wuid.dat")
    if wuid_path.is_file():
        wuid = int(wuid_path.read_text().strip())
    else:
        wuid = 0

    # get communicator
    comm = communicator.get_communicator(config)

    # Register all the results. There is  no need to ever discard found
    # processes like we do with akmc. There is no confidence to calculate.
    register_results(comm, bhstates, config)

    wuid = make_searches(comm, wuid, bhstates, config)

    Path("wuid.dat").write_text("%i\n" % wuid)

    io.save_prng_state(io.prng_state_path(config))

def make_searches(comm, wuid, bhstates, config: ConfigClass):
    num_in_buffer = comm.queued_search_count()
    logger.info("%i searches in the queue" % num_in_buffer)
    num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
    logger.info("Making %i searches" % num_to_make)

    if num_to_make == 0:
        return wuid

    searches = []

    invariants = {}

    initial_react = StringIO((Path(config.path_root) / "pos.con").read_text())

    #invariants['reactant_passed.con']=reactIO

    ini_changes = [ ('Main', 'job', 'basin_hopping') ]
    #invariants['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)
    #invariants['reactant_passed.con']  = reactIO

    searches = []
    for i in range(num_to_make):
        search = {}
        search['id'] = "%d" % wuid
        ini_changes = [ ('Main', 'random_seed', str(int(numpy.random.random()*2**32))) ]

        if config.bh_initial_state_pool_size == 0:
            reactIO = initial_react
        elif config.bh_initial_state_pool_size > 0:
            reactIO = bhstates.get_random_minimum()
            if reactIO is None:
                reactIO = initial_react
        else:
            logger.fatal("Initial state pool size negative")
            sys.exit(1)

        search['pos.con'] = reactIO
        search['config.ini'] = io.modify_config(config.config_path, ini_changes)
        searches.append(search)
        wuid += 1

    comm.submit_jobs(searches, invariants)
    logger.info( str(num_to_make) + " searches created")
    return wuid

def register_results(comm, bhstates, config):
    logger.info("Registering results")
    jobs_in = Path(config.path_jobs_in)
    if jobs_in.is_dir():
        shutil.rmtree(jobs_in)
    jobs_in.mkdir(parents=True)

    # Function used by communicator to determine whether to discard a result
    def keep_result(name):
        return True

    num_registered = 0

    for result in comm.get_results(config.path_jobs_in, keep_result):
        # The result dictionary contains the following key-value pairs:
        # product.con - an array of strings containing the reactant
        # results.dat - an array of strings containing the results
        # id - wuid

        result_info = io.parse_results(result['results.dat'])
        if 'minimum_energy' not in result_info:
            continue
        if result_info['termination_reason'] == 0:
            if bhstates.add_state(result, result_info):
                logger.info("New structure with energy %.8e",
                            result_info['minimum_energy'])

        num_registered += 1

    logger.info("%i (result) searches processed", num_registered)

def main(config: ConfigClass = None):
    if config is None:
        config = ConfigClass()
    optpar = optparse.OptionParser(usage="usage: %prog [options] config.ini")
    optpar.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False,help="only write to the log file")
    optpar.add_option("-n", "--no-submit", action="store_true", dest="no_submit", default=False,help="don't submit searches; only register finished results")
    optpar.add_option("-R", "--reset", action="store_true", dest="reset", default = False, help="reset the simulation, discarding all data")
    options, args = optpar.parse_args()

    if len(args) > 1:
        print("takes only one positional argument")
    config.init_from_cli(args)

    if config.comm_job_bundle_size != 1:
        print("error: Basin Hopping only supports a bundle size of 1")
        sys.exit(1)

    if options.no_submit:
        config.comm_job_buffer_size = 0

    if options.reset:
        res = input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
                rmdirs = [config.path_jobs_out, config.path_jobs_in, config.path_scratch,  config.path_states]
                for i in rmdirs:
                    if Path(i).is_dir():
                        io.remove_tree_and_empty_parents(i)
                log_path = Path(config.path_results) / "bh.log"
                wuid_path = Path(config.path_results) / "wuid.dat"
                prng_path = Path(io.prng_state_path(config))
                for i in [log_path, wuid_path, prng_path]:
                    i.unlink(missing_ok=True)
                print("Reset.")
                sys.exit(0)
        else:
            print("Not resetting.")
            sys.exit(1)

    # setup logging
    logging.basicConfig(level=logging.DEBUG,
            filename=str(Path(config.path_results) / "bh.log"),
            format="%(asctime)s %(levelname)s:%(name)s: %(message)s",
            datefmt="%F %T")
    logging.raiseExceptions = False

    if not options.quiet:
        rootlogger = logging.getLogger('')
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter("%(message)s")
        console.setFormatter(formatter)
        rootlogger.addHandler(console)

    lock = locking.LockFile(str(Path(config.path_results) / "lockfile"))

    if lock.aquirelock():
        if config.comm_type == 'mpi':
            from eon.mpiwait import mpiwait
            while True:
                mpiwait(config.mpi_poll_period)
                basinhopping(config)
        basinhopping(config)
    else:
        logger.warning("Couldn't get lock")
        sys.exit(1)

if __name__ == '__main__':
    main()
