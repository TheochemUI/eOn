#!/usr/bin/env python

import configparser

from io import StringIO
import logging
import logging.handlers
logger = logging.getLogger('pr')
import numpy
numpy.seterr(divide="raise", over="raise", under="print", invalid="raise")
import optparse
import shutil
import sys
from pathlib import Path

from eon import version
from eon.config import ConfigClass
from eon import atoms
from eon import communicator
from eon import fileio as io
from eon import locking
from eon import prstatelist

def parallelreplica(config: ConfigClass = None):
    if config is None:
        raise TypeError("parallelreplica requires a ConfigClass instance")
    logger.info('Eon version: %s', version)
    # First of all, does the root directory even exist?
    if not Path(config.path_root).is_dir():
        logger.critical("Root directory does not exist")
        sys.exit(1)

    # load metadata
    start_state_num, time, wuid = get_pr_metadata(config)
    logger.info("Simulation time is %e", time)
    states = get_statelist(config)
    current_state = states.get_state(start_state_num)

    # get communicator
    comm = communicator.get_communicator(config)

    # Register all the results. There is no need to ever discard processes
    # like we do with akmc. There is no confidence to calculate.
    num_registered, transition, sum_spdup = register_results(comm, current_state, states, config)

    logger.info("Time in current state is %e", current_state.get_time())

    if num_registered >= 1:
        avg_spdup = sum_spdup/num_registered
        logger.info("Total speedup is %f",avg_spdup)
    if transition:
        #current_state, previous_state = step(time, current_state, states, transition)
        time += transition['time']

    wuid = make_searches(comm, current_state, wuid, config)

    parser = configparser.RawConfigParser()
    write_pr_metadata(parser, current_state.number, time, wuid)
    io.write_info_txt(config, parser)
    io.save_prng_state(io.prng_state_path(config))

def step(current_time, current_state, states, transition, config: ConfigClass = None):
    if config is None:
        raise TypeError("step requires a ConfigClass instance")
    next_state = states.get_product_state(current_state.number, transition['process_id'])
    next_state.zero_time()
    dynamics = io.Dynamics(str(Path(config.path_results) / "dynamics.txt"))
    proc = current_state.get_process(transition['process_id'])
    dynamics.append(current_state.number, transition['process_id'],
                    next_state.number, transition['time'], transition['time']+current_time, 0, 0, current_state.get_energy())

    previous_state = current_state
    current_state = next_state

    logger.info("Currently in state %i", current_state.number)

    return current_state, previous_state

def get_statelist(config: ConfigClass = None):
    if config is None:
        raise TypeError("get_statelist requires a ConfigClass instance")
    initial_state_path = str(Path(config.path_root) / "pos.con")
    return prstatelist.PRStateList(initial_state_path, config=config)

def get_pr_metadata(config: ConfigClass = None):
    if config is None:
        raise TypeError("get_pr_metadata requires a ConfigClass instance")
    Path(config.path_results).mkdir(parents=True, exist_ok=True)
    metafile = io.info_txt_path(config)
    parser = configparser.ConfigParser()
    if Path(metafile).is_file():
        parser.read(str(metafile))
        try:
            start_state_num = parser.getint("Simulation Information",'current_state')
        except:
            start_state_num = 0
        try:
            time = parser.getfloat("Simulation Information", 'time_simulated')
        except:
            time = 0.0
        try:
            wuid = parser.getint("PR Metadata", 'wu_id')
        except:
            wuid = 0
    else:
        time = 0
        start_state_num = 0
        wuid = 0

    return start_state_num, time, wuid

def write_pr_metadata(parser, current_state_num, time, wuid):
    parser.add_section('PR Metadata')
    parser.add_section('Simulation Information')
    parser.set('PR Metadata', 'wu_id', str(wuid))
    parser.set('Simulation Information', 'time_simulated', str(time))
    parser.set('Simulation Information', 'current_state', str(current_state_num))


def end_state_table_path(config, state_number):
    """Escape-rate table for one reactant state, under ConfigClass.path_states."""
    return str(Path(config.path_states) / str(state_number) / "end_state_table")


def load_end_state_table(path):
    """Read the end-state table. The file is opened for reading from the start."""
    rows = []
    p = Path(path)
    if not p.is_file():
        return rows
    lines = p.read_text().splitlines()
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        row = {
            "state": int(parts[0]),
            "views": int(parts[1]),
            "rate": float(parts[2]),
            "time": float(parts[3]),
        }
        if len(parts) >= 5:
            row["process_id"] = int(parts[4])
        rows.append(row)
    return rows


def save_end_state_table(path, rows):
    parent = Path(path).parent
    if parent != Path("."):
        parent.mkdir(parents=True, exist_ok=True)
    with io.atomic_write(path) as g:
        g.write("state       views         rate        time        process_id\n")
        for row in rows:
            g.write(
                "%s             %s             %s             %s             %s\n"
                % (
                    row["state"],
                    row["views"],
                    row["rate"],
                    row["time"],
                    row.get("process_id", -1),
                )
            )


def find_matching_process(product, state, config):
    """Return the process id whose stored product matches *product*, or None.

    Products live at product_<xxh64>.con via State.proc_product_path, not
    product_0.con / product_1.con.
    """
    state.load_process_table()
    if not state.procs:
        return None
    for pid in state.procs:
        ppath = state.proc_product_path(pid)
        if not Path(ppath).is_file():
            continue
        product2 = io.loadcon(ppath)
        if atoms.match(
            product,
            product2,
            config.comp_eps_r,
            config.comp_neighbor_cutoff,
            True,
            check_rotation=config.comp_check_rotation,
            use_identical=config.comp_use_identical,
        ):
            return pid
    return None


def accumulate_end_state(product, state, result, config):
    """Fold one transition into the end-state table. Return the process id.

    A repeat product updates views/time/rate and does not call add_process,
    because process ids are content-addressed and a second add would collide.
    """
    table_path = end_state_table_path(config, state.number)
    rows = load_end_state_table(table_path)
    trans_time = float(result["results"]["transition_time_s"])
    time_check = max((float(r["time"]) for r in rows), default=0.0)
    matched_pid = find_matching_process(product, state, config)

    if matched_pid is not None:
        updated = False
        for row in rows:
            if row.get("process_id") == matched_pid:
                row["views"] = int(row["views"]) + 1
                row["time"] = float(row["time"]) + trans_time
                row["rate"] = (
                    float(row["views"]) / float(row["time"]) if row["time"] else 0.0
                )
                updated = True
                break
        if not updated:
            next_state = max((int(r["state"]) for r in rows), default=-1) + 1
            rows.append(
                {
                    "state": next_state,
                    "views": 1,
                    "rate": (1.0 / trans_time) if trans_time else 0.0,
                    "time": trans_time,
                    "process_id": matched_pid,
                }
            )
        process_id = matched_pid
    else:
        process_id = state.add_process(result)
        next_state = max((int(r["state"]) for r in rows), default=-1) + 1
        new_time = time_check + trans_time
        rows.append(
            {
                "state": next_state,
                "views": 1,
                "rate": (1.0 / new_time) if new_time else 0.0,
                "time": new_time,
                "process_id": process_id,
            }
        )

    save_end_state_table(table_path, rows)
    return process_id

def make_searches(comm, current_state, wuid, config: ConfigClass = None):
    if config is None:
        raise TypeError("make_searches requires a ConfigClass instance")
    reactant = current_state.get_reactant()
    num_in_buffer = comm.queued_search_count()
    logger.info("Queue contains %i searches" % num_in_buffer)
    num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
    logger.info("Making %i searches" % num_to_make)

    if num_to_make == 0:
        return wuid

    searches = []

    invariants = {}

    reactIO = StringIO()
    io.savecon(reactIO, reactant)

    searches = []
    for i in range(num_to_make):
        search = {}
        search['id'] = "%d_%d" % (current_state.number, wuid)
        search['pos.con']  = reactIO
        ini_changes = [
                        ('Main', 'job', 'parallel_replica'),
                        ('Main', 'random_seed',
                            str(int(numpy.random.random()*10**9))),
                      ]
        search['config.ini'] = io.modify_config(config.config_path, ini_changes)
        searches.append(search)
        wuid += 1

    comm.submit_jobs(searches, invariants)
    logger.info("Created " + str(num_to_make) + " searches")
    return wuid

def register_results(comm, current_state, states, config: ConfigClass = None):
    if config is None:
        raise TypeError("register_results requires a ConfigClass instance")
    logger.info("Registering results")
    jobs_in = Path(config.path_jobs_in)
    if jobs_in.is_dir():
        shutil.rmtree(jobs_in)
    jobs_in.mkdir(parents=True)
    # Function used by communicator to determine whether to discard a result
    def keep_result(name):
        return True
    transition = None
    num_registered = 0
    speedup = 0
    for result in comm.get_results(config.path_jobs_in, keep_result):
        # The result dictionary contains the following key-value pairs:
        # reactant.con - an array of strings containing the reactant
        # product.con - an array of strings containing the product
        # results.dat - an array of strings containing the results
        # id - StateNumber_WUID
        #
        # The reactant, product, and mode are passed as lines of the files because
        # the information contained in them is not needed for registering results

        state_num = int(result['name'].split("_")[0])

        state = states.get_state(state_num)

        # read in the results
        result['results'] = io.parse_results(result['results.dat'])
        speedup += result['results']['speedup']
        if result['results']['transition_found'] == 1:
            result['results']['transition_time_s'] += state.get_time()
            product = io.loadcon(result['product.con'])
            process_id = accumulate_end_state(product, state, result, config)
            time = result['results']['transition_time_s']
            logger.info("Found transition with time %.3e", time)
            if not transition and current_state.number==state.number:
                transition = {'process_id':process_id, 'time':time}
            state.zero_time()
        else:
            state.inc_time(result['results']['simulation_time_s'])
        num_registered += 1
    logger.info("Processed %i (result) searches", num_registered)
    if num_registered >=1:
        logger.info("Average speedup is  %f", speedup/num_registered)
    return num_registered, transition, speedup

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
        print("error: Parallel Replica only supports a bundle size of 1")
        sys.exit(1)

    if options.no_submit:
        config.comm_job_buffer_size = 0

    if options.reset:
        res = input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
            rmdirs = [config.path_jobs_out, config.path_jobs_in, config.path_states,
                    config.path_scratch]
            if config.debug_keep_all_results:
                rmdirs.append(Path(config.path_root) / "old_searches")
            for i in rmdirs:
                if Path(i).is_dir():
                    io.remove_tree_and_empty_parents(i)

            results = Path(config.path_results)
            for i in [
                io.info_txt_path(config),
                results / "dynamics.txt",
                results / "pr.log",
                io.prng_state_path(config),
            ]:
                p = Path(i)
                if p.is_file():
                    p.unlink()

            print("Reset")
        sys.exit(0)

    # setup logging
    logging.basicConfig(level=logging.DEBUG,
            filename=str(Path(config.path_results) / "pr.log"),
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
                parallelreplica(config)
        parallelreplica(config)
    else:
        logger.warning("Couldn't get lock")
        sys.exit(1)

if __name__ == '__main__':
    main()
