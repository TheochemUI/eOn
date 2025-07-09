#!/usr/bin/env python

import math
import sys
import configparser
import os.path
import shutil
from time import sleep
import os
import time as unix_time
import optparse
import logging
import logging.handlers
logger = logging.getLogger('akmc')
import numpy
from builtins import input
from pathlib import Path

numpy.seterr(divide="raise", over="raise", under="print", invalid="raise")

from eon.version import version
from eon.config import config as EON_CONFIG # NOTE: Since it is a global object..
from eon.config import ConfigClass
from eon import communicator
from eon import locking
from eon import akmcstatelist
from eon import explorer
from eon import fileio as io
from eon import atoms
from eon import superbasinscheme
from eon import askmc
from eon import movie

def akmc(config: ConfigClass = EON_CONFIG, steps=0):
    """Poll for status of AKMC clients and possibly make KMC steps.

    Returns the number of KMC steps in the current run at the end of
    this function. The argument "steps" gives the number of KMC steps
    done before calling this function, defaulting to zero.

    """
    #log version information
    logger.info('Eon version %s', version())
    # Here's what this does:
    # 1) Read in the state of our calculation from last time
    # 2) Initialize necessary data structures (statelist, communicator, displace)
    # 3) Get any results that have come in
    # 4) Possibly take a KMC step
    # 5) Make new work units
    # 6) Write out the state of the simulation

    # First of all, does the root directory even exist?
    if not os.path.isdir(config.path_root):
        logger.critical("Root directory does not exist, as such the " \
                        "reactant cannot exist. Exiting...")
        sys.exit(1)

    # If we are saving debug results, create the directory if it does not exist.
    if config.debug_keep_all_results:
        rp = os.path.join(config.path_root, config.debug_results_path)
        if not os.path.isdir(rp):
            os.mkdir(rp)

    # Define constants.
    kT = config.main_temperature/11604.5 #in eV

    # Load metadata, the state list, and the current state.
    start_state_num, time, previous_state_num, first_run, previous_temperature = get_akmc_metadata()

    if first_run:
        previous_temperature = config.main_temperature

    # Did the temperature change? Then the existing superbasins are invalid now.
    if abs(previous_temperature - config.main_temperature) > 1e-6:
        # Remove superbasin data.
        if os.path.isdir(config.sb_path):
            shutil.rmtree(config.sb_path)
        state_dirs = os.listdir(config.path_states)
        for state_dir in state_dirs:
            if state_dir == 'state_table':
                continue
            superbasin_file = os.path.join(config.path_states, state_dir, config.sb_state_file)
            if os.path.isfile(superbasin_file):
                os.remove(superbasin_file)
        # Keep the new temperature.
        previous_temperature = config.main_temperature

    states = get_statelist(kT)
    current_state = states.get_state(start_state_num)
    if previous_state_num == -1:
        previous_state = current_state
    else:
        previous_state = states.get_state(previous_state_num)
    # If the Novotny-based superbasining scheme is being used, initialize it.
    if config.sb_on:
        superbasining = get_superbasin_scheme(states, config)
        sb = superbasining.get_containing_superbasin(current_state)
        # If we are exploring states in a superbasin, we will always
        # explore the state with the lowest confidence.
        if sb:
            explore_state = sb.get_lowest_confidence_state()
            previous_state = explore_state # TODO: perhaps there is a better value for previous_state?
        else:
            explore_state = current_state
    else:
        superbasining = None
        sb = None
        explore_state = current_state

    state_explorer = explorer.get_minmodexplorer()(states, previous_state, explore_state, superbasin=sb)
    state_explorer.explore()

    # Take a KMC step, if it's time.
    current_state, previous_state, time, steps = kmc_step(current_state, states, time, kT, superbasining, steps)

    # Write out metadata.
    metafile = os.path.join(config.path_results, 'info.txt')
#    parser = configparser.RawConfigParser()
    parser = configparser.ConfigParser()

    if previous_state.number != current_state.number:
        previous_state_num = previous_state.number

    write_akmc_metadata(parser, current_state.number, time, previous_state_num, previous_temperature)

    parser.write(open(metafile, 'w'))

    io.save_prng_state()

    return steps

def get_akmc_metadata(config: ConfigClass = EON_CONFIG):
    if not os.path.isdir(config.path_results):
        os.makedirs(config.path_results)
    # read in metadata
    # do we want custom metadata locations?
    metafile = os.path.join(config.path_results, 'info.txt')
    parser = io.ini(metafile)
    if os.path.isfile(metafile):
        start_state_num = parser.get("Simulation Information",'current_state', 0)
        time = parser.get("Simulation Information", 'time_simulated', 0.0)
        previous_state_num = parser.get("Simulation Information", "previous_state", -1)
        previous_temperature = parser.get("Simulation Information", "previous_temperature", 0)
        first_run = parser.get("Simulation Information", "first_run", True)
    else:
        time = 0
        start_state_num = 0
        previous_state_num = -1
        first_run = True
        previous_temperature = config.main_temperature

    return start_state_num, time, previous_state_num, first_run, previous_temperature

def write_akmc_metadata(parser, current_state_num, time, previous_state_num, previous_temperature):
    parser.add_section('aKMC Metadata')
    parser.add_section('Simulation Information')
    parser.set('Simulation Information', 'time_simulated', str(time))
    parser.set('Simulation Information', 'current_state', str(current_state_num))
    parser.set('Simulation Information', 'previous_state', str(previous_state_num))
    parser.set('Simulation Information', 'first_run', str(False))
    parser.set('Simulation Information', 'previous_temperature', str(previous_temperature))

def get_statelist(kT, config: ConfigClass = EON_CONFIG):
    initial_state_path = os.path.join(config.path_root, "pos.con")
    return akmcstatelist.AKMCStateList(
        kT,
        config.akmc_thermal_window,
        config.akmc_max_thermal_window,
        initial_state_path,
        filter_hole=config.disp_moved_only,
        config=config,
    )

def get_superbasin_scheme(states, config):
    if config.sb_scheme == 'transition_counting':
        superbasining = superbasinscheme.TransitionCounting(config.sb_path, states, config.main_temperature / 11604.5, config.sb_tc_ntrans)
    elif config.sb_scheme == 'energy_level':
        superbasining = superbasinscheme.EnergyLevel(config.sb_path, states, config.main_temperature / 11604.5, config.sb_el_energy_increment)
    elif config.sb_scheme == 'rate':
        superbasining = superbasinscheme.RateThreshold(config.sb_path, states, config.main_temperature / 11604.5, config.sb_rt_rate_threshold)
    return superbasining


def kmc_step(current_state, states, time, kT, superbasining, steps=0, config: ConfigClass = EON_CONFIG):
    t1 = unix_time.time()
    previous_state = current_state
    # If the Chatterjee & Voter superbasin acceleration method is being used
    if config.askmc_on:
        pass_rec_path = None
        asKMC = askmc.ASKMC(kT, states, config.askmc_confidence, config.askmc_alpha,
                            config.askmc_gamma, config.askmc_barrier_test_on,
                            config.askmc_connections_test_on, config.sb_recycling_on,
                            config.path_root, config.akmc_thermal_window,
                            recycle_path = pass_rec_path)

    # The system might be in a superbasin
    if config.sb_on:
        sb = superbasining.get_containing_superbasin(current_state)
    else:
        sb = None

    while ((
            (not sb and current_state.get_confidence() >= config.akmc_confidence) or
            (sb and sb.get_confidence() >= config.akmc_confidence)
           ) and
           (steps < config.akmc_max_kmc_steps or config.akmc_max_kmc_steps == 0)):

        # Do a KMC step.
        steps += 1
        if config.sb_on and sb:
            mean_time, current_state, next_state, sb_proc_id_out, sb_id = sb.step(current_state, states.get_product_state)
        else:
            if config.askmc_on:
                rate_table = asKMC.get_ratetable(current_state)
            else:
                rate_table = current_state.get_ratetable()
            if len(rate_table) == 0:
                logger.error("No processes in rate table, but confidence "
                             "has been reached")

            ratesum = sum((row[1] for row in rate_table), 0.0)

            u = numpy.random.random_sample()
            p = 0.0
            nsid = 1.1 # Next state process id, will throw exception if remains unchanged.

            # If we are following another trajectory:
            if config.debug_target_trajectory != "False":
                # Get the Dynamics objects.
                owndynamics = io.Dynamics(os.path.join(config.path_results, "dynamics.txt")).get()
                targetdynamics = io.Dynamics(os.path.join(config.debug_target_trajectory, "dynamics.txt")).get()
                # Get the current_step.
                try:
                    current_step = len(owndynamics)
                except:
                    current_step = 0
                # Get the target step process id.
                if current_step > 0:
                    stateid = targetdynamics[current_step]['reactant']
                else:
                    stateid = 0
                try:
                    procid = targetdynamics[current_step]['process']
                except:
                    print("Can no longer follow target trajectory")
                    sys.exit(1)
                # Load the con file for that process saddle.
                targetSaddleCon = io.loadcon(os.path.join(config.debug_target_trajectory, "states", str(stateid), "procdata", "saddle_%d.con" % procid))
                targetProductCon = io.loadcon(os.path.join(config.debug_target_trajectory, "states", str(stateid), "procdata", "product_%d.con" % procid))
                ibox = numpy.linalg.inv(targetSaddleCon.box)
                # See if we have this process
                for i in range(len(rate_table)):
                    p1 = current_state.get_process_saddle(rate_table[i][0])
                    for dist in atoms.per_atom_norm_gen(p1.free_r() - targetSaddleCon.free_r(), targetSaddleCon.box, ibox):
                        if dist > config.comp_eps_r:
                            break
                    else:
                        p1 = current_state.get_process_product(rate_table[i][0])
                        for dist in atoms.per_atom_norm_gen(p1.free_r() - targetProductCon.free_r(), targetProductCon.box, ibox):
                            if dist > config.comp_eps_r:
                                break
                        else:
                            nsid = i
                            break
                else:
                    print("Can no longer follow target trajectory")
                    sys.exit(1)

            # We are not following another trajectory:
            else:
                for i, row in enumerate(rate_table):
                    p += row[1]/ratesum
                    if p>u:
                        nsid = i
                        break
                else:
                    logger.warning("Warning: Failed to select rate; p = " + str(p))
                    break

            next_state = states.get_product_state(current_state.number, rate_table[nsid][0])
            mean_time = 1.0/ratesum

        print("Meantime for Step "+str(steps)+": ", mean_time)
        # Accounting for time
        if config.debug_use_mean_time:
            step_time = mean_time
        else:
            #numpy.random.random_sample() uses [0,1)
            #which could produce issues with math.log()
            step_time = -mean_time*math.log(1 - numpy.random.random_sample())

        time += step_time

        # Pass transition information to extension schemes
        if config.askmc_on:
            asKMC.register_transition(current_state, next_state)
        if config.sb_on:
            superbasining.register_transition(current_state, next_state)

        if config.sb_on and sb:
            proc_id_out = -1
        else:
            proc_id_out = rate_table[nsid][0]

        # Write data to disk
        dynamics = io.Dynamics(os.path.join(config.path_results, "dynamics.txt"))
        if proc_id_out != -1:
            proc = current_state.get_process(proc_id_out)
            dynamics.append(current_state.number, proc_id_out, next_state.number, step_time, time, proc['barrier'], proc['rate'], current_state.get_energy())
            logger.info("KMC step from state %i through process %i to state %i ", current_state.number, rate_table[nsid][0], next_state.number)
        else:
            #XXX The proc_out_id was -1, which means there's a bug or this was a superbasin step.
            dynamics.append_sb(current_state.number, sb_proc_id_out, next_state.number, step_time, time, sb_id, 1.0/mean_time, current_state.get_energy())
            logger.info("SB step from state %i through process %i to state %i ", current_state.number, sb_proc_id_out, next_state.number)

        #criterion used to stop the job: currently an energy limit is used as criterion
        if current_state.get_energy() > config.debug_stop_criterion:
           sys.exit()

        previous_state = current_state
        current_state = next_state

        # The system might be in a superbasin
        if config.sb_on:
            sb = superbasining.get_containing_superbasin(current_state)
        else:
            sb = None

    if config.sb_on:
        superbasining.write_data()

    if not sb:
        logger.info("Currently in state %i with confidence %.6f", current_state.number,
                    current_state.get_confidence())
    else:
        logger.info("Currently in state %i (superbasin %i) with confidence %.6f",
                    current_state.number,
                    sb.id,
                    sb.get_confidence())
    t2 = unix_time.time()
    logger.debug("KMC finished in %.4f seconds", (t2-t1))
    logger.debug("%.2f KMC steps per second", float(steps)/(t2-t1))
    return current_state, previous_state, time, steps


def main(config: ConfigClass = EON_CONFIG):
    optpar = optparse.OptionParser(usage = "usage: %prog [options] config.ini")
    optpar.add_option("-C", "--continuous", action="store_true", dest="continuous", default=False, help="don't quit")
    optpar.add_option("-R", "--reset", action="store_true", dest="reset", default = False, help="reset the aKMC simulation, discarding all data")
    optpar.add_option("-f", "--force", action="store_true", dest="force", default = False, help="force a reset, no questions asked")
    optpar.add_option("-r", "--restart", action="store_true", dest="restart", default = False, help="restart the aKMC simulations from a clean dynamics.txt file")
    optpar.add_option("-s", "--status", action="store_true", dest="print_status", default = False, help = "print the status of the simulation and currently running jobs")
    optpar.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False,help="only write to the log file")
    optpar.add_option("-m", "--movie", action="store", dest="movie_type", default = "", help="Specify the type of movie to make [dynamics, states, fastestpath, fastestfullpath, graph, processes]. Process movies are specified like so: --movie processes,statenumber,processlimit. Where processes is the string processes, statenumber is the number of the state that you want to view, and process limit is the maximum number of processes you would like in the movie. The returned processes are reverse sorted by rate such that the fastest processes is the first in the movie.")
    optpar.add_option("-M", "--separate-movie-files", action="store_true", dest="separate_movie_files", default=False, help="Do not write the movie into a single file but use a separate POSCAR for every state. These are created in the \"movies\" directory. Only useful in conjunction with --movie.")
    optpar.add_option("-n", "--no-submit", action="store_true", dest="no_submit", default=False,help="don't submit searches; only register finished results")
    (options, args) = optpar.parse_args()
    if len(args) > 1:
        print("akmc.py takes only one positional argument")
    sys.argv = sys.argv[0:1]
    if len(args) == 1:
        sys.argv += args
        #always run from the directory where the config file is
        #os.chdir(os.path.dirname(args[0]))

    #XXX: config is ugly as it finds out where the config file is directly from
    #     sys.argv instead of being passed it.
    #import sys
    if len(sys.argv) > 1:
        config.init(sys.argv[-1])
    else:
        config.init()
    #set options.path_root to be where the config file is if given as an arg
    if config.path_root.strip() == '.' and len(args) == 1:
        config.path_root = os.path.abspath(os.path.dirname(args[0]))
        os.chdir(config.path_root)

    if options.no_submit:
        config.comm_job_buffer_size = 0

    #setup logging
    logging.basicConfig(level=logging.DEBUG,
            filename=os.path.join(config.path_results, "akmc.log"),
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

    lock = locking.LockFile(str((Path.cwd() / config.path_results / "lockfile").absolute()))

    # Some options are mutually exclusive. Let's check them now.
    exclusive_options = {}
    if len(options.movie_type) > 0:
        exclusive_options['movie_type'] = True
    else:
        exclusive_options['movie_type'] = False
    exclusive_options['print_status'] = options.print_status
    exclusive_options['reset'] = options.reset

    if sum(exclusive_options.values()) > 1:
        offending_options = [ k for k,v in list(exclusive_options.items()) if v ]
        optpar.error("Options %s are mutually exclusive" % ", ".join(offending_options))

    if len(options.movie_type) > 0:
        states = get_statelist(config.main_temperature / 11604.5)
        movie.make_movie(options.movie_type, config.path_root, states,
                         options.separate_movie_files)
        sys.exit(0)

    # From the config file: The Novotny and C&V (ASKMC) methods should not be used together.
    if config.sb_on and config.askmc_on:
        logger.error("Both superbasin methods should not be used at the same time")
        sys.exit(1)

    # From the config file: dynamics confidence can not be used with MinMode searches
    if config.saddle_method=="min_mode" and config.akmc_confidence_scheme=="dynamics":
        logger.error("Dynamics confidence scheme can not be used with MinMode saddle searches")
        sys.exit(1)

    if options.print_status:
        states = get_statelist(config.main_temperature / 11604.5)
        start_state_num, time, previous_state_num, first_run, previous_temperature =\
            get_akmc_metadata()
        current_state = states.get_state(start_state_num)
        if config.sb_on:
            sb_scheme = get_superbasin_scheme(states)
            sb = sb_scheme.get_containing_superbasin(current_state)
        else:
            sb = None
        print()
        print("General")
        print("-------")
        if not sb:
            print("Current state:", start_state_num)
        else:
            print("Current state:", start_state_num, "in superbasin", sb.id)
        print("Number of states:",states.get_num_states())
        print("Time simulated: %.3e seconds" % time)
        print()

        print("Current State")
        print("-------------")
        if not sb:
            print("Confidence: %.4f" % current_state.get_confidence())
        else:
            print("Superbasin Confidence: %.4f" % sb.get_confidence())
            non_ignored_states = set(sb._get_filtered_states())
            for s in sb.states:
                if s in non_ignored_states:
                    ignore_string = ""
                else:
                    ignore_string = " (no exit from superbasin found)"
                print("       %4i: %.4f%s" % (s.number, s.get_confidence(sb), ignore_string))
        print("Unique Saddles:", current_state.get_unique_saddle_count())
        print("Good Saddles:", current_state.get_good_saddle_count())
        print("Bad Saddles:", current_state.get_bad_saddle_count())
        print("Percentage bad saddles: %.1f" % (float(current_state.get_bad_saddle_count())/float(max(current_state.get_bad_saddle_count() + current_state.get_good_saddle_count(), 1)) * 100))
        print()

        comm = communicator.get_communicator()
        print("Saddle Searches")
        print("---------------")
        print("Searches in queue:", comm.get_queue_size())
        print()
        if config.sb_on:
            print("Superbasins")
            print("-----------")
            for i in sb_scheme.superbasins:
                print("%s: %s" % (i.id, i.state_numbers))
        sys.exit(0)
    elif options.reset:
        if options.force:
            res = 'y'
        else:
            res = input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
                def attempt_removal(thing):
                    if thing is None:
                        return
                    if os.path.isdir(thing):
                        shutil.rmtree(thing)
                        os.mkdir(thing)
                        os.removedirs(thing)
                    elif os.path.isfile(thing):
                        os.remove(thing)
                rmthings = [config.path_jobs_out,
                            config.path_jobs_in,
                            config.path_incomplete,
                            config.path_states,
                            config.path_scratch,
                            config.kdb_name,
                            config.kdb_scratch_path,
                            config.sb_path,
                            config.sb_recycling_path,
                            config.debug_results_path,
                            os.path.join(config.path_root, "searchdata"),
                            os.path.join(config.path_results, "askmc_data.txt"),
                            os.path.join(config.path_results, "searches.log"),
                            os.path.join(config.path_results, "dynamics.txt"),
                            os.path.join(config.path_results, "info.txt"),
                            os.path.join(config.path_results, "akmc.log"),
                            os.path.join(config.path_results, "jobs.tbl"),
                            os.path.join(config.path_root, "results"),
                            os.path.join(config.path_root, "prng.pkl"),
                            os.path.join(config.path_root, "explorer.pickle"),
                            os.path.join(config.path_root, "temperatures.dat"),
                            os.path.join(config.path_root, "client.log"),
                            os.path.join(config.path_root, "lockfile"),
                            ]
                for thing in rmthings:
                    attempt_removal(thing)
                if not options.quiet:
                    print("Reset")
                sys.exit(0)
        else:
            print("Not resetting")
            sys.exit(1)

    elif options.restart:
        string_sb_clear = ""

        if options.force:
            res = 'y'
        else:
            res = input("Are you sure you want to restart (remove dynamics.txt, info.txt and akmc.log)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':

            # remove akmc data that are specific for a trajectory
            dynamics_path = os.path.join(config.path_results, "dynamics.txt")
            info_path = os.path.join(config.path_results, "info.txt")
            log_path = os.path.join(config.path_results, "akmc.log")
            jobs_path = os.path.join(config.path_results, "jobs.tbl")
            for i in [info_path, dynamics_path, log_path, jobs_path]:
                if os.path.isfile(i):
                    os.remove(i)

            if config.sb_on:
                if options.force:
                    res = 'y'
                else:
                    res = input("Should the superbasins be removed? (y/N) ").lower()

                # remove superbasin data (specific for a trajectory)
                if len(res)>0 and res[0] == 'y':
                    # remove directory superbasins
                    if os.path.isdir(config.sb_path):
                        shutil.rmtree(config.sb_path)
                        #XXX: ugly way to remove all empty directories containing this one
                        os.mkdir(config.sb_path)
                        os.removedirs(config.sb_path)

                    # remove superbasins files from states dirctories
                    state_dirs = os.listdir(config.path_states)
                    for i in state_dirs:
                        if i != 'state_table':
                            superbasin_file = os.path.join(config.path_states, i)
                            superbasin_file = os.path.join(superbasin_file, config.sb_state_file)
                            if os.path.isfile(superbasin_file):
                                os.remove(superbasin_file)

                    string_sb_clear = " with directory 'superbasins' and files named '"
                    string_sb_clear += str(config.sb_state_file) + "' removed"

            if not options.quiet:
                print("Restart"+string_sb_clear+".")
            sys.exit(0)
        else:
            print("Not restarting")
            sys.exit(1)

    if lock.aquirelock():
        if options.continuous or config.comm_type == 'mpi':
            # define a wait method.
            if config.comm_type == 'mpi':
                from eon.mpiwait import mpiwait
                wait = mpiwait
            elif options.continuous:
                if config.comm_type == "local":
                    # In local, everything is synchronous, so no need to wait here.
                    wait = lambda: None
                else:
                    wait = lambda: sleep(10.0)
            else:
                raise RuntimeError("You have found a bug in EON!")
            # Run a specified number of steps or forever.
            steps = 0
            while True:
                steps = akmc(config, steps)
                if (config.akmc_max_kmc_steps > 0 and
                    steps >= config.akmc_max_kmc_steps):
                    break
                wait()
            # In MPI mode we need to signal exit to all processes.
            # TODO: This is the sledgehammer method, it would be cleaner to
            #       communicate to all clients that they should exit.
            if config.comm_type == 'mpi':
                from mpi4py import MPI
                MPI.COMM_WORLD.Abort(0)
        else:
            akmc(config)
    else:
        logger.info("Server is locked by pid %i" % lock.pid)
        sys.exit(1)

if __name__ == '__main__':
    main()
