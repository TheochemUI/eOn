#!/usr/bin/env python

##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

import ConfigParser
from cStringIO import StringIO
import logging
import logging.handlers
logger = logging.getLogger('pr')
import numpy
numpy.seterr(all='raise')
import optparse
import os
import shutil
import sys

import communicator
import config
import fileio as io
import locking
import prstatelist

def parallelreplica():
    # First of all, does the root directory even exist?
    if not os.path.isdir(config.path_root):
        logger.critical("Root directory does not exist")
        sys.exit(1)

    # load metadata
    start_state_num, time, wuid = get_pr_metadata()
    logger.info("simulation time is %e", time)
    states = get_statelist() 
    current_state = states.get_state(start_state_num)

    # get communicator
    comm = communicator.get_communicator()

    # Register all the results. There is  no need to ever discard found
    # processes like we do with akmc. There is no confidence to calculate.
    num_registered, transition = register_results(comm, current_state, states)
    logger.info("time in current state is %e", current_state.get_time())

    if transition:
        current_state, previous_state = step(time, current_state, states, transition)
        time += transition['time']

    wuid = make_searches(comm, current_state, wuid)

    # Write out metadata. XXX:ugly
    metafile = os.path.join(config.path_results, 'info.txt')
    parser = ConfigParser.RawConfigParser() 
    write_pr_metadata(parser, current_state.number, time, wuid)
    parser.write(open(metafile, 'w'))

def step(current_time, current_state, states, transition):
    next_state = states.get_product_state(current_state.number, transition['process_id'])
    next_state.zero_time()
    dynamics = io.Dynamics(os.path.join(config.path_results, "dynamics.txt"))
    proc = current_state.get_process(transition['process_id'])
    dynamics.append(current_state.number, transition['process_id'],
                    next_state.number, transition['time'], transition['time']+current_time, 0, 0)

    previous_state = current_state
    current_state = next_state

    logger.info("currently in state %i", current_state.number)

    return current_state, previous_state

def get_statelist():
    initial_state_path = os.path.join(config.path_root, 'reactant.con')
    return prstatelist.PRStateList(initial_state_path)

def get_pr_metadata():
    if not os.path.isdir(config.path_results):
        os.makedirs(config.path_results)
    metafile = os.path.join(config.path_results, 'info.txt')
    parser = ConfigParser.SafeConfigParser()
    if os.path.isfile(metafile):
        parser.read(metafile)
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

def make_searches(comm, current_state, wuid):
    reactant = current_state.get_reactant()
    #XXX:what if the user changes the bundle size?
    num_in_buffer = comm.get_queue_size()*config.comm_job_bundle_size
    logger.info("%i searches in the queue" % num_in_buffer)
    num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
    logger.info("making %i searches" % num_to_make)

    if num_to_make == 0:
        return wuid

    searches = []

    invariants = {}

    reactIO = StringIO()
    io.savecon(reactIO, reactant)
    #invariants['reactant_passed.con']=reactIO

    #Merge potential files into invariants
    #XXX: Should this be in our "science" maybe the communicator should
    #     handle this.
    invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

    searches = []
    for i in range(num_to_make):
        search = {}
        search['id'] = "%d_%d" % (current_state.number, wuid)
        search['reactant_passed.con']  = reactIO
        ini_changes = [
                        ('Main', 'job', 'parallel_replica'),
                        ('Main', 'random_seed',
                            str(int(numpy.random.random()*10**9))),
                      ]
        search['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)
        searches.append(search)
        wuid += 1

    comm.submit_jobs(searches, invariants)
    logger.info( str(num_to_make) + " searches created")
    return wuid

def register_results(comm, current_state, states):
    logger.info("registering results")
    if os.path.isdir(config.path_jobs_in):
        shutil.rmtree(config.path_jobs_in)
    os.makedirs(config.path_jobs_in)

    #Function used by communicator to determine whether to discard a result
    def keep_result(name):
        return True

    transition = None
    num_registered = 0
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
        id = int(result['name'].split("_")[1]) + result['number']

        state = states.get_state(state_num)

        #read in the results
        result['results'] = io.parse_results(result['results.dat'])
        if result['results']['transition_found'] == 1:
            result['results']['transition_time_s'] += state.get_time()
            time = result['results']['transition_time_s']
            process_id = state.add_process(result)
            logger.info("found transition with time %.3e", time)
            if not transition and current_state.number==state.number:
                transition = {'process_id':process_id, 'time':time}
            state.zero_time()
        else:
            state.inc_time(result['results']['transition_time_s'])
        num_registered += 1

    logger.info("%i (result) searches processed", num_registered)

    return num_registered, transition

def main():
    optpar = optparse.OptionParser(usage="usage: %prog [options] config.ini")
    optpar.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False,help="only write to the log file")
    optpar.add_option("-n", "--no-submit", action="store_true", dest="no_submit", default=False,help="don't submit searches; only register finished results")
    optpar.add_option("-R", "--reset", action="store_true", dest="reset", default = False, help="reset the simulation, discarding all data")
    options, args = optpar.parse_args()

    if len(args) > 1:
        print "takes only one positional argument"
    sys.argv = sys.argv[0:1]
    if len(args) == 1:
        sys.argv += args
    if len(sys.argv) > 1:
        config.init(sys.argv[-1])
    else:
        config.init()
    #set options.path_root to be where the config file is if given as an arg
    if config.path_root.strip() == '.' and len(args) == 1:
        config.path_root = os.path.abspath(os.path.dirname(args[0]))
        os.chdir(config.path_root)

    if config.comm_job_bundle_size != 1:
        print "error: Parallel Replica only supports a bundle size of 1"
        sys.exit(1)

    if options.no_submit:
        config.comm_job_buffer_size = 0

    if options.reset:
        res = raw_input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
            rmdirs = [config.path_jobs_out, config.path_jobs_in, config.path_states,
                    config.path_scratch]
            if config.debug_keep_all_results:
                rmdirs.append(os.path.join(config.path_root, "old_searches"))
            for i in rmdirs:
                if os.path.isdir(i):
                    shutil.rmtree(i)
                    #XXX: ugly way to remove all empty directories containing this one
                    os.mkdir(i)
                    os.removedirs(i)

            dynamics_path = os.path.join(config.path_results, "dynamics.txt")
            info_path = os.path.join(config.path_results, "info.txt")
            log_path = os.path.join(config.path_results, "pr.log")
            for i in [info_path, dynamics_path, log_path]:
                if os.path.isfile(i):
                    os.remove(i)

            print "Reset."
        sys.exit(0)

    #setup logging
    logging.basicConfig(level=logging.DEBUG,
            filename=os.path.join(config.path_results, "pr.log"),
            format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
            datefmt="%F %T")
    logging.raiseExceptions = False

    if not options.quiet:
        rootlogger = logging.getLogger('')
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter("%(message)s")
        console.setFormatter(formatter)
        rootlogger.addHandler(console)

    lock = locking.LockFile(os.path.join(config.path_results, "lockfile"))

    if lock.aquirelock():
        parallelreplica()
    else:
        logger.warning("couldn't get lock")
        sys.exit(1)

if __name__ == '__main__':
    main()
