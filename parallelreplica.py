#!/usr/bin/env python

import ConfigParser
from cStringIO import StringIO
import logging
import logging.handlers
logger = logging.getLogger('pr')
import optparse
import os
import shutil
import sys

import atoms
import communicator
import config
import io
import locking
import prstatelist

def parallelreplica():
    # First of all, does the root directory even exist?
    if not os.path.isdir(config.path_root):
        logger.critical("Root directory does not exist")
        sys.exit(1)

    # load metadata
    start_state_num, time, wuid = get_pr_metadata()
    logger.info("current time is %e", time)
    states = get_statelist() 
    current_state = states.get_state(start_state_num)

    # get communicator
    comm = communicator.get_communicator()

    # Register all the results. There is  no need to ever discard found
    # processes like we do with akmc. There is no confidence to calculate.
    num_registered, transition = register_results(comm, current_state, states)

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
    next_state.zero_search_count()
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
    return prstatelist.PRStateList(config.path_states, 
                               config.comp_eps_e, 
                               config.comp_eps_r, 
                               config.comp_use_identical, 
                               config.debug_list_search_results, 
                               initial_state_path)

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
    num_to_make = max(config.comm_search_buffer_size - num_in_buffer, 0)
    logger.info("making %i searches" % num_to_make)
    
    if num_to_make == 0:
        return wuid
    
    parameters_path = os.path.join(config.path_root, "parameters.dat")
    searches = []
    
    invariants = {}

    reactIO = StringIO()
    io.savecon(reactIO, reactant)
    #invariants['reactant_passed.con']=reactIO
    
    f = open(parameters_path)
    invariants['parameters_passed.dat'] = StringIO(''.join(f.readlines()))
    f.close()
    #Merge potential files into invariants
    invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

    searches = []
    for i in range(num_to_make):
        search = {}
        # The search dictionary contains the following key-value pairs:
        # id - CurrentState_WUID
        search['id'] = "%d_%d" % (current_state.number, wuid)
        search['reactant_passed.con']  = reactIO
        searches.append(search) 
        wuid += 1

    try:
        comm.submit_jobs(searches, invariants)
        logger.info( str(num_to_make) + " searches created") 
    except:
        logger.exception("Failed to submit searches.")
    return wuid

def register_results(comm, current_state, states):
    logger.info("registering results")
    if os.path.isdir(config.path_searches_in):
        shutil.rmtree(config.path_searches_in)    
    os.makedirs(config.path_searches_in)
    
    #Function used by communicator to determine whether to discard a result
    def keep_result(name):
        return True

    transition = None
    num_registered = 0
    for result in comm.get_results(config.path_searches_in, keep_result): 
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
        if result['results']['termination_reason'] == 1:
            process_id = state.add_process(result)
            time = result['results']['transition_time']
            logger.info("found transition with time %.3e", time)
            # if this is a faster transition record it
            if not transition or transition['time'] > time:
                transition = {'process_id':process_id, 'time':time}

        num_registered += 1
        state.inc_search_count()
        
    logger.info("%i (result) searches processed", num_registered)

    if transition:
        num_searches = state.get_search_count()
        logger.info("fastest time: %e", transition['time'])
        transition['time'] *= num_searches
        logger.info("real time: %e", transition['time'])

    return num_registered, transition

def main():
    optpar = optparse.OptionParser(usage="usage: %prog [options] config.ini")
    optpar.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False,help="only write to the log file")
    optpar.add_option("-n", "--no-submit", action="store_true", dest="no_submit", default=False,help="don't submit searches; only register finished results")
    optpar.add_option("-R", "--reset", action="store_true", dest="reset", default = False, help="reset the simulation, discarding all data")
    options, args = optpar.parse_args()

    config.init()

    if options.no_submit:
        config.comm_search_buffer_size = 0

    if options.reset:
        res = raw_input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
                rmdirs = [config.path_searches_out, config.path_searches_in, config.path_states,
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
