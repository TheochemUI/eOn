#!/usr/bin/env python

import math
import sys
import ConfigParser
import os.path
import shutil
import os
import time as unix_time
import optparse
import logging
import logging.handlers
import numpy
numpy.seterr(all='raise')
import pickle

import locking
import communicator
import statelist
import displace
import io
import recycling
import superbasinscheme
import kdb

#hi from maggi

def main(): 
     
    # Here's what this does:
    # 1) Read in the state of our calculation from last time
    # 2) Initialize necessary data structures (statelist, communicator, displace)
    # 3) Get any results that have come in
    # 4) Possibly take a KMC step
    # 5) Make new work units
    # 6) Write out the state of the simulation
    
    
    # Define constants. <rye> should this be a config option? </rye>    
    kT = config.akmc_temperature/11604.5 #in eV
    
    # First of all, does the root directory even exist?
    if not os.path.isdir(config.path_root):
        logger.critical("Root directory does not exist, as such the reactant cannot exist. Exiting...")
        sys.exit(1)
    
    # Load metadata, the state list, and the current state.
    start_state_num, time, wuid, searchdata = get_akmc_metadata()
    states = get_statelist(kT) 
    current_state = states.get_state(start_state_num)

    # If recycling is being used, initialize it.
    if config.recycling_on:
        recycling_start, previous_state_num = get_recycling_metadata()
        previous_state = states.get_state(previous_state_num)
    
    # If kdb is being used, initialize it.
    if config.kdb_on:
        kdber = kdb.KDB(config.kdb_path, config.kdb_querypath, config.kdb_addpath)

    # Create the communicator object.
    comm = get_communicator()

    # Handle any results returned through the communicator.
    register_results(comm, current_state, states, searchdata, kdber = kdber if config.kdb_on else None)

    # Take a KMC step, if it's time.
    current_state, previous_state, time = kmc_step(current_state, states, time) 
            
    # If we took a step, cancel old jobs and start the kdbquery.
    if current_state.number != start_state_num:
        num_cancelled = comm.cancel_state(start_state_num)
        logger.info("cancelled %i workunits from state %i", num_cancelled, start_state_num)
        if config.kdb_on:
            kdber.query(current_state, wait = True)
    
    # Create new work.
    if config.recycling_on:
        recycler = recycling.Recycling(previous_state, current_state, 
                        config.recycling_move_distance, recycling_start)
    wuid = make_searches(comm, current_state, wuid, searchdata = searchdata, recycler = recycler if config.recycling_on else None, kdber = kdber if config.kdb_on else None)
    
    # Write out metadata. XXX:ugly
    metafile = os.path.join(config.path_results, 'info.txt')
    parser = ConfigParser.RawConfigParser() 
    write_akmc_metadata(parser, current_state.number, time, wuid, searchdata)
    if config.recycling_on:
        write_recycling_metadata(parser, recycler.process_num, previous_state.number)
    parser.write(open(metafile, 'w')) 
    


def get_akmc_metadata():
    if not os.path.isdir(config.path_results):
        os.makedirs(config.path_results)
    # read in metadata
    # do we want custom metadata locations?
    metafile = os.path.join(config.path_results, 'info.txt')
    parser = ConfigParser.SafeConfigParser() 
    # give our parser some default values in case our metafile doesn't exist
    try:
        parser.read(metafile)
        start_state_num = parser.getint("Simulation Information",'current_state')
        time = parser.getfloat("Simulation Information", 'time_simulated') 
        wuid = parser.getint("aKMC Metadata", 'wu_id') 
        searchdata = eval(parser.get("aKMC Metadata", 'searchdata'))
    except:
        time = 0
        start_state_num = 0
        wuid = 0
        searchdata = {}

    if config.debug_random_seed:
        try:
            parser = ConfigParser.RawConfigParser()
            parser.read(metafile)
            seed = parser.get("aKMC Metadata", "random_state")
            from numpy import array, uint32
            numpy.random.set_state(eval(seed))
            logger.debug("Set random state from pervious run's state")
        except:
            numpy.random.seed(config.debug_random_seed)
            logger.debug("Set random state from seed")

    return start_state_num, time, wuid, searchdata



def write_akmc_metadata(parser, current_state_num, time, wuid, searchdata):
    parser.add_section('aKMC Metadata')
    parser.add_section('Simulation Information')
    parser.set('aKMC Metadata', 'wu_id', str(wuid))
    parser.set('aKMC Metadata', 'searchdata', repr(searchdata))
    parser.set('Simulation Information', 'time_simulated', str(time))
    parser.set('Simulation Information', 'current_state', str(current_state_num))
    if config.debug_random_seed:
        parser.set('aKMC Metadata', 'random_state', repr(numpy.random.get_state()))


def get_recycling_metadata():
    metafile = os.path.join(config.path_results, 'info.txt')
    parser = ConfigParser.SafeConfigParser() 
    try:
        parser.read(metafile)
        recycling_start = parser.getint('Recycling', 'start_process_id')
        previous_state_num = parser.getint('Recycling', 'previous_state')
    except:
        recycling_start = 0
        previous_state_num = 0 
    return recycling_start, previous_state_num



def write_recycling_metadata(parser, recycling_start, previous_state_num):
    parser.add_section('Recycling')
    parser.set('Recycling', 'previous_state',str(previous_state_num))  
    parser.set('Recycling', 'start_process_id',str(recycling_start)) 



def get_statelist(kT):
    initial_state_path = os.path.join(config.path_root, 'reactant.con') 
    return statelist.StateList(config.path_states, kT, config.akmc_thermal_window, config.akmc_max_thermal_window, config.comp_eps_e, config.comp_eps_r, config.comp_use_identical, initial_state_path) #might 



def get_communicator():
    if config.comm_type=='boinc':
        comm = communicator.BOINC(config.path_scratch, config.comm_boinc_project_dir, 
                config.comm_boinc_wu_template_path, config.comm_boinc_re_template_path,
                config.comm_boinc_appname, config.comm_boinc_results_path)
    elif config.comm_type=='cluster':
        comm = communicator.Cluster()
    elif config.comm_type=='local':
        comm = communicator.Local(config.path_scratch, config.comm_local_client, config.comm_local_ncpus)
    else:
        logger.error(str(config.comm_type)+" is an unknown communicator.")
        raise ValueError()
    return comm



def register_results(comm, current_state, states, searchdata, kdber = None):
    logger.info("registering results")
    t1 = unix_time.time()
    if os.path.isdir(config.path_searches_in):
        shutil.rmtree(config.path_searches_in)    
    os.makedirs(config.path_searches_in)
    results = comm.get_results(config.path_searches_in) 
    num_registered = 0
    for i in results:
        result_path = os.path.join(config.path_searches_in, i)
        if config.debug_keep_all_results:
            #XXX: We should only do these checks once to speed things up, but at the same time
            #debug options don't have to be fast
            save_path = os.path.join(config.path_root, "old_searches")
            if not os.path.isdir(save_path):
                os.mkdir(save_path)
            shutil.copytree(result_path, os.path.join(save_path, i))
        state_num = int(i.split("_")[0])
        try:
            result_data = io.parse_results_dat(os.path.join(result_path, 'results.dat'))
        except KeyError, (foo, name):
            logger.warning("Search %s's results file is missing %s" % ( i ,str(name)))
            continue
        except TypeError, (foo, name): #we need to use python 2.6
            logger.warning("%s has the wrong type in the results file of search %s" % (str(name), i))
            continue
        except IOError:
            logger.warning("Search %s did not return a results.dat" % i)
            continue

        # Store information about the search into result_data for the search_results.txt file in the state directory.
        result_data['search_id'] = int(i.split("_")[1])
        result_data['search_type'] = searchdata[i.split("_")[1] + "type"]
        
        # Remove used information from the searchdata metadata.
        del searchdata[i.split("_")[1] + "type"]
        
        if result_data['termination_reason'] == 0:
            process_id = states.get_state(state_num).add_process(result_path, result_data)
            if current_state.get_confidence() > config.akmc_confidence:
                if config.kdb_on:
                    logger.debug("Adding relevant processes to kinetic database.")
                    for process_id in states.get_state(state_num).get_process_ids():
                        output = kdber.add_process(states.get_state(state_num), process_id)
                        logger.debug("kdbaddpr.pl: %s" % output)
                break
        else:
            states.get_state(state_num).register_bad_saddle(result_path, config.debug_keep_bad_saddles, result_data)
        num_registered += 1
    t2 = unix_time.time()
    logger.info("%i results processed", num_registered)
    logger.debug("%.1f results per second", (num_registered/(t2-t1)))
    return num_registered



def get_superbasin_scheme(states):
    if config.sb_scheme == 'transition_counting':
        superbasining = superbasinscheme.TransitionCounting(config.sb_path, states, config.sb_tc_ntrans)
    return superbasining



def kmc_step(current_state, states, time):
    t1 = unix_time.time()
    previous_state = current_state 
    dynamics_file = open(os.path.join(config.path_results, "dynamics.txt"), 'a')
    start_state_num = current_state.number
    if config.sb_on:
        superbasining = get_superbasin_scheme(states)
    while current_state.get_confidence() >= config.akmc_confidence:
        if config.sb_on:
            sb = superbasining.get_containing_superbasin(current_state)

        if config.sb_on and sb:
            time, next_state = sb.step(current_state, states.get_product_state)
        else:
            rate_table = current_state.get_ratetable()
            if len(rate_table) == 0:
                logger.error("No processes in rate table, but confidence has been reached")
            ratesum = 0.0
            for i in range(len(rate_table)):
                ratesum += rate_table[i][1]
            
            u = numpy.random.random_sample()
            p = 0.0
             
            for i in range(len(rate_table)):
                p += rate_table[i][1]/ratesum
                if p>u:
                    nsid = i 
                    break
            else:
                logger.warning("Warning: failed to select rate. p = " + str(p))
                break

            next_state = states.get_product_state(current_state.number, rate_table[nsid][0])
            time -= math.log(numpy.random.random_sample())/ratesum

        if config.sb_on:
            superbasining.register_transition(current_state, next_state)    
        
        print >> dynamics_file, next_state.number, time
        logger.info("stepped from state %i to state %i", current_state.number, next_state.number)
        
        if config.recycling_on:
            recycling_start = 0
        previous_state = current_state
        current_state = next_state

    if config.sb_on:
        superbasining.write_data()

    logger.info("currently in state %i with confidence %.3f", current_state.number, current_state.get_confidence())
    dynamics_file.close()
    t2 = unix_time.time()
    logger.debug("KMC finished in " + str(t2-t1) + " seconds")
    return current_state, previous_state, time



def get_displacement(reactant):
    if config.disp_type == 'random':
        disp = displace.Random(reactant, config.disp_size, config.disp_radius)
    elif config.disp_type == 'undercoordinated':
        disp = displace.Undercoordinated(reactant, config.disp_max_coord, config.disp_size, config.disp_radius)
    elif config.disp_type == 'leastcoordinated':
        disp = displace.Leastcoordinated(reactant, config.disp_size, config.disp_radius)
    else:
        raise ValueError()
    return disp



def make_searches(comm, current_state, wuid, searchdata, kdber = None, recycler = None):
    reactant = current_state.get_reactant()
    disp = get_displacement(reactant)
    num_in_buffer = comm.get_queue_size()
    logger.info("%i jobs in the queue" % num_in_buffer)
    num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
    logger.info("making %i jobs" % num_to_make)
    parameters_path = os.path.join(config.path_root, "parameters.dat")
    
    if os.path.isdir(config.path_searches_out):
        #Clean out directory
        shutil.rmtree(config.path_searches_out)

    t1 = unix_time.time()
    os.makedirs(config.path_searches_out)
    jobpaths = []
    for i in range(num_to_make):
        job_dir = os.path.join(config.path_searches_out, str(current_state.number)+"_"+str(wuid))
        os.mkdir(job_dir)
        shutil.copy(current_state.reactant_path, os.path.join(job_dir, "reactant_passed.con"))
        shutil.copy(parameters_path, os.path.join(job_dir, "parameters_passed.dat"))
        
        # Do we want to do recycling? If yes, try. If we fail to recycle, we move to the next case
        if (config.recycling_on and current_state.number is not 0) and recycler.make_suggestion(job_dir):
            logger.debug('Recycled a saddle')
            searchdata[str(wuid)+"type"] = "recycling"
        elif config.kdb_on and kdber.make_suggestion(job_dir):
            logger.info('Made a KDB suggestion')
            searchdata[str(wuid)+"type"] = "kdb"
        else:
            disp.make_displacement(job_dir) 
            searchdata[str(wuid)+"type"] = "random"
        wuid += 1
        jobpaths.append(job_dir)

    t1 = unix_time.time()
    comm.submit_searches(jobpaths)
    t2 = unix_time.time()
    logger.info( str(num_to_make) + " jobs created") 
    logger.debug( str(num_to_make/(t2-t1)) + " jobs per second")
    
    return wuid



#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

if __name__ == '__main__':
    optpar = optparse.OptionParser(usage = "usage: %prog [options] config.ini")
    optpar.add_option("-R", "--reset", action="store_true", dest="reset", default = False, help="reset the aKMC simulation, discarding all data")
    optpar.add_option("-s", "--status", action="store_true", dest="print_status", default = False, help = "print the status of the simulation and currently running jobs")
    optpar.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False,help="only write to the log file")
    (options, args) = optpar.parse_args()
 
    if len(args) > 1:
        print "akmc.py takes only one positional argument"
    sys.argv = sys.argv[0:1]
    if len(args) == 1:
        sys.argv += args
        #always run from the directory where the config file is
        os.chdir(os.path.dirname(args[0]))

    #XXX: config is ugly as it finds out where the config file is directly from 
    #     sys.argv instead of being passed it.
    import config

    #setup logging
    logging.basicConfig(level=logging.DEBUG,
            filename=os.path.join(config.path_results, "akmc.log"),
            format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
            datefmt="%F %T")
    logging.raiseExceptions = False
    logger = logging.getLogger('akmc')

    if not options.quiet:
        rootlogger = logging.getLogger('')
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter("%(message)s")
        console.setFormatter(formatter)
        rootlogger.addHandler(console)

    lock = locking.LockFile(os.path.join(config.path_results, "lockfile"))

    #XXX: Some options are mutally exclusive. This should be handled better.
    if options.print_status:
        states = get_statelist(config.akmc_temperature/11604.5)
        start_state_num, time, wuid, searchdata = get_akmc_metadata()

        print
        print "General"
        print "-------"
        print "Current state:", start_state_num
        print "Number of states:",states.get_num_states()  
        print "Time simulated: %.3e seconds" % time
        print

        current_state = states.get_state(start_state_num)
        print "Current State"
        print "-------------"
        print "Confidence: %.3f" % current_state.get_confidence()
        print "Unique Saddles:", current_state.get_unique_saddle_count()
        print "Good Saddles:", current_state.get_good_saddle_count()
        print "Bad Saddles:", current_state.get_bad_saddle_count()
        print "Percentage bad saddles: %.1f" % (float(current_state.get_bad_saddle_count())/float(max(current_state.get_bad_saddle_count() + current_state.get_good_saddle_count(), 1)) * 100)
        print 

        comm = get_communicator()
        print "Saddle Searches"
        print "---------------" 
        #print "Searches currently running:", comm.thingy()
        print "Searches in queue:", comm.get_queue_size() 
        print

        sb = get_superbasin_scheme(states)
        print "Superbasins"
        print "-----------"
        for i in sb.superbasins:
            print i.state_numbers 
        
        sys.exit(0)
    elif options.reset:
        res = raw_input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
                rmdirs = [config.path_searches_out, config.path_searches_in, config.path_states,
                        config.path_scratch]
                if config.kdb_on:
                    rmdirs.append(config.kdb_path) #Do we want to save kdb data across resets? No.
                if config.sb_on:
                    rmdirs.append(config.sb_path)
                for i in rmdirs:
                    if os.path.isdir(i):
                        shutil.rmtree(i)
                        #XXX: ugly way to remove all empty directories containing this one
                        os.mkdir(i)
                        os.removedirs(i)
                
                dynamics_path =os.path.join(config.path_results, "dynamics.txt")  
                info_path = os.path.join(config.path_results, "info.txt") 
                for i in [info_path, dynamics_path]:
                    if os.path.isfile(i):
                        os.remove(i)
                
                print "Reset."
                sys.exit(0)
        else:
            print "Not resetting."
            sys.exit(1)

    if lock.aquirelock():
        main()
    else:
        logger.info("the server is locked by pid %i" % lock.pid)
        sys.exit(1)
