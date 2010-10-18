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
import atoms
import recycling
import superbasinscheme
import askmc
import kdb
import movie

import StringIO


def main(): 
     
    # Here's what this does:
    # 1) Read in the state of our calculation from last time
    # 2) Initialize necessary data structures (statelist, communicator, displace)
    # 3) Get any results that have come in
    # 4) Possibly take a KMC step
    # 5) Make new work units
    # 6) Write out the state of the simulation    
    
    # Define constants. 
    kT = config.akmc_temperature/11604.5 #in eV
    
    # First of all, does the root directory even exist?
    if not os.path.isdir(config.path_root):
        logger.critical("Root directory does not exist, as such the reactant cannot exist. Exiting...")
        sys.exit(1)
    
    # Load metadata, the state list, and the current state.
    start_state_num, time, wuid, searchdata, previous_state_num, first_run = get_akmc_metadata()
    states = get_statelist(kT) 
    current_state = states.get_state(start_state_num)

    # If kdb is being used, initialize it.
    if config.kdb_on:
        kdber = kdb.KDB(config.kdb_path, config.kdb_querypath, config.kdb_addpath)

    # If the Novotny-based superbasining scheme is being used, initialize it.
    if config.sb_on:
        superbasining = get_superbasin_scheme(states)
    
    # Create the communicator object.
    comm = get_communicator()

    # Handle any results returned through the communicator.
    register_results(comm, current_state, states, searchdata)
    
    # Add processes to the database if we've reached confidence. 
    if current_state.get_confidence() >= config.akmc_confidence:
        if config.kdb_on:
            logger.debug("Adding relevant processes to kinetic database.")
            for process_id in current_state.get_process_ids():
                output = kdber.add_process(current_state, process_id)
                logger.debug("kdbaddpr.pl: %s" % output)

    # Take a KMC step, if it's time.
    if config.sb_on:
        pass_superbasining = superbasining
    else:
        pass_superbasining = None
    current_state, previous_state, time = kmc_step(current_state, states, time, kT, superbasining = pass_superbasining) 
            
    # If we took a step, cancel old jobs and start the kdbquery.
    if current_state.number != start_state_num:
        num_cancelled = comm.cancel_state(start_state_num)
        logger.info("cancelled %i workunits from state %i", num_cancelled, start_state_num)
        if config.kdb_on:
            kdber.query(current_state, rhsco = config.kdb_rhsco, wait = config.kdb_wait)
    
    # If this is the first execution of akmc.py for this simulation, run kdbquery if it's on.
    if first_run:
        if config.kdb_on:
            kdber.query(current_state, rhsco = config.kdb_rhsco, wait = config.kdb_wait)
    
    # Create new work.
    recycler = None
    # If we *just* want to do simple recycling
    if config.recycling_on and not config.sb_recycling_on:
        recycler = recycling.Recycling(states, previous_state, current_state, 
                        config.recycling_move_distance, config.recycling_save_sugg)
    # Is super-basin recycling on? If so, figure out the superbasin method being used.
    if config.sb_recycling_on:
        if config.sb_on:
            sb_recycler = recycling.SB_Recycling(states, previous_state, current_state,
                          config.recycling_move_distance, config.recycling_save_sugg,
                          config.sb_recycling_path, sb_type = "novotny",
                          superbasining = superbasining)
        elif config.askmc_on:
            sb_recycler = recycling.SB_Recycling(states, previous_state, current_state,
                          config.recycling_move_distance, config.recycling_save_sugg,
                          config.sb_recycling_path, sb_type = "askmc",
                          superbasining = None)
        else:
            sb_recycler = recycling.SB_Recycling(states, previous_state, current_state,
                          config.recycling_move_distance, config.recycling_save_sugg,
                          config.sb_recycling_path, sb_type = None,
                          superbasining = None)
        # If there's nothing that the superbasin recycler can do at this point,
        # use the normal recycling process.
        if (config.recycling_on and not sb_recycler.in_progress):
            recycler = recycling.Recycling(states, previous_state, current_state, 
                            config.recycling_move_distance, config.recycling_save_sugg)
    
    if current_state.get_confidence() < config.akmc_confidence:
        if config.recycling_on and recycler:
            pass_recycler = recycler
        else:
            pass_recycler = None
        if config.kdb_on:
            pass_kdb = kdber          	
        else:
            pass_kdb = None
        if config.sb_recycling_on:
            pass_sb_recycler = sb_recycler
        else:
            pass_sb_recycler = None
        wuid = make_searches(comm, current_state, wuid, searchdata = searchdata, recycler = pass_recycler, kdber = pass_kdb, sb_recycler = pass_sb_recycler)
    
    # Write out metadata. XXX:ugly
    metafile = os.path.join(config.path_results, 'info.txt')
    parser = ConfigParser.RawConfigParser() 

    for key in searchdata.keys():
        #XXX: This may be buggy once superbasin recycling is implemented
        if int(key.split('_')[0]) < current_state.number:
            del searchdata[key]

    if previous_state.number != current_state.number:
        previous_state_num = previous_state.number

    write_akmc_metadata(parser, current_state.number, time, wuid, searchdata, previous_state_num)

    parser.write(open(metafile, 'w')) 



def get_akmc_metadata():
    if not os.path.isdir(config.path_results):
        os.makedirs(config.path_results)
    # read in metadata
    # do we want custom metadata locations?
    metafile = os.path.join(config.path_results, 'info.txt')
    parser = ConfigParser.SafeConfigParser() 
    if os.path.isfile(metafile):
        parser.read(metafile)
        try:
            start_state_num = parser.getint("Simulation Information",'current_state')
        except:
            start_state_num = 0 #Sadly, ConfigParser doesn't have a better way of specifying defaults
        try:
            time = parser.getfloat("Simulation Information", 'time_simulated') 
        except:
            time = 0.0
        try:
            wuid = parser.getint("aKMC Metadata", 'wu_id') 
        except:
            wuid = 0
        try:
            searchdata = eval(parser.get("aKMC Metadata", 'searchdata'))
        except:
            searchdata={}
        try:
            previous_state_num = parser.getint("Simulation Information", "previous_state")
        except:
            previous_state_num = -1
        try:
            first_run = parser.getboolean("Simulation Information", "first_run")
        except:
        
            first_run = True
    else:
        time = 0
        start_state_num = 0
        wuid = 0
        searchdata = {}
        previous_state_num = -1
        first_run = True

    if config.debug_random_seed:
        try:
            parser = ConfigParser.RawConfigParser()
            parser.read(metafile)
            seed = parser.get("aKMC Metadata", "random_state")
            from numpy import array, uint32
            numpy.random.set_state(eval(seed))
            logger.debug("Set random state from previous run's state")
        except:
            numpy.random.seed(config.debug_random_seed)
            logger.debug("Set random state from seed")

    return start_state_num, time, wuid, searchdata, previous_state_num, first_run



def write_akmc_metadata(parser, current_state_num, time, wuid, searchdata, previous_state_num):
    parser.add_section('aKMC Metadata')
    parser.add_section('Simulation Information')
    parser.set('aKMC Metadata', 'wu_id', str(wuid))
    parser.set('aKMC Metadata', 'searchdata', repr(searchdata))
    parser.set('Simulation Information', 'time_simulated', str(time))
    parser.set('Simulation Information', 'current_state', str(current_state_num))
    parser.set('Simulation Information', 'previous_state', str(previous_state_num))
    parser.set('Simulation Information', 'first_run', str(False))
    if config.debug_random_seed:
        parser.set('aKMC Metadata', 'random_state', repr(numpy.random.get_state()))



def get_statelist(kT):
    initial_state_path = os.path.join(config.path_root, 'reactant.con') 
    return statelist.StateList(config.path_states, 
                               kT, 
                               config.akmc_thermal_window, 
                               config.akmc_max_thermal_window, 
                               config.comp_eps_e, 
                               config.comp_eps_r, 
                               config.comp_use_identical, 
                               initial_state_path, 
                               list_search_results = config.debug_list_search_results, 
                               filter_hole = config.disp_moved_only)  



def get_communicator():
    if config.comm_type=='boinc':
        comm = communicator.BOINC(config.path_scratch, config.comm_boinc_project_dir, 
                config.comm_boinc_wu_template_path, config.comm_boinc_re_template_path,
                config.comm_boinc_appname, config.comm_boinc_results_path,
                config.comm_job_bundle_size)
    elif config.comm_type=='cluster':
        comm = communicator.Script(config.path_scratch, config.comm_job_bundle_size,
                                   comm_script_queued_jobs_cmd, comm_script_cancel_job_cmd, 
                                   comm_script_submit_job_cmd)
    elif config.comm_type=='local':
        comm = communicator.Local(config.path_scratch, config.comm_local_client, 
                                  config.comm_local_ncpus, config.comm_job_bundle_size)
    elif config.comm_type=='mpi':
        comm = communicator.MPI(config.path_scratch, config.comm_mpi_client, 
                                  config.comm_job_bundle_size, config.comm_mpi_mpicommand)
    elif config.comm_type=='arc':
        comm = communicator.ARC(config.path_scratch, config.comm_job_bundle_size, 
                                config.comm_client_path, config.comm_blacklist)
    else:
        logger.error(str(config.comm_type)+" is an unknown communicator.")
        raise ValueError()
    return comm



def register_results(comm, current_state, states, searchdata = None):
    logger.info("registering results")
    t1 = unix_time.time()
    if os.path.isdir(config.path_searches_in):
        shutil.rmtree(config.path_searches_in)    
    os.makedirs(config.path_searches_in)
    
    #Function used by communicator to determine whether to discard a result
    def keep_result(name):
        state_num = int(name.split("_")[0])
        return (config.debug_register_extra_results or \
                state_num == current_state.number or \
                states.get_state(state_num).get_confidence() < config.akmc_confidence)

    num_registered = 0
    for result in comm.get_results(config.path_searches_in, keep_result): 
        # The result dictionary contains the following key-value pairs:
        # reactant - an array of strings containing the reactant
        # saddle - an atoms object containing the saddle
        # product - an array of strings containing the product
        # mode - an array of strings conatining the mode
        # results - a dictionary containing the key-value pairs in results.dat
        # id - StateNumber_WUID
        #
        # The reactant, product, and mode are passed as lines of the files because
        # the information contained in them is not needed for registering results
        if config.debug_keep_all_results:
            #XXX: We should only do these checks once to speed things up, but at the same time
            #debug options don't have to be fast
            #save_path = os.path.join(config.path_root, "old_searches")
            #if not os.path.isdir(save_path):
            #    os.mkdir(save_path)
            #shutil.copytree(result_path, os.path.join(save_path, i))
            #XXX: This is currently broken by the new result passing scheme. Should it be 
            #     done in communicator?
            pass
        state_num = int(result['name'].split("_")[0])
        id = int(result['name'].split("_")[1]) + result['number']
        searchdata_id = "%d_%d" % (state_num, id)
        # Store information about the search into result_data for the search_results.txt file in the state directory.
        if config.debug_list_search_results:
            try:
                result['type'] = searchdata[searchdata_id + "type"]
                del searchdata[searchdata_id + "type"]
            except:
                logger.warning("Could not find search data for search %s" % searchdata_id)
            result['wuid'] = id
            # Remove used information from the searchdata metadata.
        
        #read in the results
        result['results'] = io.parse_results_dat(result['results.dat'])
        if result['results']['termination_reason'] == 0:
            process_id = states.get_state(state_num).add_process(result)
        else:
            states.get_state(state_num).register_bad_saddle(result, config.debug_keep_bad_saddles)
        num_registered += 1
        
        if current_state.get_confidence() >= config.akmc_confidence:
            if not config.debug_register_extra_results:
                break
    
    #Approximate number of searches recieved
    tot_searches = len(os.listdir(config.path_searches_in)) * config.comm_job_bundle_size
    
    t2 = unix_time.time()
    logger.info("%i (result) searches processed", num_registered)
    logger.info("Approximately %i (result) searches discarded." % (tot_searches - num_registered))
    #logger.info("%i results discarded", len(results) - num_registered + discarded * config.comm_job_bundle_size)
    if num_registered == 0:
        logger.debug("0 results per second", num_registered)
    else:
        logger.debug("%.1f results per second", (num_registered/(t2-t1)))
        
    return num_registered



def get_superbasin_scheme(states):
    if config.sb_scheme == 'transition_counting':
        superbasining = superbasinscheme.TransitionCounting(config.sb_path, states, config.akmc_temperature / 11604.5, config.sb_tc_ntrans)
    elif config.sb_scheme == 'energy_level':
        superbasining = superbasinscheme.EnergyLevel(config.sb_path, states, config.akmc_temperature / 11604.5, config.sb_el_energy_increment)
    return superbasining



def kmc_step(current_state, states, time, kT, superbasining, previous_state_num = None):
    t1 = unix_time.time()
    previous_state = current_state 
    start_state_num = current_state.number
    steps = 0
    # If the Chatterjee & Voter superbasin acceleration method is being used
    if config.askmc_on:
        if config.sb_recycling_on:
            pass_rec_path = config.sb_recycling_path
        else:
            pass_rec_path = None
        asKMC = askmc.ASKMC(kT, states, config.askmc_confidence, config.askmc_alpha,
                            config.askmc_gamma, config.askmc_barrier_test_on,
                            config.askmc_connections_test_on, config.sb_recycling_on,
                            config.path_root, config.akmc_thermal_window,
                            recycle_path = pass_rec_path)
    while current_state.get_confidence() >= config.akmc_confidence and steps < config.akmc_max_kmc_steps:
        steps += 1
        if config.sb_on:
            sb = superbasining.get_containing_superbasin(current_state)

        if config.sb_on and sb:
            mean_time, next_state = sb.step(current_state, states.get_product_state)
        else:
            if config.askmc_on:
                rate_table = asKMC.get_ratetable(current_state)
            else:
                rate_table = current_state.get_ratetable()
            if len(rate_table) == 0:
                logger.error("No processes in rate table, but confidence has been reached")
            ratesum = 0.0
            for i in range(len(rate_table)):
                ratesum += rate_table[i][1]
            
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
                    print "Can no longer follow target trajectory."
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
                    print "Can no longer follow target trajectory."
                    sys.exit(1)

            # We are not following another trajectory:
            else:
                for i in range(len(rate_table)):
                    p += rate_table[i][1]/ratesum
                    if p>u:
                        nsid = i
                        break
                else:
                    logger.warning("Warning: failed to select rate. p = " + str(p))
                    break
            next_state = states.get_product_state(current_state.number, rate_table[nsid][0])
            mean_time = 1.0/ratesum
        if config.debug_use_mean_time:
            time += mean_time
        else:
            time -= mean_time*math.log(1 - numpy.random.random_sample())# numpy.random.random_sample() uses [0,1), which could produce issues with math.log()
        if config.askmc_on:
            asKMC.register_transition(current_state, next_state)
        if config.sb_on:
            superbasining.register_transition(current_state, next_state)    
        
        if config.sb_on and sb:
            proc_id_out = -1
        else:
            proc_id_out = rate_table[nsid][0]
        dynamics = io.Dynamics(os.path.join(config.path_results, "dynamics.txt"))
        if proc_id_out != -1:            
            proc = current_state.get_process(proc_id_out)
            dynamics.append(current_state.number, proc_id_out, next_state.number, mean_time, time, proc['barrier'], proc['rate'])
        else:
            #XXX The proc_out_id was -1, which means there's a bug or this was a superbasin step.
            dynamics.append(current_state.number, proc_id_out, next_state.number, mean_time, time, 0, 0)
        logger.info("stepped from state %i to state %i", current_state.number, next_state.number)
        
        previous_state = current_state
        current_state = next_state

    if config.sb_on:
        superbasining.write_data()

    logger.info("currently in state %i with confidence %.6f", current_state.number, current_state.get_confidence())
    t2 = unix_time.time()
    logger.debug("KMC finished in " + str(t2-t1) + " seconds")
    return current_state, previous_state, time



def get_displacement(reactant, indices=None):
    if config.disp_type == 'random':
        disp = displace.Random(reactant, config.disp_size, config.disp_radius, hole_epicenters=indices)
    elif config.disp_type == 'undercoordinated':
        disp = displace.Undercoordinated(reactant, config.disp_max_coord, config.disp_size, config.disp_radius, hole_epicenters=indices)
    elif config.disp_type == 'leastcoordinated':
        disp = displace.Leastcoordinated(reactant, config.disp_size, config.disp_radius, hole_epicenters=indices)
    elif config.disp_type == 'water':
        disp = displace.Water(reactant, config.stdev_translation, config.stdev_rotation)
    else:
        raise ValueError()
    return disp



def make_searches(comm, current_state, wuid, searchdata = None, kdber = None, recycler = None, sb_recycler = None):
    
    
    reactant = current_state.get_reactant()
    num_in_buffer = comm.get_queue_size()*config.comm_job_bundle_size #XXX:what if the user changes the bundle size?
    logger.info("%i searches in the queue" % num_in_buffer)
    num_to_make = max(config.comm_search_buffer_size - num_in_buffer, 0)
    logger.info("making %i searches" % num_to_make)
    
    if num_to_make == 0:
        return wuid
    # If we plan to only displace atoms that moved getting to the current state.
    if config.disp_moved_only and current_state.number != 0:
        pass_indices = recycler.get_moved_indices()
    else:
        pass_indices = None
    disp = get_displacement(reactant, indices = pass_indices)
    parameters_path = os.path.join(config.path_root, "parameters.dat")
    searches = []
    
    invariants = {}

    reactIO = StringIO.StringIO()
    io.savecon(reactIO, reactant)
    invariants['reactant_passed.con']=reactIO
    
    f = open(parameters_path)
    invariants['parameters_passed.dat'] = StringIO.StringIO(''.join(f.readlines()))
    f.close()

    invariants = dict(invariants,  **io.load_potfiles(config.path_pot)) #Merge potential files into invariants

    t1 = unix_time.time()
    if config.recycling_on:
        nrecycled = 0
    for i in range(num_to_make):
        search = {}
        # The search dictionary contains the following key-value pairs:
        # id - CurrentState_WUID
        # displacement - an atoms object containing the point the saddle search will start at
        # mode - an Nx3 numpy array containing the initial mode 
        search['id'] = "%d_%d" % (current_state.number, wuid)
        done = False
        # Do we want to try superbasin recycling? If yes, try. If we fail to recycle the basin,
        # move to the next case
        if (config.sb_recycling_on and current_state.number is not 0):
            displacement, mode = sb_recycler.make_suggestion()
            if displacement:
                nrecycled += 1
                if config.debug_list_search_results:
                    try:
                        searchdata["%d_%d" %(current_state.number, wuid) + "type"] = "recycling"
                    except:
                        logger.warning("Failed to add searchdata for search %d_%d" % (current_state.number, wuid))
                done = True
        # Do we want to do recycling? If yes, try. If we fail to recycle, we move to the next case
        if (recycler and current_state.number is not 0):
            displacement, mode = recycler.make_suggestion()
            if displacement:
                nrecycled += 1
                if config.debug_list_search_results:                
                    try:
                        searchdata["%d_%d" % (current_state.number, wuid) + "type"] = "recycling"
                    except:
                        logger.warning("Failed to add searchdata for search %d_%d" % (current_state.number, wuid))
                done = True
        if not done and config.kdb_on:
            # Set up the path for keeping the suggestion if config.kdb_keep is set.
            keep_path = None
            if config.kdb_keep:
                if not os.path.isdir(os.path.join(current_state.path, "kdbsuggestions")):
                    os.mkdir(os.path.join(current_state.path, "kdbsuggestions"))
                keep_path = os.path.join(current_state.path, "kdbsuggestions", str(wuid))
            displacement, mode = kdber.make_suggestion(keep_path)
            if displacement:
                done = True
                logger.info('Made a KDB suggestion')
                if config.debug_list_search_results:                
                    try:
                        searchdata["%d_%d" % (current_state.number, wuid) + "type"] = "kdb"
                    except:
                        logger.warning("Failed to add searchdata for search %d_%d" % (current_state.number, wuid))
        if not done:
            displacement, mode = disp.make_displacement() 
            if config.debug_list_search_results:                
                try:
                    searchdata["%d_%d" % (current_state.number, wuid) + "type"] = "random"
                except:
                    logger.warning("Failed to add searchdata for search %d_%d" % (current_state.number, wuid))
        dispIO = StringIO.StringIO()
        io.savecon(dispIO, displacement)
        search['displacement_passed.con'] = dispIO
        
        modeIO = StringIO.StringIO()
        io.save_mode(modeIO, mode, reactant)
        search['mode_passed.dat'] = modeIO
        searches.append(search) 
        wuid += 1

    if config.recycling_on and nrecycled > 0:
        logger.debug("Recycled %i saddles" % nrecycled)

    try:
        comm.submit_jobs(searches, invariants)
        t2 = unix_time.time()
        logger.info( str(num_to_make) + " searches created") 
        logger.debug( str(num_to_make/(t2-t1)) + " searches per second")
    except:
        logger.exception("Failed to submit searches.")
    return wuid


#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

if __name__ == '__main__':
    
    optpar = optparse.OptionParser(usage = "usage: %prog [options] config.ini")
    optpar.add_option("-R", "--reset", action="store_true", dest="reset", default = False, help="reset the aKMC simulation, discarding all data")
    optpar.add_option("-s", "--status", action="store_true", dest="print_status", default = False, help = "print the status of the simulation and currently running jobs")
    optpar.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False,help="only write to the log file")
    optpar.add_option("-m", "--movie", action="store", dest="movie_type", default = "", help="Specify the type of movie to make [dynamics, states, fastestpath, fastestfullpath, graph, processes]. Process movies are specified like so: --movie processes,statenumber,processlimit. Where processes is the string processes, statenumber is the number of the state that you want to view, and process limit is the maximum number of processes you would like in the movie. The returned processes are reverse sorted by rate such that the fastest processes is the first in the movie.")
    optpar.add_option("-n", "--no-submit", action="store_true", dest="no_submit", default=False,help="don't submit searches; only register finished results")
    (options, args) = optpar.parse_args()

    if len(args) > 1:
        print "akmc.py takes only one positional argument"
    sys.argv = sys.argv[0:1]
    if len(args) == 1:
        sys.argv += args
        #always run from the directory where the config file is
        #os.chdir(os.path.dirname(args[0]))

    #XXX: config is ugly as it finds out where the config file is directly from 
    #     sys.argv instead of being passed it.
    import config
    #set options.path_root to be where the config file is if given as an arg
    if config.path_root.strip() == '.' and len(args) == 1:
        config.path_root = os.path.abspath(os.path.dirname(args[0]))
        os.chdir(config.path_root)
    if options.no_submit:
        config.comm_search_buffer_size = 0

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

    # Some options are mutually exclusive. Let's check them now.
    exclusive_options = {}
    if len(options.movie_type) > 0:
        exclusive_options['movie_type'] = True
    else:
        exclusive_options['movie_type'] = False
    exclusive_options['print_status'] = options.print_status
    exclusive_options['reset'] = options.reset

    if sum(exclusive_options.values()) > 1:
        offending_options = [ k for k,v in exclusive_options.iteritems() if v ]
        optpar.error("the options %s are mutually exclusive" % ", ".join(offending_options))

    if len(options.movie_type) > 0:
        states = get_statelist(config.akmc_temperature / 11604.5)
        movie.make_movie(options.movie_type, config.path_root, states)
        sys.exit(0)

    # From the config file: The Novotny and C&V (ASKMC) methods should not be used together.
    if config.sb_on and config.askmc_on:
        logger.error("Both superbasin methods should not be used at the same time.")
        sys.exit(1)

    if options.print_status:
        states = get_statelist(config.akmc_temperature / 11604.5)
        start_state_num, time, wuid, searchdata, previous_state_num = get_akmc_metadata()

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
        print "Confidence: %.4f" % current_state.get_confidence()
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

        if config.sb_on: 
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
                    rmdirs.append(config.kdb_path) 
                if config.sb_on:
                    rmdirs.append(config.sb_path)
                if config.sb_recycling_on:
                    rmdirs.append(config.sb_recycling_path)
                if config.debug_keep_all_results:
                    rmdirs.append(os.path.join(config.path_root, "old_searches"))
                for i in rmdirs:
                    if os.path.isdir(i):
                        shutil.rmtree(i)
                        #XXX: ugly way to remove all empty directories containing this one
                        os.mkdir(i)
                        os.removedirs(i)
                
                askmc_data_path = os.path.join(config.path_results, "askmc_data.txt")
                if os.path.isfile(askmc_data_path):
                    os.remove(askmc_data_path)
                dynamics_path = os.path.join(config.path_results, "dynamics.txt")  
                info_path = os.path.join(config.path_results, "info.txt") 
                log_path = os.path.join(config.path_results, "akmc.log") 
                for i in [info_path, dynamics_path, log_path]:
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
