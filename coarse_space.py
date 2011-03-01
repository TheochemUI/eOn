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
logger = logging.getLogger('akmc')
import numpy
numpy.seterr(all='raise')
import cPickle as pickle

import config
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

import tables
import parsers
import state

import StringIO


class CircularlyLinkedLists(list):
    """Record groups of indices by linked lists."""
    """Indices [0, 1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 11]"""
    """Data    [3, 2, 1, 4, 0, 5, 6, 11, 9, 7, 10, 8 ]"""
    """Loops in the list [0, 3, 4]; [1, 2]; [5]; [6]; [7, 11, 8, 9]; [10]"""

#    def __init__(self, listed_data=[]):
#        list.__init__(self, listed_data)
#    def __init__(self, listed_data=[]):
#        list.__init__(self, listed_data)

    def add_index(self):
        """Add a new index. The new index is equal to the length of the list (number of indices) before the addition of the new index, """\
        """so that ordered list shall always looks like: [0, 1, 2, 3, ...]. Returns the new index."""
        n = len(self)
        self.append(n)
        return n
    
    def get_linked_list(self, index):
        """List of indices belonging to the same linked list as 'index'."""
#        print index
        assert(index >= 0)
        linked_list = list()
        i = self[index]    
        linked_list.append(i)

        while i != index:
            i = self[i]
            linked_list.append(i)
        return linked_list

    def merge(self, i1, i2):
        """Merge the two linked lists (groups) to which 'i1' and 'i2' belong into a single linked list."""

        i0 = self[i1]
        self[i1] = self[i2]
        self[i2] = i0


class Cluster(list):
    """Cluster of states and coarse states."""
    #List of indices for the CoarseStates belonging to the cluster
    
    def __init__(self, linked_list, structure_cg_states):
        """'linked_list': indices for all the states making the cluster."""\
        """'cg_states': instance of CircularlyLinkedLists decribing the structure of coarse granined states."""
        list.__init__(self, linked_list)
        self._structure_cg_states = structure_cg_states
        self._active = None

    def get_contained_cg_states_nrs(self):
        """Returns a list of the nrs for the contained coarse grained states."""
        cg_states = self.get_contained_states_nrs_ordered()
        cg_state_nrs = []
        for cg_state in cg_states:
            cg_state_nrs.append(min(cg_state))
        return cg_state_nrs

    def get_contained_states_nrs_ordered(self):
        """Returns a list of lists of indices, which are the coarse grained states inside the cluster."""
        cg_states = list()
        states_in_cluster = list(self)
        
#        print states_in_cluster
        
        while len(states_in_cluster) != 0:
            j = states_in_cluster[0]
            cg_state = self._structure_cg_states.get_linked_list(j)
            cg_states.append(cg_state)
            for k in cg_state:
 #               print k
                states_in_cluster.remove(k)
        return cg_states

class SampledSpace():
    """All clusters in the sampled space."""
    #Circularly linked lists of indices for the clusters belonging to the SampledSpace
    #Indices for all states and their grouping into coarse states and clusters sampled are contained

    def __init__(self, state_list):
        self._path = config.path_coarse_states
        self._file_path = os.path.join(self._path, config.file_coarse)
        self._state_list = state_list
        self._changed = False

        if not os.path.isdir(self._path):
            os.mkdir(self._path)

        # Creates new info file for the state
        
        if not os.path.isfile(self._file_path):
            self._info = parsers.CoarseInfo(self._path)
            self._info.new()
            
            self._cg_states = CircularlyLinkedLists()
            self._clusters = CircularlyLinkedLists()
            
            self.add_state()

#            nr_states = len(state_list)
#
#            if nr_states:
#                # Add all known states
#                for state_nr in range(nr_states):
#                    self.add_state()
#
#                # Analyse the given list of states
#                for state_nr in range(nr_states):
#                    while state_nr < nr_states:
#                        next_state_nr = self.extend(state_nr, extend=0)
#                        if  state_nr == next_state_nr or next_state_nr == -1:
#                            break
#                        else:
#                            state_nr = next_state_nr
 #           else:
#                self.add_state()

        # Reads existing info file
        else:
            self._info = parsers.CoarseInfo(self._path)
            self._clusters = CircularlyLinkedLists(self._info.get_clusters())
            self._cg_states = CircularlyLinkedLists(self._info.get_cg_states())

#        print "building done"
#        print "--------------------------------"
        return

    def __del__(self):
        if self._changed:
            self._info.set_clusters(self._clusters)
            self._info.set_cg_states(self._cg_states)
            
    def __len__(self):
 #       print self._cg_states
        return len(self._cg_states)
    
    def add_state(self):
        """Add a new state, a new cluster and a new coarse grained state. Returns the index attributed to the new state."""
        nr = len(self)

        print "########### add_state ############"

        self._clusters.append(nr)
        self._cg_states.append(nr)
        cg_state_info = parsers.CoarseStateListInfo(config.path_coarse_states, nr)
        cg_state_info.new()
        cg_state_info.set_state_nrs([nr])
        del cg_state_info
        self._changed = True
        return nr

    def merge_clusters(self, i1, i2):        
        """Merge two clusters into one. 'i1' and 'i2' must be the indices of two states in two different clusters."""
        print "merge clusters"
#        print self._clusters
#        print "i1",i1
#        print "i2",i2
        print "cluster",self._clusters
        self._clusters.merge(i1, i2)
        print "cluster",self._clusters
        print "after merge clusters"

        self._changed = True
    
    def merge_cg_states(self, i1, i2):
        """Merge two coarse grained states into one. """\
        """'i1' and 'i2' must be the indices of two states in two different coase grained states."""

        print "merge cg"
        print "cg_states",self._cg_states
        index_cg_1 = min(self._cg_states.get_linked_list(i1))
        index_cg_2 = min(self._cg_states.get_linked_list(i2))

        self._cg_states.merge(i1, i2)
        indices = self._cg_states.get_linked_list(i1)

        print "cg_states",self._cg_states
        print "after merge cg"

        min_index = min(index_cg_1, index_cg_2)
        max_index = max(index_cg_1, index_cg_2)

        # The file for cg state containing the lowest state nr is expanded
        cg_state_info = parsers.CoarseStateListInfo(config.path_coarse_states, min_index)
#        cg_state_info.new()
        cg_state_info.set_state_nrs(indices)
        del cg_state_info
        # The other is removed
        cg_remove = os.path.join(config.path_coarse_states, str(max_index))
        cg_storage = os.path.join(config.path_coarse_states, "storage")
        shutil.move(cg_remove, cg_storage)
        self._changed = True

    def get_cluster(self, index):
        """Returns the cluster to which as state 'index' belongs."""
#        print self._clusters
#        print index
        linked_list = self._clusters.get_linked_list(index)
#        print linked_list
        cluster = Cluster(linked_list, self._cg_states)
        return cluster
        
    def get_cg_state_nrs(self, index):
        """Returns the indices for the states within the coarse grained state to which 'index' belongs."""
        return self._cg_states.get_linked_list(index)
    
    def get_cg_state_nr(self, index):
        """Returns the nr for the coarse state to which 'index' belongs."""
        return min(self._cg_states.get_linked_list(index))
        
    def get_exit_process(self):
        return self._exit_state_nr, self._exit_process_nr
        
        
    def get_lowest_escape_barrier_in_cluster(self, cluster):
        lowest_barrier = float('inf')
        lowest_barrier_info = None

        # find the lowest barrier in cluster
        for cluster_cg_state_nr in cluster.get_contained_cg_states_nrs():
            cluster_cg_state_nrs = self.get_cg_state_nrs(cluster_cg_state_nr)
            cluster_cg_state_list = statelist.CoarseStateList(cluster_cg_state_nr, cluster_cg_state_nrs)
            barrier_info = cluster_cg_state_list.get_lowest_escape_barrier()
            
            if barrier_info['barrier'] < lowest_barrier:
                lowest_barrier_info = barrier_info

        return lowest_barrier_info
        
    def extend(self, state_nr):#, extend=1):
        """Extend coarse grained space. \nReturns True if closed space, False if open."""

#        temperature = 1.
        time_cutoff = 10.
#        print "state_nr",state_nr
        cluster = self.get_cluster(state_nr)
        cg_state_nr = self.get_cg_state_nr(state_nr)
        cg_state_nrs = self.get_cg_state_nrs(state_nr)
        
#        print cg_state_nr
#        print cg_state_nrs
#        raw_input()
        
        cg_state_list = statelist.CoarseStateList(cg_state_nr, cg_state_nrs)
        escape_rate = cg_state_list.get_escape_rate()
        print "escape_rate",escape_rate

        raw_input()

        if escape_rate == 0.0:
            # Closed system, all done! 
            stop
        elif 1. / escape_rate > time_cutoff:
            # If mean residence time more than cutoff, create a new state in a new cluster
            lowest_barrier_info = cg_state_list.get_lowest_escape_barrier()
            
            self._exit_state_nr = lowest_barrier_info['state_nr'] 
            self._exit_process_nr = lowest_barrier_info['process_nr']
            new_state_nr = self.add_state()

        # Leaved the cluster
        else:
            lowest_barrier_info = self.get_lowest_escape_barrier_in_cluster(cluster)
            cg_state_nrs_in_cluster = cluster.get_contained_cg_states_nrs()
 
            # Get the product state for this process
            exit_state_nr = lowest_barrier_info['state_nr'] 
            exit_process_nr = lowest_barrier_info['process_nr']
            exit_state = self._state_list.get_state(exit_state_nr)
            
            product_energy = exit_state.get_process_product_energy(exit_process_nr)
            product_structure = exit_state.get_process_product(exit_process_nr)
            
            # Need to check if the product state is known
            new_state_nr = self._state_list.if_known_return_state_nr(product_energy, product_structure)

            # Arrived in a state in the same cluster: merge cg states 
            if new_state_nr in cluster:
                self.merge_cg_states(exit_state_nr, new_state_nr)
                print "known in cluster, merge cg"
#                new_state_nr = exit_state_nr
#                print "exit_state_nr",exit_state_nr
#                print "new_state_nr",new_state_nr
#                print "state_nr",state_nr

            # Arrived in a new state: add state to current cluster 
            elif new_state_nr == -1:
                print "new in cluster"
#                if extend:
                print "exit_state_nr",exit_state_nr
                print "exit_process_nr",exit_process_nr

                self._exit_state_nr = exit_state_nr 
                self._exit_process_nr = exit_process_nr

                new_state_nr = self.add_state()
                self.merge_clusters(exit_state_nr, new_state_nr)
            

            # Arrived in another cluster: merge old clusters and create a new state in a new cluster 
            else:
                "new cluster, merge old clusters"
                # Determine the lowest barrier out of the two clusters
                cluster = self.get_cluster(exit_state_nr)
                lowest_barrier_info_1 = self.get_lowest_escape_barrier_in_cluster(cluster)
                cluster = self.get_cluster(new_state_nr)
                lowest_barrier_info_2 = self.get_lowest_escape_barrier_in_cluster(cluster)

                if lowest_barrier_info_1['barrier'] < lowest_barrier_info_2['barrier']:
                    self._exit_state_nr = lowest_barrier_info_1['state_nr'] 
                    self._exit_process_nr = lowest_barrier_info_1['process_nr']
                else:
                    self._exit_state_nr = lowest_barrier_info_2['state_nr'] 
                    self._exit_process_nr = lowest_barrier_info_2['process_nr']

                self.merge_clusters(exit_state_nr, new_state_nr)
                self.add_state()
                
        self._changed = True
        return new_state_nr


def test(config):

    akmc_info = parsers.AKMCInfo(config.path_results)
    comm = communicator.get_communicator()
        
    # Loads random seed if requested
    if config.main_random_seed:
        seed = akmc_info.get_random_seed()
        if seed is None:
            numpy.random.seed(config.main_random_seed)
            logger.debug("Set random state from seed")
        else:
            numpy.random.set_state(seed)
            logger.debug("Set random state from previous run's state")

    # Create the statelist directory if it does not exist.
    if not os.path.isdir(config.path_coarse_states):
        logger.warning("coarse state path does not exist, creating %s" % config.path_coarse_states)
        storage = os.path.join(config.path_coarse_states, "storage")
        os.makedirs(config.path_coarse_states)
        os.makedirs(storage)

    # Fetch the data for the current state
    state_list = get_statelist()
    current_state_nr = akmc_info.get_current_state_nr()
    current_state = state_list.get_state(current_state_nr)        

    print "confidence:",current_state.get_confidence()
    print "current_state_nr",current_state_nr
    
    if current_state.table_complete():
#        raw_input()
        
#        print "hello"
#        print len(state_list)
        sampled = SampledSpace(state_list)
#        state_nr = current_state_nr
#        print "hello"
        state_nr = sampled.extend(current_state_nr)
        
#        print "state_nr",state_nr
        
        
#        while state_nr < len(sampled)-1 and state_nr != -1:
        while state_nr < len(sampled)-1 and state_nr != -1:
            next_state_nr = sampled.extend(state_nr) #, extend=0)
            print "next_state_nr",next_state_nr
 
#            if  state_nr == next_state_nr or next_state_nr == -1:
            if next_state_nr == -1:
                print "state_nr",state_nr
                print "next_state_nr",next_state_nr
                break
            else:
                state_nr = next_state_nr

#        sampled.add_state()
#        next_state_nr = len(sampled)
        
#        print next_state_nr
        exit_state_nr, exit_process_nr = sampled.get_exit_process()
        
        print "exit_state_nr",exit_state_nr
        print "exit_process_nr",exit_process_nr
        
        state_list.register_choosen_process(exit_state_nr, exit_process_nr)        
        akmc_info.set_current_state_nr(len(sampled)-1)
#        print len(sampled)
    else:

        try:
            sd = open(os.path.join(config.path_root, "searchdata"), "r")
            searchdata = pickle.load(sd)
            sd.close()
        except:
            searchdata={}        

        register_searches(comm, searchdata, akmc_info)
        make_searches(comm, current_state, akmc_info, searchdata = searchdata)
        
        sd = open(os.path.join(config.path_root, "searchdata"), "w")
        pickle.dump(searchdata, sd)
        sd.close()



def register_searches(comm, searchdata, akmc_info):
    logger.info("registering results")
    t1 = unix_time.time()

    # Create directories if nessesary
    if os.path.isdir(config.path_searches_in):
        if config.debug_keep_all_results:
            if not os.path.isdir(os.path.join(config.path_root, "results")):
                os.makedirs(os.path.join(config.path_root, "results"))
            for d in os.listdir(config.path_searches_in):
                shutil.move(os.path.join(config.path_searches_in, d), os.path.join(config.path_root, "results", d))
        shutil.rmtree(config.path_searches_in)  
    os.makedirs(config.path_searches_in)
    
    # Function used by communicator to determine whether to discard a result
    def keep_result(name):
        state_num = int(name.split("_")[0])
        return (config.debug_register_extra_results or \
                state_num == akmc_info.get_current_state_nr())# or 
#                not states.get_state(state_num).table_complete())

    # The new processes are registered
    nr_good = 0
    nr_bad = 0
    states_data = {}
    
    for result in comm.get_results(config.path_searches_in, keep_result):
        comment = None
        state_num = int(result['name'].split("_")[0])
        id = int(result['name'].split("_")[1]) + result['number']
        searchdata_id = "%d_%d" % (state_num, id)

        # Get data for the state corresponding to the reactant of the process
        try:
            state_data = states_data[state_num]
        except KeyError:
            state_data = {}
            state_path = os.path.join(config.path_states, str(state_num))
            state_data['info'] = parsers.StateInfo(state_path, state_num)
            state_data['process_table'] = tables.ProcessTableAKMC(state_path, state_num)
            states_data[state_num] = state_data

        # Store debug information about the search into result_data 
        # for the search_results.txt file in the state directory.
        if config.debug_list_search_results:

            try:
                result['type'] = searchdata[searchdata_id]['type']
#                del searchdata[searchdata_id]["type"]
                del searchdata[searchdata_id]
            except:
                result['type'] = "unknown"
                logger.warning("Could not find search data for search %s" % searchdata_id)
            result['wuid'] = id
        
        # Add process to table and info file for the specific state 
        result['results'] = io.parse_results(result['results.dat'])

        if result['results']['termination_reason'] == 0:
                
            # Checks if the energy of the state is known, if so use that value else set it
            if state_data['info'].get_reactant_energy() is not None:
                result['results']['potential_energy_reactant'] = state_data['info'].get_reactant_energy()
            else:
                state_data["info"].set_reactant_energy(result['results']['potential_energy_reactant'])

            process_id, comment, new_process, new_barrier = state_data['process_table'].register_good_process(result)
            state_data['info'].increase_proc_repeat_count(process_id)
            state_data['info'].increase_good_saddle_count()
            if new_process:
                state_data['info'].increase_unique_saddle_count()
                if new_barrier:
                    lowest_barrier = state_data['process_table'].get_lowest_barrier()
                    state_data['info'].set_lowest_barrier(lowest_barrier)
            nr_good += 1
        else:
            state_data['info'].increase_bad_saddle_count()
            if config.debug_keep_bad_saddles:
                path_bad_saddles = os.path.join(config.path_root, 'badprocesses')
            else:
                path_bad_saddles = None
            comment = tables.register_bad_process(result, path_bad_saddles)
            nr_bad += 1
        tables.append_search_result(result, state_num, comment)
            
    # Cleans the directory containing data about new processes 
    if os.path.isdir(config.path_searches_in):
        # Storing the results files before the clean
        if config.debug_keep_all_results:
            if not os.path.isdir(os.path.join(config.path_root, "results")):
                os.makedirs(os.path.join(config.path_root, "results"))
            for d in os.listdir(config.path_searches_in):
                shutil.move(os.path.join(config.path_searches_in, d), os.path.join(config.path_root, "results", d))
        shutil.rmtree(config.path_searches_in)
    os.makedirs(config.path_searches_in)

    # Logging the data
    t2 = unix_time.time()
    nr_searches = nr_bad + nr_good
    logger.info("%i (result) searches processed.", nr_good)
    logger.info("%i (result) searches discarded.", nr_bad)
    if nr_good == 0:
        logger.debug("0 results per second")
    else:
        logger.debug("%.1f results per second", (nr_good/(t2-t1)))
        
    return nr_good     



def make_searches(comm, current_state, akmc_info, searchdata = None, kdber = None, recycler = None, sb_recycler = None):
    reactant = current_state.get_reactant()
    #XXX:what if the user changes the bundle size?
    num_in_buffer = comm.get_queue_size()*config.comm_job_bundle_size 
    logger.info("%i searches in the queue" % num_in_buffer)
    num_to_make = max(config.comm_search_buffer_size - num_in_buffer, 0)
    logger.info("making %i searches" % num_to_make)
    
    if num_to_make == 0:
        return akmc_info.get_wuid()
    # If we plan to only displace atoms that moved getting to the current state.
    if config.disp_moved_only and current_state.get_nr() != 0:
        pass_indices = recycler.get_moved_indices()
    else:
        pass_indices = None
    disp = get_displacement(reactant, indices = pass_indices)
    searches = []
    
    invariants = {}

    reactIO = StringIO.StringIO()
    io.savecon(reactIO, reactant)
    invariants['reactant_passed.con']=reactIO
    
    ini_changes = [ ('Main', 'job', 'process_search') ]
    invariants['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)

    #Merge potential files into invariants
    invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

    t1 = unix_time.time()
    if config.recycling_on:
        nrecycled = 0
    for i in range(num_to_make):
        search = {}
        # The search dictionary contains the following key-value pairs:
        # id - CurrentState_WUID
        # displacement - an atoms object containing the point the saddle search will start at
        # mode - an Nx3 numpy array containing the initial mode 
        search['id'] = "%d_%d" % (current_state.get_nr(), akmc_info.get_wuid())
        done = False
#        # Do we want to try superbasin recycling? If yes, try. If we fail to recycle the basin,
#        # move to the next case
#        if (config.sb_recycling_on and current_state.get_nr() is not 0):
#            displacement, mode = sb_recycler.make_suggestion()
#            if displacement:
#                nrecycled += 1
#                if config.debug_list_search_results:
#                    try:
#                        searchdata["%d_%d" %(current_state.get_nr(), akmc_info.get_wuid())] = {}
#                        searchdata["%d_%d" %(current_state.get_nr(), akmc_info.get_wuid())]["type"] = "recycling"
#                    except:
#                        logger.warning("Failed to add searchdata for search %d_%d" % (current_state.get_nr(), akmc_info.get_wuid))
#                done = True
#        # Do we want to do recycling? If yes, try. If we fail to recycle, we move to the next case
#        if (recycler and current_state.get_nr() is not 0):
#            displacement, mode = recycler.make_suggestion()
#            if displacement:
#                nrecycled += 1
#                if config.debug_list_search_results:                
#                    try:
#                        searchdata["%d_%d" %(current_state.get_nr(), akmc_info.get_wuid())] = {}
#                        searchdata["%d_%d" % (current_state.get_nr(), akmc_info.get_wuid())]["type"] = "recycling"
#                    except:
#                        logger.warning("Failed to add searchdata for search %d_%d" % (current_state.get_nr(), akmc_info.get_wuid))
#                done = True
#        if not done and config.kdb_on:
#            # Set up the path for keeping the suggestion if config.kdb_keep is set.
#            keep_path = None
#            if config.kdb_keep:
#                if not os.path.isdir(os.path.join(current_state.path, "kdbsuggestions")):
#                    os.mkdir(os.path.join(current_state.path, "kdbsuggestions"))
#                keep_path = os.path.join(current_state.path, "kdbsuggestions", str(akmc_info.get_wuid()))
#            displacement, mode = kdber.make_suggestion(keep_path)
#            if displacement:
#                done = True
#                logger.info('Made a KDB suggestion')
#                if config.debug_list_search_results:                
#                    try:
#                        searchdata["%d_%d" %(current_state.get_nr(), akmc_info.get_wuid())] = {}
#                        searchdata["%d_%d" % (current_state.get_nr(), akmc_info.get_wuid())]["type"] = "kdb"
#                    except:
#                        logger.warning("Failed to add searchdata for search %d_%d" % (current_state.get_nr(), akmc_info.get_wuid()))
        if not done:
            displacement, mode = disp.make_displacement() 
            if config.debug_list_search_results:                
                try:
                    searchdata["%d_%d" %(current_state.get_nr(), akmc_info.get_wuid())] = {}
                    searchdata["%d_%d" % (current_state.get_nr(), akmc_info.get_wuid())]["type"] = "random"
                except:
                    logger.warning("Failed to add searchdata for search %d_%d" % (current_state.get_nr(), akmc_info.get_wuid()))
        dispIO = StringIO.StringIO()
        io.savecon(dispIO, displacement)
        search['displacement_passed.con'] = dispIO
        
        modeIO = StringIO.StringIO()
        io.save_mode(modeIO, mode, reactant)
        search['mode_passed.dat'] = modeIO
        searches.append(search) 
        akmc_info.increase_wuid()

    if config.recycling_on and nrecycled > 0:
        logger.debug("Recycled %i saddles" % nrecycled)

    try:
        comm.submit_jobs(searches, invariants)
        t2 = unix_time.time()
        logger.info( str(num_to_make) + " searches created") 
        logger.debug( str(num_to_make/(t2-t1)) + " searches per second")
    except:
        logger.exception("Failed to submit searches.")
    return akmc_info.get_wuid

def get_statelist():
    initial_state_path = os.path.join(config.path_root, 'reactant.con') 
    return statelist.AKMCStateList(initial_state_path)#config.kT, 
#                               config.akmc_thermal_window, 
#                               config.akmc_max_thermal_window, 
#                               initial_state_path, 
#                               filter_hole = config.disp_moved_only)  


    
def get_displacement(reactant, indices=None):
    if config.disp_type == 'random':
        disp = displace.Random(reactant, config.disp_magnitude, config.disp_radius, hole_epicenters=indices)
    elif config.disp_type == 'under_coordinated':
        disp = displace.Undercoordinated(reactant, config.disp_max_coord, config.disp_magnitude, config.disp_radius, hole_epicenters=indices, cutoff=config.comp_neighbor_cutoff, use_covalent=config.comp_use_covalent, covalent_scale=config.comp_covalent_scale)
    elif config.disp_type == 'least_coordinated':
        disp = displace.Leastcoordinated(reactant, config.disp_magnitude, config.disp_radius, hole_epicenters=indices, cutoff=config.comp_neighbor_cutoff, use_covalent=config.comp_use_covalent, covalent_scale=config.comp_covalent_scale)
    elif config.disp_type == 'not_FCC_HCP_coordinated':
        disp = displace.NotFCCorHCP(reactant, config.disp_magnitude, config.disp_radius, hole_epicenters=indices, cutoff=config.comp_neighbor_cutoff, use_covalent=config.comp_use_covalent, covalent_scale=config.comp_covalent_scale)
    elif config.disp_type == 'listed_atoms':
        disp = displace.ListedAtoms(reactant, config.disp_magnitude, config.disp_radius, hole_epicenters=indices, cutoff=config.comp_neighbor_cutoff, use_covalent=config.comp_use_covalent, covalent_scale=config.comp_covalent_scale)
    elif config.disp_type == 'water':
        disp = displace.Water(reactant, config.stdev_translation, config.stdev_rotation, config.molecule_list, config.disp_at_random)
    else:
        raise ValueError()
    return disp







def main():
    optpar = optparse.OptionParser(usage = "usage: %prog [options] config.ini")
#    optpar.add_option("-R", "--reset", action="store_true", dest="reset", default = False, help="reset the aKMC simulation, discarding all data")
#    optpar.add_option("-f", "--force", action="store_true", dest="force", default = False, help="force a reset, no questions asked")
#    optpar.add_option("-r", "--restart", action="store_true", dest="restart", default = False, help="restart the aKMC simulations from a clean dynamics.txt file")
#    optpar.add_option("-s", "--status", action="store_true", dest="print_status", default = False, help = "print the status of the simulation and currently running jobs")
    optpar.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False,help="only write to the log file")
#    optpar.add_option("-m", "--movie", action="store", dest="movie_type", default = "", help="Specify the type of movie to make [dynamics, states, fastestpath, fastestfullpath, graph, processes]. Process movies are specified like so: --movie processes,statenumber,processlimit. Where processes is the string processes, statenumber is the number of the state that you want to view, and process limit is the maximum number of processes you would like in the movie. The returned processes are reverse sorted by rate such that the fastest processes is the first in the movie.")
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
        config.comm_search_buffer_size = 0

    #setup logging
    logging.basicConfig(level=logging.DEBUG,
            filename=os.path.join(config.path_results, "test.log"),
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
        test(config)
    else:
        logger.info("the server is locked by pid %i" % lock.pid)
        sys.exit(1)

if __name__ == '__main__':

# Create the statelist directory if it does not exist.
    if not os.path.isdir(config.path_coarse_states):
        storage = os.path.join(config.path_coarse_states, "storage")
        logger.warning("coarse state path does not exist, creating %s" % config.path_coarse_states)
        os.makedirs(config.path_coarse_states)
        os.makedirs(storage)


    sampled = SampledSpace()
    extend_sampled_space(sampled, 0)
