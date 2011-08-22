import logging
logger = logging.getLogger('explorer')
from time import time
import shutil
import StringIO
import os
import numpy

import communicator
import config
import displace
import fileio as io
import kdbing
import recycling

def get_minmodexplorer():
    if config.akmc_server_side:
        return ServerMinModeExplorer
    else:
        return ClientMinModeExplorer

class Explorer:
    def __init__(self):
        self.wuid_path = os.path.join(config.path_scratch, "wuid")
        self.load_wuid()

    def load_wuid(self):
        try:
            f = open(self.wuid_path)
            self.wuid = int(f.read())
            f.close()
        except IOError:
            self.wuid = 0

    def save_wuid(self):
        f = open(self.wuid_path, 'w')
        f.write("%i\n" % self.wuid)
        f.close()



class MinModeExplorer(Explorer):
    def __init__(self, states, previous_state, state):
        Explorer.__init__(self)
        self.states = states
        self.state = state
        self.previous_state = previous_state
        self.comm = communicator.get_communicator()

        if config.recycling_on: 
            self.nrecycled = 0
            self.recycler = recycling.Recycling(self.states,
                                           self.previous_state, 
                                           self.state, 
                                           config.recycling_move_distance,
                                           config.recycling_save_sugg)

        # If we plan to only displace atoms that moved getting to the current state.
        if config.disp_moved_only and self.state.number != 0:
            moved_atoms = self.recycler.process_atoms
        else:
            moved_atoms = None


        if config.kdb_on:
            self.kdber = kdbing.KDB()
            if len(self.state.get_ratetable()) <= 1:
                self.kdber.query(self.state, wait = config.kdb_wait)


        self.reactant = self.state.get_reactant()
        self.displace = displace.DisplacementManager(self.reactant, moved_atoms)

    def explore(self):
        self.register_results()
        if self.state.get_confidence() < config.akmc_confidence:
            self.make_searches()
        else:
            num_cancelled = self.comm.cancel_state(self.state.number)
            logger.info("cancelled %i workunits from state %i", 
                        num_cancelled, self.state.number)
            #XXX: Do we ever call explore a completed state twice?
            if config.kdb_on:
                logger.debug("Adding relevant processes to kinetic database.")
                for process_id in self.state.get_process_ids():
                    output = self.kdber.add_process(self.state, process_id)
                    logger.debug("kdbaddpr.pl: %s" % output)

    def generate_displacement(self):
        if config.recycling_on and self.state.number is not 0:
            displacement, mode = self.recycler.make_suggestion()
            if displacement:
                self.nrecycled += 1
                return displacement, mode, 'recycling'

        if config.kdb_on:
            displacement, mode = self.kdber.make_suggestion()
            if displacement:
                return displacement, mode, 'kdb'
                logger.info('Made a KDB suggestion')

        displacement, mode = self.displace.make_displacement() 
        return displacement, mode, 'random'

class ServerMinModeExplorer(MinModeExplorer):
    def __init__(self, states, previous_state, state):
        MinModeExplorer.__init__(self, states, previous_state, state)
        job_table_path = os.path.join(config.path_root, "jobs.tbl")
        job_table_columns = [ 'state', 'search_id', 'wuid', 'job_type',
                              'displace_type', 'running', 'fcs', 'finished',
                              'exit_status' ]
        self.job_table = io.Table(job_table_path, job_table_columns)

        if len(self.job_table) > 0:
            self.search_id = self.job_table.max_value('search_id')+1
            logger.info("read in search id of %i" % (self.search_id))
        else:
            self.search_id = 0

    def register_saddle_search(self, result):
        if "results.dat" not in result:
            logger.fatal("client crashed")
        termination_code = result['results']['termination_reason']
        reason_strings = [ "good", "", "no_convex", "high_energy",
                           "max_concave_iterations", 
                           "max_iterations" ]
        search_finished = False
        if termination_code == 0:
            search_finished = True

        reason_string = reason_strings[termination_code]
        self.job_table.get_row('wuid', result['wuid'])['exit_status'] = reason_string
        logger.info("saddle search status: %s" % reason_string)
        self.save_incomplete_result(result)

    def save_incomplete_result(self, result):
        dir_path = os.path.join(config.path_incomplete, result['name'])
        os.makedirs(dir_path)
        for k in result:
            if hasattr(result[k], 'getvalue'):
                fn = os.path.join(dir_path, k)
                f = open(fn, "w")
                f.write(result[k].getvalue())
                f.close()

    def register_results(self):
        logger.info("registering results")
        t1 = time()
        if os.path.isdir(config.path_jobs_in):
            shutil.rmtree(config.path_jobs_in)  
        os.makedirs(config.path_jobs_in)

        if not os.path.isdir(config.path_incomplete):
            os.makedirs(config.path_incomplete)
        
        #Function used by communicator to determine whether to discard a result
        def keep_result(name):
            state_num = int(name.split("_")[0])
            return (config.debug_register_extra_results or \
                    state_num == self.state.number or \
                    self.state.get_confidence() < config.akmc_confidence)

        num_registered = 0
        for result in self.comm.get_results(config.path_jobs_in, keep_result): 
            state_num = int(result['name'].split("_")[0])
            id = int(result['name'].split("_")[1]) + result['number']
            searchdata_id = "%d_%d" % (state_num, id)

            job_type = self.job_table.get_row('wuid', id)['job_type']
            self.job_table.get_row('wuid', id)['finished'] = "true"
            self.job_table.get_row('wuid', id)['running'] = "false"
            result['wuid'] = id
            result['results'] = io.parse_results(result['results.dat'])
            self.job_table.get_row('wuid', id)['fcs'] = result['results']['total_force_calls']
            num_registered += 1

            if job_type == "saddle_search":
                self.register_saddle_search(result)
            
            if self.state.get_confidence() >= config.akmc_confidence:
                if not config.debug_register_extra_results:
                    break

            self.job_table.write()
        
        #Approximate number of searches recieved
        tot_searches = len(os.listdir(config.path_jobs_in)) * config.comm_job_bundle_size
        
        t2 = time()
        logger.info("%i (result) searches processed", num_registered)
        logger.info("Approximately %i (result) searches discarded." % (tot_searches - num_registered))
        if num_registered == 0:
            logger.debug("0 results per second", num_registered)
        else:
            logger.debug("%.1f results per second", (num_registered/(t2-t1)))
            
        return num_registered

    def make_searches(self):
        num_in_buffer = self.comm.get_queue_size()*config.comm_job_bundle_size 
        logger.info("%i searches in the queue" % num_in_buffer)
        num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
        logger.info("making %i searches" % num_to_make)
        
        if num_to_make == 0:
            return
        
        jobs = []
        
        invariants = {}

        #Merge potential files into invariants
        invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

        t1 = time()
        for i in range(num_to_make):
            job = {}
            search_status = {}
            for row in self.job_table.get_rows("state", self.state.number):
                if row['job_type'] == "saddle_search":
                    if row['exit_status'] == "good":
                        search_status[row['search_id']] = "complete"
                    elif row['running'] == 'true':
                        search_status[row['search_id']] = "inprogress"
                    else:
                        search_status[row['search_id']] = "incomplete"

            new_job = True
            print "search_status", search_status
            for search_id, status in search_status.iteritems():
                if status != "incomplete": continue
                #need to result search from last one
                rows = self.job_table.get_rows("search_id", search_id)
                rows = [ r for r in rows if r["job_type"] == "saddle_search" ]
                row = rows[-1]
                resume_id = "%i_%i" % (row["state"], row["wuid"]) 
                resume_path = os.path.join(config.path_incomplete, resume_id)
                reactIO = StringIO.StringIO()
                io.savecon(reactIO, self.reactant)
                dispIO = StringIO.StringIO()
                displacement = io.loadcon(os.path.join(resume_path, "saddle.con"))
                io.savecon(dispIO, displacement)
                modeIO = StringIO.StringIO(file(os.path.join(resume_path, "mode.dat")).read())
                new_job = False
                self.job_table.add_row( {'state':self.state.number,
                                         'search_id':search_id,
                                         'wuid':self.wuid,
                                         'job_type':"saddle_search",
                                         'displace_type':row['displace_type'],
                                         'running':"true",
                                         'fcs':0,
                                         'finished':"false",
                                         'exit_status':"-"} )

            if new_job:
                displacement, mode, disp_type = self.generate_displacement()
                dispIO = StringIO.StringIO()
                io.savecon(dispIO, displacement)
                modeIO = StringIO.StringIO()
                io.save_mode(modeIO, mode)
                self.job_table.add_row( {'state':self.state.number,
                                         'search_id':self.search_id,
                                         'wuid':self.wuid,
                                         'job_type':"saddle_search",
                                         'displace_type':disp_type,
                                         'running':"true",
                                         'fcs':0,
                                         'finished':"false",
                                         'exit_status':"-"} )

                 
            job['id'] = "%d_%d" % (self.state.number, self.wuid)
            reactIO = StringIO.StringIO()
            job['reactant_passed.con'] = reactIO
            ini_changes = [ ('Main', 'job', 'saddle_search') ]
            job['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)
            io.savecon(reactIO, self.reactant)

            job['displacement_passed.con'] = dispIO
            job['mode_passed.dat'] = modeIO

            self.search_id += 1

            jobs.append(job) 
            self.wuid += 1
            # eager write
            self.save_wuid()

        if config.recycling_on and self.nrecycled > 0:
            logger.info("recycled %i saddles" % self.nrecycled)

        try:
            self.comm.submit_jobs(jobs, invariants)
            t2 = time()
            logger.info( str(len(jobs)) + " searches created") 
            logger.debug( str(num_to_make/(t2-t1)) + " searches per second")
        except:
            logger.exception("Failed to submit searches.")


class ClientMinModeExplorer(MinModeExplorer):
    def __init__(self, states, previous_state, state):
        MinModeExplorer.__init__(self, states, previous_state, state)
        job_table_path = os.path.join(config.path_root, "jobs.tbl")
        job_table_columns = [ 'state', 'wuid', 'type']
        self.job_table = io.Table(job_table_path, job_table_columns)

    def make_searches(self):
        #XXX:what if the user changes the bundle size?
        num_in_buffer = self.comm.get_queue_size()*config.comm_job_bundle_size 
        logger.info("%i searches in the queue" % num_in_buffer)
        num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
        logger.info("making %i searches" % num_to_make)
        
        if num_to_make == 0:
            return
        
        searches = []
        
        invariants = {}

        reactIO = StringIO.StringIO()
        io.savecon(reactIO, self.reactant)
        invariants['reactant_passed.con']=reactIO
        
        ini_changes = [ ('Main', 'job', 'process_search') ]
        invariants['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)

        #Merge potential files into invariants
        invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

        t1 = time()
        for i in range(num_to_make):
            search = {}
            # The search dictionary contains the following key-value pairs:
            # id - CurrentState_WUID
            # displacement - an atoms object containing the point the saddle search will start at
            # mode - an Nx3 numpy array containing the initial mode 
            search['id'] = "%d_%d" % (self.state.number, self.wuid)
            displacement, mode, disp_type = self.generate_displacement()
            self.job_table.add_row( {'state':self.state.number,
                                     'wuid':self.wuid,
                                     'type':disp_type } )

            if displacement:
                dispIO = StringIO.StringIO()
                io.savecon(dispIO, displacement)
                search['displacement_passed.con'] = dispIO
                modeIO = StringIO.StringIO()
                io.save_mode(modeIO, mode)
                search['mode_passed.dat'] = modeIO
                searches.append(search) 
                self.wuid += 1
                # eager write
                self.save_wuid()

        if config.recycling_on and self.nrecycled > 0:
            logger.info("recycled %i saddles" % self.nrecycled)

        try:
            self.comm.submit_jobs(searches, invariants)
            t2 = time()
            logger.info( str(len(searches)) + " searches created") 
            logger.debug( str(num_to_make/(t2-t1)) + " searches per second")
        except:
            logger.exception("Failed to submit searches.")
        self.job_table.write()

    def register_results(self):
        logger.info("registering results")
        t1 = time()
        if os.path.isdir(config.path_jobs_in):
            shutil.rmtree(config.path_jobs_in)  
        os.makedirs(config.path_jobs_in)
        
        #Function used by communicator to determine whether to discard a result
        def keep_result(name):
            state_num = int(name.split("_")[0])
            return (config.debug_register_extra_results or \
                    state_num == self.state.number or \
                    self.state.get_confidence() < config.akmc_confidence)

        num_registered = 0
        for result in self.comm.get_results(config.path_jobs_in, keep_result): 
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
                #XXX: We should only do these checks once to speed things up, 
                #     but at the same time debug options don't have to be fast
                #save_path = os.path.join(config.path_root, "old_searches")
                #if not os.path.isdir(save_path):
                #    os.mkdir(save_path)
                #shutil.copytree(result_path, os.path.join(save_path, i))
                #XXX: This is currently broken by the new result passing 
                #      scheme. Should it be done in communicator?
                pass
            state_num = int(result['name'].split("_")[0])
            id = int(result['name'].split("_")[1]) + result['number']
            searchdata_id = "%d_%d" % (state_num, id)
            # Store information about the search into result_data for the 
            # search_results.txt file in the state directory.
            job_type = self.job_table.get_row('wuid', id)['type']
            result['type'] = job_type
            if job_type == None:
                logger.warning("Could not find search data for search %s" 
                               % searchdata_id)
            else:
                self.job_table.delete_row('wuid', id)
            result['wuid'] = id
            
            #read in the results
            result['results'] = io.parse_results(result['results.dat'])
            if result['results']['termination_reason'] == 0:
                self.state.add_process(result)
            else:
                self.state.register_bad_saddle(result, config.debug_keep_bad_saddles)
            num_registered += 1
            
            if self.state.get_confidence() >= config.akmc_confidence:
                if not config.debug_register_extra_results:
                    break
        
        #Approximate number of searches recieved
        tot_searches = len(os.listdir(config.path_jobs_in)) * config.comm_job_bundle_size
        
        t2 = time()
        logger.info("%i (result) searches processed", num_registered)
        logger.info("Approximately %i (result) searches discarded." % (tot_searches - num_registered))
        if num_registered == 0:
            logger.debug("0 results per second", num_registered)
        else:
            logger.debug("%.1f results per second", (num_registered/(t2-t1)))
            
        self.job_table.write()
        return num_registered
