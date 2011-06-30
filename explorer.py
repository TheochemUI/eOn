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
import io
import kdbing
import recycling

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

        job_table_path = os.path.join(config.path_root, "jobs.tbl")
        job_table_columns = [ 'state', 'wuid', 'type']
        self.job_table = io.Table(job_table_path, job_table_columns)

        #clean jobs from old states
        #self.job_table.delete_row_func('state', lambda s: s != self.state.number)

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
        if config.disp_type == 'random':
            self.displace = displace.Random(self.reactant, config.disp_magnitude, 
                                   config.disp_radius, hole_epicenters=moved_atoms)
        elif config.disp_type == 'under_coordinated':
            self.displace = displace.Undercoordinated(self.reactant, config.disp_max_coord, 
                                             config.disp_magnitude, 
                                             config.disp_radius, 
                                             hole_epicenters=moved_atoms, 
                                             cutoff=config.comp_neighbor_cutoff,
                                             use_covalent=config.comp_use_covalent,
                                             covalent_scale=config.comp_covalent_scale)
        elif config.disp_type == 'least_coordinated':
            self.displace = displace.Leastcoordinated(self.reactant, config.disp_magnitude,
                                             config.disp_radius, 
                                             hole_epicenters=moved_atoms,
                                             cutoff=config.comp_neighbor_cutoff,
                                             use_covalent=config.comp_use_covalent, 
                                             covalent_scale=config.comp_covalent_scale)
        elif config.disp_type == 'not_FCC_HCP_coordinated':
            self.displace = displace.NotFCCorHCP(self.reactant, config.disp_magnitude,
                                        config.disp_radius, hole_epicenters=moved_atoms,
                                        cutoff=config.comp_neighbor_cutoff,
                                        use_covalent=config.comp_use_covalent,
                                        covalent_scale=config.comp_covalent_scale)
        elif config.disp_type == 'listed_atoms':
            self.displace = displace.ListedAtoms(self.reactant, config.disp_magnitude, 
                                        config.disp_radius, hole_epicenters=moved_atoms,
                                        cutoff=config.comp_neighbor_cutoff,
                                        use_covalent=config.comp_use_covalent,
                                        covalent_scale=config.comp_covalent_scale)
        elif config.disp_type == 'water':
            self.displace = displace.Water(self.reactant, config.stdev_translation,
                                  config.stdev_rotation, config.molecule_list, 
                                  config.disp_at_random)
        else:
            logger.error("Unknown displacement type: %s", config.disp_type)
            raise ValueError()


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
                io.save_mode(modeIO, mode, self.reactant)
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
            
        return num_registered
