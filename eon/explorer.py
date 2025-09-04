
import logging
logger = logging.getLogger('explorer')
from time import time
import shutil
import io
import os
import sys
import pickle as pickle
from copy import copy
import numpy

from eon import atoms
from eon import communicator
from eon import displace
from eon import fileio as io
#import kdb
from eon import recycling
from eon import eon_kdb as kdb

from eon.config import config as EON_CONFIG
from eon.config import ConfigClass # Typing

def get_minmodexplorer(config: ConfigClass = EON_CONFIG):
    if config.akmc_server_side_process_search:
        return ServerMinModeExplorer
    else:
        return ClientMinModeExplorer

class Explorer:
    def __init__(self, superbasin=None, config: ConfigClass = EON_CONFIG):
        self.config = config
        self.wuid_path = os.path.join(self.config.path_scratch, "wuid")
        self.superbasin = superbasin
        self.load_wuid()

    def load_wuid(self):
        try:
            f = open(self.wuid_path)
            self.wuid = int(f.read())
            f.close()
        except IOError:
            if not self.superbasin:
                self.wuid = 0
            else:
                # In a superbasin search, WUIDs would be repeated. So
                # we assume there will never be more than 1e6 searches
                # per state, so we start at superbasin number times
                # one million.
                self.wuid = self.superbasin.id * 1000000

    def save_wuid(self):
        f = open(self.wuid_path, 'w')
        f.write("%i\n" % self.wuid)
        f.close()


class MinModeExplorer(Explorer):
    def __init__(self, states, previous_state, state, superbasin=None, config: ConfigClass = EON_CONFIG):
        """Init MinModeExplorer.

        If superbasin is passed, all states in the superbasin will be
        searched until all of them have the required confidence. It
        will also lead to a different confidence definition being
        used, which takes only processes exiting the SB into account.

        """
        Explorer.__init__(self, superbasin, config = config)
        self.states = states
        self.state = state
        self.previous_state = previous_state
        self.comm = communicator.get_communicator(self.config)

        if self.config.recycling_on:
            self.nrecycled = 0
            self.recycler = recycling.Recycling(
                self.states,
                self.previous_state,
                self.state,
                self.config.recycling_move_distance,
                self.config.recycling_save_sugg,
                config=self.config,
            )

        # If we plan to only displace atoms that moved getting to the current state.
        if self.config.disp_moved_only and self.state.number != 0:
            moved_atoms = self.recycler.process_atoms
        else:
            moved_atoms = None

        if self.config.kdb_on:
            if not os.path.isdir(self.config.kdb_scratch_path):
                os.makedirs(self.config.kdb_scratch_path)
            try:
                queried = [int(q) for q in open(os.path.join(self.config.kdb_scratch_path, "queried"), 'r').readlines()]
            except:
                queried = []
            if self.state.number not in queried:
                queried.append(self.state.number)
                f = open(os.path.join(self.config.kdb_scratch_path, "queried"), 'w')
                for q in queried:
                    f.write("%d\n" % q)
                f.close()
                kdb.query(self.state)

        self.reactant = self.state.get_reactant()
        self.displace = displace.DisplacementManager(self.reactant, moved_atoms, config = self.config)

    def explore(self):
        self.register_results()
        if self.state.get_confidence(self.superbasin) < self.config.akmc_confidence:
            self.make_jobs()
        else:
            num_cancelled = self.comm.cancel_state(self.state.number)
            logger.info("Cancelled %i workunits from state %i",
                        num_cancelled, self.state.number)
            #XXX: Do we ever call explore on a completed state twice?
            if self.config.kdb_on:
                logger.info("Adding relevant processes to kinetic database")
                for process_id in self.state.get_process_ids():
                    output = kdb.insert(self.state, process_id)
                    logger.debug("kdb insert: %s", output)

    def generate_displacement(self):
        if self.config.recycling_on and self.state.number != 0:
            displacement, mode = self.recycler.make_suggestion()
            if displacement:
                self.nrecycled += 1
                return displacement, mode, 'recycling'

        if self.config.kdb_on:
            displacement, mode = kdb.make_suggestion()
            if displacement:
                displacement.mass = self.reactant.mass
                logger.info('Made a KDB suggestion')
                return displacement, mode, 'kdb'

        if self.config.saddle_method == 'dynamics':
            return None, None, 'dynamics'

        displacement, mode = self.displace.make_displacement()
        return displacement, mode, 'random'


class ClientMinModeExplorer(MinModeExplorer):
    def __init__(self, states, previous_state, state, superbasin=None):
        MinModeExplorer.__init__(self, states, previous_state, state, superbasin)
        job_table_path = os.path.join(self.config.path_root, "jobs.tbl")
        job_table_columns = [ 'state', 'wuid', 'type']
        self.job_table = io.Table(job_table_path, job_table_columns)

        if self.superbasin:
            self.job_table.delete_row_func('state', lambda s: s not in self.superbasin.state_dict)
        else:
            self.job_table.delete_row_func('state', lambda s: s != state.number)

    def make_jobs(self):
        #XXX:what if the user changes the bundle size?
        num_in_buffer = self.comm.get_queue_size()*self.config.comm_job_bundle_size
        logger.info("Queue contains %i searches" % num_in_buffer)
        num_to_make = max(self.config.comm_job_buffer_size - num_in_buffer, 0)
        logger.info("Making %i process searches" % num_to_make)

        if num_to_make == 0:
            return

        searches = []

        invariants = {}

        reactIO = io.StringIO()
        io.savecon(reactIO, self.reactant)
        file_permission = os.stat("pos.con").st_mode
        invariants['pos.con'] = (reactIO, file_permission)

        t1 = time()
        if self.config.saddle_method == 'dynamics' and \
                self.config.recycling_on and \
                self.config.disp_moved_only and \
                self.state.number != 0:

            moved_atoms = self.recycler.process_atoms
            mass_weights = self.reactant.mass.copy()
            mass_weights *= self.config.recycling_mass_weight_factor

            for i in range(len(self.reactant)):
                if i in moved_atoms:
                    mass_weights[i] = self.reactant.mass[i]

            weightsIO = io.StringIO()
            numpy.savetxt(weightsIO, mass_weights)
#            file_permission = os.stat("masses.dat").st_mode
            invariants['masses.dat'] = (weightsIO, file_permission)

        # Merge potential files into invariants
        invariants = dict(invariants, **io.load_potfiles(self.config.path_pot))

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

            ini_changes = [ ('Main', 'job', 'process_search'),
                            ('Main', 'random_seed',
                                str(int(numpy.random.random()*10**9))),
                          ]
            # if we are recycling a saddle, but using "dynamics saddle search" we need
            # to switch to min_mode searches
            if self.config.saddle_method == 'dynamics' and disp_type != 'dynamics':
                ini_changes.append( ('Saddle Search', 'method', 'min_mode') )

            search['config.ini'] = io.modify_config(self.config.config_path, ini_changes)

            if displacement:
                dispIO = io.StringIO()
                io.savecon(dispIO, displacement)
                search['displacement.con'] = dispIO
                modeIO = io.StringIO()
                io.save_mode(modeIO, mode)
                search['direction.dat'] = modeIO

            searches.append(search)
            self.wuid += 1
            # eager write
            self.save_wuid()

        if self.config.recycling_on and self.nrecycled > 0:
            logger.info("Recycled %i saddles" % self.nrecycled)

        try:
            self.comm.submit_jobs(searches, invariants)
            t2 = time()
            logger.info( "Created " + str(len(searches)) + " searches")
            #logger.debug( "Created " + str(num_to_make/(t2-t1)) + " searches per second")
            logger.debug( "Created %.2f searches per second", num_to_make/(t2-t1))
        except:
            logger.exception("Failed to submit searches")
        self.job_table.write()

    def register_results(self):
        logger.info("Registering results")
        t1 = time()
        if os.path.isdir(self.config.path_jobs_in):
            try:
                shutil.rmtree(self.config.path_jobs_in)
            except (OSError, IOError):
                pass
        if not os.path.isdir(self.config.path_jobs_in):
            os.makedirs(self.config.path_jobs_in)

        # Function used by communicator to determine whether to discard a result
        def keep_result(name):
            # note that all processes are assigned to the current state
            state_num = int(name.split("_")[0])
            if self.superbasin:
                return (state_num in self.superbasin.state_dict and
                        self.superbasin.get_confidence() < self.config.akmc_confidence)
            else:
                return (state_num == self.state.number and
                        self.state.get_confidence() < self.config.akmc_confidence)

        num_registered = 0
        for result in self.comm.get_results(self.config.path_jobs_in, keep_result):
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
            if self.config.debug_keep_all_results:
                #XXX: We should only do these checks once to speed things up,
                #     but at the same time debug options don't have to be fast
                # save_path = os.path.join(self.config.path_root, "old_searches")
                # if not os.path.isdir(save_path):
                #    os.mkdir(save_path)
                # shutil.copytree(result_path, os.path.join(save_path, i))
                # XXX: This is currently broken by the new result passing
                #      scheme. Should it be done in communicator?
                pass
            if len(result) == 0: continue
            state_num = int(result['name'].split("_")[0])
            id = int(result['name'].split("_")[1]) + result['number']
            searchdata_id = "%d_%d" % (state_num, id)
            # Store information about the search into result_data for the
            # search_results.txt file in the state directory.
            try:
                job_type = self.job_table.get_row('wuid', id)['type']
            except TypeError:
                logger.warning("Could not find job type for search %s"
                               % searchdata_id)
                continue
            result['type'] = job_type
            if job_type is None:
                logger.warning("Could not find search data for search %s"
                               % searchdata_id)
            else:
                self.job_table.delete_row('wuid', id)
            result['wuid'] = id

            # If we are doing a search for a superbasin the results
            # could be for a different state.
            if self.superbasin:
                try:
                    state = self.superbasin.state_dict[state_num]
                except KeyError:
                    logger.warning("State of job %s is not part of "
                                   "the superbasin" % result['name'])
                    continue
            else:
                state = self.state

            # read in the results
            result['results'] = io.parse_results(result['results.dat'])
            if result['results']['termination_reason'] == 0:
                state.add_process(result, self.superbasin)
            else:
                state.register_bad_saddle(result, self.config.debug_keep_bad_saddles, superbasin=self.superbasin)
            num_registered += 1

            if ((self.superbasin and self.superbasin.get_confidence() >= self.config.akmc_confidence) or
                (not self.superbasin and self.state.get_confidence() >= self.config.akmc_confidence)):
                if not self.config.debug_register_extra_results:
                    break

        # Approximate number of searches received
        tot_searches = len(os.listdir(self.config.path_jobs_in)) * self.config.comm_job_bundle_size

        t2 = time()
        logger.info("Processed %i results", num_registered)
        if tot_searches != num_registered:
            logger.info("Discarded approximately %i results" % (tot_searches - num_registered))
        logger.debug("Registered %.1f results per second", (num_registered/(t2-t1)))

        self.job_table.write()
        return num_registered


class ServerMinModeExplorer(MinModeExplorer):

    def __init__(self, states, previous_state, state, superbasin=None):
        #XXX: need to init somehow
        self.search_id = 0

        self.wuid_to_search_id = {}
        self.process_searches = {}
        self.job_info = {}

        if os.path.isfile("explorer.pickle"):
            f = open("explorer.pickle", "rb")
            tmp_dict = pickle.load(f)
            f.close()
            self.__dict__.update(tmp_dict)

        MinModeExplorer.__init__(self, states, previous_state, state, superbasin)

    def save(self):
        f = open("explorer.pickle", "w")
        d = self.__dict__.copy()
        del d['states']
        del d['previous_state']
        del d['state']
        del d['comm']
        pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
        f.close()

        f = open("searches.log", 'w')
        f.write("%9s %13s %11s %10s\n" % ("search_id", "job_type", "status", "job_name"))
        for search_id in self.job_info:
            fmt = "%9i %13s %11s %10s\n"
            lines = {'saddle_search':None, 'min1':None, 'min2':None}
            for name, job in list(self.job_info[search_id].items()):
                lines[job['type']] = fmt % (search_id, job['type'], job['status'], name)

            if lines['saddle_search']:
                f.write(lines['saddle_search'])
            if lines['min1']:
                f.write(lines['min1'])
            else:
                f.write(fmt % (search_id, 'min1', 'not_running', ''))
            if lines['min2']:
                f.write(lines['min2'])
            else:
                f.write(fmt % (search_id, 'min2', 'not_running', ''))
        f.close()

    def explore(self):
        if not os.path.isdir(self.config.path_jobs_in): #XXX: does this condition ever happen?
            os.makedirs(self.config.path_jobs_in)
            if self.state.get_confidence(self.superbasin) >= self.config.akmc_confidence:
                self.process_searches = {}
                self.save()

        MinModeExplorer.explore(self)

    def register_results(self):
        logger.info("Registering results")
        t1 = time()
        if os.path.isdir(self.config.path_jobs_in):
            try:
                shutil.rmtree(self.config.path_jobs_in)
            except OSError as msg:
                logger.error("Error cleaning up %s: %s", self.config.path_jobs_in, msg)
            else:
                os.makedirs(self.config.path_jobs_in)

        if not os.path.isdir(self.config.path_incomplete):
            os.makedirs(self.config.path_incomplete)

        # Function used by communicator to determine whether to keep a result
        def keep_result(name):
            # note that all processes are assigned to the current state
            state_num = int(name.split("_")[0])
            return (state_num == self.state.number and \
                    self.state.get_confidence(self.superbasin) < self.config.akmc_confidence)

        num_registered = 0
        for result in self.comm.get_results(self.config.path_jobs_in, keep_result):
            state_num = int(result['name'].split("_")[0])
            # XXX: doesn't this doesn't give the correct id wrt bundling
            id = int(result['name'].split("_")[1]) + result['number']
            searchdata_id = "%d_%d" % (state_num, id)

            search_id = self.wuid_to_search_id[id]
            if search_id not in self.process_searches:
                continue
            self.job_info[search_id][searchdata_id]['status'] = 'complete'
            #logger.info("got result for search_id %i" % search_id)
            final_result = self.process_searches[search_id].process_result(result)
            if final_result:
                results_dict = io.parse_results(final_result['results.dat'])
                reason = results_dict['termination_reason']
                if reason == 0:
                    self.state.add_process(final_result)
                else:
                    final_result['wuid'] = id
                    self.state.register_bad_saddle(final_result, self.config.debug_keep_bad_saddles)
            else:
                ps = self.process_searches[search_id]
                saddle = ps.get_saddle()
                if saddle:
                    barrier = ps.data['barrier_reactant_to_product']
                    if self.state.find_repeat(ps.get_saddle_file(), barrier):
                        self.state.add_process(ps.build_result())
                        del self.process_searches[search_id]
            num_registered += 1
            if self.state.get_confidence(self.superbasin) >= self.config.akmc_confidence:
                if not self.config.debug_register_extra_results:
                    break

        # Approximate number of searches received
        tot_searches = len(os.listdir(self.config.path_jobs_in)) * self.config.comm_job_bundle_size

        t2 = time()
        logger.info("Processed %i results", num_registered)
        if tot_searches != num_registered:
            logger.info("Discarded approximately %i results" % (tot_searches - num_registered))
        logger.debug("Registered %.1f results per second", (num_registered/(t2-t1)))

        self.save()
        return num_registered

    def make_jobs(self):
        num_unsent = self.comm.get_queue_size()*self.config.comm_job_bundle_size
        logger.info("Queued %i jobs" % num_unsent)
        num_in_progress = self.comm.get_number_in_progress()*self.config.comm_job_bundle_size
        logger.info("Running %i jobs" % num_in_progress)
        num_total = num_unsent + num_in_progress
        num_to_make = max(self.config.comm_job_buffer_size - num_unsent, 0)
        if self.config.comm_job_max_size != 0:
            if num_total + num_to_make>= self.config.comm_job_max_size:
                num_to_make = max(0, self.config.comm_job_max_size - num_total)
                logger.info("Reached max_jobs")
        logger.info("Making %i jobs" % num_to_make)

        if num_to_make == 0:
            return

        jobs = []

        invariants = {}

        # Merge potential files into invariants
        invariants = dict(invariants, **io.load_potfiles(self.config.path_pot))

        t1 = time()

        # start new searches
        for i in range(num_to_make):
            job = None
            for ps in list(self.process_searches.values()):
                job, job_type = ps.get_job(self.state.number)
                if job:
                    self.wuid_to_search_id[self.wuid] = ps.search_id
                    sid = ps.search_id
                    break

            if not job:
                displacement, mode, disp_type = self.generate_displacement()
                reactant = self.state.get_reactant()
                process_search = ProcessSearch(reactant, displacement, mode,
                                               disp_type, self.search_id, self.state.number)
                self.process_searches[self.search_id] = process_search
                self.wuid_to_search_id[self.wuid] = self.search_id
                job, job_type = process_search.get_job(self.state.number)
                sid = self.search_id
                self.search_id += 1

            job['id'] = "%i_%i" % (self.state.number, self.wuid)
            job_entry = {'type':job_type, 'status':'running'}

            if sid not in self.job_info:
                self.job_info[sid] = {}
            self.job_info[sid][job['id']] = job_entry

            self.wuid += 1
            self.save_wuid()
            jobs.append(job)

        self.save()
        if self.config.recycling_on and self.nrecycled > 0:
            logger.info("Recycled %i saddles" % self.nrecycled)

        try:
            self.comm.submit_jobs(jobs, invariants)
            t2 = time()
            logger.info( "Created " + str(len(jobs)) + " searches")
            logger.debug( "Created " + str(num_to_make/(t2-t1)) + " searches per second")
        except:
            logger.exception("Failed to submit searches")


class ProcessSearch:
    def __init__ (self, reactant, displacement, mode, disp_type, search_id, state_number):
        self.reactant = reactant
        self.displacement = displacement
        self.mode = mode
        self.search_id = search_id
        self.displacement_type = disp_type
        self.state_number = state_number

        # valid statuses are 'not_started', 'running', 'complete', 'unneeded', and 'error'
        self.job_statuses = {
                             'saddle_search':'not_started',
                             'min1':'not_started',
                             'min2':'not_started'
                            }

        unknown = "unknown_exit_code"
        self.job_termination_reasons = {
                'saddle_search':[ "good", unknown, "no_convex", "high_energy",
                                  "max_concave_iterations",
                                  "max_iterations", unknown, unknown, unknown,
                                  unknown, unknown, "potential_failed",
                                  "nonnegative_abort", "nonlocal abort"],
                'minimization':[ "good", "max_iterations", "potential_failed", ]}

        self.finished_jobs = []

        self.finished_saddle_name = None
        self.finished_min1_name = None
        self.finished_min2_name = None
        self.finished_reactant_name = None
        self.finished_product_name = None

        self.data = {
                      'termination_reason':0,
                      'potential_energy_saddle':None,
                      'potential_energy_reactant':None,
                      'potential_energy_product':None,
                      'barrier_reactant_to_product':None,
                      'barrier_product_to_reactant':None,
                      'prefactor_reactant_to_product':self.config.process_search_default_prefactor,
                      'prefactor_product_to_reactant':self.config.process_search_default_prefactor,
                      'displacement_saddle_distance':0.0,
                      'force_calls_saddle':0,
                      'force_calls_minimization':0,
                      'force_calls_prefactors':0,
                    }

    def get_job(self, state_number):
        if state_number != self.state_number:
            return None, None

        if True in [ s == 'error' for s in list(self.job_statuses.values()) ]:
            return None, None

        if self.job_statuses['saddle_search'] == 'not_started':
            self.job_statuses['saddle_search'] = 'running'
            return self.start_search(), 'saddle_search'

        if self.job_statuses['min1'] == 'not_started' and \
           self.job_statuses['saddle_search'] == 'complete':
            self.job_statuses['min1'] = 'running'
            return self.start_minimization('min1'), 'min1'

        if self.job_statuses['min2'] == 'not_started' and \
           self.job_statuses['saddle_search'] == 'complete':
            self.job_statuses['min2'] = 'running'
            return self.start_minimization('min2'), 'min2'

        return None, None

    def process_result(self, result):
        results_dat = io.parse_results(result['results.dat'])
        #XXX: can remove this line now
        result['results.dat'].seek(0)
        job_type = results_dat['job_type']
        termination_code = results_dat['termination_reason']

        self.save_result(result)
        self.finished_jobs.append(result['name'])

        if job_type == 'saddle_search':
            self.data['termination_reason'] = termination_code
            logger.info("Search_id: %i saddle search complete" % self.search_id)
            if termination_code == 0:
                self.job_statuses[job_type] = 'complete'
            else:
                self.job_statuses[job_type] = 'error'
            self.finished_saddle_name = result['name']
            self.finish_search(result)

        elif job_type == 'minimization':
            if self.job_statuses['min1'] == 'running':
                min_name = 'min1'
                min_number = 1
                self.finished_min1_name = result['name']
            else:
                min_name = 'min2'
                min_number = 2
                self.finished_min2_name = result['name']

            self.job_statuses[min_name] = 'complete'
            logger.info("Search_id: %i minimization %i complete" % (self.search_id, min_number))

            if min_number == 2:
                self.finish_minimization(result)

        done = False not in  [ s == 'complete' for s in list(self.job_statuses.values()) ]
        if not done:
            done = True in [ s == 'error' for s in list(self.job_statuses.values()) ]

        if done:
            return self.build_result()

    def build_result(self):
        result = {}
        saddle_result = self.load_result(self.finished_saddle_name)

        if self.finished_reactant_name and self.finished_product_name:
            reactant_result = self.load_result(self.finished_reactant_name)
            result['reactant.con'] = reactant_result['min.con']
            product_result = self.load_result(self.finished_product_name)
            result['product.con'] = product_result['min.con']

        result['saddle.con'] = saddle_result['saddle.con']
        result['mode.dat'] = saddle_result['mode.dat']
        result['results'] = self.data

        results_string = '\n'.join([ "%s %s" % (v,k) for k,v in list(self.data.items()) ])
        result['results'] = self.data
        result['results.dat'] = io.StringIO(results_string)
        result['type'] = self.displacement_type
        result['search_id'] = self.search_id

        return result

    def get_saddle(self):
        if self.finished_saddle_name:
            saddle_result = self.load_result(self.finished_saddle_name)
            saddle = io.loadcon(saddle_result['saddle.con'])
        else:
            saddle = None
        return saddle

    def get_saddle_file(self):
        if self.finished_saddle_name:
            saddle_result = self.load_result(self.finished_saddle_name)
            saddle = saddle_result['saddle.con']
        else:
            saddle = None
        return saddle

    def start_minimization(self, which_min):
        job = {}

        saddle_path = os.path.join(self.config.path_incomplete, self.finished_saddle_name)

        mode_file = open(os.path.join(saddle_path, "mode.dat"))
        mode = io.load_mode(mode_file)
        mode_file.close()

        reactant_file = open(os.path.join(saddle_path, "saddle.con"))
        reactant = io.loadcon(reactant_file)
        reactant_file.close()

        if which_min == "min2":
            mode = -mode

        reactant.r += self.config.process_search_minimization_offset*mode

        reactIO = io.StringIO()
        io.savecon(reactIO, reactant)
        job['pos.con'] = reactIO

        ini_changes = [ ('Main', 'job', 'minimization') ]
        job['config.ini'] = io.modify_config(self.config.config_path, ini_changes)

        return job

    def finish_minimization(self, result):
        result1 = self.load_result(self.finished_min1_name)
        result2 = result

        atoms1 = io.loadcon(result1['min.con'])
        atoms2 = io.loadcon(result2['min.con'])

        results_dat1 = io.parse_results(result1['results.dat'])
        results_dat2 = io.parse_results(result2['results.dat'])
        self.data['force_calls_minimization'] += results_dat1['total_force_calls']
        self.data['force_calls_minimization'] += results_dat2['total_force_calls']

        is_reactant = lambda a: atoms.match(a, self.reactant,
                                            self.config.comp_eps_r,
                                            self.config.comp_neighbor_cutoff, False)

        tc1 = io.parse_results(result1['results.dat'])['termination_reason']
        tc2 = io.parse_results(result2['results.dat'])['termination_reason']

        termination_reason1 = self.job_termination_reasons['minimization'][tc1]
        termination_reason2 = self.job_termination_reasons['minimization'][tc2]
        if termination_reason1 == 'max_iterations' or termination_reason2 == 'max_iterations':
            self.data['termination_reason'] = 9
            self.data['potential_energy_saddle'] = 0.0
            self.data['potential_energy_reactant'] = 0.0
            self.data['potential_energy_product'] = 0.0
            self.data['barrier_reactant_to_product'] = 0.0
            self.data['barrier_product_to_reactant'] = 0.0
            return

        # Check the connectivity of the process
        if (not is_reactant(atoms1) and not is_reactant(atoms2)) or \
           (is_reactant(atoms1) and is_reactant(atoms2)):
            # Not connected
            self.data['termination_reason'] = 6
            self.data['potential_energy_saddle'] = 0.0
            self.data['potential_energy_reactant'] = 0.0
            self.data['potential_energy_product'] = 0.0
            self.data['barrier_reactant_to_product'] = 0.0
            self.data['barrier_product_to_reactant'] = 0.0
            return
        elif is_reactant(atoms1):
            reactant_results_dat = results_dat1
            product_results_dat = results_dat2
            self.finished_reactant_name = self.finished_min1_name
            self.finished_product_name = self.finished_min2_name
        elif is_reactant(atoms2):
            reactant_results_dat = results_dat2
            product_results_dat = results_dat1
            self.finished_reactant_name = self.finished_min2_name
            self.finished_product_name = self.finished_min1_name

        self.data['potential_energy_reactant'] = reactant_results_dat['potential_energy']
        self.data['potential_energy_product'] = product_results_dat['potential_energy']

        self.data['barrier_reactant_to_product'] = self.data['potential_energy_saddle'] - \
                self.data['potential_energy_reactant']
        self.data['barrier_product_to_reactant'] = self.data['potential_energy_saddle'] - \
                self.data['potential_energy_product']

    def start_search(self):
        job = {}

        dispIO = io.StringIO()
        io.savecon(dispIO, self.displacement)
        job['displacement.con'] = dispIO

        modeIO = io.StringIO()
        io.save_mode(modeIO, self.mode)
        job['direction.dat'] = modeIO

        reactIO = io.StringIO()
        io.savecon(reactIO, self.reactant)
        job['pos.con'] = reactIO

        ini_changes = [ ('Main', 'job', 'saddle_search') ]
        job['config.ini'] = io.modify_config(self.config.config_path, ini_changes)

        return job

    def finish_search(self, result):
        results_dat = io.parse_results(result['results.dat'])
        self.data.update(results_dat)
        reactant_energy = results_dat['potential_energy_reactant']
        barrier = results_dat['potential_energy_saddle'] - reactant_energy
        self.data['potential_energy_reactant'] = reactant_energy
        self.data['barrier_reactant_to_product'] = barrier

    def save_result(self, result):
        dir_path = os.path.join(self.config.path_incomplete, result['name'])
        os.makedirs(dir_path)
        for k in result:
            if hasattr(result[k], 'getvalue'):
                fn = os.path.join(dir_path, k)
                f = open(fn, "w")
                f.write(result[k].getvalue())
                f.close()

    def load_result(self, result_name):
        dir_path = os.path.join(self.config.path_incomplete, result_name)
        result = {}
        for file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, file)
            f = open(file_path)
            result[file] = io.StringIO(f.read())
            f.close()
        return result
