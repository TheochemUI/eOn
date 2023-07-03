#!/usr/bin/env python

from io import StringIO
import logging
import logging.handlers
logger = logging.getLogger('basinhopping')
import numpy
numpy.seterr(divide="raise", over="raise", under="print", invalid="raise")
import optparse
import os
import shutil
import sys
import random

from eon import atoms
from eon import communicator
from eon.config import config
from eon import fileio as io
from eon import locking
from eon.version import version

#class RandomStructure:
#    def __init__(self, structure):
#        self.structure = structure
#        self.radii = numpy.array([ atoms.elements[name]['radius']
#                                   for name in structure.names ])
#        self.box_center = numpy.diagonal(structure.box)/2.0
#        self.p = 0.1
#
#    def generate(self):
#        indexes = range(len(self.structure))
#        random.shuffle(indexes)
#
#        rs = atoms.Atoms(0)
#        first = indexes[0]
#        rs.append(numpy.zeros(3), True, self.structure.names[first],
#                  self.structure.mass[first])
#
#        failures = 0
#        for i in indexes[1:]:
#            rs.append(numpy.zeros(3), True,
#                      self.structure.names[i], self.structure.mass[i])
#            rs.box = self.box_p(rs)
#            valid = False
#            while not valid:
#                rs.r[-1] = numpy.random.uniform(-rs.box[0][0]/2.0,rs.box[0][0]/2.0,3)
#                distances = numpy.zeros(len(rs)-1)
#                for j in range(len(rs)-1):
#                    distances[j] = numpy.linalg.norm(rs.r[-1]-rs.r[j])
#                    bond_length = self.radii[i]+self.radii[j]
#
#                if min(distances) < 0.8*bond_length:
#                    failures += 1
#                    valid = False
#                elif min(distances) > 1.3*bond_length:
#                    failures += 1
#                    valid = False
#                else:
#                    #print '%i/%i' % (i+1,len(self.structure))
#                    valid = True
#        rs.box = self.structure.box
#        #XXX: shift to center of box, only correct for cubic cells
#        rs.r += rs.box.diagonal()/2.0
#        return rs
#
#    def box_p(self, rs):
#        V_atoms = sum(4./3*3.14159*self.radii[0:len(rs)])
#        a = (V_atoms/self.p)**(1/3.)
#        return numpy.array( ((a,0,0),(0,a,0),(0,0,a)) )

class BHStates:
    def __init__(self):
        if not os.path.isdir(config.path_states):
            os.mkdir(config.path_states)

        state_table_path = os.path.join(config.path_states, 'state_table')
        self.energy_table = io.Table(state_table_path, ['state', 'energy', 'repeats'])

    def get_random_minimum(self):
        if len(self.energy_table.rows) == 0: return None
        N = config.bh_initial_state_pool_size
        self.energy_table.rows.sort(key=lambda r:-r['energy'])
        lowest_N = self.energy_table.rows[:N]
        i = random.choice(lowest_N)['state']
        f = open(os.path.join(config.path_states, str(i), 'minimum.con'))
        return StringIO(f.read())

    def add_state(self, result_files, result_info):
        energy = result_info['minimum_energy']

        energetically_close = []
        added = True

        if len(self.energy_table) != 0:
            for row in self.energy_table:
                if abs(energy-row['energy']) < config.comp_eps_e:
                    energetically_close.append(row['state'])
            if len(energetically_close) != 0:
                a1 = io.loadcon(result_files['min.con'])
                for state_number in energetically_close:
                    state_con_path = os.path.join(config.path_states,
                                                  str(state_number),
                                                  'minimum.con')
                    a2 = io.loadcon(state_con_path)
                    if atoms.match(a1, a2, config.comp_eps_r, config.comp_neighbor_cutoff, True):
                        logger.info("Found a repeat of state %i", state_number)
                        added = False
                        for row in self.energy_table.rows:
                            if row['state'] == state_number:
                                row['repeats'] += 1
                                self.energy_table.write()
                                break

        if added:
            state_number = len(self.energy_table)

            row = { 'state':state_number, 'energy':energy, 'repeats':0 }
            self.energy_table.add_row(row)
            self.energy_table.rows.sort(key=lambda r:-r['energy'])
            self.energy_table.write()

            state_path = os.path.join(config.path_states, str(state_number))
            os.mkdir(state_path)

            result_files['minimum.con'] = result_files['min.con']
            del result_files['min.con']

            for fn, fh in list(result_files.items()):
                if hasattr(fh, 'getvalue') == False:
                    continue
                p = os.path.join(state_path, fn)
                f = open(p, 'w')
                f.write(fh.getvalue())
                f.close()

        return added

def basinhopping():
    logger.info('Eon version %s', version())
    # First of all, does the root directory even exist?
    if not os.path.isdir(config.path_root):
        logger.critical("Root directory does not exist")
        sys.exit(1)

    # load metadata
    bhstates = BHStates()

    if os.path.isfile("wuid.dat"):
        wuid_file = open("wuid.dat")
        wuid = int(wuid_file.readline().strip())
        wuid_file.close()
    else:
        wuid = 0

    # get communicator
    comm = communicator.get_communicator()

    # Register all the results. There is  no need to ever discard found
    # processes like we do with akmc. There is no confidence to calculate.
    register_results(comm, bhstates)

    wuid = make_searches(comm, wuid, bhstates)

    wuid_file = open("wuid.dat","w")
    wuid_file.write("%i\n" % wuid)
    wuid_file.close()

    io.save_prng_state()

def make_searches(comm, wuid, bhstates):
    num_in_buffer = comm.get_queue_size()*config.comm_job_bundle_size
    logger.info("%i searches in the queue" % num_in_buffer)
    num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
    logger.info("Making %i searches" % num_to_make)

    if num_to_make == 0:
        return wuid

    searches = []

    invariants = {}

    f = open(os.path.join(config.path_root, 'pos.con'))
    initial_react = StringIO(f.read())
    f.close()

    #invariants['reactant_passed.con']=reactIO

    ini_changes = [ ('Main', 'job', 'basin_hopping') ]
    #invariants['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)
    #invariants['reactant_passed.con']  = reactIO

    #Merge potential files into invariants
    invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

    searches = []
    for i in range(num_to_make):
        search = {}
        search['id'] = "%d" % wuid
        ini_changes = [ ('Main', 'random_seed', str(int(numpy.random.random()*2**32))) ]

        #if config.bh_random_structure:
        #    reactIO = StringIO()
        #    io.savecon(reactIO, rs.generate())
        #else:
        #    reactIO = initial_react

        if config.bh_initial_state_pool_size == 0:
            reactIO = initial_react
        elif config.bh_initial_state_pool_size > 0:
            reactIO = bhstates.get_random_minimum()
            if reactIO is None:
                reactIO = initial_react
        else:
            logger.fatal("Initial state pool size negative")
            sys.exit(1)

        search['pos.con'] = reactIO
        search['config.ini'] = io.modify_config(config.config_path, ini_changes)
        searches.append(search)
        wuid += 1

    comm.submit_jobs(searches, invariants)
    logger.info( str(num_to_make) + " searches created")
    return wuid

def register_results(comm, bhstates):
    logger.info("Registering results")
    if os.path.isdir(config.path_jobs_in):
        shutil.rmtree(config.path_jobs_in)
    os.makedirs(config.path_jobs_in)

    # Function used by communicator to determine whether to discard a result
    def keep_result(name):
        return True

    num_registered = 0

    for result in comm.get_results(config.path_jobs_in, keep_result):
        # The result dictionary contains the following key-value pairs:
        # product.con - an array of strings containing the reactant
        # results.dat - an array of strings containing the results
        # id - wuid

        #result_id = result['id']
        #del result['id']
        result_info = io.parse_results(result['results.dat'])
        if 'minimum_energy' not in result_info:
            continue
        if result_info['termination_reason'] == 0:
            if bhstates.add_state(result, result_info):
                logger.info("New structure with energy %.8e",
                            result_info['minimum_energy'])

            #logger.info("found new structure with energy %.3e", fe)
            #if bhminima.add_minimum(fe, result['product.con'].getvalue()):
            #    logger.info("found new low energy structure with energy %.3e", fe)

        num_registered += 1

    logger.info("%i (result) searches processed", num_registered)

def main():
    optpar = optparse.OptionParser(usage="usage: %prog [options] config.ini")
    optpar.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False,help="only write to the log file")
    optpar.add_option("-n", "--no-submit", action="store_true", dest="no_submit", default=False,help="don't submit searches; only register finished results")
    optpar.add_option("-R", "--reset", action="store_true", dest="reset", default = False, help="reset the simulation, discarding all data")
    options, args = optpar.parse_args()

    if len(args) > 1:
        print("takes only one positional argument")
    sys.argv = sys.argv[0:1]
    if len(args) == 1:
        sys.argv += args
    if len(sys.argv) > 1:
        config.init(sys.argv[-1])
    else:
        config.init()
    # set options.path_root to be where the config file is if given as an arg
    if config.path_root.strip() == '.' and len(args) == 1:
        config.path_root = os.path.abspath(os.path.dirname(args[0]))
        os.chdir(config.path_root)

    if config.comm_job_bundle_size != 1:
        print("error: Basin Hopping only supports a bundle size of 1")
        sys.exit(1)

    if options.no_submit:
        config.comm_job_buffer_size = 0

    if options.reset:
        res = input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
                rmdirs = [config.path_jobs_out, config.path_jobs_in, config.path_scratch,  config.path_states]
                for i in rmdirs:
                    if os.path.isdir(i):
                        shutil.rmtree(i)
                        #XXX: ugly way to remove all empty directories containing this one
                        os.mkdir(i)
                        os.removedirs(i)
                log_path = os.path.join(config.path_results, "bh.log")
                wuid_path = os.path.join(config.path_results, "wuid.dat")
                prng_path = os.path.join(config.path_results, "prng.pkl")
                for i in [log_path, wuid_path, prng_path ]:
                    if os.path.isfile(i):
                        os.remove(i)
                print("Reset.")
                sys.exit(0)
        else:
            print("Not resetting.")
            sys.exit(1)

    # setup logging
    logging.basicConfig(level=logging.DEBUG,
            filename=os.path.join(config.path_results, "bh.log"),
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

    lock = locking.LockFile(os.path.join(config.path_results, "lockfile"))

    if lock.aquirelock():
        if config.comm_type == 'mpi':
            from eon.mpiwait import mpiwait
            while True:
                mpiwait()
                basinhopping()
        basinhopping()
    else:
        logger.warning("Couldn't get lock")
        sys.exit(1)

if __name__ == '__main__':
    main()
