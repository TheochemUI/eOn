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

from cStringIO import StringIO
import logging
import logging.handlers
logger = logging.getLogger('basinhopping')
import numpy
numpy.seterr(all='raise')
import optparse
import os
import shutil
import sys

import atoms
import communicator
import config
import io
import locking

class BHMinima:
    def __init__(self):
        if not os.path.isdir(config.path_bh_minima):
            os.mkdir(config.path_bh_minima)

        self.minima = []

        for fn in os.listdir(config.path_bh_minima):
            if 'dat' in fn:
                energy = float(open(os.path.join(config.path_bh_minima,fn)).readline())
                number = int(fn.split('_')[-1].split('.')[0])
                structure_file = fn[:-3]+'con'
                structure_data = open(os.path.join(config.path_bh_minima,structure_file)).read()
                self.minima.append({'number':number, 'energy':energy, 'structure':structure_data})

    def get_energies(self):
        energies = [ m['energy'] for m in self.minima ]
        return numpy.array(energies)

    def add_minimum(self, energy, structure_data):
        if len(self.minima) == 0:
            max_energy = 0
        else:
            max_energy = self.minima[-1]['energy']
        added = False
        candidate = {
            'number':len(self.minima), 
            'energy':energy, 
            'structure':structure_data
        }

        if len(self.minima) < config.bh_number_of_minima:
            added = True
        elif energy < max_energy:
            ediff = self.get_energies()-energy

            if ediff.min() > config.comp_eps_e:
                self.minima[ediff.argmin()] = candidate
                added = True

        if added:
            self.minima.append(candidate)
            self.minima.sort(cmp=lambda a,b: cmp(a['energy'],b['energy']))
            self.minima = self.minima[:config.bh_number_of_minima]
            for i,m in enumerate(self.minima):
                m['number'] = i
            self.write_all()
        return added

    def write_all(self):
        prefix = "minimum"
        for m in self.minima:
            dat_path = os.path.join(config.path_bh_minima, "%s_%i.dat" % (prefix, m['number']))
            con_path = os.path.join(config.path_bh_minima, "%s_%i.con" % (prefix, m['number']))

            f = open(dat_path, 'w')
            f.write("%.5f\n"%m['energy'])
            f.close()

            f = open(con_path, 'w')
            f.write(m['structure'])
            f.close()

def basinhopping():
    # First of all, does the root directory even exist?
    if not os.path.isdir(config.path_root):
        logger.critical("Root directory does not exist")
        sys.exit(1)

    # load metadata
    bhminima = BHMinima()

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
    register_results(comm, bhminima)

    wuid = make_searches(comm, wuid, bhminima)

    wuid_file = open("wuid.dat","w")
    wuid_file.write("%i\n" % wuid)
    wuid_file.close()
    
def make_searches(comm, wuid, bhminima):
    num_in_buffer = comm.get_queue_size()*config.comm_job_bundle_size
    logger.info("%i searches in the queue" % num_in_buffer)
    num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
    logger.info("making %i searches" % num_to_make)
    
    if num_to_make == 0:
        return wuid
    
    searches = []
    
    invariants = {}

    f = open(os.path.join(config.path_root, 'reactant.con'))
    initial_react = StringIO(f.read())
    f.close()

    #invariants['reactant_passed.con']=reactIO
    
    ini_changes = [ ('Main', 'job', 'basin_hopping') ]
    #invariants['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)
    #invariants['reactant_passed.con']  = reactIO

    #Merge potential files into invariants
    invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

    searches = []
    number_random = 0
    number_minima = 0
    for i in range(num_to_make):
        search = {}
        search['id'] = "%d" % wuid
        ini_changes = [ ('Main', 'random_seed', str(int(numpy.random.random()*10**9))) ]
        reactIO = StringIO()
        if len(bhminima.minima) == 0 or \
           numpy.random.random() < config.bh_md_probability:
            number_random += 1
            reactIO = initial_react
            ini_changes.append( ('Basin Hopping', 'md_first', 'true') )
        else:
            number_minima += 1
            i = numpy.random.randint(0,len(bhminima.minima)-1)
            reactIO = StringIO(bhminima.minima[i]['structure'])
            ini_changes.append( ('Basin Hopping', 'md_first', 'false') )
        search['reactant_passed.con'] = reactIO
        search['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)
        searches.append(search)
        wuid += 1

    logger.info("%i from random structures %i from previous minima", number_random, number_minima)
    comm.submit_jobs(searches, invariants)
    logger.info( str(num_to_make) + " searches created") 
    return wuid

def register_results(comm, bhminima):
    logger.info("registering results")
    if os.path.isdir(config.path_jobs_in):
        shutil.rmtree(config.path_jobs_in)    
    os.makedirs(config.path_jobs_in)
    
    #Function used by communicator to determine whether to discard a result
    def keep_result(name):
        return True

    num_registered = 0

    for result in comm.get_results(config.path_jobs_in, keep_result): 
        # The result dictionary contains the following key-value pairs:
        # product.con - an array of strings containing the reactant
        # results.dat - an array of strings containing the results
        # id - wuid

        #id = result['name']

        #read in the results
        result['results'] = io.parse_results(result['results.dat'])
        if result['results']['termination_reason'] == 0:
            fe = result['results']['minimum_energy']
            #logger.info("found new structure with energy %.3e", fe)
            if bhminima.add_minimum(fe, result['product.con'].getvalue()):
                logger.info("found new low energy structure with energy %.3e", fe)
        num_registered += 1
        
    logger.info("%i (result) searches processed", num_registered)

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
        print "error: Basin Hopping only supports a bundle size of 1"
        sys.exit(1)

    if options.no_submit:
        config.comm_job_buffer_size= 0

    if options.reset:
        res = raw_input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
                rmdirs = [config.path_jobs_out, config.path_jobs_in, config.path_scratch,  config.path_bh_minima]
                for i in rmdirs:
                    if os.path.isdir(i):
                        shutil.rmtree(i)
                        #XXX: ugly way to remove all empty directories containing this one
                        os.mkdir(i)
                        os.removedirs(i)
                log_path = os.path.join(config.path_results, "bh.log") 
                wuid_path = os.path.join(config.path_results, "wuid.dat") 
                min_path = os.path.join(config.path_results, "min_energy.con") 
                mine_path = os.path.join(config.path_results, "min_energy.dat") 
                for i in [log_path, wuid_path, min_path, mine_path]:
                    if os.path.isfile(i):
                        os.remove(i)
                print "Reset."
                sys.exit(0)
        else:
            print "Not resetting."
            sys.exit(1)


    #setup logging
    logging.basicConfig(level=logging.DEBUG,
            filename=os.path.join(config.path_results, "bh.log"),
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
        basinhopping()
    else:
        logger.warning("couldn't get lock")
        sys.exit(1)

if __name__ == '__main__':
    main()
