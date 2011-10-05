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
logger = logging.getLogger('Distributedreplica')
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

def distributedreplica():
    # First of all, does the root directory even exist?
    if not os.path.isdir(config.path_root):
        logger.critical("Root directory does not exist")
        sys.exit(1)

    # load metadata
    min_energy_dat_path = "min_energy.dat"
    min_energy_con_path = "min_energy.con"
    if os.path.isfile(min_energy_dat_path):
        min_energy_file = open(min_energy_dat_path, "r")
        min_energy = float(min_energy_file.readline().strip())
        min_energy_file.close()
    else:
        min_energy = 1e100

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
    num_registered, new_min_energy_structure = register_results(comm, min_energy)

    if new_min_energy_structure:
        min_energy_file = open(min_energy_dat_path, "w")
        min_energy_file.write("%.3e\n"%new_min_energy_structure[0])
        min_energy_file.close()
        min_energy_file = open(min_energy_con_path, "w")
        min_energy_file.write("".join(new_min_energy_structure[1]))
        min_energy_file.close()
    elif not os.path.isfile(min_energy_con_path):
        #this should only happen on the first run
        shutil.copy("reactant.con", min_energy_con_path)

    wuid = make_searches(comm, wuid, min_energy_con_path)

    wuid_file = open("wuid.dat","w")
    wuid_file.write("%i\n" % wuid)
    wuid_file.close()

def make_searches(comm, wuid, min_energy_con_path):
    num_in_buffer = comm.get_queue_size()*config.comm_job_bundle_size
    logger.info("%i searches in the queue" % num_in_buffer)
    num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
    logger.info("making %i searches" % num_to_make)

    if num_to_make == 0:
        return wuid

    searches = []

    invariants = {}

    f = open(min_energy_con_path)
    reactIO = StringIO(''.join(f.readlines()))
    f.close()
    #invariants['reactant_passed.con']=reactIO

    ini_changes = [ ('Main', 'job', 'distributed_replica') ]
    invariants['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)
    invariants['reactant_passed.con']  = reactIO

    #Merge potential files into invariants
    invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

    searches = []
    for i in range(num_to_make):
        search = {}
        search['id'] = "%d" % wuid
        param_ini_str = "[Main]\nrandom_seed=%i\n" % int(numpy.random.random()*10**9)
        paramIO = StringIO(param_ini_str)
        search['config_passed.ini'] = paramIO
        searches.append(search)
        wuid += 1

    comm.submit_jobs(searches, invariants)
    logger.info( str(num_to_make) + " searches created") 
    return wuid

def register_results(comm, min_energy):
    logger.info("registering results")
    if os.path.isdir(config.path_jobs_in):
        shutil.rmtree(config.path_jobs_in)
    os.makedirs(config.path_jobs_in)

    #Function used by communicator to determine whether to discard a result
    def keep_result(name):
        return True

    new_min_energy_structure = None
    num_registered = 0

    for result in comm.get_results(config.path_jobs_in, keep_result): 
        # The result dictionary contains the following key-value pairs:
        # product.con - an array of strings containing the reactant
        # results.dat - an array of strings containing the results
        # id - wuid

        id = result['name']

        #read in the results
        result['results'] = io.parse_results(result['results.dat'])
##        if result['results']['termination_reason'] == 0:
        fe = result['results']['minimum_energy']
        logger.info("found new structure with energy %.3e", fe)
        if fe < min_energy:
              new_min_energy_structure = ( fe, result['product.con'] )
              min_energy = fe
        num_registered += 1

    logger.info("%i (result) searches processed", num_registered)

    return num_registered, new_min_energy_structure

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
        config.comm_job_buffer_size = 0

    if options.reset:
        res = raw_input("Are you sure you want to reset (all data files will be lost)? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
                rmdirs = [config.path_jobs_out, config.path_jobs_in, config.path_scratch]
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
        distributedreplica()
    else:
        logger.warning("couldn't get lock")
        sys.exit(1)

if __name__ == '__main__':
    main()
