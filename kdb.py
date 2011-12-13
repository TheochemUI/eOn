##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

import subprocess
import os
import shutil
import signal
import glob
import numpy
import logging
logger = logging.getLogger('kdb')    

import config
import fileio as io

try:
    from tsase import kdb
    import tsase
except:
    logger.debug('tsase module not found, kdb unavailable')

def insert(state, process_id):
    try:
        from tsase import kdb
        import tsase
    except:
        logger.error('tsase module not found, kdb unavailable')
        return
    logger.debug("Inserting process.")
    reactant = tsase.io.read_con(state.proc_reactant_path(process_id))
    saddle = tsase.io.read_con(state.proc_saddle_path(process_id))
    product = tsase.io.read_con(state.proc_product_path(process_id))
    mode = state.get_process_mode(process_id)
    process = state.get_process(process_id)
    barrier1 = process["saddle_energy"] - state.get_energy()
    barrier2 = process["saddle_energy"] - process["product_energy"]
    prefactor1 = process["prefactor"] 
    prefactor2 = process["product_prefactor"] 
    kdb.insert(reactant, saddle, product, mode, config.kdb_path)

def query(state):
    try:
        from tsase import kdb
        import tsase
    except:
        logger.error('tsase module not found, kdb unavailable')
        return
    if os.path.isdir(os.path.join(config.kdb_scratch_path, "kdbmatches")):
        shutil.rmtree(os.path.join(config.kdb_scratch_path, "kdbmatches"))        
    os.makedirs(os.path.join(config.kdb_scratch_path, "kdbmatches"))
    kdbpath = os.path.abspath(os.path.join(config.path_root, config.kdb_path))
    if state.number == 0:
        reactant = tsase.io.read_con(os.path.abspath(os.path.join(config.path_root, "reactant.con")))
    else:
        reactant = tsase.io.read_con(os.path.abspath(state.reactant_path))
    kdb.query(reactant, kdbpath, os.path.join(config.kdb_scratch_path, "kdbmatches"), nodupes = config.kdb_nodupes)

def make_suggestion():
    try:
        from tsase import kdb
        import tsase
    except:
        logger.error('tsase module not found, kdb unavailable')
        return None, None
    if os.path.isdir(os.path.join(config.kdb_scratch_path, "kdbmatches")):
        dones = glob.glob(os.path.join(config.kdb_scratch_path, "kdbmatches",".done_*"))
        if len(dones) > 0:
            number = dones[0].split("_")[1]
            displacement = io.loadcon(os.path.join(config.kdb_scratch_path, "kdbmatches", "SADDLE_%s" % number))
            mode = [[float(i) for i in l.strip().split()] for l in
                    open(os.path.join(config.kdb_scratch_path, "kdbmatches",
                    "MODE_%s" % number), 'r').readlines()[:]]
            os.remove(os.path.join(config.kdb_scratch_path, "kdbmatches", ".done_%s" % number))
            os.remove(os.path.join(config.kdb_scratch_path, "kdbmatches", "SADDLE_%s" % number))
            os.remove(os.path.join(config.kdb_scratch_path, "kdbmatches", "MODE_%s" % number))
            return displacement, mode
    return None, None


