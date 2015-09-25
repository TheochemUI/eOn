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

def insert(state, process_id):
    try:
        import eon.kdb as kdb
    except:
        logger.error('Python module kdb not found, kdb will not be used.')
        return
    logger.debug("KDB inserting process")
    reactant = kdb.aselite.read_any(state.proc_reactant_path(process_id))
    saddle   = kdb.aselite.read_any(state.proc_saddle_path(process_id))
    product  = kdb.aselite.read_any(state.proc_product_path(process_id))
    mode     = state.get_process_mode(process_id)
    db = kdb.local_db.LocalDB(config.kdb_name)
    params = db.get_params()
    insert_sub_class = kdb.local_insert.LocalInsert()
    insert_sub_class.insert(reactant, saddle, product, mode=mode, nf=params['nf'], 
               dc=params['dc'], mac=params['mac'],kdbname=config.kdb_name)

def query(state):
    try:
        import eon.kdb as kdb
    except:
        logger.error('Python module kdb not found, kdb will not be used.')
        return
    if os.path.isdir(os.path.join(config.kdb_scratch_path, "kdbmatches")):
        shutil.rmtree(os.path.join(config.kdb_scratch_path, "kdbmatches"))
    os.makedirs(os.path.join(config.kdb_scratch_path, "kdbmatches"))
    kdbpath = os.path.abspath(os.path.join(config.path_root, config.kdb_path))
    if state.number == 0:
        reactant = kdb.aselite.read_any(os.path.abspath(os.path.join(config.path_root, "pos.con")))
    else:
        reactant = kdb.aselite.read_any(os.path.abspath(state.reactant_path))
    db = kdb.local_db.LocalDB(config.kdb_name)
    params = db.get_params()
    query_sub_class = kdb.local_query.LocalQuery()
    query_sub_class.query(reactant, os.path.join(config.kdb_scratch_path, "kdbmatches"), 
                          nodupes = config.kdb_nodupes, kdbname=config.kdb_name, 
                          dc=params['dc'], nf=params['nf'])

def make_suggestion():
    try:
        import eon.kdb as kdb
    except:
        logger.error('Python module kdb not found, kdb will not be used.')
        return None, None
    if os.path.isdir(os.path.join(config.kdb_scratch_path, "kdbmatches")):
        dones = glob.glob(os.path.join(config.kdb_scratch_path, "kdbmatches",".done_*"))
        if len(dones) > 0:
            number = dones[0].split("_")[1]
            try:
                displacement = io.loadcon(os.path.join(config.kdb_scratch_path, "kdbmatches", "SADDLE_%s" % number))
            except FloatingPointError:
                displacement = io.loadposcar(os.path.join(config.kdb_scratch_path, "kdbmatches", "SADDLE_%s" % number))
            except ValueError:
                displacement = io.loadposcar(os.path.join(config.kdb_scratch_path, "kdbmatches", "SADDLE_%s" % number))
            mode = [[float(i) for i in l.strip().split()] for l in
                    open(os.path.join(config.kdb_scratch_path, "kdbmatches",
                    "MODE_%s" % number), 'r').readlines()[:]]
            os.remove(os.path.join(config.kdb_scratch_path, "kdbmatches", ".done_%s" % number))
            os.remove(os.path.join(config.kdb_scratch_path, "kdbmatches", "SADDLE_%s" % number))
            os.remove(os.path.join(config.kdb_scratch_path, "kdbmatches", "MODE_%s" % number))
            return displacement, mode
    return None, None

