
import subprocess
import os
import shutil
import signal
import glob
import numpy
import logging
logger = logging.getLogger('kdb')

from eon.config import config
from eon import fileio as io

def insert(state, process_id):
    try:
        from kdb import aselite
        from kdb import local_db
        from kdb import local_insert
    except:
        logger.error('Python module kdb not found, kdb will not be used.')
        return
    logger.debug("KDB inserting process")
    reactant = aselite.read_any(state.proc_reactant_path(process_id))
    saddle   = aselite.read_any(state.proc_saddle_path(process_id))
    product  = aselite.read_any(state.proc_product_path(process_id))
    mode     = state.get_process_mode(process_id)
    db = local_db.LocalDB(config.kdb_name)
    params = db.get_params()
    insert_sub_class = local_insert.LocalInsert()
    insert_sub_class.insert(reactant, saddle, product, mode=mode, nf=params['nf'],
               dc=params['dc'], mac=params['mac'], kdbname=config.kdb_name)

def query(state):
    try:
        from kdb import aselite
        from kdb import local_db
        from kdb import local_query
    except:
        logger.error('Python module kdb not found, kdb will not be used.')
        return
    if os.path.isdir(os.path.join(config.kdb_scratch_path, "kdbmatches")):
        shutil.rmtree(os.path.join(config.kdb_scratch_path, "kdbmatches"))
    os.makedirs(os.path.join(config.kdb_scratch_path, "kdbmatches"))
    kdbpath = os.path.abspath(os.path.join(config.path_root, config.kdb_path))
    if state.number == 0:
        reactant = aselite.read_any(os.path.abspath(os.path.join(config.path_root, "pos.con")))
    else:
        reactant = aselite.read_any(os.path.abspath(state.reactant_path))
    db = local_db.LocalDB(config.kdb_name)
    params = db.get_params()
    query_sub_class = local_query.LocalQuery()
    query_sub_class.query(reactant, os.path.join(config.kdb_scratch_path, "kdbmatches"),
                          nodupes = config.kdb_nodupes, kdbname=config.kdb_name,
                          dc=params['dc'], nf=params['nf'])

def make_suggestion():
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
