import logging
import shutil
from pathlib import Path

import numpy

logger = logging.getLogger('kdb')

from eon import fileio as io

def insert(state, process_id, config):
    try:
        from kdb import aselite
        from kdb import local_db
        from kdb import local_insert
    except Exception:
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

def query(state, config):
    try:
        from kdb import aselite
        from kdb import local_db
        from kdb import local_query
    except Exception:
        logger.error('Python module kdb not found, kdb will not be used.')
        return
    matches = Path(config.kdb_scratch_path) / "kdbmatches"
    if matches.is_dir():
        shutil.rmtree(matches)
    matches.mkdir(parents=True)
    if state.number == 0:
        reactant = aselite.read_any(str((Path(config.path_root) / "pos.con").resolve()))
    else:
        reactant = aselite.read_any(str(Path(state.reactant_path).resolve()))
    db = local_db.LocalDB(config.kdb_name)
    params = db.get_params()
    query_sub_class = local_query.LocalQuery()
    query_sub_class.query(reactant, str(matches),
                          nodupes = config.kdb_nodupes, kdbname=config.kdb_name,
                          dc=params['dc'], nf=params['nf'])

def make_suggestion(config):
    matches = Path(config.kdb_scratch_path) / "kdbmatches"
    if matches.is_dir():
        dones = sorted(matches.glob(".done_*"))
        if dones:
            number = dones[0].name.split("_")[1]
            saddle = matches / ("SADDLE_%s" % number)
            try:
                displacement = io.loadcon(str(saddle))
            except OSError:
                # readcon reports a parse failure as an OSError subclass; the
                # match may be a POSCAR instead.
                displacement = io.loadposcar(str(saddle))
            mode = io.load_mode(str(matches / ("MODE_%s" % number)))
            dones[0].unlink()
            saddle.unlink()
            (matches / ("MODE_%s" % number)).unlink()
            return displacement, mode
    return None, None
