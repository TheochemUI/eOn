#!/usr/bin/env python

import os
import glob
from io import StringIO

from eon import akmc
from eon import basinhopping
from eon import parallelreplica
from eon import escaperate
from eon import fileio as io
from eon.config import config

def server():
    config.init()

    # Should we have some kind of sanity-check module/function somewhere?
    fnames = [os.path.basename(f) for f in glob.glob(os.path.join(config.path_pot, '*'))]
    if 'pos.con' in fnames:
        print("WARNING: pos.con found in potfiles path. Are you sure you want this? It will overwrite the pos.con in the calculation directory when your jobs are being run.")

    job = config.main_job.lower()
    if job == 'akmc':
        akmc.main()
    elif job == 'parallel_replica' or job == 'unbiased_parallel_replica':
        parallelreplica.main()
    elif job == 'basin_hopping':
        basinhopping.main()
    elif job == 'escape_rate':
        escaperate.main()
    else:
        from eon import communicator
        import shutil
        config.path_scratch = config.path_root
        comm = communicator.get_communicator()

        invariants = {}

        # Merge potential files into invariants
        invariants = dict(invariants, **io.load_potfiles(config.path_pot))

        job = {}
        files = [ f for f in os.listdir(".") if os.path.isfile(f) ]
        for f in files:
            fh = open(f)
            if(len(f.split('.')) > 1):
                #f_passed = f.split('.')[0] + "_passed." + f.split('.')[1]
                f_passed = f.split('.')[0] + "." + f.split('.')[1]
                job[f_passed] = StringIO(fh.read())
            fh.close()
        job["id"] = "output"
        if os.path.isdir("output_old"):
            shutil.rmtree("output_old")
        if os.path.isdir("output"):
            shutil.move("output", "output_old")
        comm.submit_jobs([job], invariants)

if __name__ == '__main__':
    server()
