#!/usr/bin/env python
import config

import akmc
import basinhopping
import parallelreplica

def main():
    config.init()
    job = config.main_job.lower()
    if job == 'akmc':
        akmc.main()
    elif job == 'parallel_replica':
        parallelreplica.main()
    elif job == 'basin_hopping':
        basinhopping.main()
    else:
        #TODO: Work on running the client directly for unknown job types
        import communicator
        comm = communicator.Local(".", eonclient, 1, 1)

if __name__ == '__main__':
    main()
