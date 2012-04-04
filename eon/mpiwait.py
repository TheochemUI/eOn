import numpy
import os
import logging
import logging.handlers
logger = logging.getLogger('mpiwait')
import config

def mpiwait():
    from mpi4py import MPI
    from array import array
    stopcar_path = os.path.join(config.path_root, "STOPCAR")
    if os.path.isfile(stopcar_path):
        os.unlink(stopcar_path)
    client_ranks = [ int(r) for r in os.environ['EON_CLIENT_RANKS'].split(":") ]

    if os.path.isfile(stopcar_path):
        logging.info("stopping due to STOPCAR")
        for i in range(MPI.COMM_WORLD.Get_size()):
            buf = array('c', 'STOPCAR\0')
            MPI.COMM_WORLD.Isend(buf, i)
        MPI.COMM_WORLD.Abort()

    MPI.COMM_WORLD.Probe(MPI.ANY_SOURCE, MPI.ANY_TAG)
    #we now need to clear out any other mpi_send to us
    for r in client_ranks:
        if MPI.COMM_WORLD.Iprobe(source=r, tag=1):
            tmp = numpy.empty(1, dtype='i')
            MPI.COMM_WORLD.Recv(tmp, source=r, tag=1)
