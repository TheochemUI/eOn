##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

import logging
import logging.handlers
logger = logging.getLogger('mpiwait')
import config
from time import sleep
from sys import exit
import signal

QUIT = False
def signal_handler(signum, frame):
    global QUIT
    QUIT = True
    logger.info('Signal handler caught signal %i', signum)


def mpiwait():
    from mpi4py import MPI

    signal.signal(signal.SIGINT,  signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGUSR1, signal_handler)
    signal.signal(signal.SIGUSR2, signal_handler)

    while True:
        sleep(config.mpi_poll_period)
        if QUIT:
            exit(0)

        if MPI.COMM_WORLD.Iprobe(MPI.ANY_SOURCE, MPI.ANY_TAG):
            break

