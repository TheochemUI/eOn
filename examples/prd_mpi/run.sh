#!/bin/bash
export EON_NUMBER_OF_CLIENTS=7
export EON_SERVER_PATH=~/code/eon/eon/parallelreplica.py
mpirun -n 8 ~/code/eon/client/eonclientmpi
