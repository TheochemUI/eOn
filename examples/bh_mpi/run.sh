#!/bin/bash
export EON_NUMBER_OF_CLIENTS=2
export EON_SERVER_PATH=~/code/eon/eon/basinhopping.py
mpirun -n 3 ~/code/eon/client/eonclientmpi
