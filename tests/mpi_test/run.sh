#!/bin/sh
../../akmc.py --reset -f
mpirun -n 4 ../../tools/gpaw_sp.py : \
       -n 1 ../../akmc.py : \
       -n 2 ../../client/client
