#!/bin/sh
../../akmc.py --reset -f
mpirun -n 1 ../../tools/emt_sp.py : \
       -n 1 ../../akmc.py : \
       -n 1 ../../client/client
