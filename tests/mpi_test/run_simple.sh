#!/bin/sh
../../akmc.py --reset -f
mpirun -n 1 ../../akmc.py : \
       -n 2 ../../client/client
