#!/bin/bash
prefix=/home/chill/cross
setups=${prefix}/share/gpaw-setups-0.8.7929
pythonpath=${prefix}/lib/python2.6/site-packages
gpaw_python=${prefix}/bin/gpaw-python
lib_dir=${prefix}/lib

rm lockfile
cobalt-mpiexec -verbose 2 -label -mode vn \
    : -n 256 -wdir . -env PYTHONPATH=$pythonpath \
                     -env LD_LIBRARY_PATH=$lib_dir \
                     -env GPAW_SETUP_PATH=$setups \
                     $gpaw_python gpaw_sp.py 4 \
    : -n 256 -wdir . bgpeonclient ~/eon/akmc.py 4
