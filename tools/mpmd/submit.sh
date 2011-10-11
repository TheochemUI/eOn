#!/bin/bash

set -e
set -x

qsub -V -b y -j y -S /bin/bash -cwd \
     -pe mpi8 `./run.sh -n` \
     -l h_rt=1:00:00 \
     -N IrH \
     -o "\$JOB_NAME.\$JOB_ID" \
     ./run.sh
