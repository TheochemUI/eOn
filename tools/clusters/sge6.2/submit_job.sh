#!/bin/sh
qsub -V -b y -cwd -j y -o ll_out -S /bin/bash \
    -N "akmc_$1" \
    -wd "$2" \
     ~/bin/eonclient | awk '{print $3}'
