#!/bin/sh
qsub -V -b y -j y -o ll_out -S /bin/bash -N "$1" -wd "$2" ~/code/eon/client/client | awk '{print $3}'
