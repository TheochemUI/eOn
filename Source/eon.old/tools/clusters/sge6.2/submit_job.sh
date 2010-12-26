#!/bin/sh
qsub -V -b y -j y -o ll_out -S /bin/bash -N "$1" -wd "$2" ~/bin/eonclient | awk '{print $3}'
