#!/bin/sh
# queued_jobs.sh: print the ids of this user's pending and running jobs,
# one per line. A failing squeue must fail this script: eOn would otherwise
# read an empty queue as "every job finished" and harvest running jobs.
set -eu
squeue -h -u "${USER:-$(id -un)}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED -o %i
