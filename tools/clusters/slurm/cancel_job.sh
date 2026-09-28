#!/bin/sh
# cancel_job.sh JOBID: cancel one eOn job.
set -eu
scancel "$1"
