#!/bin/sh
# submit_job.sh JOBNAME JOBPATH: submit one eOn job to Slurm, print its id.
#
# Site options come from the environment the eOn server runs in:
#   EON_SBATCH_ARGS  extra sbatch options, e.g.
#                    "-A myaccount -p cpu -N 1 --ntasks=48 -t 01:00:00"
#   EON_CLIENT       the command the job runs (default: eonclient). A
#                    potential that needs several ranks is launched by the
#                    ExtPot wrapper inside the job (srun -n ...), so the
#                    client itself stays one process.
set -eu
jobname=$1
jobpath=$2
# shellcheck disable=SC2086
id=$(sbatch --parsable ${EON_SBATCH_ARGS:-} -J "$jobname" -D "$jobpath" \
    -o "$jobpath/slurm-%j.out" --wrap="${EON_CLIENT:-eonclient}")
# --parsable prints "jobid" or "jobid;cluster".
echo "${id%%;*}"
