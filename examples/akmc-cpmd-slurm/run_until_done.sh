#!/bin/bash
# Run the eOn server until the AKMC run stops making progress or reaches
# max_kmc_steps. Each pass submits and harvests Slurm jobs; start it inside a
# tmux or screen session on a login node.
set -euo pipefail
while :; do
  eon
  if [ -f dynamics.txt ] && [ "$(wc -l < dynamics.txt)" -ge 3 ]; then
    break
  fi
  sleep 30
done
