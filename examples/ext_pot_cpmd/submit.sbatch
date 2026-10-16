#!/bin/bash
#SBATCH --job-name=eon-cpmd
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=eon-cpmd-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
exec bash ./run.sh "$@"
