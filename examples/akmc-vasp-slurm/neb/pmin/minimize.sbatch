#!/bin/bash
#SBATCH --job-name=eon_pmin
#SBATCH --output=min.log
#SBATCH --qos=debug
#SBATCH -N 1
#SBATCH --time=0:5:00
#SBATCH -C haswell
#SBATCH --account=user_account
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user_email

date
pwd

module load gcc
module load vasp/5.4.1_vtst-gcc
module load ase

eon
