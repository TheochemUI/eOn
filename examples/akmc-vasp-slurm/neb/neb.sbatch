#!/bin/bash
#SBATCH --job-name=eon_neb_test
#SBATCH --output=neb.log
#SBATCH --qos=debug
#SBATCH -N 1
#SBATCH --time=0:20:00
#SBATCH -C haswell
#SBATCH --account=m3179
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michael.j.waters@northwestern.edu


date
pwd

module load gcc
module load vasp/5.4.1_vtst-gcc
module load ase

eonclient
