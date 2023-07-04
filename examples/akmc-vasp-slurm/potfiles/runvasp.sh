#!/bin/sh
srun -n 32 vasp_std
#srun -N 1 -n 8 -c 1 --cpu_bind=cores vasp_std
#vasp_std
