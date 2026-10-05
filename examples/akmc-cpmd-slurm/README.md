# AKMC with CPMD on a Slurm cluster

A silicon vacancy (seven atoms in an eight-atom diamond cell, a = 5.43 A)
explored by adaptive kinetic Monte Carlo. Every process search is one Slurm
job; inside it, eOn's `ext_pot` potential calls `potfiles/cpmd_extpot.py`,
which runs a BLYP single point with `cpmd.x` on several MPI ranks.

## What the wrapper does

- Writes a CPMD deck from `from_eon_to_extpot`: `OPTIMIZE WAVEFUNCTION`,
  `PRINT FORCES ON`, and `CELL VECTORS` from eOn's cell.
- Runs `cpmd.x` under the launcher in `CPMD_LAUNCH`.
- Reads the energy from `TOTAL ENERGY` and the forces from `GEOMETRY`, where
  columns 4 to 6 hold the ionic force in Ha/Bohr at full precision. The force
  table in the CPMD output carries the same numbers to three significant
  digits under a `GRADIENTS (-FORCES)` header; the wrapper uses it as a cross
  check. The table holds forces, not gradients.
- Converts to eV and eV/A and restores eOn's atom order; CPMD groups atoms by
  species.
- Fails when the wavefunction optimisation reports `BUT NO CONVERGENCE`.
- Keeps `RESTART.1` in the ExtPot exchange directory, so every call after
  the first in a job starts from the previous wavefunction.

## Start from a minimum

`pos.con` is the vacancy cell already relaxed with this wrapper at 20 Ry
(eOn `job = minimization`, L-BFGS to 0.01 eV/A; E = -737.740 eV). AKMC
compares every search's endpoints with the reactant, so an unrelaxed
reactant makes each saddle look "not connected to initial state". For a
new system, relax first:

```ini
[Main]
job = minimization

[Potential]
potential = ext_pot
ext_pot_path = ./ext_pot

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
```

and use its `min.con` as `pos.con`.

## Run it

```bash
export CPMD_BIN=/path/to/cpmd.x
export CPMD_LAUNCH="mpirun -np 4 --bind-to none"  # or "srun -n 4 --overlap"
export PP_LIBRARY_PATH=/path/to/pseudopotentials/
export CPMD_PP="14:Si_MT_BLYP.psp:LMAX=P"    # Z:file:options, comma separated
export CPMD_CUTOFF=20                        # Ry; raise for production
export EON_SBATCH_ARGS="-p cpu -N 1 --ntasks=4 --cpus-per-task=1 -t 02:00:00"
export EON_CLIENT=$(command -v eonclient)
# Set script_path in config.ini to eOn's tools/clusters/slurm.
./run_until_done.sh
```

The server runs on the login node, where it submits and harvests jobs; start
it inside `tmux`. `sbatch` passes the exported variables into each job, where
the wrapper reads them. Each job asks for four Slurm tasks: `sbatch --wrap`
starts the client once, and `mpirun` (or `srun --overlap`) starts the four
CPMD ranks on those task slots. One task with four CPUs gives `mpirun` a
single slot, and OpenMPI 5 then refuses `-np 4` ("not enough slots"); the
search fails and the server makes new ones without end. `--bind-to none`
keeps OpenMPI from refusing to bind four ranks to cores when Slurm hands
out hardware threads (four CPUs on two physical cores).

## Checked

With this `pos.con`, the first process search on terra returned a
process with a 0.029 eV barrier to a product 3.86 eV lower, from one
Slurm job running the client and a four-rank `cpmd.x`.


OpenCPMD was built against OpenBLAS and OpenMPI on a 32-core node. On a
displaced Si8 cell at 20 Ry and four ranks, the wrapper gives
Fx = -0.63414 eV/A; a central finite difference gives -0.63416 eV/A. A
second call from `RESTART.1` reproduces the energy to 1e-6 eV.
