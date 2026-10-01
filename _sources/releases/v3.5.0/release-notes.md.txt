# eOn 3.5.0

A minor release on `v3.4.0`. See monorepo `CHANGELOG.md` section 3.5.0
for the full list.

## Tunnelling rates

`job = instanton` finds the ring-polymer instanton between two minima
for a tunnelling splitting. With `mode = rate` it finds the closed-ring
instanton through a saddle, for the thermal rate below the crossover
temperature (Richardson and Althorpe, doi:10.1063/1.3267318).

- The rate search is an index-1 Newton step on the ring Hessian. Each
  solve, determinant and inertia count goes through a block LU of the
  open bead chain plus a low-rank Woodbury correction, at O(N f^3) with
  no dense ring matrix.
- Without a band, the ring is seeded from a steepest-descent path out of
  the saddle by the period condition. A supplied band (`initial_path`)
  seeds it the same way and adds a one-dimensional WKB rate.
- A ring that converges onto another saddle is refused. On the LJ13 test
  pair, cooling a cosine ring from the crossover slides onto a
  neighbouring saddle that leads to a different minimum.
- An even ring runs on one half, from one turning point to the other. A
  half ring that has converged, or stopped stationary at the wrong index,
  is checked for an unstable mode odd under the mirror. That mode is two
  copies of the instanton on one ring, and the search finishes on the
  whole ring.
- The rigid motions of the whole ring are rebuilt from the current beads
  at every step and kept out of the step and the classification.
- Above the crossover, the parabolic-barrier factor multiplies quantum
  harmonic TST, so `ln k` is continuous at `T_c`.
- `temperatures` takes a list and writes one row per temperature to
  `rate_instanton.dat`. `instanton_centroid.con` writes the centroid with
  the readcon `spreads` section, the per-atom spread of the beads.

On the LJ13 pair at 198 K, 0.6 of the 330 K crossover, with 16 beads the
search converges in 32 Newton iterations and 2369 force calls. The ring
has one negative mode, and ln(k s) = -2.49 against -24.82 from harmonic
TST.

## Path-integral QTST and the ring-polymer rate

`[Instanton] pi_planes` samples a ring polymer with its centroid held on
parallel planes from behind the reactant to the saddle. The centroid
mean force integrates to a free-energy profile, which gives the
path-integral QTST rate (Voth, Chandler and Miller,
doi:10.1063/1.457242). `pi_recrossing_parents` adds the ring-polymer
recrossing factor. Parent rings on the top plane launch thermostat-free
RPMD children in momentum-reversed pairs, which give the
Bennett-Chandler plateau kappa, and `k_RPMD = kappa k_QTST`. On a two-dimensional curved
barrier the classical kappa is 0.779 +- 0.005, against 0.777 +- 0.003
from independent trajectories.

## Path-integral dynamics

The dynamics job runs path-integral trajectories with a normal-mode PILE
or a GLE read from a matrix file, a separate centroid Langevin
thermostat, and economised ring-polymer springs. Bead forces are one
batch. The normal-mode transform is one matrix product, about eight
times faster per step at 32 beads. The economised spring fit iterates to
its least-squares minimum and raises an error when it does not converge.

## Solid-state NEB

Solid-state NEB relaxes the lower-triangular cell of each interior image
together with the atoms, with the tangent and spring of
doi:10.1063/1.3684549. A calculator's Cauchy stress (xTB, the i-PI
virial, LAMMPS pressure) is read on the force call.

## In-process CPMD

- `[RgpotPot] ranks_per_image` splits the MPI world into cpmdc calculator
  groups, one CPMD session each, and a NEB spreads its images over the
  groups through `forceBatch`. The finite-difference Hessian spreads its
  columns the same way.
- `params_path` loads a CPMDParams file, and `input_block` passes CPMD
  `&SECTION` text ahead of the generated sections. The Python server
  accepts the CPMD key spellings the client reads.
- A failed engine call in one group raises on every rank, and the run
  aborts the MPI world instead of hanging.
- The rgpot subproject pins 3.4.0, and the build requires readcon-core
  0.16.0.

## Other changes

- The Hessian job writes `modes.con`, one frame per normal mode, and
  accepts the trivial-mode count of the structure's symmetries.
- `neb.con` frames carry the mass-weighted reaction coordinate and the
  band's one-dimensional tunnelling estimate.
- aKMC with `use_kdb = true` stores processes through `amsel.KdbStore`
  and the run's readcon-db corpus. The tsase `kdb` path is removed.
- `with_mpi` is a feature option and stays off unless `enabled`.
  `-Dwith_parallel_neb` has no effect; the NEB image pool uses
  `std::thread`.
- The GP dimer links a `gpr_optim` checkout in `subprojects/gpr_optim`
  under `-Dwith_gprd=auto`.
