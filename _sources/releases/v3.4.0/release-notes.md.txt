# eOn 3.4.0

A minor release on `v3.3.1`. See monorepo `CHANGELOG.md` section 3.4.0
for the full list.

## AKMC on a cluster with CPMD

`examples/akmc-cpmd-slurm` explores a silicon vacancy with CPMD on four
MPI ranks per search. Each process search is one Slurm job: the client
runs as one process and its `ext_pot` wrapper starts `cpmd.x` under
`mpirun`. On a test node the wrapper's force matched a central finite
difference to 2e-5 eV/A, and the first search returned a process with a
0.029 eV barrier. The ExtPot guide explains how CPMD reports forces:
the table under `GRADIENTS (-FORCES)` holds forces, and `GEOMETRY`
carries them at full precision.

Fixes on this path:

- A dimer or process search on a per-image potential such as `ext_pot`
  evaluated both images on one instance at the same time, with
  `parallel = true` (the default). Two external programs then shared one
  exchange directory. Each image now keeps its own instance, and an
  error in the second thread stops the job instead of calling
  `std::terminate`.
- The cluster communicator cancels only the jobs still queued when a
  state reaches its confidence, keeps finished results for harvest, and
  logs a failed cancel. It used to cancel every job and delete the
  scratch directory, results included.
- The MPI communicator no longer stops idle clients when a state reaches
  its confidence.
- The Slurm scripts take the account, partition and job shape from
  `EON_SBATCH_ARGS` and the client from `EON_CLIENT`, and list the
  current user's jobs.
- The local communicator kills each client's whole process group, so
  wrappers and their MPI launchers do not outlive the server.
- The MPI potential sends the working directory as `int`, as the
  `MPI_INT` it declares, and sleeps `mpi_poll_period` seconds between
  polls.

## pyeonclient

- Linux wheels are `manylinux_2_28` and vendor the Fortran and OpenMP
  runtimes, so they install from PyPI and import on a minimal system.
  The bundled rgpot keeps the private SONAME `libeon_rgpot.so.3`, so a
  pip `rgpot` wheel imports beside pyeonclient in either order.
- In-process jobs record termination, energy, force calls and job type
  on `job_result`. `LocalInProcess` takes a `Structure` or `ConFrame`
  and returns saddle and product as frames. In-process jobs poll a
  cancel token between jobs.

## Nudged elastic band

- Zoom-NEB packs images onto a window around the climbing image and
  keeps OCI-NEB available on that band.
- The climbing image applies only when an interior image is above both
  endpoints.

## Structures and numerics

- Structure matching, rotational matching and space groups go through
  readcon-ops. `identical` keeps a one-to-one atom map.
- Finite-difference Hessians gain a fourth-order stencil and
  cutoff-colored columns, which take fewer force calls on sparse
  systems.
- Per-atom force norms and periodic wrapping use Highway where it is
  available. NaN forces no longer read as converged.

## Build

- The direct in-process rgpot arm is linked on every non-Windows build;
  `-Dwith_rgpot` is deprecated and ignored. eOn builds against
  rgpot 3.3.0.
- `Parameters` keeps its state behind a private `Impl`.
- The pixi environment carries `nlohmann_json`, so a host copy cannot
  pull `/usr/include` into a conda build.
