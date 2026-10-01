

<!-- towncrier release notes start -->

## [3.5.0](https://github.com/TheochemUI/eOn/tree/3.5.0) - 2026-10-01

### Removed

- The tsase `kdb` path for aKMC (`eon/eon_kdb.py` and its import of the PyPI `kdb` package) is removed. `use_kdb = true` stores and suggests processes through `amsel.KdbStore` and the run's readcon-db corpus.

### Deprecated

- `-Dwith_parallel_neb` is deprecated and has no effect. The NEB image pool uses `std::thread`, so `[Main] parallel = true` needs no TBB.

### Added

- Above the crossover, `mode = rate` writes `parabolic_factor` and `rate_parabolic`, the factor `(pi*T_c/T)/sin(pi*T_c/T)` times the harmonic TST rate from the reactant and saddle Hessians. The ring search still runs only below the crossover. At the crossover the factor diverges and that temperature records no rate.
- An even bead count evaluates the rate-instanton potential from one turning point to the other and copies that half onto the closed ring. The rate is the same as a full evaluation, and each bead still contributes its Hessian.
- Path-integral quantum transition-state theory after the rate instanton: `[Instanton] pi_planes` samples a ring polymer with its centroid held on parallel planes from behind the reactant to the saddle, integrates the centroid mean force to a free-energy profile with block-average errors, and writes the quantum free-energy barrier and the PI-QTST rate beside the classical and instanton barriers, with `rate_piqtst.dat` and one centroid-plus-spread frame per plane.
- Path-integral trajectories on the dynamics job, with normal-mode PILE, a normal-mode GLE read from a matrix file and a separate centroid Langevin thermostat, and economised ring-polymer springs. Bead forces are one batch. A centroid hyperplane reports its mean force. Economised springs are refused with the normal-mode GLE and with the instanton.
- Ring-polymer MD transmission factor for PI-QTST. `[Instanton] pi_recrossing_parents` samples parent rings on the top plane and launches momentum-reversed pairs of thermostat-free RPMD children from each. The job reports the Bennett-Chandler plateau kappa with a jackknife error and `k_RPMD = kappa k_PI-QTST` in `results.dat` and `rate_piqtst.dat`. `kappa_piqtst.dat` holds the kappa(t) curve.
- Rings of at most 4096 active coordinates take an index-1 Newton step with a Bofill update. Below three quarters of the crossover that search starts at 0.85 of the crossover. An even ring runs from one turning point to the other. Larger rings follow the minimum mode and copy one half of an even ring. The rate uses the bead Hessians.
- Setting `[amsel] discover_decide = true` runs on the current state's process table while `use_mcamc` stays off. A barrier below `e_min_init` stays in the basin (0.20 eV on Si6N8 isomer 1, against a 0.30 eV flip out). The exit time and channel come from the mean-rate method or first-passage-time analysis. The `[amsel]` keys are `AmselConfig`.
- Solid-state NEB relaxes the lower-triangular cell of each interior image together with the atoms. The tangent and spring use the Jacobian of doi:10.1063/1.3684549. A potential without a stress tensor is differentiated on the cell. OCINEB, zoom, and the action-based springs stay off when `solid_state` is set.
- The Hessian job writes `modes.con`: one frame per normal mode, with the unit Cartesian mode as a readcon `displacements` section and `mode_eigenvalue`, `hbar_omega` and `wavenumber` in the frame metadata. `[Hessian] write_modes = false` turns it off.
- The finite-difference Hessian evaluates its displaced structures through `forceBatch` when the potential batches, so under `[RgpotPot] ranks_per_image` its columns spread over the CPMD calculator groups. Prefactor and instanton Hessians gain the same. A column checkpoint keeps the one-column path.
- The instanton rate on a long ring comes from the cyclic block fluctuation product. A small ring still uses the dense product, and the two have to agree.
- The rate instanton's Newton step, ring determinant and inertia count run through a block LU of the open bead chain plus a low-rank Woodbury correction (closure, cycle, rigid modes, eigenvector-following flips), O(N f^3) with no dense ring matrix; bead curvature blocks start from the saddle Hessian and take Bofill updates (`[Instanton] initial_hessians`), and the Eckart flux is checked against the exact transmission.
- The search accepts a list of temperatures and starts at the highest. Each temperature writes one row of rate_instanton.dat. A supplied band seeds the ring and reports a one-dimensional semiclassical rate.
- `[RgpotPot] params_path` loads a CPMDParams file, and `ranks_per_image` splits the launch into calculator groups. The user guide renders the model and a 7-image band on 28 ranks.
- `[RgpotPot] ranks_per_image` splits the MPI world into cpmdc calculator groups, one CPMD session per group, and a NEB spreads its images over the groups through `forceBatch`. Each rank receives every image's energy and forces, so all ranks keep the same band; single force calls run on group 0 and are shared the same way. Needs rgpot built with `-Drgpot:with_mpi=enabled` and libcpmdc with `cpmdc_bind_calculator`.
- `job = instanton` computes the ring-polymer instanton between two minima and its tunnelling splitting, beyond one-dimensional WKB along a band. It writes the beads to `instanton.con` and the splitting, action and diagnostics to `results.dat`. The beads are evaluated in one batch per iteration, so CPMD calculator groups run them in parallel.
- `job = instanton` with `mode = rate` computes the thermal rate through a saddle below the crossover temperature, from the closed ring-polymer instanton of Richardson and Althorpe (doi:10.1063/1.3267318). It reads the reactant and `saddle.con`, takes `temperature` in kelvin, and writes `rate_instanton` to `results.dat`. Above the crossover it writes the parabolic barrier factor times the harmonic TST rate. At the crossover the factor diverges and the job stops. The reactant Hessian shows which rotations to omit for a cluster inside a periodic cell. The cubic-well check is the decay of Caldeira and Leggett (doi:10.1016/0003-4916(83)90202-6).
- `job = instanton` writes `instanton_centroid.con`: the beads' centroid with the readcon `spreads` section (per-atom root-mean-square spread in Å), the delocalised configuration as centroid plus spread in both the splitting and the rate mode.
- `neb.con` frames carry `reaction_coordinate_mw`, the mass-weighted arc
  length along the band, and the first frame carries the band's
  one-dimensional tunnelling estimate: `hbar_omega_reactant`,
  `hbar_omega_product`, `tunnel_action`, `tunnel_splitting` (WKB,
  Landau and Lifshitz prefactor), `tls_energy` and `tunnel_deep_wells`, the
  inputs a two-level-system screen of glass minima needs. On a quartic double
  well sampled by 21 images the splitting is within 5 percent of the exact
  two-level gap once the barrier stands three quanta high.

### Developer

- The metatomic `setupeon` pixi task and `scripts/run_coverage_cpp.sh` pass `-Dwith_gprd=enabled` and `-Dwith_gprd=disabled`. `with_gprd` is a feature option, and meson rejects `True` and `false` for it at `meson setup`.

### Changed

- A library force call writes no RESTART.1, LATEST, GEOMETRY, or GEOMETRY.xyz. The orbitals for the next call stay in memory.
- Configuration models render as a field list. The PDF title page carries the logo, and chapters open on the next page.
- On a potential that batches (CPMD calculator groups), the NEB evaluates its two endpoints in one call, the improved dimer evaluates its two ends together, a band update routes image `i` to the same calculator whether or not the images before it moved, and an OCI restore puts the saved forces back instead of asking the potential again. L-BFGS seeds its first step from the min-mode curvature when the objective knows one, so the MMF dimer inside OCI-NEB spends no force call on the finite-difference probe. `NEBObjectiveFunction::getEnergy` goes through the batched band update.
- Ring-polymer dynamics transforms the beads to normal modes as one matrix product and holds a centroid hyperplane without transforming the ring, about eight times faster per step at 32 beads.
- The build requires readcon-core 0.16.0 or newer, from pkg-config or the subproject, and the wrap pins v0.16.0. `ConFileIO` writes spreads with `set_spreads_from_flat`, which 0.15.x does not have, so a 0.15 host package used to pass configure and fail to compile. Packagers raise the host readcon-core floor to `>=0.16`.
- The rgpot subproject pins the 3.4.0 release (calculator groups, MPI finalised at exit, CPMD session reuse and the abort after a failed engine call).
- `-Dwith_mpi=enabled` requires an embedded Python (`python3-embed`): the server rank of the MPI client runs `EON_SERVER_PATH` through `Py_Main`. Configure stops when the embed dependency is not found.
- `[cpmd]` holds the scalar CPMD message when `[RgpotPot] params_path` is empty. `cutOffRy` there overrides `cutoff_ry`. A loaded CPMDParams file keeps its sections, and `input_block` is appended to that file's blocks.
- `python -m eon.server --reset` also removes the `[Paths] kdb` catalog directory, wherever it points. A shared or curated catalog needs a copy outside the run. The readcon-db corpus `readcon.db` is not removed. The kinetic database guide describes both.
- `with_mpi` is a feature option and stays disabled unless set to `enabled`. `auto` does not link MPI. `-Dwith_mpi=enabled` builds the client/server program; calculator groups stay `-Drgpot:with_mpi=enabled`.
- aKMC with `use_kdb = true` writes each good process into `amsel.KdbStore` when it is registered: the barrier in eV, the prefactor, the mode, and the readcon-db frame keys. The reactant, saddle, and product frames go into the run's `readcon.db`. The next search of a matching state refines from the stored saddle, then uses a random displacement when no suggestion remains. The `kdb` module is not imported. `kdb_nf` (fraction, default 0.2), `kdb_dc` (angstroms, default 0.3), and `kdb_mac` (minimum mode cosine, default 0.7) are read from the ini. `Paths.kdb` is the catalog directory.

### Fixed

- A calculator-group run whose rgpot engine call failed (rgpot requested `MPI_Abort` at exit) now aborts the MPI world on exit. The driver skips the stop broadcast and eOn's exit handler calls `MPI_Abort` instead of `MPI_Finalize`, so ranks left inside a CPMD collective no longer hold the job until the walltime kill.
- A calculator-group run with the rgpot wrap (`-Drgpot:with_mpi=enabled`) constructs its groups again. rgpot 3.4.0 reports MPI only after `MPI_Init`, and the construction agreement that calls `MPI_Init` was gated on that report, so every `ranks_per_image` greater than 0 stopped with "ranks_per_image needs rgpot built with MPI". The agreement runs whenever the process was started as one of several ranks.
- A cpmdc job under mpirun calls cpmdc_finalize on every rank, then MPI_Finalize, then _Exit, so CPMD and MPI library destructors do not run on a finalized world. Rank 0 prints the engine error from another calculator and exits non-zero. A rank that cannot read params_path stops every rank before MPI_Comm_split. eOn's -Dwith_mpi build is the client/server program; calculator groups are mpirun -np N eonclient with -Drgpot:with_mpi=enabled.
- A failed evaluation in one `[RgpotPot] ranks_per_image` calculator raises on every rank after all results are shared. Other ranks no longer wait in a broadcast the failed calculator never joins, and workers keep serving the driver.
- A half-ring instanton that has converged, or that is stationary at the wrong index, is probed for an unstable mode odd under the ring mirror. That mode marks two copies of the instanton on one ring, and the search finishes on the whole ring. A cooling schedule probes only its last temperature.
- A process whose product column is -1 is an amsel exit when its barrier is at or above `e_min_init`. The absorbing label is a 32-bit id. The hop creates the product state from the process id.
- In the MPI communicator a client whose job fails (a potential error, a bad `config.ini`, an unknown job, a missing job directory) logs the error, stages its logs into the job directory and hands the directory back without `results.dat`. The server skips that result and the rank takes the next job, so one failed CPMD call no longer ends the allocation.
- JSON configuration (the Python binding and serve mode) reads every `[Hessian]` key that `to_json` writes, including `write_modes`, with `atom_list` accepted for `phva_atoms` as in the ini. JSON loading also derives the internal `path_pile_tau` and Andersen collision period from the femtosecond inputs when the keys are absent, so a `pile` or `piglet` run loaded from JSON uses the same damping time as one loaded from config.ini instead of stopping on a zero time.
- NEB image forces with `[Main] parallel = true` run on at most `std::thread::hardware_concurrency()` threads, each taking the next image, instead of one thread per image. An exception from one image's potential is rethrown after every thread joins; it used to reach `std::terminate`. GCC and Clang builds no longer need TBB; nvc++ keeps `-Dstdpar=cpu|gpu`.
- On the Morse Pt cell, `gprdimer` and `dimer` both stop at -1462.008706 eV. `gprdimer` uses 17 force calls and `dimer` uses 38. A Linux build links that method under `-Dwith_gprd=auto` when a checkout of the private gpr_optim repository at `b55c89e2` sits in `subprojects/gpr_optim`.
- The Hessian accepts the trivial-mode count the structure's symmetries give:
  6 for a free cluster, 5 for a linear one, 3 for a periodic cell and none once
  any atom is fixed. It used to report every periodic bulk cell as an error
  for having 3 instead of 6.
- The MPI communicator reads the job path a client returns with `ndarray.tobytes`. `tostring` is gone in NumPy 2, so the server stopped with an `AttributeError` at the first returned job.
- The Python server accepts a potential name in any case. An exact listed spelling wins.
- The Python server accepts the CPMD cutoff and functional spellings the client reads: `cutOffRy`, `cutoff_ry` and `cpmd_cut_off_ry`, and `cpmd_functional`, on `[cpmd]` and on `[RgpotPot]`. A config.ini with `[cpmd] cutoff_ry` or `[RgpotPot] cutOffRy` used to stop `python -m eon` with "unknown option". The eon-schema `Cpmd` and `RgpotPot` models accept the aliases and leave them out of `model_dump`, so a written file carries one cutoff key.
- The amsel superbasin gate returns its configured fallback when an amsel call fails, instead of raising `NameError` on an unbound `config`; the policy is `[amsel] on_error` (`fallback_single`, `unavailable_mcamc` or `raise`). Processes whose product is not yet linked to a state (-1) are left out of the graph handed to amsel, and the discover cutoff stays at `[amsel] e_min_init` instead of being lifted above every known barrier, which left the basin without an absorbing border.
- The calculator-group safety code in the RGPOT potential (the agreement before `MPI_Comm_split` when a rank cannot read `params_path`, and the cpmdc_finalize, MPI_Finalize, _Exit sequence at exit) compiles whenever rgpot links MPI: `-Drgpot:with_mpi=enabled` with eOn's `with_mpi` off, the documented calculator-group build, or an installed rgpot whose `rgpot.pc` defines `RGPOT_HAS_MPI`. It used to compile only with eOn's `-Dwith_mpi=enabled`.
- The economised ring-polymer spring fit iterates to the least-squares minimum, with the 10000-iteration cap of Zeng and Manolopoulos, and raises an error when it does not converge. The fit stopped after 500 Newton steps, partway along a flat valley of the objective, and returned frequencies that differed between compilers and platforms by up to 5% on the high-frequency modes (48 beads at a maximum reduced frequency of 20).
- The installed `eon` server package ships `cancel.py` and `process_id.py`.
  An installed AKMC server no longer stops with `ModuleNotFoundError: No module
  named 'eon.process_id'` on the first registered saddle.
- The installed `eonclient` finds `libeonclib` with a multiarch or other non-default `libdir`. Its run path is `$ORIGIN` plus the relative path from `bindir` to `libdir`, and the installed libraries carry `$ORIGIN`, where the potential modules are installed. Meson's own libdir entry is not added once `install_rpath` is set, so such an install used to exit 127.
- The minimum-mode rate-instanton search on half of an even ring checks a converged ring for an unstable mode odd under the mirror and, when it finds one, finishes on the whole ring. The half-ring search no longer returns two copies of the instanton on one ring with two negative modes.
- The rgpot potential links the Message Passing Interface (MPI) when `with_mpi` is enabled. A system `mpic++` on the default path no longer adds its include directory to a build that left MPI off.
- The server catalog lists the path-integral keys on `[Dynamics]` and `springs` on `[Instanton]`, with the same names and defaults as the client. JSON load reads those keys, including an instanton section.
- The solid-state band keeps the Cauchy stress a calculator already computed. xTB, the i-PI virial, and LAMMPS pressure are read on the force call. A cell difference is used only when that tensor is absent.
- With `[amsel] discover_decide` on, a kinetic Monte Carlo step no longer waits for the old repeat-count confidence. A lone barrier at or above `e_min_init` leaves through the mean-rate method or first-passage-time analysis, and a faster edge stays in the basin. A missing amsel package logs `status=unavailable` and does not step early.
- Without `initial_path`, the rate instanton seeds its ring from a steepest-descent path out of the saddle. A ring that converged onto another saddle is refused. The climb follows the last climb by overlap. A converged ring is classified with finite-difference bead Hessians. The rigid motions of the whole ring come from the current beads.
- `[RgpotPot] input_block` reaches the cpmdc backend as CPMD `&SECTION` text ahead of the generated sections, so an in-process CPMD run can be periodic instead of the isolated cold deck. A new `permanent_dir` key sets the CPMD `FILEPATH` for `RESTART` files.
- `eon.atoms` loads without minimage again. The `Cell.wrap_many` binding for
  readcon-ops ran at import time, so the aKMC server stopped with
  `ModuleNotFoundError: No module named 'minimage'` wherever minimage, which
  is not on PyPI or conda-forge, was not installed by hand. The binding now
  runs just before the readcon-ops calls that need it, and is skipped when
  minimage is absent.
- `pip install eon` declares PyYAML, which `eon.config` imports at start-up.
  Without it `python -m eon` failed with `ModuleNotFoundError: No module named
  'yaml'` in any environment other than the pixi and conda ones, which already
  listed it.
- `with_gprd=auto` leaves the GP dimer off when the `gpr_optim` fetch fails. `with_gprd=enabled` still stops configuration.

  A build directory that stored `with_gprd` as `true` or `false` rejects `meson setup --reconfigure` (`Option "with_gprd" value auto is not boolean`). `python scripts/migrate_with_gprd_option.py <builddir>` maps `true` to `enabled` and `false` to `disabled`, and the reconfigure then runs.
- pyeonclient 0.4.1 carries the eOn 3.4.0 client. PyPI already held a 0.4.0 built from July sources, and the 0.4.0 publish of the 3.4.0 client left those files in place; the wheel publish now refuses a version the index already has.
- pyeonclient builds with the default GP dimer: the min-mode bindings include `AtomicGPDimer.h` when `WITH_GPRD` is set.


## [3.4.0](https://github.com/TheochemUI/eOn/tree/3.4.0) - 2026-09-29

### Added

- Cutoff-colored finite differences on a four-atom line take 19 central force calls, against 25 for one column per coordinate, and match the serial matrix to 1e-8. The fourth-order stencil stays within 1e-10 of the complex-step derivative of z^4 + 0.3 z^2 at 0.8, where the central stencil is more than 1e-3 away.
- Hessian `fd_scheme = fourth` builds the assembled Hessian, and Lanczos and Davidson products, with a real fourth-order central stencil.
- In-process jobs poll a cancel token before each job and before each dispatch. ``cancel_state`` returns 1 and stops the next job only while a batch is running; an idle call returns 0 and does not stick. A compiled relax still finishes the current call.
- New example `examples/akmc-cpmd-slurm`: adaptive kinetic Monte Carlo on a silicon vacancy with CPMD on several MPI ranks behind `ext_pot`, one Slurm job per search. The ExtPot guide explains how the wrapper reads CPMD forces from `GEOMETRY`.
- The in-process record stores termination, energy, force calls, and job type on job_result. results.dat is written from that dict.
- The production neighbor list stays on vesin. ``neighbor_list_linkcell`` compares that cutoff list with linkcell k-nearest pairs.
- Wire Geometry carries forces, atom ids, and a per-axis fixed bitmask at ordinals 8, 9, and 10. Empty lists mean those fields are absent.
- Zoom-NEB packs NEB images onto a window around the climbing image and keeps OCINEB available on that band.
- `eon_schema.jobs` writes an `eon.trajectory.v1` manifest of exact frame geometry digests and the rgpot potential identity. Landfold embedding and FES inputs are read from that manifest, not from `results.dat`.
- `opt_method = xtsci` can borrow an eindir objective when eindir-core is installed. The adapter checks the ABI stamp and the eV/angstrom descriptor, and minimization results record the optimizer name.

### Developer

- Docs CI minimizes the LJ13 cluster and regenerates the plt-min tutorial figures.
- The readcon-core wrap pins v0.14.11. The v0.14.10 tag recorded its version as 0.14.9 in Cargo.toml and meson.build, so pkg-config reported 0.14.9.

### Changed

- Finite-difference and saddle-search jobs log through the client logger instead of printf, and the finite-difference step list is a fixed array rather than a sentinel. Monte Carlo writes its final structure once.
- Indistinguishable structure matches, rotational matches, and crystal space groups go through readcon-ops. The space-group result includes the Hall number. ``identical`` is the readcon-ops one-to-one map, and it binds ``Cell.wrap_many`` when the installed minimage only has ``displacement``.
- Lanczos, Hessian, prefactor, bundling, IDPP, and the BGSD and basin-hopping saddle searches use standard library parsing and drop unused helpers. Basin-hopping refuses a band that has no interior tangent, and both searches record the status they return.
- LocalInProcess takes a Structure or ConFrame on the job and returns saddle and product as ConFrames. In-process job dicts no longer carry .con text, and those frames are not written to disk.
- Per-atom force convergence (`Matter::maxForce` and the NEB `max_atom` metric) reduces row norms with Highway when the library is available. Fully fixed atoms stay out of the maximum.
- Periodic wrapping uses Highway for the minimum-image floor (`x - floor(x + 0.5)`) and the legacy unit-interval wrap. Builds without Highway keep the same scalar formulas.
- The direct in-process RGPOT arm is linked on every non-Windows build, and Cap'n Proto is required there. `-Dwith_rgpot` is deprecated and ignored, so existing invocations still configure. Windows builds still omit the NWChem/CPMD frontends.
- ``Parameters`` stores load state and option groups in a private ``Impl``.
  ``sizeof(Parameters)`` is that pointer. Const accessors and
  ``ParametersLoadAccess`` are the read and write surface. Option-group
  types stay in ``ParametersOptions.h``. ``Matter`` still exposes Eigen.
- eOn builds against rgpot 3.3.0, and the pyeonclient wheel check installs `rgpot>=3.3.0`, the first rgpot wheel that imports without a system OpenBLAS or OpenMP runtime.

### Fixed

- A LAMMPS worker ``waitpid`` interrupted by a signal is retried. A failed ``fork`` closes the pipes it just opened.
- A process search or dimer on `ext_pot` with `parallel = true` evaluated its two images on one instance at the same time, so two external programs shared one exchange directory; CPMD then read a corrupt `RESTART.1`. The images now keep separate instances, cloned from the job's potential, and a force call that fails in the second thread stops the job with its error instead of `std::terminate`.
- A signal during ``waitpid`` no longer makes a live VASP job look dead. The child closes the extra ``vaspout`` descriptor after redirecting it.
- GLE `apply` leaves velocities unchanged when any mass is non-positive, so a zero-mass DOF cannot write inf.
- LAMMPSPot locks the fixed-atom mask and the in-process force call on one mutex, including the Windows and MPI paths.
- LBFGS, conjugate gradients, and steepest descent allocate zero history vectors when the optimizer is constructed.
- Least-coordinated AKMC displacements accept per-coordinate free masks. Atoms with at least one free coordinate remain eligible, and fully frozen atoms are excluded.
- NEB skips the initial L-BFGS finite-difference curvature probe while preserving automatic scaling from successive steps. Climbing-image forces apply only when an interior image exceeds both endpoint energies; paths whose highest energy lies at an endpoint retain their spring forces. The convergence tests exercise the default climbing-image setting.
- The AKMC server honours `remove_translation` the way the client does: with no fixed atoms, a saddle or state that differs only by a rigid drift of the whole periodic cell matches its earlier copy, so repeated processes count as repeats and the state confidence rises. Every repeat was counted as a new process before.
- The Hessian job opens pos.con through getRelevantFile, so pos_cp.con and pos_in.con are honoured like the other jobs.
- The MPI potential sends the working directory as `int` character codes, matching the `MPI_INT` it declares, and sleeps `mpi_poll_period` seconds between polls instead of spinning. A run with no potential ranks, or a rank count that does not divide by the number of clients, stops with a message.
- The MSVC xtb import library aliases lowercase MinGW exports to the
  camel-case names in the xtb header.
- The NEB maximum-image tangent, the dynamics saddle search modes, and the Lanczos and Davidson eigenvectors no longer turn a collapsed vector into NaN. They stay zero.
- The Slurm scripts in `tools/clusters/slurm` read the account, partition and job shape from `EON_SBATCH_ARGS` and the client command from `EON_CLIENT`. `queued_jobs.sh` lists the current user's jobs and fails when `squeue` fails.
- The cluster communicator cancels only the jobs still queued when a state reaches its confidence. Finished results stay for harvest, and a failed cancel is logged instead of stopping the server. The MPI communicator no longer stops idle clients at that point.
- The local communicator starts each client in its own session and kills the whole process group on exit, so ExtPot wrappers and the MPI launchers they start do not outlive the server.
- The rgpot wrap builds ExprPot and the client defines RGPOT_HAS_EXPR, so PotType.EXPR works in wrap builds instead of stopping with the generic error.
- The saddle-search job opens `pos.con` through `getRelevantFile`.
- Windows MSVC builds of the native xtb potential link a generated xtb.lib from the conda-forge libxtb DLL instead of the MinGW import library.
- ``identical`` keeps a one-to-one map. A crossed pair within the tolerance still matches. Two atoms cannot both match the same partner.
- ``identical`` no longer maps two atoms onto one site. ``internal_motion`` translates atom 0 onto the reference atom. ``point_energy_match`` forwards ``use_identical`` into the match.
- ``identical`` rematches a close index pair whose elements differ. ``internal_motion`` turns a reversed bond around and removes the rotation about that bond.
- ``internal_motion`` places atom 0 on the reference atom before it removes the rotation.
- ``lammps_logging`` keeps ``client_lammps-N.log``. The worker child writes that file. The parent copies new lines into the process log after each force call. A new LAMMPS open truncates that file, and the copy starts over. The child does not call the process logger.
- `eonclient -c` compares a copy so `Matter::compare` cannot translate the first structure.
- `resolveMobileAtoms` / `freeAtomIndices` throw on a null Matter.
- pyeonclient Linux wheels are manylinux_2_28 again and import on a system without gfortran or OpenMP: the private librgpot SONAME step now rewrites RECORD and runs before auditwheel, which vendors those runtimes. A failed repair fails the build instead of shipping a `linux_x86_64` wheel.
- pyeonclient wheels bundle ``librgpot`` as ``libeon_rgpot.so.3`` so a
  side-by-side ``rgpot`` install does not share that SONAME.


## [3.3.1](https://github.com/TheochemUI/eOn/tree/3.3.1) - 2026-09-27

### Removed

- Drop leftover AMS debug cout dumps.
- Drop leftover Eclipse Created-on banners from Obs, CuH2,
  AtomsConfiguration, and GPRDimer unit tests.
- Drop leftover Eclipse auto-generated ctor/dtor TODOs in
  Obs, CuH2, AtomsConfiguration, and GPRDimer unit tests.
- Drop leftover MPI Recv and Table.write debug prints.
- Drop leftover Meantime debug stdout from each AKMC step.
- Drop leftover Python 2 `from builtins import input` in AKMC.
  Live reset/restart input() prompts stay.
- Drop leftover canary `os.system` tee next to the live redirect.
  Keep the redirect and pass/fail stdout.
- Drop leftover commented AMS debug readFile and validate_order. Live force
  and real TODO notes stay.
- Drop leftover commented BasinHopping force-call and recentRatio debug dumps.
  Live basin-hopping accept and displacement adjust stay.
- Drop leftover commented Broken AMS `test_one_pt_ams_dimer`.
  Live Morse `test_one_pt` stays.
- Drop leftover commented BytesIO alias, pylab import, and Python 2
  canary energy prints. Live StringIO I/O and pass/fail prints stay.
- Drop leftover commented DEBUGGER recieveFromSystem and old AMS force
  bodies. Live AMS force path stays.
- Drop leftover commented EONC_LOG_TRACE in GPSurrogateJob. Live traces stay.
- Drop leftover commented Matter operator==/!= and unused force/max
  variance stubs. Live compare() and getEnergyVariance stay.
- Drop leftover commented Matter::setPotentialEnergy stub. Live getPotentialEnergy stays.
- Drop leftover commented MinimizationJob include, unused
  status/returnFiles/epot_hop/earr, and fSPDLOG/fprintf debug dumps
  in GlobalOptimizationJob.cpp. Live hop and escape path stays.
- Drop leftover commented PR search-result and superbasin basin calls.
  Live allocate_process_id, make_basin_from_sets, and _get_filtered_states stay.
- Drop leftover commented Python 2 debug prints and dead max-rate
  blocks in AKMC state connect and Water displace.
- Drop leftover commented Python 2 prints in eon-state-stats.
  Live table output stays.
- Drop leftover commented dead branches in the eon server.
  Live saddle, basin-hopping, superbasin, match, and water-displace paths stay.
- Drop leftover commented dead statements and the rescaleVelocity
  stub in GlobalOptimizationJob.cpp. Live hop and escape stay.
- Drop leftover commented debug prints in ASE_NWCHEM client.py and
  blah.py. Live socket force and energy path stays.
- Drop leftover commented debug prints in disconnectivity_graph,
  gpaw_sp, and mk_cuh2_vid. XXX notes stay.
- Drop leftover commented energy and force printf dumps in
  MPIPot.cpp. Live MPI pot send/recv path stays.
- Drop leftover commented map_out_pes and plot helpers from the
  extpot 2D PES test fixture. Live _calculate_landscape stays.
- Drop leftover commented print_help/sys.exit in eon-minimize and
  the unused compareStru body in disconnectivity_graph.
- Drop leftover commented try/except around AKMC search-result append
  and Main random_seed. Live append_search_result write and getint stay.
- Drop leftover commented try/except around live dg.draw_minima
  in disconnectivity_graph. Live draw_minima stays.
- Drop leftover commented unit_system include and ERGS_PER_ANGSTROM2
  notes from Water CCL and SPC/E. Unit-system fallbacks stay.
- Drop leftover commented unused members in GlobalOptimizationJob, ReplicaExchangeJob, Hessian, and BondBoost. Live members stay.
- Drop leftover debug stdout of cwd in akmc-gui and bh-gui.
  GTK unused imports, pathfix, and os.fork daemonize stay.
- Drop leftover pot rank / my_client_rank debug stdout in
  tools/emt-sp.py. Commented GPAW constructor and live MPI rank
  split stay.
- Drop leftover register_process debug prints and the commented rate = cur_rate line.
  Live rate assignment, eq-rate clamp, and process-table writes stay.
- Drop leftover unused PoissonSolver from the GPAW import in tools/emt-sp.py. The commented GPAW constructor stays.
- Drop leftover unused algorithm includes in Obs, CuH2,
  AtomsConfiguration, and GPRDimer unit tests.
- Drop leftover unused fileio cell reexport and unused imports in server, state, and superbasinscheme.
- Drop leftover unused imports and debug prints in explorer
  and communicator. MPI Recv and harvest stay.
- Drop leftover unused imports in ASE_NWCHEM blah.py, client.py, and
  server.py. Live socket energy and force path stays.
- Drop leftover unused imports in analyze, basinhopping, config, eon.__init__, and eon_kdb.
- Drop leftover unused imports in tools disconnectivity_graph, emt-sp, evt, eon-ini, eon-config-docs, and mk_cuh2_vid. pathfix side-effect imports and the commented GPAW constructor imports stay.
- Drop leftover unused imports in unit tests. Keep pytest where mark, raises, approx, fixture, or importorskip is used.
- Drop leftover unused numpy in tools/modedot.py and unused PoissonSolver plus ase.utils.devnull in tools/gpaw_sp.py.
- Drop leftover unused sys, atoms, and numpy imports from tools.
- Drop the leftover commented eonclient path in eon-minimize.
  The live call is already `eonclient` on PATH.

### Added

- `[Xtsci] method = lbfgs` honours `[LBFGS] lbfgs_step` and `lbfgs_precon`. The pair matrix stays in eOn; xtsci applies it as \(H_0\) or as a Newton/RFO step.
- ``[Xtsci] qn_step`` and ``precon`` are first-class (Cap'n Proto, INI, JSON, YAML, Pydantic). They drive ``xts_solver_set_qn_step`` and the host pair Hessian. Native ``[LBFGS]`` is unchanged; ``lbfgs_step`` / ``lbfgs_precon`` remain a fallback when the Xtsci fields are left at their defaults.
- `ci_mmf_restore_unhelpful` (default false) restores the climbing image after an alignment reject or force increase. The Frontiers OCI-NEB article (doi:10.3389/fchem.2026.1807063, Algorithm 1) restores only on positive curvature; that is the default.
- `opt_method = xtsci` selects the xtsci-optimize engine. `[Xtsci] method` picks the solver inside it (L-BFGS, BFGS, SR1/SR2, Newton, RFO, NLCG conjugacies, Adam, PSO). YAML, JSON, pydantic, and SSOT defaults carry the same catalog.

### Developer

- The NEB regression tests run between two real LJ13 minima joined by a saddle. The climbing-image case converges at -O2 with NDEBUG as well as at -O3.
- The manylinux wheel build retries the Cap'n Proto download and uses the GitHub tag archive when capnproto.org drops the TLS handshake.
- Unit tests alias eonc types, including EigenmodeStrategy, from the shared test header.

### Changed

- AKMC displacement sampling, eon-minimize, and GPAW
  single-point leftover path I/O uses pathlib.
- AKMC leftover path I/O uses pathlib. Reset and restart share
  remove_tree_and_empty_parents and info_txt_path.
- AKMC movie leftover path I/O uses pathlib for POSCAR, graph.dot,
  and dynamics.txt.
- AKMC, BH, PR, and escape-rate pass the config path into
  `ConfigClass.init_from_cli` instead of rewriting sys.argv. Dead
  commented RandomStructure in basin hopping is gone.
- AKMC, PR, and escape-rate write `info.txt` through `fileio.write_info_txt`.
  The unused CatLearnPot `variance` member is gone.
- AKMC/BH helper leftover path I/O uses pathlib for search-stats,
  state-neighbors, fillkdb, and GTK glade loaders.
- AKMCState leftover path I/O uses pathlib for badprocdata,
  jobs.tbl, search results, and saddle/mode files.
- ARTnSaddleSearch accepts an injected IARTnResource. Production still uses the
  process-default ARTnResource singleton.
- AS-KMC leftover path I/O uses pathlib for askmc_data.txt,
  askmc_processtable, and recycle_path.
- Amsel superbasin split cache leftover path I/O uses pathlib.
- Basin hopping leftover path I/O uses pathlib for states, wuid.dat,
  bh.log, and the lockfile.
- CI concurrency groups use the PR number or git ref so develop pushes cancel superseded runs.
- CNA helpers drop the leftover unused brute argument. Importing
  atoms no longer mutates sys.getrecursionlimit.
- Communicator leftover path I/O uses pathlib: harvest, unbundle,
  scratch, script submit, and client lookup.
- Config leftover path I/O uses pathlib for config.yaml, config.ini,
  and the PRNG state file.
- Develop package version is 3.3.1.dev0 after the 3.3.0 release, not 3.2.1.
- Explorer leftover path I/O uses pathlib for debug results, KDB
  scratch, jobs.tbl, incomplete saddles, and explorer.pickle.
- GPRHelpers remaps Matter Z to dense 0..n-1 because gpr_optim pairtype
  is indexed that way. Stale Fortran meson and statelist TODOs are gone.
- Helper tools leftover path I/O uses pathlib for pathfix, INI
  dump, GTK glade paths, CNA walks, and state-stats.
- IRACompare accepts an injected IIRAResource. Production still uses the
  process-default IRAResource singleton.
- In-process communicator leftover path I/O uses pathlib for
  config.ini basename.
- LockFile leftover path I/O uses pathlib. Stale and missing
  locks unlink with missing_ok.
- MetatomicDynPot accepts an injected IMetatomicLoader. Production still uses
  the process-default MetatomicLoader singleton.
- Parallel replica and escape-rate leftover path I/O uses pathlib.
  info.txt goes through info_txt_path.
- Path-based con reads and writes keep a copy in a readcon-db corpus next to the run. readcon still parses the file, and the con file remains the copy the client reads.
- Potential accepts an injected IPotRegistry. Production still uses the
  process-default PotRegistry::get() singleton.
- Queue estimates use the communicator's construction-time bundle size.
  KDB insert on a finished state writes a marker so a second explore does
  not insert again. Table float precision is a constructor argument.
- Replace leftover Python 2 ``file()`` in config docs and ``xrange`` in
  mastereqn. Keep the live ``open()`` write; do not convert tools/toykmc.
- Replace leftover Python 2 print and SafeConfigParser in tools
  and dynamics analyze with configparser.ConfigParser and print().
- Repo tests, get_version, and Sphinx cache leftover path I/O
  use pathlib.
- RgpotAdapter accepts an injected IPluginLoader. Production still uses the
  process-default PluginLoader singleton.
- Server leftover path I/O uses pathlib for potfiles, cwd listing,
  and output/output_old.
- State leftover path I/O uses pathlib for procdata, reactant,
  processtable, and staging.
- Superbasin and KDB path I/O uses pathlib. Dead commented
  get_superbasins is gone.
- Superbasin merge always archives basins to storage. The MPI client
  stop path is STOPCAR. ASE-NWChem warns when a non-identity cell is
  ignored.
- Superbasin recycling leftover path I/O uses pathlib for
  recycling_data.txt, current_sb_states, and saddle_suggestions.
- The server package imports numpy as np. Optional annotations in modules that already postpone evaluation are written as X | None, and the config loader's section local is section_name.
- `load_potfiles`, `make_bundles`, and the leftover reset/KDB path
  helpers use pathlib instead of `os.path`.
- `make_bundles` loads potfiles. `load_potfiles` skips subdirectories
  under the pot directory, not same-named entries in the CWD.
- clang-format wrap leftover after the GlobalOptimization debug drop.
- eon package tests leftover path I/O uses pathlib for canaries
  and ConfigClass injection tests.
- eon-schema MainConfig leftover path I/O uses pathlib for
  jobs, states, potfiles, kdb, and superbasin defaults.
- eonc::Runtime is a move-only composition root. ClientEON constructs one and moves it into Job; Job builds the Potential from runtime.pots() rather than PotRegistry::get(). Catch2 uses a local Runtime. NEB image-force std::execution::par is unchanged.
- fileio leftover path I/O uses pathlib for atomic_write, savecon,
  Dynamics, and Table.
- gitversion and wheel-repair leftover path I/O uses pathlib.
  os.path.relpath stays for the meson dist printout.
- gitversion leftover meson-dist printout uses Path.relative_to.
- pyeonclient Session uniquely owns Runtime. make_job(params, session) borrows it with nanobind keep_alive so Jobs cannot outlive the Session. write_potcall_summary uses session.pots(); CLI Job borrows a stack Runtime. Stock nanobind, not jaxlib nb_class_ptr.

### Fixed

- NEB now rejects reactant and product structures with different atom counts before the Eigen path interpolation, instead of aborting inside the subtract. ([#540](https://github.com/TheochemUI/eOn/issues/540))
- A blank or non-integer `nperp_limitation` token returns an ARTn error and calls `artn_destroy` instead of throwing out of the search.
- AKMC movie fastest_path(full=True) sorts graph nodes with a
  Python 3 key instead of a Python 2 cmp.
- ARTn sends each lattice row to pARTn as a box column. A non-orthogonal cell is no longer transposed, so the periodic minimum image uses the cell eOn stored.
- ARTnSaddleSearch::run returns an error when has_sad is set but tau_sad was not retrieved, so the pre-convergence geometry is not recorded as the saddle.
- AtomicGPDimer restores floating-point traps if execute throws, so a saddle search that catches the failure does not leave later force calls running with traps masked.
- AtomicGPDimer::compute copies the Matter it is given before it builds the GP midpoint. A saddle search that moves the geometry after the solver is constructed is the geometry that is searched and written back.
- Basin hopping `total_normal_displacement_steps` subtracts only quench displacements that ran. A `stop_energy` break no longer counts configured `quenching_steps` that never executed, which could make that total negative.
- Basin hopping includes the last interior NEB image when it chooses the climb and when it writes ``neb_initial_band.con``. An image count below two no longer reads the bead before the reactant.
- Basin hopping keeps the minimized energy as the Metropolis reference when a jump does not relax. With significant_structure off, the raw jump energy is no longer stored and that geometry is not saved as the global minimum.
- Basin hopping linear and quadratic displacement keeps the configured step when every atom sits at the center, instead of dividing by a zero radius.
- Basin hopping skips displacement adjustment when adjust_period is zero, instead of taking the Monte Carlo step modulo that period.
- Basin hopping takes the dimer direction from the minimum-image of (r_3 - r_1) / 2. A bead that crosses the cell no longer starts the climb along a full box jump.
- BasinHoppingJob::getElements records atomic number 118 instead of writing past the element table.
- BasinHoppingJob::run applies exp(-dE/(kB T)) only when the hop is uphill and temperature is positive. Temperature at or below zero rejects that hop instead of dividing by temperature.
- BasinHoppingSaddleSearch::describeStatus labels a Metropolis rejection as Basin hop rejected instead of Initialized.
- BasinHoppingSaddleSearch::run applies `exp(-dE/(kB T))` only when the quenched hop is uphill and temperature is positive. Temperature at or below zero rejects that hop instead of dividing by temperature before the test.
- BasinHoppingSaddleSearch::run returns the min-mode climb status, so an unconverged or aborted climb is no longer reported as a good saddle.
- Bond-boost forces follow the minimum-image bond from distance(), so the stored bias force is the gradient of the bias energy in a non-orthogonal cell.
- Classic dimer no longer applies `rotation_angle` after the torque check accepts a direction, so `Dimer::compute` returns that orientation with the curvature measured there.
- Classic dimer rotation rejects non-finite `forceBatch` energies and forces in `Dimer::calcRotationalForceReturnCurvature` before `setComputedPotential`. A NaN torque no longer skips every exit check in that loop.
- ClientEON starts the wall clock at each job, so `time_seconds` is that job's duration instead of the time since the process began.
- Default L-BFGS takes the two-loop step (ASE). Energy and Grippo accept are opt-in.
- Doubly nudged projection now leaves the perpendicular spring remainder unscaled unless doubly_nudged_switching is set.
- Drop leftover dead rotation comments and the TShacked marker
  in atoms.py. point_energy_match returns False on mismatch.
- Dynamics saddle search keeps the detecting configuration when refineTransition has no MD snapshots, instead of indexing an empty list.
- Dynamics saddle search no longer takes a modulus of a state-check interval that rounds to zero steps. The interval is at least one dynamics step, so the state check still runs.
- Dynamics saddle search with bias_potential set to bond_boost now calls setBiasPotential and BondBoost::advance once per MD step, so the trajectory includes the bond-boost force.
- DynamicsSaddleSearch::refineTransition returns the first snapshot that does not minimize to the reactant. The reported transition time stays at or after the preceding reactant frame.
- Energy-weighted NEB now sets E_ref to the higher endpoint energy, so springs stay at k_min until an image rises above that minimum.
- FiniteDifferenceJob steps only free axes from Matter::getFixed(atom, axis). A column-4 mask that freezes one axis no longer moves that axis or includes it in the normalized curvature direction.
- Fix leftover TabError in tools/akmc-search-stats.py mixed tabs and spaces.
- GPR dimer check_derivatives is passed to gpr_optim as the strings "true" and "false". Assigning the bool stored one character, so derivative checks never turned on.
- GPSurrogateJob takes sub_job and linear_path_always from gp_surrogate_options and builds the surrogate with helpers::create::makeSurrogatePotential.
  Later bands turn climbing images off with neb_options().climbing_image.enabled and restart from copied Matter images.
- LBFGS compiles again: leftover OCINEB helpers::maxAtomMotionAppliedV calls now use eonc::geometry:: like FIRE, CG, SD, and Quickmin.
- Langevin dynamics no longer moves a coordinate frozen on only one axis.
  `Dynamics::langevinVerlet` applies noise and the position step per free axis, not only when the whole atom is fixed.
- MSVC Cap'n Proto builds force-include a guard that undefines
  windows.h `interface` after Meson's cl.exe sanity check.
- Matter passes a zero cell into ``force()`` when periodic boundaries are
  off. GFN2 treats a stored .con box as PBC and dies on multipoles.
- MonteCarloJob lists out.con once. Bundling can rename that file instead of failing on a second copy that is already gone.
- NEB L-BFGS no longer runs the auto_scale finite-difference H0 probe. That probe treats the projected band force as a potential gradient and never built memory. The LJ13 job fixtures also pin climbing_image off, matching the unit tests.
- NEB `minimize_endpoints` defaults to false, matching SVN and the
  pyeonclient NEB API. The old true default added endpoint relaxations
  to `total_force_calls` on fixtures that never asked for them.
- NEB file initialization keeps the endpoint structures passed to NudgedElasticBand. Frames listed in initial_path_in still fill the interior images, and minimize_endpoints_for_ipath is no longer overwritten by those files.
- NEB peak modes on the first and last spline intervals now use the geometric endpoint tangent. setup_mmf_peaks no longer writes the neighboring interior tangent or a zero mode for those peaks.
- Nose-Hoover (`noseHooverVerlet`) counts each unfixed Cartesian axis in `nFreeCoords`. A partly fixed atom is no longer thermostatted as three degrees of freedom.
- OCINEB restores a climbing image that finishes a min-mode walk at positive curvature before a band force under force_tolerance can mark the NEB converged.
- OCINEB scores MMF help on the walked image so a CI-index hop cannot rewind a finished saddle step.
- OH-TST `planes_used` now counts every sampled plane, including the plane where the progression converges.
- OH-TST keeps the unwrapped free coordinates in `samplePlane` when the next hyperplane is loaded. A reload after `setPositionsFreeV` shifted a coordinate by a lattice vector, and the projection then left the constrained restart.
- OH-TST symmetry products now minimum-image and remove rigid drift with the same fixed point as the primary guideline. `symmetryReflect` compares those rays in one frame.
- OHTSTJob::symmetryReflect measures distance to each product half-line, so a configuration behind the reactant uses the distance back to that endpoint instead of the perpendicular distance to the infinite line.
- On Apple x86_64, the SIGFPE handler now sets MXCSR exception masks in the saved context, and disableFPE writes those masks back instead of restoring the unmasked environment from enableFPE.
- PMF scan places OH-TST planes on a uniform grid from the reactant to the product. The scan no longer starts at s_init and then steps backward.
- Parallel replica dephase keeps the bond-boost pointer when it copy-assigns the trajectory back to the saved state. Later dynamics steps still apply that bias while the boosted clock advances.
- Parallel replica keeps its bond boost across the copies it makes of the trajectory: the dephase reset and, when ``stop_after_transition`` is false, the transition structure. ``Matter::assignKeepingBias`` does the copy.
- Parallel replica reads `stop_after_transition`. When it is false, the dynamics keep the full step count after the first new state instead of stopping.
- PotentialBase::computePt includes a platinum Lennard-Jones pair when one atom is fixed, so the fixed atom still contributes energy and force on its movable neighbor.
- Same-size ``Matter::resize`` keeps .con column-5 atom ids. A leftover
  duplicate resize body was resetting them to ``0..N-1``.
- Session supports Python weak references, so a Job that borrows it keeps the Session alive after the Python name is dropped. The installed ARTn saddle header no longer names the feature build flag.
- Shared plugin TUs include Parameters.h so accessors compile when
  Potential.h only forward-declares Parameters. Potential constructors
  that take Parameters live in eoncbase so plugin .so files do not need
  eonclib. Duplicate leftover field-style option blocks are gone.
  readcon is resolved before plugin subdirs so SocketNWChem links rkr_*.
- The MPI build links the MPI potential against MPI.
- The MPI client copies `client_quill.log` and `client_traceback.log` into the job directory before `return_files.dat` names them. Harvest keeps those logs instead of dropping names that exist only in the launch directory.
- The MPI client stops when `enterJobDirectory` cannot enter the job directory. It no longer runs that job in the launch directory or sends the path back as a finished job.
- The MPI client waits for the job-path send to finish before that buffer is released.
- The NEB no longer puts a projected force on fixed atoms. When the endpoints placed a fixed atom differently, the spring force along the tangent leaked onto it, and the norm convergence metric counted it.
- The Windows floating-point exception filter stores the exception-mask bits in ContextRecord->MxCsr before EXCEPTION_CONTINUE_EXECUTION. _controlfp_s alone does not survive that restore, so the faulting instruction trapped again.
- The Zhang-Xu LBFGS quadratic test writes optimizer options through ParametersLoadAccess, matching the private Parameters layout.
- The release dry-run accepts a development version when pyproject.toml and pixi.toml match. A .devN tip does not need a changelog section yet.
- Wheel and sdist builds pin nanobind to ``>=2.2,<3``, matching pixi.
  Unconstrained ``pip install -U nanobind`` pulled 3.0, which dropped
  ``nb::detail::keep_alive`` used by pyeonclient ``tie_lifetime``.
- Wheel repair searches for ``librgpot.so*`` next to ``libreadcon_core``,
  so pyeonclient can import without a host ``librgpot.so.3``.
- When NEB `max_iterations` is not positive, `DynamicsSaddleSearch::run` passes the initial-band tangent to `MinModeSaddleSearch` instead of an empty mode.
- With no equilibration samples (hyperdynamics rmd_time of zero), BondBoost::advance and BondBoost::boost measure the current tagged-bond lengths before BondSelect. The bias is no longer stuck at dvmax, so a stretched bond can turn the hypertime factor off.
- `AtomicGPDimer` passes a zero cell into the GP force box when
  `Matter::getPeriodic` is false. XTB no longer treats a stored .con box
  as periodic boundaries on those force calls.
- `ConjugateGradients::line_search` clamps a negative max move by its length, so a forward secant step stays forward when bowl breakout passes that value.
- `Lanczos::compute` diagonalizes the current Krylov block when the residual vanishes, instead of returning the Ritz pair from the previous subspace.
- `NudgedElasticBand::updateForces` passes a zero cell into `forceBatch` when `Matter::getPeriodic` is false. Interior images no longer look periodic to a batch evaluator that infers boundaries from the stored box.
- `Parameters::load` accepts an AMS or AMS_IO file whose engine is DFTB with resources, or FORCEFIELD with only the engine.
- `Parameters::load`, `Parameters::load(FILE*)`, and `Parameters::load_ini_text` return failure when the INI parser reports an error, including a missing '=' or a line past the maximum line length. A partial config is not applied.
- ``ARTnSaddleSearch::run`` holds ``library_mutex`` from ``artn_create`` through ``artn_destroy``, including force calls between ``artn_step`` invocations, so another search cannot touch the process-global Fortran state mid-run.
- ``structure_to_matter`` treats ``Structure.periodic`` as a 3-vector.
  ``bool(array)`` is a ValueError; use ``np.any``.
- `eonclient -p` no longer starts serve mode unless a serve flag is set. A structure file on that command is read as a normal client run.
- `geometry::identical` searches for a one-to-one pairing inside distanceDifference instead of locking an index-aligned pair that blocks the only valid permutation.
- `matter2xyz` append inserts a newline when the existing file does not end with one, so the next frame atom count starts on its own line.
- `opt_method = xtsci` holds an `xts_solver_t` session. Each `step()` is one outer iteration, so L-BFGS pairs and NLCG conjugacy survive the host loop.
- `pushApart` steps along the minimum-image separation. A coincident pair takes a finite step instead of dividing by zero.
- `unbundle` restores compressed names that `bundle` writes, so `results_3.con.gz` is copied back to `results.con.gz` instead of being skipped.
- findSplineExtrema stores one stationary point when the cubic coefficient is zero or the cubic root is repeated, instead of writing that point twice.
- fpe_signal_handler no longer returns into an x86 integer divide. FPE_INTDIV and FPE_INTOVF step past the faulting instruction instead of only masking floating-point traps, so sigreturn does not re-execute that divide.
- get_param and get_runparam now take a void** out pointer. pARTn allocates the value and writes that pointer through the argument, so a void* out parameter stored it at the wrong address.
- pyeonclient wheels carry the vendored rgpot as libeon_rgpot.so.3. Under the shared librgpot.so.3 SONAME, the pip rgpot wheel bound to that copy. Importing rgpot after pyeonclient then failed on a missing rgpot::D3Pot symbol, and the metatomic engine was not found.
- sortedR compares neighbor shells by the number of distances in the shell, not by the atom count. An identical geometry, including a homonuclear diatomic compared with itself, now matches.
- xtsci rematch now uses the native per-atom `max_move` clip, rigid-mode projection, cautious pair filter, and extra-updates instead of a Euclidean cap on the whole 3N step.
- xtsci sessions evaluate energy and forces in one potential call. ``[Xtsci] accept`` (`none` / `energy` / `nonmonotone`) replaces a mandatory energy backtrack that spent extra force calls on every Newton / pair step.


## [3.3.0](https://github.com/TheochemUI/eOn/tree/3.3.0) - 2026-09-15

### Removed

- Drop the unused `with_pybind11` meson option and the dead `EONMPIBGP` include in ClientEON. Feature-gated public headers remain.
- OCI-NEB MMF triggers on the adaptive threshold only. The floor is `2 * force_tolerance`; `ci_mmf_after` is no longer an extra OR.

### Added

- AKMC CI builds with `-Dwith_pyeonclient=true`. XTB CI includes windows-2022. Lanczos docs name stiff EAM/Tersoff. `load_ini_text` and AMSEl `on_error` / `fallback_single` are in.
- Document Potential/Matter thread-safe interface. `matter_to_ase` copies atom ids. Testing guide names Baker 05/17/19/20 as the reduced parameter-check set.
- External C++ includes `eon/api.h` for Matter, Parameters, NEB, IDPP/SIDPP, and dimer. Job headers stay internal.
- Hash-pin GitHub Actions (zizmor), add `ci_zizmor.yml`, and run zizmor from prek.
- Hessian, Prefactor, and Dynamics jobs write results.dat from JobResultEnvelope. LocalInProcess builds CON text only when a caller reads min.con.
- IRACompare.cpp is always compiled so NEB can call match_endpoints when
  IRA is off (the methods stub). `matchArrays` plus `pyeonclient.ira_match`
  let `eon.atoms.rot_match` use IRA when the wheel has it, else Kabsch.
- In-process `process_search` puts `saddle.con` on the result record so explorer can register the saddle without an `eonclient` subprocess.
- JobResult Cap'n Proto names historical `termination_reason` integers as `TerminationCode` (0–21). `statusCode` stays the results.dat integer.
- JobResult has a `body` union (`minimization` / `neb` / `processSearch`). Shared scalars stay on the outer struct.
- Jobs register themselves into a table; `makeJob` looks up the factory. AMSEl persist_split writes a sidecar so a resume does not re-split the basin.
- LAMMPSPot accepts an injected ILammpsLoader. Production still uses the
  process-default LammpsLoader singleton.
- LocalInProcess also dispatches dynamics, Monte Carlo, basin hopping, Hessian, prefactor, and finite-difference. Those C++ jobs write results.dat from JobResultEnvelope.
- LocalInProcess dispatches minimization, point, process_search, and saddle_search. PointJob writes results.dat from JobResultEnvelope. JobType.OH_TST is bound.
- Meson `-Dstdpar=cpu|gpu` turns on nvc++ `-stdpar=multicore|gpu` for the NEB image-force `std::execution::par` path without TBB. `-Dstdpar=gpu` also passes `-gpu=cc80,cc90` (override with `-Dstdpar_gpu_cc=`). GCC/Clang still use `-Dwith_parallel_neb=true`. Morse atom loops stay on the host. `QMC_GPU` is a QMCPACK flag, not an eOn option. eOn skips the cmake Highway wrap on nvc++ so Meson does not import `hwy_list_targets` (PGICompiler has no `get_pie_args`).
- Metatomic `forceBatch` tries a single `model.forward` over N systems and falls back to sequential `force()` if that path throws.
- MetatomicPotential can clone the already-loaded Torch module for per-image NEB instances without a second disk load.
- MinimizationJob fills an in-memory JobResultEnvelope and writes results.dat from it. The Cap'n Proto schema is unchanged.
- NEB writes `peakNN_pos.con` and `peakNN_mode.dat` for every interior spline maximum (`setup_mmf_peaks`, default on). Those files are the dimer seeds for a follow-up saddle search.
- NEB, dimer, and process search query pot thread-safety through `potAllowsSharedInstance` in `eon/PotCapabilities.h`. Virtual `isThreadSafe` stays the override point.
- Optional ASE ``batch_calculate`` sets ``supportsBatchEvaluation``.
  Without the hook the pot stays on sequential ``force()``. From #409.
- Parameters can load INI from a memory string (`load_ini_text`). LocalInProcess uses that instead of a temp file. Docs floors drop the invalid rgpycrumbs `analysis` extra.
- RgpotAdapter rejects a kernel whose `caps().reentrancy` is not SharedInstance, PerInstance, or ProcessSerial. OCI-NEB hybrid-dimer defaults stay `ci_stability_count=5` and `max_steps=1000`.
- SaddleSearchJob and ProcessSearchJob accept an in-memory Matter via `runFromMatter`, so pyeonclient does not have to write pos.con first.
- With `write_movies = true`, saddle search writes `mode_000.dat`, `mode_001.dat`, … so dimer mode-evolution plots have per-iteration eigenvectors.
- `Potential::layoutFlags()` reports in-process vs cwd vs subprocess execution. ExtPot sets Subprocess and NeedsWorkingDirectory.
- `RGPOT_NWCHEM_ENGINE` is accepted as an alias of `RGPOT_NWCHEMC_ENGINE`. Devdocs list which pots live in rgpot vs in-tree.
- `[Nudged Elastic Band] match_endpoints = true` rigid-aligns and permutes the reactant onto the product with IRA before interpolation.
- `[Potential] thread_safe = false` keeps a shared Potential serial. EAM and EMT refuse shared-instance threading (cell lists / ASAP state).
- ``eon.geometry.pbc`` uses minimage ``wrap_many`` for packed displacements.
  ``neighbor_list_pairs`` returns vesin ``ijS`` rows without unique-index
  MIC reduction.
- `eon_schema.jobs` encodes and decodes JobResult dicts (`job_result_dumps` / `job_result_loads`). pycapnp is used when installed; otherwise JSON of the wire field names.
- `job = test` now constructs TestJob. The implementation is linked into eonclib and writes results.dat (OK/FAIL/SKIP per pot). Missing `pos_test.con` is SKIP, not a link error.
- pyeonclient exports ProcessSearchJob with min1/min2/saddle/prefactors. LocalInProcess uses that job when present. NEB writes results.dat from JobResultEnvelope plus image keys.

### Changed

- AMS and AMS_IO look up symbols with `readcon::z_to_symbol`. GPR dropped
  an unused private element table.
- ASE and CatLearn pot headers no longer inject `using namespace
  pybind11::literals`. The `_a` literal is function-local in the TUs.
- ASE-NWChem `basis` and `memory` come from `[ASE_NWCHEM]` instead of
  hard-coded `3-21G` / `2 gb`. Defaults are unchanged.
- ConjugateGradients, FIRE, LBFGS, Quickmin, and SteepestDescent headers no longer inject `using eonc::…`. Implementations live in `namespace eonc`.
- Dimer `opt_method` is `OptType` (CG/LBFGS/SD) like the other optimizers. INI still accepts `cg`/`lbfgs`/`sd`.
- Dimer, ImprovedDimer, Lanczos, Davidson, LOR, LowestEigenmode, MinModeSaddleSearch, and EigenmodeStrategy headers no longer inject `using eonc::…`.
- IRA CShDA marks unassigned pairs with `huge()` instead of 999.9.
  A 38-atom core-plus-outlier pair that used to return a non-bijective
  permutation now assigns every index.
- Matter, Parameters, NEB, and the RPC server TUs live in `namespace eonc`. `NudgedElasticBand.h`, `Optimizer.h`, and `ServeRpcServer.h` no longer inject `using eonc::…`.
- Minimization, Point, Dynamics, Hessian, FiniteDifference, and Test job headers no longer inject `using eonc::…` into the global namespace. The job factory qualifies `eonc::` at the register site.
- More public headers (dynamics, Hessian, IDPP, IRA, saddle methods) no longer inject `using eonc::…`. Implementations live in `namespace eonc`.
- NEB image forces use `std::execution::par` when built with `-Dwith_parallel_neb=true` (TBB), instead of one raw thread per bead.
- Option structs live in ``ParametersOptions.h``. ``Parameters.h`` is the
  loader aggregate and keeps ``Parameters::neb_options_t`` aliases.
- Parameters option groups now expose only const accessors. INI/JSON loaders and
  bindings write through ``ParametersLoadAccess``; MPI comm/rank use dedicated
  setters so ``ParametersMpi.h`` no longer takes the address of a temporary.
- Pot headers inherit `eonc::Potential` directly. `Potential.h` no longer injects `using eonc::Potential`.
- Production client sources no longer use `using namespace eonc::helpers`. Call sites are qualified.
- Public headers no longer include mpi.h or gate EigenmodeStrategy on
  WITH_GPRD. Potential and Matter forward-declare Parameters. Optional
  MPI helpers live in ParametersMpi.h.
- Public headers no longer inject `using eonc::Matter`, `Parameters`, or the
  BaseStructures enums. Helper rng/geometry names are called as `eonc::rng`
  and `eonc::geometry`. `Potential::get_ef` is defined on `eonc::Potential`.
- Python CNA now builds unique-index adjacency from vesin ``ijS`` pair
  arrays. The local element table is a Z-keyed radius/color overlay;
  symbol and Z lookups stay on readcon.
- Python CNA uses the same 0/1/2 labels as the C++ client (fcc/hcp/other).
  Common-neighbor lookup is a set membership test that keeps neighbor-list order.
- Remaining job headers no longer inject `using eonc::…`. Implementations live in `namespace eonc`. `Job.h` still has the `using eonc::Job` alias.
- ServeMode, BondBoost, Rgpot, and AMS use `std::ranges::transform` for case folding instead of `std::transform`.
- SocketNWChem looks up symbols with `readcon::z_to_symbol` instead of a
  private element table. The last owned file-scope `using` directives are gone.
- The ASE potential maps incoming C arrays with Eigen::Map, matching
  ASE-ORCA, instead of building pybind11 shape vectors by hand.
- The client Fortran compile no longer passes `-w`. Named `-Wno-*`
  flags stay. Water displacements are no longer a fake `Displace` subclass.
- The command-line parser qualifies Argum types instead of a file-scope
  using-directive. Helper file helpers use `std::string` / `std::ifstream`.
- The eon-schema `rgpycrumbs` extra matches the tree pin (`>=1.10.4`). Dimer docs mark LOR as added in 2.17.2.
- Water and Water_Pt potential TUs wrap implementations in
  ``namespace forcefields`` instead of file-scope ``using namespace``.
- `Job.h` and `ServeMode.h` no longer inject `using eonc::…` into the global namespace.
- ``Parameters::load``, ``load_ini_text``, and ``load_json`` take
  ``std::string_view`` instead of owning ``std::string``.
- ``Parameters`` now has a private load-state ``Impl``
  (``last_load_source`` / ``last_load_error``). Option-group layout stays
  in the installed header and is not ABI-stable. ``Matter`` still exposes
  Eigen.
- ``Parameters`` option groups are private. Callers use ``main_options()``
  and the other accessors. ``load`` / ``to_json`` still round-trip.
- ``Potential::force`` has a ``std::span`` overload that checks sizes
  before the raw C-array virtual. Matter, ``get_ef``, and surrogate
  ``get_ef_var`` use it. Fortran/FFI loaders keep the pointer API.
- ``eon.atoms.atomic_number`` / ``symbol_for_z`` call readcon's Python
  ``symbol_to_atomic_number`` / ``atomic_number_to_symbol`` (0.14.9+). The
  ctypes probe of ``rkr_symbol_to_z`` is gone.
- `eon.atoms.atomic_number` / `symbol_for_z` prefer readcon for Z/symbol. The local table keeps radius and colour for viz.
- `eon.geometry.neighbor_list` takes vesin's per-axis `periodic` flag.
  `Structure.periodic` defaults to all-periodic and is used when the
  argument is omitted.
- `eon.geometry.pbc` accepts a bool or length-3 mask so a free axis is
  not minimum-image wrapped. Partial masks use the numpy wrap.
- `import eon` no longer imports the AKMC server. `eon.server` loads on first access.

### Fixed

- NEB `minimize_endpoints` defaults to false, matching SVN and the
  pyeonclient NEB API. The old true default added endpoint relaxations
  to `total_force_calls` on fixtures that never asked for them.
- A failed `hessian.dat` write now fails the Hessian calculation instead of returning a successful eigen solve.
- AKMCState stores the search count in MetaData and increments it on each append, instead of rereading the search log.
- AtomicGPDimer throws on a null Matter and returns a zero mode if the GP orientation size does not match the free-atom count.
- BGSD uses `safe_normalize` on the trial mode so an all-fixed Matter does not write NaN.
- Bond-boost returns zero bias when QRR or the selected bond count is not positive, and skips a bond whose length or equilibrium is zero.
- Collective IDPP `lastMaxForce` uses free-atom residuals only, so frozen pair forces no longer pin SIDPP at `max_iterations`.
- Collective IDPP leaves a collapsed-image tangent as zero instead of `normalize()` to NaN.
- Cookbook tests skip unless ``EON_PET_NEB_ROOT`` / ``EON_PET_MAD_MODEL`` are set. ASE NWChem examples take ``NWCHEM_COMMAND``.
- Dephase and INI dynamics step counts no longer divide by a zero `time_step`. A non-positive dt throws in the dephase paths and yields 0 steps in INI load.
- Document that XTB NEB is serial: `isThreadSafe` and `needsPerImageInstance` are both false because of Fortran restart units.
- Documentation PRs no longer get `contents: write`. The build job is read-only; deploy is a separate job that runs only on push.
- Dynamics throws on a null Matter. Kinetic temperature, thermal velocities, and Andersen redraws skip when there are no free atoms or a mass is not positive.
- EpiCenters throw on a null Matter instead of dereferencing it.
- Finite-difference displacements use `structure_comparison.neighbor_cutoff` and reject a zero displacement instead of normalizing it to NaN.
- Fix Matter size/cache bugs, dimer rotate/PBC/fixed-atom mode, NEB/Helper/min-mode signed max-component, FIRE zero-force/max-move, Hessian failed-cache, CG extra step, Dynamics extra MD step, BondBoost 0/dt and listed-atom range, Prefactor frozen atoms in the Hessian filter, GP surrogate target layout, OCINEB/GP-dimer eigenvector size, min-mode status/forcecalls/confine loop, ARTn direction.dat load, process-search pot/mode/displacement, fCallsMin double-count, replica-exchange config swap and per-thread ran2, basin-hopping accept/compare/swap, Nose-Hoover G2/Q2 plus post-Verlet KE, EpiCenters empty-set throw, MonteCarlo Metropolis sign, Quickmin zero-force, HessianJob failed-freq status, TAD/SafeHyper clock clamp and refine bounds, SafeHyper bias binding, BGSD NaN/iteration cap, OCINEB converged_only plus SIDPP collapse, TAD empty product, ParallelReplica first-hop stop, NEB path size, ASE position alias, SD two-point PBC, Prefactor empty moved set, ASE mixed-PBC throw, ASE unknown-constraint throw, readonly forces view, non-finite force reject, LAMMPS extract null and missing in.lammps, setForces mask persist, Structure-Matter atom_id copy, job catch of runtime_error, Lanczos/Davidson max_iterations floor, LBFGS negative s·y skip, LAMMPS neighbor rebuild and restricted cell, and Structure id uniqueness. Morse Pt saddle energy is unchanged; the dimer-rotate fix shifts the FD curvature from -1.014995 to -1.010564 and process-search force calls from 67 to 69.
- GPRD CI no longer reports Success when `SUBMODULE_PRIVATE` is missing. A `gprd-gate` job records availability; `gprd-integration` is Skipped on forks instead of a fake pass.
- IRA compare returns an error on a zero-atom Matter instead of indexing `cand1[0]`.
- LAMMPS `setforce` follows Matter's per-axis freeze mask. An atom frozen
  only in z is no longer free to move in z inside LAMMPS.
- LAMMPS applies `fix setforce 0 0 0` to atoms with `Atom.fixed` set on all three axes before `run 1`, so buffer atoms stay put during post-saddle minimization.
- MPI `stop_clients` sends `STOPCAR` to ready client ranks. AKMC uses that instead of `Abort`.
- Matter atom accessors throw `out_of_range` on a bad atom index or axis instead of indexing Eigen unchecked.
- Meson accepts readcon-core 0.14.9. The v0.14.10 wrap tag still reports crate version 0.14.9, so a cargo-c install failed the old `>=0.14.10` check.
- Metatomic seeds the Torch RNG from `main.randomSeed` when rotation averaging is on. NC force blocks must have `nAtoms*3` values before reshape.
- Min-mode `bowl_breakout` does not index atom 0 on an empty Matter and caps the bowl set at `nAtoms`.
- Minima-hopping Boltzmann accept/reject compares the hop to the current minimum and uses `exp(-dE/(kB T))` with `params.constants.kB`.
- NEB `isUncertain` now includes the product endpoint (`path[numImages+1]`).
- NEB spline extrema leave a collapsed endpoint tangent as zero instead of `normalize()` to NaN.
- Process search opens `pos.con` through `getRelevantFile`, throws on a failed load, and compares a copy so `Matter::compare` cannot translate `initial`.
- Removed the empty `CUH2_POT` ifdef from `Potential.cpp`.
- Replica exchange rejects `replicas < 1` and applies the same temperature fallback as the Python bind (`low` from main T or 300 K, `high` at least `1.5 * low`).
- SIDPP throws if adjacent images collapse. `climbing_image_converged_only` no longer reports GOOD while the rest of the band exceeds `climbing_image_band_slack` times the force tolerance.
- Same-size `Matter::resize` and displacement.con loads keep the reactant's `.con` atom ids instead of restamping 1..N.
- SocketNWChem drops a dead i-PI connection and accepts a new NWChem client once, so process_search / AKMC can continue after `task scf optimize` exits.
- Superbasin `make()` also calls `connect_states` on the new basin so intra-basin reverse processes register.
- TAD and SafeHyper rewind to `mdBuffer[0]` when refine returns 0. SafeHyper skips the boost exponential when T or kB is not positive.
- The SIGFPE continue path now demotes ARM FPCR trap-enable bits (Linux fpsimd_context and Apple neon state). enableFPE on Apple aarch64 sets `__fpcr`, not `__fpsr`. x86 MXCSR masking is unchanged.
- The client calls `getBundleSize()` again and enables bundling when a numbered bundle is present. A thrown `int` from a job is a failed exit, not a silent continue.
- The point job opens `pos.con` through `getRelevantFile`, so `pos_cp.con` / `pos_in.con` are used when present.
- The prefactor job writes `prefactor_reactant_to_product` and `prefactor_product_to_reactant`, honors `getPrefactors` failure, and does not assert on an empty moved set.
- Unbundled results name each slot `{state}_{first_wuid+slot}` so explorer
  does not add `result['number']` to recover the wuid. Reset paths share
  `fileio.remove_tree_and_empty_parents`.
- When a superbasin is active, `previous_state` stays the last KMC hop. Only `explore_state` is the lowest-confidence basin member.
- XTB sets periodicity from the box diagonal.
- `Matter::pbc` is the raw difference when periodic boundaries are off. OCINEB does not resample a non-periodic climbing image, so a 3-atom vacuum cell cannot fold to a 60 A reaction coordinate.
- `Parameters::load(FILE*)` returns 1 on a null FILE, a failed seek, or a negative `ftell` size.
- `Prefactor::getPrefactors` returns -1 if any of min1, saddle, or min2 is null.
- `WeightedSpring::compute` throws if the image index is 0 or past the spring table.
- `[ASE_ORCA] charge` is read from the INI and passed to `ase.Atoms`.
- ``SaddleSearchConfig.displace_atom_list`` accepts the scalar int ``-1``
  used by ``examples/akmc-al/run_akmc_al.py``, as well as ``[-1]`` and
  ``"-1"``.
- ``SaddleSearchJob`` and ``ProcessSearchJob`` apply
  ``client_displace_type = listed_atoms`` (and ``random`` / ``last_atom`` /
  ``least_coordinated``) through ``listedAtomEpiCenter`` plus a
  ``displace_radius`` / ``displace_magnitude`` kick instead of copying the
  reactant unchanged.
- ``listedAtomEpiCenter`` now picks with ``randomDouble(size)`` so the last
  listed (or last free, for lone ``-1``) atom can be chosen. The previous
  ``size-1`` interval always dropped that last index.
- `accuratePES` compares max `|pred-true|` and restores each image's incoming potential. The unused `sqrt(pred^2-true^2)` leftover is gone.
- `debug_keep_all_results` writes `results.dat` and con payloads from the result dict into `debug_results_path`.
- `eon_matter_to_atmconf` throws on a null Matter or when no atoms move.
- `eonclient -m` does not write a `.con` when no output path is given.
- `geometry::identical` will not assign two left atoms to the same right partner.
- `getMidSlice` keeps the CatLearn order (endpoints, then a two-thirds interior image) and throws if the path has fewer than three images. It is not `n/2`.
- `getRelevantFile` no longer throws on a name without a `.`. `_cp` / `_in` suffixes still win when those files exist.
- `gh release edit` failures are no longer ignored. A stale notes update now fails the release job.
- `job = dynamics` writes a `results.dat` envelope (`termination_reason`, `job_type`, `potential_energy`, `total_force_calls`) and still returns `final.con`.
- `job = structure_comparison` loads `matter1.con` and `matter2.con`, compares a copy, and writes `results.dat` (`match`, `distance`, `per_atom_norm`, energies).
- `min_mode_method=gprdimer` constructs `AtomicGPDimer` in place in the eigenmode variant. A WITH_GPRD Catch2 case covers that path.
- `removeNetForce` no longer runs on a one-atom Matter. Subtracting the mean force was identically zero and made a 1-atom NEB report immediate GOOD.
- `rgpot_pot` now compiles with the same `_args` as the other pot plugins (`WITH_RGPOT`). Dead Metatomic `n_avg` and the empty `.cargo/config.toml` are gone.
- `rotm`, `get_rotation_matrix`, and `internal_motion` call `numpy.cos` /
  `sin` / `arccos`. They used to call bare `cos`/`sin`/`acos` and raise
  NameError. `get_mappings` uses vector PBC distances.
- `unbundle` drops only those bundle slots that lack a `results*` file.
- `unbundle` only copies `*_N.ext` when `N` is all digits and matches the bundle number. Non-numeric suffixes such as `pos_final.con` are left alone.
- eonclib compiles after the Parameters accessor split: pot TUs include
  Parameters.h, write-holes use ParametersLoadAccess, and Dimer/ConFileIO
  qualify helpers and Matter.
- pyeonclient ConFrame export uses mkstemps (exclusive) and releases the GIL while writing. The guessable world-writable temp name is gone.


## [3.2.1](https://github.com/TheochemUI/eOn/tree/3.2.1) - 2026-09-13

### Added

- Exposed `opt_method`, `neb_opt_method`, `refine_opt_method`, and `refine_threshold` on pyeonclient `Parameters` so dimers and NEBs can pick FIRE/LBFGS/CG without a minimization-only path. ([#406](https://github.com/TheochemUI/eOn/issues/406))
- Added `PotType.EXPR` so rgpot 3.1 ExprPot can compose named kernels (`0.5*lj + d3`) from `[ExprPot]` expression and terms.
- Agreement tests for ``minimage`` and ``linkcell`` against eOn's
  numpy/C++ wrap oracles, LAMMPS ``minimum_image`` (ortho and
  restricted triclinic), and GROMACS ``pbc_dx``.
- Exposed rgpot 3.1 Grimme DFT-D3 and DFT-D4 as `PotType.DFTD3` / `DFTD4` (`potential = dftd3` / `dftd4`, `[D3Pot]` / `[D4Pot]`).
- Exposed rgpot 3.2 MOPACPot as `PotType.MOPAC` (`potential = mopac`, `[MOPACPot]`, Expr term `mopac`). Default model is AM1.
- Standalone doxyYoda C++ API HTML is published at `/api-cpp/` next to the Sphinx book.
- The C++ API landing page links to the Sphinx book. `pixi run -e docs makedocs` writes `/api-cpp/` next to the book.
- ``eon.geometry.pbc`` uses ``minimage`` when installed. ``neighbor_list_linkcell``
  compares the vesin neighbor list to ``linkcell.knearest`` on the same
  ``Structure``.

### Developer

- Pinned potentials-schema wrap to e825c207 so it matches rgpot 3.2.0's vendored Potentials.capnp.
- Updated wrap pins to current releases: rgpot v3.1.2, readcon-core v0.14.10, Highway 1.4.0, potentials-schema v1.15.1, xtb v6.7.1, plus current IRA and ARTn heads.

### Changed

- PotentialConfig now uses parseable pot tokens (socketnwchem, ase_pot, catlearn) and exposes emt_rasmussen and potentials_path. ([#414](https://github.com/TheochemUI/eOn/issues/414))
- Docs footer credits antics and loads only `antics.js`.
- Raised Python floors to readcon 0.14.9, vesin 0.6.1, and rgpot 3.1.2.

### Fixed

- The C++ API mainpage no longer prints undefined Makefile Doxygen aliases, and the theme header version matches 3.2.0. ([#401](https://github.com/TheochemUI/eOn/issues/401))
- Fixed IDPP, collective IDPP, and SIDPP path initialization moving frozen atoms when a mover sat next to a constraint. ([#410](https://github.com/TheochemUI/eOn/issues/410))
- Fixed ASE from_ase wrapping coordinates before PBC was applied, and stopped treating every get_indices constraint as FixAtoms. Positions views are not writeable in place. to_ase attaches SinglePoint energy/forces. NEB keeps the GIL when the pot is not thread-safe. ([#414](https://github.com/TheochemUI/eOn/issues/414))
- Write ``hessian.dat`` even when ``[Main] quiet = true``. Quiet still
  suppresses the log line; the Hessian job's artifact is no longer gated
  on it.
- ``LammpsLoader::require_loaded`` now says whether ``liblammps.so`` is
  missing, failed ``dlopen`` (glibc), or opened without
  ``lammps_open_no_mpi`` (the eOn plugin is not LAMMPS).
- ``displace_atom_list`` is CON file-order. ``ListedAtoms`` remaps those
  rows through the ``atom_id`` sort when the raw list is all frozen, so
  a movable-first active-volume ``.con`` no longer raises
  "Listed atoms are all frozen".



## [3.2.0](https://github.com/TheochemUI/eOn/tree/3.2.0) - 2026-08-16

### Added

- The basic-eon workflow runs the tests/ suite on Unix after the
  client install, including the `sh` package those process tests import.
  Job modules import `version` from the `eon` package so a source-tree
  checkout without generated `version.py` still dispatches.
- User guide for the AiiDA plugin (`pip install aiida-eon`) at
  `docs/source/user_guide/aiida.md`. The communicator page now links
  there.
- `RgpotAdapter` reports `caps().batched` and forwards a band of images to
  `rgpot::forceBatchImpl`, so a kernel that evaluates several systems in
  one call sees the whole band. The wrap tracks the rgpot branch that
  adds `ForceBatch` until a tag carries it.

### Changed

- Appending a frame to a `.con` no longer reads and rewrites the whole movie. `matter2con(filename, append=true)` serializes the new frame on its own and concatenates it, so writing an N-frame trajectory costs N frame writes instead of N(N+1)/2, and per-step movie writers such as `dynamics.con`, `movie.con`, and the basin hopping trial movies stop growing quadratically with step count. Output bytes are unchanged: a movie built by repeated appends is identical to the same frames written in one call. A target eOn did not write, or whose size or modification time moved since eOn wrote it, is still parsed once before anything is added, so an unparseable file yields an append error with its bytes intact. Gzip and zstd targets keep the read-and-rewrite path because a compressed member cannot be extended in place.
- C++ wrap and Python pins move to readcon 0.14.5 together, so the client
  can read the spec-3 files the 0.14 writer emits. x-only constraints
  round-trip on the client and the Python Structure writer (readcon-core #25).
  0.14.5 ships a win_amd64 wheel, so Windows CI does not build the sdist.
- Rewrite the Sphinx user guide, tutorials, install notes, developer
  docs, and older release pages so they no longer trip the house
  prose rules (robust, comprehensive, below is, out of the box,
  highlighting, first-class, ship-as-verb, and the rest of that list).
  Drop remaining X-not-Y contrast frames: state what each path is.

### Fixed

- ASE-NWChem and ASE-ORCA calculators each get a private work directory instead
  of `directory='.'`, so concurrent LocalInProcess jobs do not clobber each
  other's scratch files. Calculator errors throw rather than abort the process.
- Bond-boost hyperdynamics advances the equilibration counter once per MD step.
  `boost()` only evaluates the current bias, so ParallelReplica (which also
  installs the bias potential on the trajectory) and SafeHyper share the same
  `rmd_time` schedule.
- Constructing a `VASP` potential no longer deletes the working directory's contents. The constructor removed fifteen files, among them `WAVECAR`, `CHGCAR`, `TMPCAR` and `OUTCAR`, so a client restarted in a directory lost the wavefunction and charge density an `ISTART`/`ICHARG` restart reads, along with the previous run's record; and a potential built after the first force call, as `AtomicGPDimer` does, deleted `FU` and `NEWCAR` out from under a VASP process that was still running. A run now clears only the three handshake files it owns, `FU`, `NEWCAR` and `STOPCAR`, once, immediately before it starts VASP. `STOPCAR` was never removed although eOn writes it at shutdown, so a client restarted after a clean shutdown handed the new VASP an abort instruction on its first ionic step. Discarding the results and restart files is available as `VASP::removeStaleFiles()` for a caller that means to start from an empty directory; the `examples/akmc-vasp-slurm` scripts already do the same removals in shell.
- Copying a Matter leaves biasPotential null. Copy-assign used to skip
  that pointer, so getBiasForces read an indeterminate BondBoost and
  the Python suite died on main.
- Library code no longer calls `std::exit` on a bad `convergence_metric`, FIRE time-step collapse, or basin-hopping / global-optimization enum typo. Those errors throw, so a Python caller inside `gil_scoped_release` gets an exception instead of a silent interpreter kill. Unknown enumerated strings are rejected when `Parameters` is loaded.
- On Windows, ExtPot launches a suffix-less ext_pot wrapper with python.
  cmd.exe does not honor a shebang, so a quoted path to the script was
  not a runnable command.
- Potential loaders probe plugin libraries on disk before `dlopen`, and `LammpsLoader` no longer loads `liblammps` just to ask whether it is present. Availability checks therefore skip the banner-printing static initializers those libraries run on load.
- Setting ``params.write_con_forces`` on a pyeonclient Parameters object now writes force sections for that Matter. A ``ConFrameMetadata.write_con_forces`` value overrides the process-wide flag so two writers can disagree.
- The external-program potential backends (`AMS`, `AMS_IO`, `ExtPot`, `VASP`) now check every file open, read and shell command instead of trusting them. A missing, truncated or stale result file previously left the force and energy arrays holding whatever they held before, and those values went straight into the optimizer; they now raise an error naming the file and how many atoms were read. `AMS_IO` also truncates `ams_output` rather than appending to it, so a failed run cannot hand back the previous run's forces, and `VASP` rejects a structure whose species are not in contiguous runs rather than writing a `POSCAR` whose counts disagree with its coordinates.
- The readcon wrap overlay still cargo-builds both crate types. A generated
  .c that depends on that custom_target is the ninja order edge, so eonclib
  does not link before the outputs exist and does not record NEEDED
  libreadcon_core.so or pass the DLL to MSVC.
- `ExtPot` and `VASP` no longer claim to be safe to call from several threads on one instance. Both exchange structures and forces through files at fixed names, so the parallel image evaluation in `NudgedElasticBand`, `ImprovedDimer` and `ProcessSearchJob` had every thread writing one input file and reading one result file, and images took each other's forces. `ExtPot` also asks for a potential instance per image, which gives each thread an exchange directory of its own and keeps the evaluation parallel; `VASP` cannot, since every instance drives the same VASP process through the same files, so its image evaluation runs sequentially. The other external-program backends (`ASE`, `ASE_NWCHEM`, `ASE_ORCA`, `LAMMPS`) already declared both.
- `ExtPot` runs the external program in a private exchange directory named `extpot_<pid>_<n>` instead of the client's working directory. The exchange files keep their names, `from_eon_to_extpot` and `from_extpot_to_eon`, so a wrapper that opens them by relative name needs no change; two clients started in one directory can no longer read each other's structures and forces. A relative `ext_pot_path` that names an existing file, such as the default `./ext_pot`, is resolved before the command runs; a command line with arguments needs an absolute path to the script, and `EON_EXTPOT_RUN_DIR` names the directory eOn runs in for a wrapper that has to reach it. On Windows a resolved Python wrapper (shebang or `.py`) is invoked through `python.exe`, because `cmd.exe` will not run an extensionless file. The previous result file is removed before each call, so an external program that exits successfully without writing one raises an error rather than handing back the previous call's forces.
- `PotRegistry` is a process-lifetime heap singleton, so a `Potential` destroyed during interpreter finalization no longer calls into a destroyed registry or locks a destroyed mutex.
- `loadposcar` accepts both VASP 4 (integer counts after the cell, as the kdb tool writes) and VASP 5 (species names then counts, as `saveposcar` and the server movie files write). eOn can read back `movie.poscar`, `dynamics.poscar`, and the other `movie.py` outputs it wrote.
- `matter2xyz` writes extended XYZ: the comment carries `Lattice="..."` so a reader can reconstruct the cell, and coordinates use 17 significant digits to match the CON path. Appending a frame whose atom count differs from the last frame in the file is rejected and leaves the file unchanged.
- `min_mode_method = gprdimer` without `-Dwith_gprd=true` is an error.
  The search no longer falls through to ImprovedDimer and reports a
  nonnegative-mode abort.
- `savecon(..., w="a")` appends one serialized `.con` frame instead of reading the whole movie back and rewriting it. Gzip and zstd targets still rewrite, because a compressed member cannot be extended in place.
- pyeonclient Structure/Matter conversion looks up every element (H–Og) without importing the server package, and raises on an unknown atomic number instead of writing a fabricated ``Z79`` symbol.


## [3.1.0](https://github.com/TheochemUI/eOn/tree/3.1.0) - 2026-08-08

### Changed

- pyeonclient 0.4.0: the potentials it exposes come from ``librgpot`` now
  rather than eOn's in-tree kernels, so the wheels bundle that library and
  no longer carry per-pot plugin objects. The Python API is unchanged --
  the same ``potential`` names select the same physics.

  Its wheel workflow gained a ``target`` input so an upload can be
  rehearsed against TestPyPI before the irreversible one; a
  ``pyeonclient-v*`` tag still goes straight to PyPI.

### Fixed

- Fat ``MetatomicPotential`` loads exported PET-MAD (and similar) models on
  metatensor-torch 0.10.3: scripted modules without ``_mts_buffer_names`` no
  longer trip the mixed-dict ``.to()`` walk.


## [3.0.0](https://github.com/TheochemUI/eOn/tree/3.0.0) - 2026-07-26

### Added

- Continuous ASV dashboard: results history on the `asv-results` orphan branch, asv-tachyon UI under `gh-pages/bench/`, optional Netlify deploy for bench.eondocs.org.
  Windows CI: robust MSVC activation (VS 18 runners) and strip GNU link.exe from PATH for meson test. ([#378](https://github.com/TheochemUI/eOn/issues/378))
- PHVA mobile/active sets for dense Hessian and matrix-free min-mode via
  ``phva_atoms`` (default ``All`` = all free atoms). Free/fixed stays the
  optimizer mask; Krylov dimension is ``3 * N_active``. Shared
  ``resolveMobileAtoms`` drives ``HessianJob``, Lanczos, and Davidson.
  pyeonclient: ``Lanczos``/``Davidson.compute(..., atoms=)``,
  ``resolve_mobile_atoms``, ``free_atom_indices``, and
  ``Parameters.hessian_phva_atoms`` / ``lanczos_phva_atoms`` /
  ``davidson_phva_atoms``. INI keys: ``[Hessian|Lanczos|Davidson] phva_atoms``. ([#379](https://github.com/TheochemUI/eOn/issues/379))
- Cap'n Proto ``schema/eon_job_result.capnp`` defines typed ``JobRequest`` /
  ``JobResult`` / flat ``Geometry`` envelopes for the in-process control plane
  (kill-file-IPC). ``eon_schema.jobs`` provides ``results.dat`` adapters for
  legacy paths. ([#381](https://github.com/TheochemUI/eOn/issues/381))
- Vendored vesin updates to 0.6.0 with three local extensions carried as
  upstream candidates: a CPU implementation of the declared-but-missing
  ``VesinBruteForce`` algorithm (nearest-image MIC pair search), lazy
  thread-pool worker spawn (serial consumers stop paying
  ``hardware_concurrency()`` thread creations per process), and a fused pair
  visitation API (``vesin_neighbors_visit`` plus the header-only
  ``vesin_visit.hpp``).

  Every Fortran-backed potential draws its neighbours from vesin rather than
  from pot-local scaffolding: the EDIP and Lenosky cell/ghost machinery is
  gone, Tersoff trades its O(N^3) sweeps for per-atom lists, and SW's silent
  ``MAXNEI`` overflow and FeHe's 800-neighbour gather ceiling are hard errors
  instead of quiet truncation. Those kernels moved to rgpot in the same
  release (see below), so the vesin Fortran interface is vendored there
  rather than here; eOn's vendored copy is the C++ translation unit only. ([#389](https://github.com/TheochemUI/eOn/issues/389))
- Documentation systems tutorials (Morse Pt NEB, LJ minimization, Pt saddle) with
  built-in potentials and current `rgpycrumbs` / `plt-neb` / `plt-min` conventions
  (1:1 reaction-valley landscapes, full structure strips, one min landscape per
  endpoint).
- Optimized hyperplanar TST (OH-TST) job (Johannesson and Jonsson, J. Chem. Phys.
  115, 9644 (2001)): a `job = oh_tst` mode that progresses a hyperplanar dividing
  surface by reversible work with thermostatted, plane-constrained sampling, and
  reports the free-energy barrier and crossing rate. Supports Andersen or GLE
  colored-noise (Ceriotti-Bussi-Parrinello) thermostats and symmetry-restricted
  sampling of equivalent product minima. Configured under `[OH_TST]`.

### Developer

- Golden masters for CPython server atoms helpers vs pre-#368 scalar code; eon.fileio load/save/round-trip via live readcon fixtures; document chemfiles as optional and unused on server .con path. ([#370](https://github.com/TheochemUI/eOn/issues/370))
- Windows metatomic CI follows the metatomic/metatensor torch workflow:
  windows-2022, setup-python, pip CPU torch (PIP_EXTRA_INDEX_URL), and MSVC,
  with pixi for C++ deps and flang. In-tree Fortran pots and CuH2 stay enabled
  via feedstock-style flang_rt LIBPATH and MSVC AR=lib. Catch2 runs on the basic
  multi-OS matrix including windows-2022; ConFileIO tests close temp streams
  before remove (Windows file locks); inih example tests strip CR for MSVC
  text-mode stdout. ([#377](https://github.com/TheochemUI/eOn/issues/377))
- The Windows CI step that used to load ``eon_sw.dll`` and assert it exports a
  bare ``sw_`` now audits ``librgpot`` instead: the kernels ship inside it with
  hidden linkage, so the check is that *no* legacy Fortran name reaches the
  export table, which is the PE counterpart of the Linux ``FortranSymbolAudit``.
  It accepts a ``--default-library=static`` build, where the kernels land in an
  archive that has no export table at all. ([#390](https://github.com/TheochemUI/eOn/issues/390))

### Changed

- Client public headers live under ``include/eon/`` (numpy/fmt layout). Sources
  stay in ``client/``. Use ``#include "eon/Potential.h"`` (and
  ``eon/fpe_handler.h``, ``eon/potentials/...``) with ``-I$prefix/include``.
  Relative ``../`` includes are gone; headers install via ``install_subdir``. ([#379](https://github.com/TheochemUI/eOn/issues/379))
- Classical C++ pair pots (LJ, Morse, LJCluster, QSC) use the shared
  ``eonc::VesinNeighbors`` wrapper for neighbor lists instead of pot-local
  Verlet / O(N^2) loops. Vesin is always linked into ``eoncbase`` so Metatomic
  and classical pots share one NL backend. Client wall time continues to be
  tracked by the existing ASV suite (``TimePointMorsePt``,
  ``TimeMinimizationLJCluster``, Morse saddle/NEB); see ``benchmarks/README.md``. ([#386](https://github.com/TheochemUI/eOn/issues/386))
- The Fortran-backed potentials (Stillinger-Weber, EDIP, Lenosky, Tersoff,
  EAM aluminium, FeHe, CuH2, and TIP4P-H) are Fortran 2018 kernels inside
  ``librgpot`` and evaluate through ``RgpotAdapter`` like the classical
  pots. Each kernel was rewritten rather than wrapped: modules with
  ``implicit none``, kinds from ``iso_fortran_env`` checked against the C
  types at compile time, derived-type parameters in place of COMMON
  blocks, ``intent`` on every argument, ``pure`` kernels, structured
  control flow, and status returns instead of ``stop``. Neighbours come
  from vesin, and the pair sums are restated as gathers so the atom loops
  run under ``do concurrent``.

  eOn no longer builds, installs, or dlopens Fortran: the ``eon_*.so``
  plugin modules, their Windows ``.def`` export files, and the flang
  runtime handling are gone. ``FortranPotLoader`` becomes ``PluginLoader``,
  which still finds engine plugins (the rgpot metatomic and xtb backends)
  across ``EON_POTENTIALS_PATH`` and ``[Potential] potentials_path``.

  Configuration is unchanged: the same ``potential`` names select the same
  physics, pinned by the existing reference energies in ``SiPotTest``,
  ``EAMAlTest``, ``FeHeTest``, and ``cuh2Test``. ([#390](https://github.com/TheochemUI/eOn/issues/390))
- ``[Hessian] phva_atoms`` replaces ``atom_list`` for the PHVA mobile set
  (default still ``All``). The C++ client still accepts the legacy
  ``atom_list`` key when ``phva_atoms`` is absent. Lanczos and Davidson
  gain the same ``phva_atoms`` key in INI, schema, and ``config.yaml``.

### Fixed

- eOn no longer requires a Fortran compiler. It compiles no Fortran since the
  kernels moved into ``librgpot``, but ``client/meson.build`` still asked for
  the language with ``required: true`` whenever ``with_fortran`` or
  ``with_cuh2`` was set -- both default on -- so every consumer had to supply
  an unused toolchain. That bit hardest on builds against an *installed*
  rgpot, which need none at all. The language is detected rather than
  required now; a wrap build still finds it for rgpot's kernels, and the
  Windows flang runtime discovery is gated on actually having one.

  ``with_fortran`` and ``with_cuh2`` keep their meaning as potential
  selectors: they gate the ``RgpotAdapter`` arms, not any compilation.
- Local communicator writes client stderr to ``stderr.dat`` (no undrained
  ``PIPE`` deadlock). The FPE continue handler masks the fault class in the
  restored MXCSR so a single divide-by-zero cannot re-storm. The LAMMPS pot
  worker demotes floating-point traps after fork and ``forceLocal`` uses
  ``eat_fpe`` (same external-pot contract as ASE/Metatomic); SafeMath guards
  remain on CG/dimer/min-mode bare divisions. ([#379](https://github.com/TheochemUI/eOn/issues/379))
- Server-side process search and superbasin amsel gating receive a ``ConfigClass``
  (no bare ``config`` / missing ``self.config``). Superbasin recycling passes
  config into ``Recycling``. Restored ``atoms.identical`` for
  indistinguishable-atom matching. ``LocalInProcess`` unpacks
  ``Matter.relax`` as ``(Matter, converged)``. LAMMPS ``forceLocal`` restores
  FE traps after ``eat_fpe``. Tip4p implements the full ``Potential::force``
  signature (variance parameter). ([#380](https://github.com/TheochemUI/eOn/issues/380))
- ClientEON and pyeonclient ``append_results_timing`` write the
  ``results.dat`` timing footer as ``<value> <key>`` (same contract as
  job writers and ``parse_results``), so ``time_seconds`` /
  ``user_time`` / ``system_time`` parse as floats under those keys. ([#382](https://github.com/TheochemUI/eOn/issues/382))
- AMS_IO writes the XC functional name into the run script: the
  ``fprintf("xc %s\\n")`` call was missing its ``xc`` argument (undefined
  behavior). ([#384](https://github.com/TheochemUI/eOn/issues/384))
- ASE potential import and force failures throw ``std::runtime_error``
  instead of ``exit(1)``, so pyeonclient / in-process callers can recover
  instead of killing the whole process. ([#385](https://github.com/TheochemUI/eOn/issues/385))
- LAMMPS worker path survives a bad geometry: non-finite or failed evaluations
  reject the structure without ending the client, with bounded respawns,
  SIGPIPE ignored on the pipe, serialised concurrent exchanges, and process
  search reporting of how far a rejected endpoint landed from the reactant. ([#387](https://github.com/TheochemUI/eOn/issues/387))
- Classical pair pots (LJ, Morse, LJCluster, QSC) use ``eonc::PairListCache``:
  a process-global pool of Verlet-skin cached pair lists. The candidate list
  builds at ``cutoff + skin`` once; force evaluations on geometries whose atoms
  have moved less than ``skin/2`` since the build re-use the cached pairs and
  derive exact vectors from current positions, with the true-cutoff filter
  keeping results identical to a fresh build. On a cache miss in the MIC regime
  the pair kernel inlines into the build's single brute-force scan, so one-shot
  evaluations (point jobs) pay one pair sweep like the pre-list code did.
  Proximity-matched pool slots keep NEB's per-image force path both race-free
  and cache-warm, whether images are evaluated serially or in parallel — the
  pool survives NEB's per-iteration worker threads. Fixes the ASV regressions
  from the #386 vesin port (point / min / saddle / NEB wall times); minimization,
  saddle-search, and NEB fixtures run faster than the pre-#386 baseline. ([#389](https://github.com/TheochemUI/eOn/issues/389))
- A build that resolves rgpot through the subproject wrap now takes vesin
  from rgpot instead of compiling eOn's vendored copy beside it. Both trees
  carry the same upstream release, but rgpot's carries local patches its
  Fortran interface binds to, so linking both left those objects calling a
  symbol eOn's copy did not define. The vendored translation unit still
  serves builds against an installed rgpot.

  pyeonclient wheels need no shared object beside ``librgpot``: the Cap'n
  Proto schema moved inside it, so importing ``pyeonclient._core`` no longer
  depends on finding a separate ``libptlrpc.so`` at run time. ([#390](https://github.com/TheochemUI/eOn/issues/390))
- Process registration no longer freezes the process-id counter when the
  ``processtable`` contains duplicate ids. New process ids are
  content-addressed with xxHash (``allocate_process_id`` / xxh64 over the
  saddle payload and barrier), not ``len(procs)`` or ``max(id)+1``, so
  procdata files are not overwritten. Duplicate appends raise.


## [2.17.10](https://github.com/TheochemUI/eOn/tree/2.17.10) - 2026-07-20

### Fixed

- Optimizer file logs (``_lbfgs.log`` etc.) use a durable temp directory so
  fixture workdir cleanup no longer races Quill file sinks.
- Packaging: Catch2 ``--allow-running-no-tests`` for the optional rgpot embed
  suite so all-SKIP without nwchemc/cpmdc is a clean success.


## [2.17.9](https://github.com/TheochemUI/eOn/tree/2.17.9) - 2026-07-20

### Fixed

- EDIP OpenMP: zero shared energy/force accumulators under ``!$omp single``
  before partial reduction (fixes wrong PointJob energies with multi-thread OMP).
- Basin-hopping force-call reference updated for default LBFGS auto_scale path
  (1692); energy and acceptance assertions unchanged.


## [2.17.8](https://github.com/TheochemUI/eOn/tree/2.17.8) - 2026-07-20

### Fixed

- Evaluate EDIP serially under OpenMP (Fortran energy reduction race under
  OMP_NUM_THREADS>1 gave wrong PointJob energies with plausible forces).
- Basin-hopping integration test pins LBFGS auto_scale off so force-call
  budget matches the SVN reference path.


## [2.17.7](https://github.com/TheochemUI/eOn/tree/2.17.7) - 2026-07-20

### Fixed

- Un-nest SW CG minimization TEST_CASE in SiPotTest (was inside LBFGS case after
  packaging SKIP refactor, breaking with_tests packaging builds).


## [2.17.6](https://github.com/TheochemUI/eOn/tree/2.17.6) - 2026-07-20

### Fixed

- Keep JobIntegrationFixture::runJob a class method after SKIP macro refactor
  (stray brace broke TEST_CASE_METHOD inheritance under packaging builds).


## [2.17.5](https://github.com/TheochemUI/eOn/tree/2.17.5) - 2026-07-20

### Fixed

- Catch2 SKIP only via macros in job integration tests (helpers cannot call SKIP).


## [2.17.4](https://github.com/TheochemUI/eOn/tree/2.17.4) - 2026-07-20

### Fixed

- Fix Catch2 SKIP usage in unit-test helpers (macros / fixture require path)
so packaging builds compile with ``-Dwith_tests=true``.


## [2.17.3](https://github.com/TheochemUI/eOn/tree/2.17.3) - 2026-07-20

### Fixed

- Pin metatomic builds to vesin>=0.6 and rgpot>=2.5.2 so RGPOT
  ``libmetatomic_engine`` and fat Metatomic share a matching VesinOptions ABI
  (skin/n_threads). Bump the rgpot wrap to the 2.5.2 fix commit.
- Unit tests SKIP when packaging fixtures or optional engines are missing
  (``EON_TEST_SYSTEMS_DIR`` / ``EON_POTENTIALS_PATH`` / nwchemc·cpmdc), so
  ``meson test`` succeeds for installable client builds without local test data.


## [2.17.2](https://github.com/TheochemUI/eOn/tree/2.17.2) - 2026-07-17

### Added

- pyeonclient 0.3.3: in-memory ``path_frames`` / ``to_conframes`` for NEB (same stamps as writePathCon), ``Matter.relax(retain_frames=…)`` movie frames, and ``MinModeSaddleSearch.run_retain_frames`` climb frames. ([#pyeonclient-path-frames](https://github.com/TheochemUI/eOn/issues/pyeonclient-path-frames))
- Dimer supports `rotation_backend = classical|lanczos|davidson|lor` (Leng et al. JCP 2013 LOR for softest-mode with force translation; Lanczos/Davidson FD min-mode; classical constrained rotation). Climb remains on the dimer path.
- Optional Superbasin gate via amsel.discover_decide_status ([amsel] / amsel_discover_decide) before MCAMC; rejected_no_metastable_basin falls back to AKMC.

### Fixed

- HessianJob FD eigen solve: ColMajor copy for SelfAdjointEigenSolver (avoids
  RowMajor segfaults on partial VTST blocks), finite-force/index guards, optional
  `[Hessian] atom_list` intersected with free atoms, and a stable mobile VectorXi
  copy so moving-atom Hessians match PHVA-class active sets. ([#357](https://github.com/TheochemUI/eOn/issues/357))
- Release CI installs quill/capnp so ``eon-akmc`` sdist and pyeonclient wheels configure cleanly on GitHub runners. ([#418](https://github.com/TheochemUI/eOn/issues/418))
- Republish pyeonclient as **0.3.4** with the full wheel matrix (abi3 + freethreading) after incomplete 0.3.3 tag CI.


## [2.17.1](https://github.com/TheochemUI/eOn/tree/2.17.1) - 2026-07-17

### Fixed

- Restore ``with_pyeonclient`` and ``install_eon_server`` meson options dropped
  from ``meson_options.txt`` in the v2.17.0 merge (configure failed with
  ``Option with_pyeonclient does not exist``).
- RGPotEngine accepts rgpot >= 2.5 ``operator()`` return type
  ``std::tuple<double, AtomMatrix, double>`` (energy, forces, variance).


## [2.17.0](https://github.com/TheochemUI/eOn/tree/2.17.0) - 2026-07-17

### Added

- Distribution packaging: prefer installed ``rgpot`` via pkg-config (wrap pin
  ``v2.5.0``) and installed ``nlohmann_json`` (module with wrap fallback, same
  pattern as readcon). RGPOT ``backend=metatomic`` / ``backend=xtb`` dlopen
  ``libmetatomic_engine.so`` / ``libxtb_engine.so`` so packaging keeps
  ``-Dwith_xtb=false`` (native XTBPot deprecation warning) and avoids fat
  metatomic links into eOn.
- SafeMath accepts Eigen 5 ``EIGEN_CORE_MODULE_H`` as well as Eigen < 5
  ``EIGEN_CORE_H``.

- Add ``eon_schema.config`` INI helpers (``write_ini``, ``hydrate_ini``,
  ``unknown_ini_keys``, ``write_models_ini``) so tooling can author validated
  ``config.ini`` from L0/L1 without importing eon-akmc. ([#eon-schema-ini](https://github.com/TheochemUI/eOn/issues/eon-schema-ini))
- Parameter field graph for Main, Potential, Optimizer, Structure Comparison, and Process Search is authored in ``schema/eon_params.capnp`` (Cap'n Proto L0 SSoT) with INI/JSON adapters; parity tests gate ``config.yaml`` / ``eon.schema`` against the catalog. ([#params-ssot](https://github.com/TheochemUI/eOn/issues/params-ssot))
- Public step composition for Matter-first NEB/min: ``write_neb_results``,
  ``pot_registry_total_force_calls`` at package root; ``NEB.find_extrema``;
  GIL release on ``forces_free`` / ``max_force``. NebSpec covers EW/CI/OCI-MMF. ([#pyeonclient-step-compose](https://github.com/TheochemUI/eOn/issues/pyeonclient-step-compose))
- Add packages/eon-schema (PyPI eon-schema): vendored Cap'n Proto SSoT copy, optional pydantic API models; full-tree/conda-forge eon path unchanged. ([#eon-schema-0.1](https://github.com/TheochemUI/eOn/issues/eon-schema-0.1))
- Add pyeonclient nanobind module (Matter/Parameters/Potential); stable ABI abi3 on CPython 3.12+, free-threaded when Py_GIL_DISABLED; no pybind11. ([#372](https://github.com/TheochemUI/eOn/issues/372))
- pyeonclient Matter end-to-end: LocalInProcess communicator (comm_type local_lib/inprocess); NbGuard without pybind11; ASE embed polarity documented. ([#373](https://github.com/TheochemUI/eOn/issues/373))
- Expand pyeonclient: full enums, Job/make_job/run_job_in_directory, Potential.get_ef, ASE to_ase/from_ase and Structure helpers. ([#374](https://github.com/TheochemUI/eOn/issues/374))
- Standalone PyPI package pyeonclient (pyproject-pyeonclient.toml): abi3 + cp313t + cp314t wheels via pyeonclient-wheels.yml. ([#375](https://github.com/TheochemUI/eOn/issues/375))
- RGPOT backend ``metatomic`` loads ``libmetatomic_engine.so`` (thin host path).
  Native ``potential = Metatomic`` is unchanged for conda-forge packaging. ([#377](https://github.com/TheochemUI/eOn/issues/377))
- Document fat vs RGPOT-dlopen vs ASE metatomic backends; benchmark figure is
  generated from committed JSON at ``sphinx-build`` time (not checked in).
- Add ``pixi`` environment ``mta-bench`` and ``mta-backend-bench`` task to
  reproducibly build fat/ASE pyeonclient, refresh the metatomic backend compare
  JSON, and regenerate the docs figure.

### Developer

- Multi-flag Codecov coverage for Python, C++, and Fortran (OIDC uploads). ([#367](https://github.com/TheochemUI/eOn/issues/367))
- CPython server entry points fixed (`eon-server` → `eon.server:main`); vectorized atoms neighbor/free-atom hot paths. ([#368](https://github.com/TheochemUI/eOn/issues/368))

### Changed

- Move full job-config pydantic models into ``eon-schema`` (``eon_schema.config``).
  ``eon.schema`` and ``pyeonclient.models`` re-export the shared package so both
  eon-akmc and pyeonclient share one schema surface. ([#eon-schema-l1](https://github.com/TheochemUI/eOn/issues/eon-schema-l1))
- Structure is readcon-backed (ConFrame bridge); geometry PBC/NL via vesin; process-atom selection uses vesin shells. Retires aselite-era atoms container as storage source of truth. ([#371](https://github.com/TheochemUI/eOn/issues/371))

### Fixed

- Release the GIL around Matter energy/force accessors and NEB construction,
  ``update_forces``, and per-image energy reads so metatomic/torch autograd can
  run from C++ without deadlocking. Optional ``torch/cuda.h`` and ``torch/mps.h``
  includes (``__has_include``) so CPU-only pip torch builds compile without
  CUDA/MPS headers. ([#gil-torch-autograd](https://github.com/TheochemUI/eOn/issues/gil-torch-autograd))
- get_process_atoms returns plain Python ints so recycling metadata repr/eval round-trips (no np.int64). ([#369](https://github.com/TheochemUI/eOn/issues/369))
- Standalone pyeonclient wheels omit eon server package and eonclient binary (install_eon_server=false). ([#376](https://github.com/TheochemUI/eOn/issues/376))


## [2.16.0](https://github.com/TheochemUI/eOn/tree/2.16.0) - 2026-07-03

### Added

- Matter exposes `PbcConvention` (Legacy vs MinimumImage) for position wrapping, selectable without merging the historical `v3c_pbcs` branch (issue #176). ([#176](https://github.com/TheochemUI/eOn/issues/176))
- ASE NWChem calculator takes `mpi_launcher` (`mpirun` default or `srun`) for Slurm-friendly execution (issue #193). ([#193](https://github.com/TheochemUI/eOn/issues/193))
- Metatomic accepts explicit `energy_output` / `energy_uncertainty_output` keys when models use non-default output names (issue #215). ([#215](https://github.com/TheochemUI/eOn/issues/215))
- Metatomic can apply a random SO(3) rotation per evaluation with
  `random_rotation` (issue #287), rotating forces back to the lab frame. ([#287](https://github.com/TheochemUI/eOn/issues/287))
- Metatomic can average energy and forces over `n_symmetry_rotations` random
  rotations for approximate O(3) symmetrization (issue #292). ([#292](https://github.com/TheochemUI/eOn/issues/292))
- Metatomic supports non-conservative forces via `non_conservative` and
  `variant_force` / `force_output` (issue #296), matching metatomic-ase. ([#296](https://github.com/TheochemUI/eOn/issues/296))
- In-process RgpotPot for NWChem/CPMD via rgpot (optional build); user guide and examples for BLYP/DFT XC blocks.
- Metatomic exposes `deterministic` and `deterministic_strict` knobs for
  reproducible PyTorch evaluation (JIT profiling off, deterministic algorithms).
- Minimum-mode finding supports `min_mode_method = davidson`: Ritz subspace with FD Hessian-vector products and residual correction (alternative to dimer rotation minimization and Lanczos).
- The `rgpot` potential is selectable from the Python configuration layer: `eon/config.yaml` gains an `[RgpotPot]` section and `eon.schema` a matching `RgpotPot` model, mirroring the C++ INI keys (backend, basis, theory, scf_type, functional, cutoff_ry, charge, multiplicity, engine paths, title, memory_mb, scratch_dir, input_block).

### Changed

- `Matter::getForcesFree` / `getForcesFreeV` are const, and `setPositionsFree`
  applies PBC wrapping like `setPositions` so free-atom optimizers stay consistent
  for TIP4P/SPCE and other PBC pots (issue #171). ([#171](https://github.com/TheochemUI/eOn/issues/171))
- Metatomic potential requests per-atom outputs via ``sample_kind == "atom"``
  instead of the removed ``get_per_atom`` / ``set_per_atom`` API. ([#362](https://github.com/TheochemUI/eOn/issues/362))
- Build against metatomic-torch >=0.1.15 and metatensor-torch >=0.10, including
  the renamed pip layout (`metatensor_torch/`) and `sample_kind` on ModelOutput.
- Morse pair force hot path inlines the potential, uses inv-r scaling, and accumulates energy in a register (cachegrind: ~60% Ir in Morse::force on NEB Morse workloads).
- `ConFileIO` targets readcon-core **0.13** with a nanobind-ready `IoStatus` enum
  (replacing `bool` returns), bulk/zero-copy builder maps
  (`positions_data` / `set_*_from_flat` / `set_atom_velocity`), and
  `writeNebPath` using `ConFrameBuilder::clone()` so NEB bands write in one pass
  without re-reading the movie file per image. Matter exposes atom-index and CON
  header accessors for bindings.

### Fixed

- Min-mode saddle search no longer reports success (`STATUS_GOOD`) when the climb
  objective is still unconverged (issue #20). ([#20](https://github.com/TheochemUI/eOn/issues/20))
- Saddle and process search synthesize a missing `displacement.con` from `pos.con` plus `saddle_search.displace_magnitude` times a unit `direction.dat` mode (issue #79). ([#79](https://github.com/TheochemUI/eOn/issues/79))
- Molecular QM potentials (NWChem socket, ASE-ORCA, ASE-NWChem) auto-disable PBC on
  `Matter` attach and hard-fail if PBC is re-enabled, preventing wraps from tearing
  non-centered molecules (issue #188). ([#188](https://github.com/TheochemUI/eOn/issues/188))
- MCAMC `energy_level` superbasin scheme tracks a true run-wide minimum energy for increment scaling and no longer passes a generator into `min()` (issue #212). ([#212](https://github.com/TheochemUI/eOn/issues/212))
- NEB with `initializer = file` and `initial_path_in` loads endpoints from the first and last path frames and no longer requires (or crashes on) separate `reactant.con` / `product.con` (issue #278). ([#278](https://github.com/TheochemUI/eOn/issues/278))
- .con outputs are readable by classic con parsers again: "Forces of Component" sections are opt-in via `[Main] write_con_forces` (default off) instead of always written. ASE's eon reader rejects frames carrying force sections, which broke every ASE-based consumer of eOn-written con files (including the atomistic-cookbook eon-pet-neb example). Enable the knob to keep force+energy co-loading for warm restarts; reading force-bearing frames works regardless.
- RgpotPot engine-path environment overrides are backend-scoped: `NWCHEMC_LIBRARY` / `RGPOT_NWCHEMC_ENGINE` apply only to the NWChem backend and `CPMDC_LIBRARY` / `RGPOT_CPMDC_ENGINE` only to the CPMD backend, so setting both no longer sends the NWChem engine library into a CPMD configure.
- The `examples/rgpot_cpmd_blyp` example and the rgpot-vs-SocketNWChem comparison scripts now use the shipped in-process `[RgpotPot]` configuration (backend / functional / cutoff_ry / engine_library) instead of the abandoned potserv host/port keys, which eonclient never read.


## [2.15.0](https://github.com/TheochemUI/eOn/tree/2.15.0) - 2026-06-24

### Added

- MPI build revived on the C API (`with_mpi`) after retiring the obsolete C++ MPI bindings path. ([#339](https://github.com/TheochemUI/eOn/issues/339))
- Runtime-loadable Fortran potentials via `dlopen` and `[Potential] potentials_path` (no rebuild to swap pot modules). ([#342](https://github.com/TheochemUI/eOn/issues/342))

### Developer

- Documented and automated the full maintainer release path (cog/towncrier lockstep assert, GitHub tarball release, PyPI, conda-forge/eon-feedstock checklist, incomplete-release recovery). ([#343](https://github.com/TheochemUI/eOn/issues/343))

### Fixed

- Windows/conda-forge builds enable in-tree Fortran including CuH2 (m2w64 gfortran + MSVC C++, 16 MiB stack, xtb MinGW import lib) per feedstock#15. ([#15](https://github.com/TheochemUI/eOn/issues/15))
- NEB isolates per-image LAMMPS instances on private MPI communicators so parallel image evaluation does not corrupt neighbor state. ([#340](https://github.com/TheochemUI/eOn/issues/340))
- Process search restores fixed-atom rows after loading `displacement.con`, avoiding drift of constrained atoms. ([#341](https://github.com/TheochemUI/eOn/issues/341))


## [2.14.0](https://github.com/TheochemUI/eOn/tree/2.14.0) - 2026-04-24

### Added

- Structured per-iteration metadata is now embedded directly in minimization, saddle-search, and NEB movie `.con` files via `readcon-core`. The legacy minimization and saddle-search sidecar tables remain available temporarily behind `write_deprecated_outs = true`. ([#trajectory-output](https://github.com/TheochemUI/eOn/issues/trajectory-output))
- Added the `filin` option under `[ARTn]` to expose pARTn's input-file name. Empty (the default) leaves pARTn with its post-`artn_create` sentinel so no file is read; setting it to a path wires that file into pARTn setup and aborts early with a clear error if it is missing. ([#334](https://github.com/TheochemUI/eOn/issues/334))

### Fixed

- Fixed a regression in threaded force-evaluation paths by keeping legacy Fortran-backed potentials out of shared-instance parallel execution. ([#332](https://github.com/TheochemUI/eOn/issues/332))


## [2.13.0](https://github.com/TheochemUI/eOn/tree/2.13.0) - 2026-04-18

### Removed

- Removed 9 dead legacy pointer+size math functions from HelperFunctions (``dot``, ``length``, ``add``, ``subtract``, ``multiplyScalar``, ``divideScalar``, ``copyRightIntoLeft``, ``normalize``, ``makeProjection``), all superseded by Eigen. ([#dead-code](https://github.com/TheochemUI/eOn/issues/dead-code))
- Removed dead code: INIFile.cpp/h (superseded by inih), legacy CMakeLists.txt
  (36 files), Makefile/buildRules.mk, old unittest framework (unittests/).
  Total: ~2600 lines of dead code removed. ([#dead-code-removal](https://github.com/TheochemUI/eOn/issues/dead-code-removal))
- Removed dead ``TADJob::saddleSearch()`` function and associated ``dimerSearch`` member (dead since 2013). ([#dead-tad](https://github.com/TheochemUI/eOn/issues/dead-tad))

### Added

- Release builds now use ``-O3 -flto=auto``. Optional ``native_arch`` and ``fast_math`` meson options for ``-march=native`` and ``-ffast-math``. ([#compiler-flags](https://github.com/TheochemUI/eOn/issues/compiler-flags))
- Added ARTn (Activation-Relaxation Technique nouveau) saddle search method via
  the pARTn Fortran library, as a complementary explorer for pushing away from
  minima. Supports two modes: ``method = artn`` (standalone, pARTn drives push +
  internal eigenmode search + relaxation) and ``min_mode_method = artn``
  (drop-in, eOn's displacement seeds the mode). Configurable via ``[ARTn]``
  section with ``push_step_size``, ``force_threshold``, ``max_iterations``,
  ``ninit``, ``nperp_limitation``, ``lanczos_min_size``, ``nsmooth``, and
  ``nnewchance`` parameters. Requires ``-Dwith_artn=true`` at build time.
  For automated saddle point refinement in production workflows, prefer OCINEB
  (``ci_mmf = true``). ([#feat-artn-integration](https://github.com/TheochemUI/eOn/issues/feat-artn-integration))
- Added IRA (Iterative Rotations and Assignments) structure comparison via the
  libira Fortran library. Provides ``IRACompare::match()`` (CShDA + SVD
  alignment), ``matchPBC()`` (periodic boundary conditions), and
  ``findSymmetry()`` (SOFI point group detection). Requires ``-Dwith_ira=true``
  at build time. ([#feat-ira-integration](https://github.com/TheochemUI/eOn/issues/feat-ira-integration))
- Added Highway SIMD as optional subproject (cmake wrap). When available,
  ``-DWITH_HIGHWAY`` is set so future potential kernels can opt in.
  Hand-written SIMD kernels for Morse/LJ/EAM are staged for a follow-up
  release (see the separate ``highway-simd-potentials`` developer note);
  today the main benefit of enabling Highway is compile-time availability
  of the subproject for those downstream kernels. ([#highway-simd](https://github.com/TheochemUI/eOn/issues/highway-simd))
- Highway SIMD subproject builds and is detected at configure time. SIMD
  kernels for Morse/LJ/EAM force loops are planned but not yet implemented
  (core algorithm files already benefit from Eigen's auto-vectorization). ([#highway-simd-potentials](https://github.com/TheochemUI/eOn/issues/highway-simd-potentials))
- IDPP (Image Dependent Pair Potential) path initialization for NEB. ([#neb-idpp](https://github.com/TheochemUI/eOn/issues/neb-idpp))
- NEB decomposed into modular strategy-pattern components: tangent, projection, spring force, OCINEB controller, spline extrema, initial paths, and objective function. ([#neb-modularize](https://github.com/TheochemUI/eOn/issues/neb-modularize))
- OCINEB (Off-Path Climbing Image NEB): recommended hybrid CI-NEB + Min-Mode Following with Hessian eigenmode alignment for automated saddle point refinement (``ci_mmf = true``). See Goswami, Gunde, Jónsson, *Enhanced Climbing Image Nudged Elastic Band Method with Hessian Eigenmode Alignment*, 2026, [arXiv:2601.12630](https://arxiv.org/abs/2601.12630). ([#neb-ocineb](https://github.com/TheochemUI/eOn/issues/neb-ocineb))
- Onsager-Machlup action-based NEB for minimum action paths (``onsager_machlup = true``). ([#neb-om](https://github.com/TheochemUI/eOn/issues/neb-om))
- Parallel image force evaluation for NEB (requires TBB, ``-Dwith_parallel_neb=true``). ([#neb-parallel](https://github.com/TheochemUI/eOn/issues/neb-parallel))
- Parallel improved dimer gradient evaluation via ``std::thread``. The two dimer replicas evaluate forces concurrently when the potential is thread-safe. ``std::thread`` rather than ``std::jthread`` for Apple Clang libc++ compatibility. ([#parallel-dimer-stdthread](https://github.com/TheochemUI/eOn/issues/parallel-dimer-stdthread))
- Parallel NEB force evaluation via ``std::thread``. Each image evaluates its potential concurrently when the potential is thread-safe. Achieves 2.5x speedup over SVN baseline on 5-image NEB with Morse potential. No external dependencies (replaces TBB-based ``std::execution::par``). Apple Clang's libc++ does not yet ship ``std::jthread``, so the implementation uses ``std::thread`` with explicit exception-safe joins. ([#parallel-neb-stdthread](https://github.com/TheochemUI/eOn/issues/parallel-neb-stdthread))
- ReplicaExchangeJob now runs replica MD steps in parallel via ``std::thread``
  when ``parallel=true`` and the potential supports it. Each replica gets its
  own potential instance when ``needsPerImageInstance()`` is true (e.g. ML
  potentials). ([#parallel-replica-exchange](https://github.com/TheochemUI/eOn/issues/parallel-replica-exchange))
- JSON serialization for Parameters via `nlohmann/json <https://github.com/nlohmann/json>`_. New ``load_json()`` and ``to_json()`` methods enable programmatic configuration for library usage and RPC transport. ([#params-json](https://github.com/TheochemUI/eOn/issues/params-json))
- Added readcon-core v0.7.1 as a Meson subproject for .con/.convel file I/O. The Rust FFI library replaces the hand-written C FILE*-based parser with an mmap-based reader and a type-safe ConFrameBuilder/ConFrameWriter for output. Requires Rust >= 1.88 and cbindgen >= 0.29 at build time.

### Developer

- Planned: approval tests (snapshot-based regression) and fuzz testing
  for parser robustness (INI, .con format, command line arguments). ([#approval-fuzz-testing](https://github.com/TheochemUI/eOn/issues/approval-fuzz-testing))
- Competitive NEB benchmarking framework (private repo HaoZeke/eon_benchmarks).
  Compare eOn against ASE, ORCA, OPTIM, ARTn, NWChem, Sella on Baker test
  set and Pt surfaces. Snakemake workflow with reproducible environments. ([#competitive-benchmarks](https://github.com/TheochemUI/eOn/issues/competitive-benchmarks))
- Added developer design docs for Parameters decomposition and NEB modularization. Updated testing inventory, user guides (minimization, dynamics, NEB, saddle search), and added algorithm selection guide with rgpycrumbs examples and atomistic-cookbook links. ([#docs-design](https://github.com/TheochemUI/eOn/issues/docs-design))
- Replace remaining new[]/delete[] in EMT/Asap (NeighborList.cpp, EMT.cpp)
  with std::vector. 5 allocations in Asap library internals. ([#emt-asap-raii](https://github.com/TheochemUI/eOn/issues/emt-asap-raii))
- Added Fortran column-major layout conversion helpers in ``Eigen.h``:
  ``AtomMatrixF`` type alias, zero-copy ``to/from_fortran_layout()``, and
  ``map_from_flat_colmajor/rowmajor()`` for interfacing with Fortran libraries. ([#feat-eigen-fortran-helpers](https://github.com/TheochemUI/eOn/issues/feat-eigen-fortran-helpers))
- Modernize Fortran source files (SW, Tersoff, EDIP, Lenosky, Aluminum,
  CuH2). Convert to free-form F90, add intent declarations, replace
  common blocks with modules. Requires SVN regression verification. ([#fortran-modernization](https://github.com/TheochemUI/eOn/issues/fortran-modernization))
- Highway SIMD kernels for Morse, LJ, and EAM pair-potential force loops.
  Subproject builds and -DWITH_HIGHWAY is set; needs actual vectorized
  inner loop implementations using HWY_DYNAMIC_DISPATCH. ([#highway-potential-kernels](https://github.com/TheochemUI/eOn/issues/highway-potential-kernels))
- SafeHyperJob integration test with SVN reference data. Needs Morse Pt
  system with element-specific BondBoost parameters (SIGFPE on generic LJ). ([#safe-hyper-test](https://github.com/TheochemUI/eOn/issues/safe-hyper-test))
- Add SVN-verified integration tests for remaining unused reference data:
  neb_morse_pt, global_optimization_lj, minimization_eam_fire,
  minimization_sw_cg, min_lj_sd, replica_exchange_lj, bh50. ([#svn-reference-coverage](https://github.com/TheochemUI/eOn/issues/svn-reference-coverage))
- Migrated test suite from GoogleTest to Catch2. Added 5 new tests (DimerTest, SaddleSearchTest, OptimizerTest, ConFileIOTest, HessianTest) and revived 4 (MatterTest, PotTest, ImpDimerTest, StringHelpersTest). Total: 16 tests. ([#tests-new](https://github.com/TheochemUI/eOn/issues/tests-new))
- Removed abandoned v3c TOML migration artifacts from 11 documentation files. ([#docs-v3c-cleanup](https://github.com/TheochemUI/eOn/issues/docs-v3c-cleanup))

### Changed

- Eigenmode methods (Dimer, ImprovedDimer, Lanczos, GPRDimer) now use ``std::variant`` instead of an abstract base class, eliminating virtual dispatch overhead. ([#eigenmode-variant](https://github.com/TheochemUI/eOn/issues/eigenmode-variant))
- LAMMPS potential now uses runtime dynamic loading (``dlopen``/``LoadLibrary``)
  instead of compile-time linking. A single eOn binary can use LAMMPS
  potentials if ``liblammps`` is installed, without requiring LAMMPS at build
  time. Install via ``conda install -c conda-forge lammps``. ([#lammps-runtime](https://github.com/TheochemUI/eOn/issues/lammps-runtime))
- Extracted ``RandomNumbers``, ``GeometryAnalysis``, and ``ConFileIO`` modules from Matter and HelperFunctions. Matter.cpp reduced from 1253 to ~500 lines, HelperFunctions.cpp from 838 to ~300 lines. ([#module-extract](https://github.com/TheochemUI/eOn/issues/module-extract))
- Removed all ``using namespace std;`` from the entire client codebase
  (60+ files). All standard library symbols now explicitly qualified with
  ``std::``, preventing ADL-related bugs and improving code clarity. ([#namespace-cleanup](https://github.com/TheochemUI/eOn/issues/namespace-cleanup))
- INI parsing extracted from ``Parameters.cpp`` into ``ParametersINI.cpp`` with a ``validate_and_link()`` function for cross-group dependency resolution. ([#params-ini-extract](https://github.com/TheochemUI/eOn/issues/params-ini-extract))
- Replaced vendored 724-line ``CIniFile`` INI parser with `inih <https://github.com/benhoyt/inih>`_ (r62) via meson wrap. ([#params-inih](https://github.com/TheochemUI/eOn/issues/params-inih))
- Narrowed parameter passing: Matter, Potential base, Optimizer hierarchy (``OptimizerConfig``), and Dynamics (``DynamicsConfig``) no longer require the full ``Parameters`` object. NEB no longer mutates its Parameters copy. ([#params-narrow](https://github.com/TheochemUI/eOn/issues/params-narrow))
- All 33 Parameters option-group structs now use C++20 NSDMI (Non-Static Data Member Initialization) for defaults. The constructor shrunk from 481 to 3 lines. ([#params-nsdmi](https://github.com/TheochemUI/eOn/issues/params-nsdmi))
- Fixed pass-by-value of ``VectorXd`` in ``ObjectiveFunction``, ``LBFGS``, and optimizer interfaces. Dynamic Eigen types now passed by ``const`` reference per Eigen documentation, eliminating unnecessary 40-160KB copies per optimizer step. ([#perf-eigen-passbyref](https://github.com/TheochemUI/eOn/issues/perf-eigen-passbyref))
- ``Matter::getForces()`` now returns ``const AtomMatrix&`` to a cached masked-force result instead of copying and zeroing fixed atoms on every call. ([#perf-masked-forces](https://github.com/TheochemUI/eOn/issues/perf-masked-forces))
- NEB tangent and projection strategies cached as class members (built once in constructor). Spring strategy still rebuilt per iteration as it depends on per-step energy data. SIMD-optimized ``Eigen::Map<VectorXd>.dot()`` replaces ``(a.array() * b.array()).sum()`` in all NEB force projections. ([#perf-neb-strategy-cache](https://github.com/TheochemUI/eOn/issues/perf-neb-strategy-cache))
- Vectorized PBC wrapping: replaced scalar ``fmod`` loop with ``floor``-based Eigen array operation. Single x86 ``vroundsd`` instruction instead of expensive ``fmod`` library call per element. ([#perf-pbc-vectorize](https://github.com/TheochemUI/eOn/issues/perf-pbc-vectorize))
- Eliminated ``std::pow()`` with integer exponents from all hot-path force
  loops: LJ (pow(x,6) -> x*x*x), EAM (pow(r,5/6), simplified Morse pair
  from 3 exp to 1), IDPP (pow(r,4/5)), NEB spline (Horner's method),
  Water_Pt (11 pow calls replaced with explicit multiplies). ([#perf-pow-elimination](https://github.com/TheochemUI/eOn/issues/perf-pow-elimination))
- Extracted ``ReplicaDynamicsJob`` base class from TADJob and SafeHyperJob, deduplicating ~190 lines of shared code (checkState, refine, dephase, saveData). ([#replica-base](https://github.com/TheochemUI/eOn/issues/replica-base))
- All job classes now use ``std::shared_ptr<Matter>`` and RAII (``ForceCallTimer``, ``std::ofstream``). No more raw ``new/delete`` for Matter objects or unchecked ``fopen`` calls. ([#smart-ptrs](https://github.com/TheochemUI/eOn/issues/smart-ptrs))
- Converted commented-out SPDLOG_LOGGER_DEBUG calls to active ``QUILL_LOG_TRACE_L1`` (compiled out in release builds). ([#spdlog-quill](https://github.com/TheochemUI/eOn/issues/spdlog-quill))
- Replaced C-style headers (``math.h``, ``string.h``, ``time.h``, etc.) with C++ equivalents (``cmath``, ``cstring``, ``ctime``) in 10 client files. Added ``using enum JobType`` in Job.cpp. ([#cpp20-headers](https://github.com/TheochemUI/eOn/issues/cpp20-headers))
- Python `eon/fileio.py` now uses the `readcon` package (PyPI) for loading and saving .con files, replacing ~60 lines of hand-written parsing. The `loadcon`, `loadcons`, and `savecon` functions delegate to `readcon.read_con()` / `readcon.write_con()`.
- Replaced all FILE*-based con/convel I/O in the C++ client with readcon-core. Reading uses `readcon::read_first_frame()` (mmap). Writing uses `ConFrameBuilder` + `ConFrameWriter` with 17-digit precision for positions. All FILE* overloads removed; callers now pass filenames directly.

### Fixed

- Fixed bare ``abs()`` calls on ``double`` values in BondBoost, LBFGS, and
  Hessian that resolved to the C integer ``abs(int)`` overload on non-MSVC
  compilers, silently truncating floating-point values. ([#bugfix-bare-abs](https://github.com/TheochemUI/eOn/issues/bugfix-bare-abs))
- ``cellInverse`` now recomputed in ``Matter::setCell()`` (was stale after cell changes). ([#bugfix-cellinverse](https://github.com/TheochemUI/eOn/issues/bugfix-cellinverse))
- Fixed ``confine_positive`` Eigen indexing (flat ``3*i+0`` replaced with proper ``(i,k)`` row-column access) and replaced raw ``new[]/delete[]`` with ``std::vector``. ([#bugfix-confine-positive](https://github.com/TheochemUI/eOn/issues/bugfix-confine-positive))
- Fixed ``convergenceForce()`` loop variable mutation that could skip images. ([#bugfix-convergenceforce](https://github.com/TheochemUI/eOn/issues/bugfix-convergenceforce))
- Fixed LBFGS aborting on degenerate curvature updates. The ``abs()`` ->
  ``std::abs()`` fix exposed that the ``s0.y0 < LBFGS_EPS`` check was
  previously disabled by integer truncation. Now resets L-BFGS memory
  (standard restart strategy) instead of aborting the optimization. ([#bugfix-lbfgs-curvature](https://github.com/TheochemUI/eOn/issues/bugfix-lbfgs-curvature))
- Fixed uninitialized ``cuttOffU`` member in LJ potential (caused garbage energy with ``MALLOC_PERTURB_``). ([#bugfix-lj-cutoffu](https://github.com/TheochemUI/eOn/issues/bugfix-lj-cutoffu))
- Fixed ``maxAtomMotionV`` out-of-bounds read when vector size < 3 elements
  (e.g. 2D optimizer objectives). Valgrind caught ``segment<3>`` reading past
  buffer, causing silent data corruption with ``MALLOC_PERTURB_`` enabled. ([#bugfix-maxatommotionv](https://github.com/TheochemUI/eOn/issues/bugfix-maxatommotionv))
- ``maxEnergyImage`` now default-initialized to 0 (was uninitialized). ([#bugfix-maxenergyimage](https://github.com/TheochemUI/eOn/issues/bugfix-maxenergyimage))
- Added null ``FILE*`` guard in NEB write_movies to prevent crash when movie writing is disabled. ([#bugfix-neb-write-movies](https://github.com/TheochemUI/eOn/issues/bugfix-neb-write-movies))
- Added bounds check for ``numExtrema`` to prevent out-of-range access in NEB spline extrema. ([#bugfix-numextrema](https://github.com/TheochemUI/eOn/issues/bugfix-numextrema))
- Fixed Prefactor Hessian size validation that checked min1 vs saddle on both
  sides of the OR condition, never validating min2 frequency array size. ([#bugfix-prefactor-hessian](https://github.com/TheochemUI/eOn/issues/bugfix-prefactor-hessian))
- Fixed DynamicsSaddleSearch MD snapshot recording using shared_ptr aliasing
  instead of deep copy, causing all snapshots to point to the same live object
  and making transition time refinement unreliable. ([#bugfix-snapshot-aliasing](https://github.com/TheochemUI/eOn/issues/bugfix-snapshot-aliasing))


## [2.12.0](https://github.com/TheochemUI/eOn/tree/2.12.0) - 2026-03-08

### Removed

- Remove `_potcalls.log` text logger and `QUILL_LOG_TRACE_L3` from `Potential::get_ef()`. Remove static counters `Potential::fcalls`, `fcallsTotal`, `wu_fcallsTotal`, `totalUserTime`. Output is now `_potcalls.json` with structured per-instance records.

### Added

- Add `PotRegistry` singleton for thread-safe, enum-indexed force call tracking. Replaces per-instance `FileScoped` text loggers and dead static counters (`Potential::fcalls` et al.). Tracks per-instance lifecycle (created_at, destroyed_at, force_calls, unique ID) with JSON output (`_potcalls.json`). Restore force call delta tracking in SafeHyperJob, TADJob, ParallelReplicaJob, ReplicaExchangeJob, NudgedElasticBandJob, BasinHoppingJob, HessianJob, and PrefactorJob.
- Add `SafeMath.h` utility header with guarded arithmetic functions (`safe_div`, `safe_recip`, `safe_acos`, `safe_sqrt`, `safe_atan_ratio`) and an Eigen-aware `safe_normalized` template. These prevent floating-point exceptions (SIGFPE) from division-by-zero and domain errors in numerical code without changing results for valid inputs.
- Drop spdlog and fmt for quill and cpp20

### Changed

- Replace cxxopts with argum for command line parsing. This change updates the
  CLI argument handling library and requires C++20 support. ([#320](https://github.com/TheochemUI/eOn/issues/320))
- Replace per-site `scoped_interpreter` guards with lazy singleton `eonc::ensure_interpreter()` in `PyGuard.h`. Python interpreter is only started when a Python-based potential is actually used. ([#324](https://github.com/TheochemUI/eOn/issues/324))
- Migrated from spdlog/fmt to quill logging library with std::format for modern C++20 logging infrastructure. Quill provides lock-free asynchronous logging with lower latency and better performance characteristics. All LOG_* macros now use quill backend with configurable formatters and sinks. ([#327](https://github.com/TheochemUI/eOn/issues/327))
- Migrate to modern C++20 logging API (EonLogger.h).

  Replaced verbose quill logger initialization throughout codebase with new `eonc::log::Scoped` RAII helper and `eonc::log::get_file()` utility. Eliminates 14-line boilerplate for file loggers and manual initialization for default loggers. Net reduction of 74 lines while improving code clarity and safety.
- Optimize quill backend for improved logging performance.

  Configure `BackendOptions` with reduced sleep duration (100us to 10us), larger initial transit buffer (256 to 2048), zero timestamp ordering grace period (single-threaded, SPSC guarantees ordering), faster flush interval (200ms to 100ms), and disabled printable char checking (numeric data only).
- Switched all logging call sites from bare `LOG_*` to `QUILL_LOG_*` prefixed macros and enabled `QUILL_DISABLE_NON_PREFIXED_MACROS` to prevent macro collisions when eOn is compiled alongside other libraries.
- Updated metatomic ecosystem dependencies: torch 2.10, metatomic-torch 0.1.9+, metatensor-torch 0.8.4+, vesin/vesin-torch 0.5.2+, metatrain 2026.2.1+. Ensures compatibility with latest machine learning potential infrastructure.
- Updated vesin to v0.5.2: adapted to new API with per-dimension periodicity (`bool[3]` instead of single `bool`) and VesinDevice struct syntax. Ensures compatibility with latest metatomic/metatensor ecosystem.
- Wrap all client classes, enums, and helper namespaces under `namespace eonc`. Rename `helper_functions` namespace to `helpers`. `BaseStructures.h` enums (`PotType`, `JobType`, `AtomState`) are now scoped under `eonc`. Backward-compatible `using` aliases are provided for all classes and enums at global scope. Remove `using namespace std` from all client headers to prevent symbol leakage into downstream translation units.

### Fixed

- Fix ASE_POT compilation errors (wrong constructor, PotType, FPE calls, undeclared variable) and rename `-DASE_POT` to `-DWITH_ASE_POT` to avoid macro-enum collision. ([#321](https://github.com/TheochemUI/eOn/issues/321))
- Fix `[IDimerRot]` column misalignment: widen force placeholder from 10 to 18
  dashes, change angle precision from `{:6.2f}` to `{:6.3f}`, and add missing
  "Align" column specifier to the `[Dimer]` header. ([#322](https://github.com/TheochemUI/eOn/issues/322))
- Suppress FPE trapping during libtorch operations in `MetatomicPotential`
  constructor and `force()`, preventing SIGFPE from benign NaN/Inf produced by
  SiLU (sleef) and autograd internals. Follows the existing `FPEHandler` pattern
  from ASE_ORCA, ASE_NWCHEM, and AtomicGPDimer. ([#323](https://github.com/TheochemUI/eOn/issues/323))
- Change `uncertainty_threshold` default from `0.1` to `-1` (disabled) in both
  C++ and Python. Most models lack uncertainty outputs, so the previous default
  triggered a noisy exception+catch in the metatomic constructor for no benefit. ([#325](https://github.com/TheochemUI/eOn/issues/325))
- Fixed quill migration test failures: added logger initialization to all test fixtures (XTBTest, ASEPotTest, ServeSpecParseTest, EpiCentersTest, MetatomicTest) to prevent segfaults from uninitialized quill backend. Restored correct ConfigParser defaults in config.yaml (49 path interpolations accidentally replaced during migration). Added `-DNOMINMAX` for Windows builds to fix MSVC compilation errors in quill headers. ([#327](https://github.com/TheochemUI/eOn/issues/327))
- Added `safe_normalize_inplace` to SafeMath.h and guarded remaining unprotected `.normalize()` / `.normalized()` calls in Dimer, ImprovedDimer, ConjugateGradients, and LBFGS that could trigger FPE on zero vectors.
- Guard unprotected floating-point divisions and domain-error-prone math across 12 source files using `eonc::safemath` utilities. Eliminates spurious SIGFPE signals during saddle search (Dimer, ImprovedDimer, Lanczos), optimization (LBFGS, CG, SteepestDescent), and infrastructure (Matter, Hessian, HelperFunctions, NEB, ReplicaExchange). Fallback values preserve existing branch/skip/reset behavior so valid inputs produce identical results.
- Make the POSIX FPE signal handler async-signal-safe by replacing `std::cerr` (undefined behavior in signal context) with `write(STDERR_FILENO, ...)`. Windows SEH handler switched from `std::cerr` to `fprintf(stderr, ...)` for consistency.
- Use `-isystem` instead of `-I` for pip-installed metatomic and vesin include paths to suppress third-party compiler warnings when building with `-Wall -Wextra`.


## [2.12.0](https://github.com/theochemui/eongit/tree/2.12.0) - 2026-03-04

### Changed

- Replace cxxopts with argum for command line parsing. This change updates the
  CLI argument handling library and requires C++20 support. ([#320](https://github.com/theochemui/eongit/issues/320))
- Replace per-site `scoped_interpreter` guards with lazy singleton `eonc::ensure_interpreter()` in `PyGuard.h`. Python interpreter is only started when a Python-based potential is actually used. ([#324](https://github.com/theochemui/eongit/issues/324))

### Fixed

- Fix ASE_POT compilation errors (wrong constructor, PotType, FPE calls, undeclared variable) and rename `-DASE_POT` to `-DWITH_ASE_POT` to avoid macro-enum collision. ([#321](https://github.com/theochemui/eongit/issues/321))
- Fix `[IDimerRot]` column misalignment: widen force placeholder from 10 to 18
  dashes, change angle precision from `{:6.2f}` to `{:6.3f}`, and add missing
  "Align" column specifier to the `[Dimer]` header. ([#322](https://github.com/theochemui/eongit/issues/322))
- Suppress FPE trapping during libtorch operations in `MetatomicPotential`
  constructor and `force()`, preventing SIGFPE from benign NaN/Inf produced by
  SiLU (sleef) and autograd internals. Follows the existing `FPEHandler` pattern
  from ASE_ORCA, ASE_NWCHEM, and AtomicGPDimer. ([#323](https://github.com/theochemui/eongit/issues/323))
- Change `uncertainty_threshold` default from `0.1` to `-1` (disabled) in both
  C++ and Python. Most models lack uncertainty outputs, so the previous default
  triggered a noisy exception+catch in the metatomic constructor for no benefit. ([#325](https://github.com/theochemui/eongit/issues/325))


## [2.11.1](https://github.com/theochemui/eongit/tree/2.11.1) - 2026-03-01

### Added

- External potential (`ext_pot`) documentation with protocol spec, DeePMD and ASE wrapper examples, and conda-forge availability badges on all potential pages. ([#318](https://github.com/theochemui/eongit/issues/318))

### Developer

- Add `ExtPotTest` unit test verifying the file-based ext_pot protocol with a harmonic spring calculator. ([#318](https://github.com/theochemui/eongit/issues/318))

### Fixed

- Rename `PotType::EXT` to `EXT_POT` so `magic_enum` matches the `ext_pot` config string. Previously `potential = ext_pot` was silently mapped to `UNKNOWN`. ([#318](https://github.com/theochemui/eongit/issues/318))


## [2.11.0](https://github.com/theochemui/eongit/tree/2.11.0) - 2026-02-24

### Added

- Add `eonclient --serve` mode that wraps any eOn potential as an
  rgpot-compatible RPC server over Cap'n Proto. Supports four serving modes:
  single-potential (`--serve-port`), multi-model (`--serve "lj:12345,eam_al:12346"`),
  replicated (`--replicas N` on sequential ports), and gateway (single port with
  round-robin pool via `--gateway`). All options are also available through a
  `[Serve]` INI config section. Requires `-Dwith_serve=true` at build time. ([#316](https://github.com/theochemui/eongit/issues/316))
- Add dictionary-style configuration examples using `rgpycrumbs` to the user
  guide, demonstrating programmatic config generation alongside INI files. ([#317](https://github.com/theochemui/eongit/issues/317))

### Developer

- Switch benchmark PR comment workflow from hand-rolled scripts to the
  `asv-perch` GitHub Action, and parallelize benchmark execution with a matrix
  strategy for main and PR HEAD. ([#315](https://github.com/theochemui/eongit/issues/315))
- Add rgpot subproject wrap, `with_serve` meson option, `serve` pixi environment,
  CI workflow for serve mode builds, and Catch2 unit tests for serve spec parsing. ([#316](https://github.com/theochemui/eongit/issues/316))

### Fixed

- Skip `torch_global_deps` on Windows where the conda-forge libtorch package does not ship it. ([#314](https://github.com/theochemui/eongit/issues/314))
- Fixed serve mode segfault caused by `AtomMatrix` type collision between eOn's
  Eigen-based type and rgpot's custom type. Replaced the `rgpot::PotentialBase`
  virtual interface with a flat-array `ForceCallback`, eliminating the name
  collision entirely. The serve code now only links the capnp schema dependency
  (`ptlrpc_dep`) from rgpot, not the full library. ([#316](https://github.com/theochemui/eongit/issues/316))


## [2.10.2](https://github.com/theochemui/eongit/tree/2.10.2) - 2026-02-22

### Fixed

- Fixed a significant performance regression in NEB calculations caused by incorrect Eigen matrix storage order mapping. Added a regression test and updated CI to automatically mark PRs as draft if benchmark regressions exceed 10x. ([#310](https://github.com/theochemui/eongit/issues/310))
- Absorbed conda-forge Windows patches upstream: replace C99 VLA in XTBPot with `std::vector`, guard empty-string indexing in INIFile, decouple xtb from Fortran requirement, add Windows library search paths for libtorch/metatensor/vesin, guard POSIX headers, and replace shell commands in IMD with `std::filesystem`. ([#312](https://github.com/theochemui/eongit/issues/312))


## [2.10.1](https://github.com/theochemui/eongit/tree/2.10.1) - 2026-02-18

### Developer

- Added a CI-NEB XTB regression test (`CINEBXTBTest.cpp`) that runs a 10-image
  climbing-image NEB with GFN2-xTB on a 9-atom molecule.  The test completes in
  under 2 seconds and guards against storage-order regressions that corrupt force
  projections.

### Fixed

- Replaced the `EIGEN_DEFAULT_TO_ROW_MAJOR` preprocessor macro with explicit
  row-major type aliases in `client/Eigen.h`.  The macro made eOn's Eigen types
  binary-incompatible with other Eigen-based libraries; removing it without
  updating bare `MatrixXd` types caused NEB force projections to silently corrupt
  and the optimizer to diverge from the first step.
- Use `datetime.timezone.utc` instead of `datetime.UTC` in `get_version.py` for
  Python 3.10 compatibility (the `datetime.UTC` alias was added in 3.11).


## [v2.10.0](https://github.com/theochemui/eongit/tree/v2.10.0) - 2026-02-15

### Added

- Added ASV benchmark CI workflow with asv-spyglass for PR performance comparison
- Added adsorbate_region.py example script for identifying adsorbate atoms and nearby surface atoms by element or z-coordinate
- Added displacement scripts tutorial with worked examples for vacancy diffusion (PTM) and adsorbate-on-surface scenarios
- Added displacement strategies prose section to saddle search docs explaining epicenters, weight-based selection, and dynamic atom lists
- Expose gprd_linalg_backend option for selecting GPR-dimer linear algebra backend (eigen, cusolver, kokkos, stdpar)

### Developer

- Added macOS arm64 to metatomic CI matrix using Homebrew gfortran (conda-forge gfortran_osx-arm64 wrapper is broken)
- Cleanup to build on windows
- Expanded ASV benchmark suite with point evaluation, LJ minimization, and NEB workloads
- Use internal pick output helper
- bld(meson): reduce build times by linking to xtb by default

### Changed

- Eliminated unnecessary Eigen matrix copies in Matter, Potential, and HelperFunctions hot paths
- Replace per-typedef `Eigen::RowMajor` with a single `eOnStorageOrder` constant in `client/Eigen.h`
- Enriched schema descriptions for displace_atom_kmc_state_script, displace_all_listed, displace_atom_list, and client_displace_type
- Refactored MetatomicPotential variant resolution to use upstream metatomic_torch::pick_output
- Updated pinned gpr_optim commit with new linear algebra backends and performance improvements

### Fixed

- Fix Windows `STATUS_STACK_OVERFLOW` crash caused by large Fortran local arrays in the EAM Al potential (`gagafeDblexp.f`) exceeding the 1 MB default stack; request 16 MB via linker flags
- Fix Windows silent client failure by using non-color spdlog sink when stdout is redirected
- Use Goswami & Jonsson 2025 for removing rotations through projections


## [v2.9.0](https://github.com/theochemui/eongit/tree/v2.9.0) - 2026-01-27

### Added

- Add support for 'charge' and 'uhf' (multiplicity) parameters in the xTB potential
- Introduce custom Catch2 Eigen matchers and add comprehensive regression tests for GFN2-xTB forces
- Setup Collective-IDPP path generation for NEB runs
- Setup IDPP path generation for NEB runs
- Setup sequential Collective-IDPP path generation for NEB runs
- feat(mtapot): handle variants for energy and energy uncertainty within Metatomic models
- feat(neb): add a zbl+sidpp penalty for initial path generation
- feat(neb): implement the OCI-NEB/RONEB/enhanced CI via MMF
- feat(neb): implement the onsager machlup action logic
- feat(neb): write out peaks and modes for subsequent dimer runs

### Changed

- Optimize xTB potential performance by persisting internal state and using coordinate updates between force calls
- Update installation guide to recommend Pixi and clarify dependency management

### Fixed

- bug(ewneb): do not turn on if cineb threshold is not met!
- fix(mtapot): stop double counting mta calls


## [v2.8.2](https://github.com/theochemui/eongit/tree/v2.8.2) - 2025-12-01

### Added

- Metatomic is now uncertainty aware
- Metatomic variance reports per-atom uncertainty mean

### Changed

- Reworked metatomic to use torch 2.9

### Developer

- Update to use `metatensor_torch::Module`
- Use `metatomic::pick_device` correctly


## [v2.8.1](https://github.com/theochemui/eongit/tree/v2.8.1) - 2025-11-03

### Added

- Enable minimization for given initial paths

### Changed

- Reworked metatomic to use torch 2.8

### Fixed

- Generate neb.dat correctly without clobbering neb_000.dat


## [2.8.0](https://github.com/theochemui/eongit/tree/2.8.0) - 2025-09-04

### Added

* **Potentials & Interfaces**
    * Expanded potential interfaces to a variety of new quantum chemistry and ML
      potentials via an embedded Python interpreter:
        * **NWChem**: A high-performance, socket-based interface.
          ([#244](https://github.com/theochemui/eongit/issues/244))
        * **ORCA**: Interface to the ORCA quantum chemistry program via ASE.
        * **AMS**: Interface for the Amsterdam Modeling Suite.
        * **XTB**: Interface for semi-empirical GFN-xTB methods.
        * **ASE**: A general-purpose interface to any calculator supported by
          the Atomic Simulation Environment.
    * Added the Ziegler-Biersack-Littmark (ZBL) universal screening potential,
      useful for collision cascade simulations.
      ([#241](https://github.com/theochemui/eongit/issues/241))
    * Integrated support for `metatomic` machine-learned potentials via the
      `vesin` library, enabling high-performance simulations with models from
      the metatensor ecosystem.
      ([#201](https://github.com/theochemui/eongit/issues/201))

* **Nudged Elastic Band (NEB)**
    * NEB calculations can now pre-optimize the initial and final states,
      improving path quality and convergence. This feature is fully compatible
      with restarts. ([#221](https://github.com/theochemui/eongit/issues/221))
    * NEB calculations can now be initialized from a user-provided sequence of
      structures, offering greater control over the initial reaction pathway.
    * Introduced energy-weighted springs to improve the stability and quality of
      paths with high energy barriers.
    * Enabled the use of dual optimizers (e.g., a starting with QuickMin and
      switching to LBFGS after a convergence threshold).
    * Implement the novel RO-NEB-CI (Rohit's Optimal NEB with MMF CI steps)
      method ([#239](https://github.com/theochemui/eongit/issues/239))


### Developer

- Consistent formatting and counting
- Support for M1 MacOS machines

#### Build & Tooling

* **Build System**
    * Overhauled the build system, migrating from legacy Makefiles/CMake to
      **Meson** for a faster, more reliable, and truly cross-platform build
      experience. This change also lays the groundwork for a future pure Python
      `eon-server` package.
      ([#124](https://github.com/theochemui/eongit/issues/124))
* **Dependency Management**
    * Adopted `pixi` and `conda-lock` for robust, reproducible dependency
      management across all platforms.
* **Cross-Platform Support & CI**
    * Established a full Continuous Integration (CI) pipeline, testing on Linux,
      Windows, and macOS (Intel & Apple Silicon).
    * The Command Line Interface (CLI) is now fully compatible with Windows
      environments.


#### Code Quality & Refactoring

* **C++ Modernization**
    * Modernized the C++ backend to the C++17 standard, improving code clarity
      and performance.
    * Enhanced memory safety by replacing raw pointers with smart pointers
      (`std::unique_ptr`, `std::shared_ptr`).
    * Adopted the `<filesystem>` library for platform-independent file I/O.
* **Logging**
    * Replaced the internal logging system with `spdlog` for high-performance,
      asynchronous, and more informative configurable output.
* **Code Style**
    * Enforced a consistent code style and formatting across the entire C++ and
      Python codebase.

#### Documentation

* **Configuration & Schema**
    * Implemented a comprehensive **Pydantic schema** for all configuration
      files, providing automatic input validation and clear error messages. This
      forms the foundation for automated API documentation.
* **User Guides**
    * Added detailed user documentation for the Nudged Elastic Band (NEB)
      module, covering theory, keywords, and practical examples.
