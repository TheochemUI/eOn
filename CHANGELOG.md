

<!-- towncrier release notes start -->

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
