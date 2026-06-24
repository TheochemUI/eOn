

<!-- towncrier release notes start -->

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
