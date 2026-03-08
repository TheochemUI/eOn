

<!-- towncrier release notes start -->

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
