

<!-- towncrier release notes start -->

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
