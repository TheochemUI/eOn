## [v2.8.2](https://github.com/theochemui/eongit/tree/v2.8.2) - 2025-12-01

### Added

- Metatomic is now uncertainity aware
- Metatomic variance reports per-atom uncertainity mean

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
