---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.8.0, covering new features, enhancements, build system improvements, and documentation updates."
    "keywords": "eOn v2.8.0 release notes, new features, enhancements, bug fixes"
---

# Release notes

## [2.8.0] - 2025-09-04

Five years of work: new transition-state methods, more potential
interfaces, and a Meson build that runs on Linux, macOS, and Windows.

### ✨ Major New Features

  * **RO-NEB-CI**: climbing-image NEB with a rotating-orbit refinement for
    locating complex transition states.
    ([\#239](https://github.com/theochemui/eongit/issues/239))
  * **Machine-Learned Potentials with `metatomic`**: Integrated full support for
    `metatomic` potentials, enabling high-performance simulations with
    current models from the [metatensor
    ecosystem](https://docs.metatensor.org/latest/index.html).
    ([\#201](https://github.com/theochemui/eongit/issues/201))
  * **Expanded Potential Interfaces**: Drastically increased interoperability by
    adding direct interfaces for a wide range of popular computational chemistry
    packages:
      * **NWChem**: High-performance, socket-based interface for large-scale
        calculations. ([\#244](https://github.com/theochemui/eongit/issues/244))
      * **ORCA**: Direct support for the versatile and efficient ORCA quantum
        chemistry program.
      * **AMS**: Interface for the Amsterdam Modeling Suite (ADF, BAND, DFTB).
      * **XTB**: Fast and reliable calculations using the GFN-xTB semi-empirical
        tight-binding methods.
      * **ASE**: General-purpose interface to any calculator supported by the
        Atomic Simulation Environment (ASE).
  * **Input Validation and Schema**: A Pydantic schema for all configuration
    files. Automatic input validation, clearer error messages, and a
    base for auto-generated documentation.
  * **ZBL Universal Potential**: Added support for the Ziegler-Biersack-Littmark
    (ZBL) universal screening potential, ideal for simulating high-energy
    collision cascades and ion-implantation effects.
    ([\#241](https://github.com/theochemui/eongit/issues/241))

### 🚀 Enhancements & Improvements

  * **Nudged Elastic Band (NEB) Enhancements**: The NEB module has received
    significant upgrades for flexibility and performance:
      * **Endpoint Minimization**: Added an option to pre-optimize the initial
        and final states of a NEB path, improving convergence and accuracy. This
        feature is fully compatible with restarts.
        ([\#221](https://github.com/theochemui/eongit/issues/221))
      * **Custom Initial Paths**: Users can now provide a custom series of
        intermediate structures to initialize a NEB calculation, offering
        greater control over the reaction pathway.
      * **Dual Optimizers**: Enabled the use of different optimizers for the
        climbing image versus the regular images in a CI-NEB calculation.
      * **Energy-Weighted Springs**: Implemented energy-weighted dynamic
        springs, improving stability and performance for reaction paths with
        high energy barriers.

### 🛠️ Build, Installation & Developer Experience

  * **Modernized Build System with Meson**: The entire project has been ported
    to the Meson build system. One build on Linux, macOS, and Windows. The
    same tree is the base for a later pure-Python `eon-server` package.
    ([\#124](https://github.com/theochemui/eongit/issues/124))
  * **Cross-Platform CI & Support**:
      * Continuous Integration (CI) pipeline, automatically
        testing builds and features across Linux, Windows, and macOS (Intel &
        Apple Silicon ARM).
      * Official support for Apple Silicon (M1/M2/M3) machines.
      * The command-line interface is now fully compatible with Windows
        environments.
  * **Dependency Management**: Added official support for Conda and
    Mamba, allowing for simple, one-command environment setup.
  * **High-Performance Logging**: Replaced the internal logging system with
    `spdlog` for asynchronous, configurable, and more informative output with
    minimal performance overhead.
  * **Code Modernization**:
      * The C++ backend has been upgraded to the C++17 standard.
      * Adopted modern C++ features like smart pointers and the `<filesystem>`
        library for safer memory management and platform-independent file I/O.
      * Enforced consistent code style and formatting across the entire
        codebase.

### 📚 Documentation

  * Added documentation for the Nudged Elastic Band (NEB) module,
    covering theory, keywords, and practical examples.
