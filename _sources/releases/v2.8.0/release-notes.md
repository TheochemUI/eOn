---
myst:
  html_meta:
    "description": "Detailed release notes for EON v2.8.0, covering new features, enhancements, build system improvements, and documentation updates."
    "keywords": "EON v2.8.0 release notes, new features, enhancements, bug fixes"
---

# Release notes

## [2.8.0] - 2025-09-04

This is a landmark release for EON, representing over five years of dedicated
development. The entire framework has been modernized to be more powerful,
easier to use, and more accessible to researchers across all major platforms.
This release introduces cutting-edge simulation methods, a vast expansion of
supported computational chemistry codes, and a complete overhaul of the build
and installation process.

### ✨ Major New Features

  * **Advanced Transition State Finding with RO-NEB-CI**: Implemented the novel
    RO-NEB-CI method, providing a powerful new tool for accurately locating
    complex transition states.
    ([\#239](https://github.com/theochemui/eongit/issues/239))
  * **Machine-Learned Potentials with `metatomic`**: Integrated full support for
    `metatomic` potentials, enabling high-performance simulations with
    state-of-the-art models from the [metatensor
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
  * **Input Validation and Schema**: Implemented a comprehensive Pydantic schema
    for all configuration files. This provides automatic input validation,
    clearer error messages for users, and a robust foundation for auto-generated
    documentation.
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
    to the Meson build system. This provides a significantly faster, more
    reliable, and truly cross-platform installation experience on Linux, macOS,
    and Windows. This work also lays the foundation for a future pure Python
    `eon-server` package.
    ([\#124](https://github.com/theochemui/eongit/issues/124))
  * **Cross-Platform CI & Support**:
      * Established a robust Continuous Integration (CI) pipeline, automatically
        testing builds and features across Linux, Windows, and macOS (Intel &
        Apple Silicon ARM).
      * Official support for Apple Silicon (M1/M2/M3) machines.
      * The command-line interface is now fully compatible with Windows
        environments.
  * **Streamlined Dependency Management**: Added official support for Conda and
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

  * Added comprehensive documentation for the Nudged Elastic Band (NEB) module,
    covering theory, keywords, and practical examples.
