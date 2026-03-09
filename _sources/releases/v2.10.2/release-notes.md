---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.10.2, covering Windows MSVC fixes, Eigen row-major mapping, and version harmonization."
    "keywords": "eOn release notes, Windows, MSVC, VLA, INIFile, xtb, libtorch, metatensor, Eigen row-major"
---

# Release notes

## [v2.10.2] - 2026-02-22

### Fixed

#### Eigen row-major storage for atom data mapping

`Eigen::Map` calls that wrap raw `double*` atom data now specify
`Eigen::RowMajor` explicitly, ensuring consistent `[x0, y0, z0, x1, ...]`
layout regardless of compile-time defaults.  This completes the storage-order
cleanup started in v2.10.1.

#### Windows MSVC compatibility (conda-forge upstream absorption)

Several fixes originally carried as conda-forge recipe patches have been
absorbed into the source tree:

- **XTBPot VLA removal**: replaced the C99 variable-length array in
  `XTBPot.cpp` with `std::vector`, since MSVC does not support VLAs.
- **INIFile empty-string guard**: added a bounds check before indexing the
  last character of a string, preventing undefined behavior on empty lines
  under MSVC.
- **POSIX header guards**: wrapped `<unistd.h>` and `<sys/wait.h>` includes
  in `#ifdef` guards and replaced `system()` shell calls in IMD with
  `std::filesystem` equivalents.
- **Meson build-system fixes**: decoupled xtb from the Fortran requirement,
  added Windows library search paths for libtorch, metatensor, and vesin.

### Changed

#### Version handling harmonization

All version consumers (C++ client, Python package, Sphinx docs) now derive
from `pyproject.toml` as the single source of truth.

### Developer

#### Benchmark CI improvements

The benchmark comparison step was moved into the commenter workflow and now
uses `uvx --quiet` for cleaner PR comment output.
