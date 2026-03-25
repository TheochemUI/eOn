---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.11.1, covering ext_pot enum fix, external potential documentation, and conda-forge Windows serve mode fixes."
    "keywords": "eOn release notes, ext_pot, external potential, conda-forge, Windows, MSVC, Cap'n Proto"
---

# Release notes

## [v2.11.1] - 2026-03-01

### Added

#### External potential documentation

The user guide now includes a full specification of the file-based external
potential (`ext_pot`) protocol, covering the `FORCE_INPUT` and `FORCE_OUTPUT`
file formats, atomic number and Cartesian coordinate conventions, and expected
energy/force units.  Wrapper examples for DeePMD-kit and ASE are provided.
Conda-forge availability badges have been added to all potential documentation
pages.

### Developer

#### `ExtPotTest` unit test

A Catch2 unit test exercises the full `ext_pot` round-trip: writes atomic
positions, invokes a harmonic spring calculator, reads back forces, and verifies
energy and force values against analytic expectations.

### Fixed

#### `ext_pot` enum mapping

`PotType::EXT` was renamed to `PotType::EXT_POT` so that `magic_enum` correctly
maps the `potential = ext_pot` configuration string.  Previously the string
`ext_pot` did not match the enum name `EXT`, causing silent fallback to
`UNKNOWN` and a confusing error about unsupported potential types.

#### Conda-forge Windows serve mode (packaging)

Three fixes were applied to the rgpot subproject for Windows MSVC builds on
conda-forge (eon-feedstock PR #23).  These do not affect eOn source code
directly but are required for the conda-forge package:

- **`capnpc` to `capnp compile`**: Windows `capnpc.EXE` does not support the
  `-o` flag; the portable `capnp compile` subcommand is used instead.
- **`.c++` to `.cpp` rename**: MSVC `cl.exe` does not recognize the `.c++`
  extension that Cap'n Proto generates; a Python wrapper renames the output.
- **`ws2_32` linkage**: Cap'n Proto's `kj-async` uses Winsock2 on Windows but
  the conda-forge package does not export it as a transitive dependency;
  `ws2_32` is now linked conditionally through `ptlrpc_dep`.
