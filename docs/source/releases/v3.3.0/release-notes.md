# eOn 3.3.0

A minor on `v3.2.1`. See monorepo `CHANGELOG.md` section 3.3.0 for the
fragment list after towncrier consumes `docs/newsfragments/`.

## Parameters

Option groups are private. Callers read through const accessors
(`main_options()`, `potential_options()`, …). INI/JSON loaders and
nanobind setters write through `ParametersLoadAccess`. MPI comm/rank
use dedicated setters. `load`, `load_ini_text`, and `load_json` take
`std::string_view`.

## Potential

`Potential::force(std::span…)` checks sizes and forwards to the raw
C-array virtual. Matter, `get_ef`, and surrogate `get_ef_var` use it.
Fortran and FFI loaders keep pointers.

## Python atoms

CNA builds unique-index adjacency from vesin `ijS`. The local table is
a Z-keyed radius/color overlay. Symbol and Z go through readcon
`symbol_to_atomic_number` / `atomic_number_to_symbol` (0.14.9+).
`rot_match` prefers `pyeonclient._core.ira_match`.

## Channels

PyPI project is `eon-akmc`. conda-forge package is `eon` from the
GitHub fat tarball `eon-v3.3.0.tar.xz`.
