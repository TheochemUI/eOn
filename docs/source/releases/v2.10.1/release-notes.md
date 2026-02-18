---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.10.1, covering the Eigen row-major fix and CI-NEB regression test."
    "keywords": "eOn release notes, Eigen row-major, NEB regression, storage order"
---

# Release notes

## [v2.10.1] - 2026-02-18

This is a bug-fix release.  There are no API or configuration changes.

### Bug Fixes

#### Eigen storage-order regression

The v2.10.0 commit `6e8461c3` replaced the `EIGEN_DEFAULT_TO_ROW_MAJOR`
preprocessor macro with a `constexpr` constant and explicit row-major aliases
for `AtomMatrix` and `RotationMatrix`.  However, ~135 bare `Eigen::MatrixXd`
uses across the codebase were not updated, causing them to revert to Eigen's
default column-major layout.

Row-major layout is required so that `.data()` yields
`[x0, y0, z0, x1, y1, z1, ...]`, which is what the Fortran potentials and
`VectorXd::Map` round-trips expect.  The column-major matrices stored
`[x0, x1, ..., y0, y1, ..., z0, z1, ...]` instead, silently corrupting force
projections and causing the NEB to diverge from the first step.

The fix removes `using namespace Eigen;` from `client/Eigen.h` (which would
shadow the explicit row-major aliases with Eigen's column-major defaults) and
replaces it with selective `using` declarations for vector types that are
unaffected by storage order (`VectorXd`, `VectorXi`, `Vector3d`).  All bare
`Eigen::MatrixXd`, `Eigen::Matrix3d`, and `Eigen::Matrix4d` uses are replaced
with eOn's row-major aliases.

See {doc}`../../devdocs/design/client/types` for the full type-system
documentation.

#### Python 3.10 compatibility

`client/get_version.py` used `datetime.UTC`, an alias added in Python 3.11.
Replaced with `datetime.timezone.utc` for compatibility with the conda-forge
Python 3.10 builds.

### Developer

#### CI-NEB XTB regression test

A new Catch2 test (`CINEBXTBTest.cpp`) runs a 10-image climbing-image NEB with
GFN2-xTB on a 9-atom C2H4ON2 molecule.  The test converges in ~64 NEB steps
(under 2 seconds wall time) and will fail immediately if a storage-order
regression reappears.
