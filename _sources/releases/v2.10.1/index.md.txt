---
myst:
  html_meta:
    "description": "Overview of the eOn v2.10.1 patch release, fixing the Eigen row-major storage-order regression."
    "keywords": "eOn v2.10.1, release, bug fix"
---

## [v2.10.1] - 2026-02-18

Patch release that fixes a critical NEB convergence regression introduced in
v2.10.0, where removing the `EIGEN_DEFAULT_TO_ROW_MAJOR` macro left bare
`MatrixXd` types as column-major, silently corrupting force projections.  Also
fixes a Python 3.10 compatibility issue and adds a CI-NEB regression test to
prevent similar regressions in the future.

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
```
