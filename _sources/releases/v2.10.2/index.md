---
myst:
  html_meta:
    "description": "Overview of the eOn v2.10.2 release, absorbing conda-forge Windows patches upstream and fixing Eigen storage-order mapping."
    "keywords": "eOn v2.10.2, release, Windows, MSVC, Eigen, conda-forge"
---

## [v2.10.2] - 2026-02-22

A patch release that absorbs the conda-forge Windows patches upstream,
fixing MSVC compilation across several potentials and build-system modules.
Also fixes a remaining Eigen row-major storage-order issue in `Eigen::Map`
calls for atom data, and harmonizes version handling so all consumers derive
from `pyproject.toml` as the single source of truth.

See {doc}`../v2.10.1/index` for the Eigen storage-order regression fix
introduced in v2.10.1.

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
```
