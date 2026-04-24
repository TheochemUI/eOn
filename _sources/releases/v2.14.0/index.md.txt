---
myst:
  html_meta:
    "description": "Overview of eOn v2.14.0: structured per-iteration metadata in .con movies via readcon-core, ARTn filin exposure, batched parallel force evaluation for NEB, and threading fixes for legacy Fortran potentials."
    "keywords": "eOn v2.14.0, structured trajectory output, readcon-core, ARTn filin, batched force evaluation, parallel NEB"
---

## [v2.14.0] - 2026-XX-XX

Trajectory output gets a first-class metadata format and NEB gets a real
batched force path. `readcon-core` now ships a per-frame metadata
channel so minimization, saddle-search, and NEB movies carry the
per-iteration state they previously scattered across sidecar `.dat`
files; the legacy tables remain behind an opt-in flag. pARTn's `filin`
namelist input is finally exposed to the `[ARTn]` section so users can
point at an `artn.in` file without patching sources. Parallel image
evaluation learns to batch potential calls across images instead of
forking one thread per image, and a threading regression that had crept
into the shared-instance path for legacy Fortran potentials is fixed.

```{toctree}
:maxdepth: 2
:caption: Release notes

release-notes
```
