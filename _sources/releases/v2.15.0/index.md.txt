---
myst:
  html_meta:
    "description": "Overview of the planned eOn v2.15.0 cut: runtime Fortran potentials, process-search fixed-atom fix, LAMMPS/MPI NEB isolation, release process (GH/PyPI/feedstock) automation."
    "keywords": "eOn v2.15.0, runtime potentials, process search, LAMMPS MPI, release process, PyPI"
---

## [v2.15.0] - pending cut

Post-`v2.14.0` work on `main` plus maintainer release-process hardening, consolidated
here for the **next** semver (supersede incomplete multi-channel `2.14.0` per
{doc}`/devdocs/release` §5 rather than re-tagging). Run `cog bump --version 2.15.0`
only after pre-release gates pass; this page is authored **before** the version
chore so the cut is not blocked on writing notes.

Highlights: runtime-loaded Fortran potentials (`potentials_path` / `dlopen`),
process-search fixed-atom restore, per-image LAMMPS MPI isolation, MPI C API
build path, and the full release/PyPI/feedstock/Nickel GHA process.

```{toctree}
:maxdepth: 2
:caption: Release notes

release-notes
```
