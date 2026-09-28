---
myst:
  html_meta:
    "description": "eOn v3.4.0: adaptive kinetic Monte Carlo on a Slurm cluster with multi-rank codes such as CPMD, portable pyeonclient wheels, zoom-NEB and the in-process job API."
    "keywords": "eOn v3.4.0, AKMC, CPMD, ext_pot, Slurm, pyeonclient, zoom-NEB, release"
---

## [v3.4.0] - unreleased

Minor release on `v3.3.1`. AKMC now runs on a Slurm cluster with a
multi-rank external code such as CPMD, one job per search. The
pyeonclient wheels are manylinux again and import without the build
host's libraries. The release also adds zoom-NEB and the in-process
job API, and it requires rgpot 3.3.0.

```{toctree}
:maxdepth: 2
:caption: Release notes

release-notes
```
