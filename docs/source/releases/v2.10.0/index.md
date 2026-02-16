---
myst:
  html_meta:
    "description": "Overview of the eOn v2.10.0 release, featuring displacement script documentation, performance improvements, and developer tooling."
    "keywords": "eOn v2.10.0, release"
---

## [v2.10.0] - 2026-02-15

This release focuses on documentation, developer tooling, and performance. It
adds comprehensive user-facing documentation for the displacement atom list and
script features (targeted saddle search displacements), exposes new GPU and
parallel linear algebra backends for the GPR-dimer method, introduces ASV
benchmark infrastructure, and includes C++ performance improvements that
eliminate unnecessary Eigen matrix copies in hot paths.

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
```
