---
myst:
  html_meta:
    "description": "Overview of the eOn v2.16.0 cut: in-process rgpot potential (NWChem/CPMD), Davidson min-mode solver, readcon 0.13 ConFileIO rework, Metatomic output controls, PbcConvention."
    "keywords": "eOn v2.16.0, rgpot, RgpotPot, Davidson, readcon, Metatomic, PbcConvention"
---

## [v2.16.0] - pending cut

Post-`v2.15.0` work on `main`, cut with the standard release flow
({doc}`/devdocs/release`). This page is authored **before** the version chore so
the cut is not blocked on writing notes.

Highlights: the in-process **RgpotPot** potential (NWChem / CPMD through
[rgpot](https://github.com/OmniPotentRPC/rgpot) with `dlopen`ed engine
libraries), a **Davidson** minimum-mode solver as a dimer/Lanczos alternative,
the readcon-core 0.13 `ConFileIO` rework with force+energy co-loading,
Metatomic output/determinism/rotation controls on the current
metatomic-torch stack, `PbcConvention` for position wrapping, and a Morse
force hot-path rework.

```{toctree}
:maxdepth: 2
:caption: Release notes

release-notes
```
