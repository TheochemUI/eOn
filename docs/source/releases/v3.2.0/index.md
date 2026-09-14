---
myst:
  html_meta:
    "description": "eOn v3.2.0: CON spec 3 output through readcon 0.14.5, per-axis atom constraints, streaming .con append, external-potential IO checks, batched rgpot force evaluation, and pytest on Unix CI."
    "keywords": "eOn v3.2.0, readcon 0.14.5, CON spec 3, per-axis constraints, ExtPot, VASP, rgpot forceBatch, pytest CI"
---

## [v3.2.0] - 2026-08-16

Post-`v3.1.0` work on `main`, cut with the standard release flow
({doc}`/devdocs/release`). This page is authored **before** the version chore so
the cut is not blocked on writing notes.

The release is dominated by an audit of the client and server IO paths: 170 of
the 331 commits since `v3.1.0` carry a `fix` prefix. It moves both the C++ wrap
and the Python package to **readcon 0.14.5**, which changes the CON spec version
eOn writes; see the compatibility note below before upgrading a pipeline. It
also stores atom constraints **per Cartesian axis**, appends `.con` frames
without rewriting the movie, checks every file and command the
external-program potentials rely on, adds a batched force entry point on the
rgpot adapter, and runs `tests/` under Unix CI for the first time.

```{toctree}
:maxdepth: 2
:caption: Release notes

release-notes
```
