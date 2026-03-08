---
myst:
  html_meta:
    "description": "Overview of the eOn v2.11.2 release, fixing metatomic FPE trapping, uncertainty default, and improved dimer log alignment."
    "keywords": "eOn v2.11.2, release, metatomic, SIGFPE, FPE, uncertainty, dimer"
---

## [v2.11.2] - 2026-03-02

A patch release that fixes SIGFPE crashes when using the built-in `metatomic`
potential (libtorch triggers benign FPE), changes the `uncertainty_threshold`
default to disabled (`-1`) so models without uncertainty outputs no longer
produce noisy warnings, and corrects column alignment in `[IDimerRot]` and
`[Dimer]` log output.

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
```
