---
myst:
  html_meta:
    "description": "Overview of the eOn v2.11.0 release, adding rgpot-compatible RPC serve mode, rgpycrumbs config examples, and asv-perch benchmark CI."
    "keywords": "eOn v2.11.0, release, serve mode, rgpot, RPC, asv-perch"
---

## [v2.11.0] - 2026-02-24

This release adds a new serve mode that exposes any eOn potential over Cap'n
Proto RPC via the rgpot protocol, enabling integration with external tools such
as ChemGP.  Four serving modes are supported: single-potential, multi-model,
replicated, and gateway (round-robin pool behind a single port).  The feature
is gated behind the `-Dwith_serve=true` meson option and documented in the
{doc}`../../user_guide/serve_mode` guide.

Also included: dictionary-style configuration examples using rgpycrumbs, a
switch to asv-perch for benchmark PR comments, and a Windows torch_global_deps
fix.

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
```
