---
myst:
  html_meta:
    "description": "Overview of the eOn v2.11.1 release, fixing the ext_pot enum mapping, adding ext_pot documentation, and improving conda-forge packaging."
    "keywords": "eOn v2.11.1, release, ext_pot, external potential, conda-forge"
---

## [v2.11.1] - 2026-03-01

A patch release that fixes the `ext_pot` potential type so it can be selected
from configuration files, adds documentation and protocol specification for the
external potential interface, and includes a unit test for the file-based
protocol.  This release also ships the conda-forge packaging fixes for
Windows serve mode builds (Cap'n Proto codegen, MSVC `.c++` extension,
`ws2_32` linkage).

See {doc}`../v2.11.0/index` for the serve mode feature introduced in v2.11.0.

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
```
