---
myst:
  html_meta:
    "description": "Overview of the eOn v2.12.0 release, replacing spdlog/fmt with quill, adding namespace eonc, PotRegistry, SafeMath, and FPE hardening."
    "keywords": "eOn v2.12.0, release, quill, namespace, PotRegistry, SafeMath"
---

## [v2.12.0] - 2026-03-08

This release is a major internal modernization. The logging stack moves from
spdlog/fmt to quill (lock-free, async, lower latency). All client code is
wrapped in `namespace eonc` with backward-compatible `using` aliases. A new
`PotRegistry` singleton replaces the dead static force call counters and the
per-call `_potcalls.log` text logger with structured JSON output and
per-instance lifecycle tracking. `SafeMath.h` guards all division-by-zero and
domain-error-prone math to eliminate spurious SIGFPE crashes.

Other highlights: FPE hardening across the codebase, metatomic ecosystem
updates (torch 2.10, metatomic-torch 0.1.9+), vesin v0.5.2 per-dimension
periodicity API, and argum replacing cxxopts for CLI parsing.

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
```
