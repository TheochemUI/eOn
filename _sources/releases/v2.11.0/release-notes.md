---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.11.0, covering rgpot serve mode, rgpycrumbs examples, asv-perch CI, and Windows fixes."
    "keywords": "eOn release notes, serve mode, rgpot, RPC, gateway, asv-perch, rgpycrumbs"
---

# Release notes

## [v2.11.0] - 2026-02-24

### Added

#### Serve mode

eOn can now serve any potential over Cap'n Proto RPC using the rgpot protocol.
Run `eonclient -p lj --serve-port 12345` for the simplest case, or use
`--serve "lj:12345,eam_al:12346"` for multi-model serving.  Replicated mode
(`--replicas N`) starts N independent servers on sequential ports; gateway mode
(`--gateway`) places a round-robin dispatcher in front of the pool so clients
only need one address.

All options are also available through a `[Serve]` INI config section for
fully config-driven operation.  The feature requires `-Dwith_serve=true` at
build time and a Cap'n Proto installation (provided by the `serve` pixi
environment).

See {doc}`../../user_guide/serve_mode` for full documentation.

#### rgpycrumbs configuration examples

The user guide now includes dictionary-style configuration examples using
rgpycrumbs, demonstrating programmatic config generation alongside INI files.

### Developer

#### asv-perch benchmark CI

The benchmark PR comment workflow now uses the asv-perch GitHub Action instead
of hand-rolled shell/JS scripts.  Benchmark execution for main and PR HEAD is
parallelized via a matrix strategy.

#### Serve mode build infrastructure

Added rgpot subproject wrap, `with_serve` meson option, `serve` pixi
environment, a CI workflow for serve mode builds, and Catch2 unit tests for
serve spec parsing.

### Fixed

#### Windows torch_global_deps

Skip `torch_global_deps` on Windows where the conda-forge libtorch package
does not ship it.

#### Serve mode AtomMatrix collision

Fixed a segfault caused by a type collision between eOn's Eigen-based
`AtomMatrix` and rgpot's custom `AtomMatrix`.  The `rgpot::PotentialBase`
virtual interface was replaced with a flat-array `ForceCallback`, and the serve
code now only links the capnp schema dependency from rgpot.
