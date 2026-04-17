---
myst:
  html_meta:
    "description": "Overview of eOn v2.13.0: ARTn saddle search via pARTn, OCINEB climbing-image + min-mode hybrid, IRA structure comparison, NEB strategy decomposition, Highway SIMD subproject, parallel dimer/NEB via std::thread."
    "keywords": "eOn v2.13.0, ARTn, pARTn, OCINEB, IRA, NEB modularization, Highway SIMD, parallel dimer"
---

## [v2.13.0] - 2026-XX-XX

Saddle-search and NEB get new capabilities. ARTn (Activation-Relaxation
Technique nouveau) is wired in via the [pARTn](https://gitlab.com/mammasmias/artn-plugin)
Fortran library, both as a standalone method (`method = artn`) and as a
drop-in eigenmode backend (`min_mode_method = artn`). OCINEB (Off-Path
Climbing Image NEB) is the recommended hybrid method for automated
saddle-point refinement in production workflows
({cite:t}`neb-goswamiEnhancedClimbingImage2026`). IRA
(Iterative Rotations and Assignments) provides structure comparison and
SOFI point-group detection.

NEB is decomposed into a strategy-pattern: tangent, projection, spring
force, OCINEB controller, spline extrema, initial paths, and objective
function are now pluggable components. Parallel image evaluation via
C++20 threads lands in both NEB (per-image) and improved dimer (two
replicas). Highway SIMD is added as an optional subproject; hand-written
kernels for Morse/LJ/EAM pair loops are staged for a follow-up release.

The `artn-plugin` subproject tracks upstream
`mammasmias/artn-plugin` directly; consumers no longer need to pin the
HaoZeke personal fork.

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
```
