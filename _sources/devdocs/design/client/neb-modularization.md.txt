---
myst:
  html_meta:
    "description": "Design rationale for NEB modularization using strategy patterns, parallel evaluation, and eigenmode variants."
    "keywords": "eOn NEB, strategy pattern, modularization, parallel, eigenmode, variant"
---

# NEB modularization

The Nudged Elastic Band implementation was decomposed from a monolithic 1,225-line
file into modular, testable components using the Strategy pattern.

## Module structure

| Module | Lines | Purpose |
|---|---|---|
| `NudgedElasticBand.cpp` | ~530 | Orchestration: image loop, convergence, CI activation |
| `NEBTangent.cpp` | ~80 | Tangent estimation strategies |
| `NEBProjection.cpp` | ~80 | Force projection strategies |
| `NEBSpringForce.cpp` | ~140 | Spring force strategies (constant, energy-weighted, OM) |
| `NEBForceProjection.cpp` | ~120 | Perpendicular force, CI force, zero-translation |
| `NEBOcinebController.cpp` | ~180 | OCINEB hybrid dimer controller |
| `NEBObjectiveFunction.cpp` | ~60 | Optimizer adapter |
| `NEBSplineExtrema.cpp` | ~80 | Cubic spline peak detection |
| `NEBInitialPaths.cpp` | ~200 | Path initialization (linear, IDPP, file) |

## Strategy pattern

Three interchangeable strategies are built at NEB construction time via factory
functions:

### Tangent strategies (`NEBTangent.cpp`)

- **Improved tangent** (default): uses energy-weighted bisection at maxima/minima,
  standard interpolation elsewhere. Better behavior near transition states.
- **Old tangent**: simple finite-difference `tau = R[i+1] - R[i-1]`. Retained
  for backward compatibility.

### Projection strategies (`NEBProjection.cpp`)

- **Standard NEB**: projects out parallel spring force, keeps perpendicular true force.
- **Elastic band**: no projection (spring + true force).
- **Doubly-nudged elastic band (DNEB)**: adds perpendicular spring force component.
- **Onsager-Machlup**: modifies spring constant per-image based on force magnitude.

### Spring force strategies (`NEBSpringForce.cpp`)

- **Constant spring**: uniform `k` for all images.
- **Energy-weighted**: adjusts `k` based on image energy relative to endpoints.
  Higher `k` near the barrier, lower in basins.
- **Onsager-Machlup**: adaptive `k` based on local curvature for minimum-action paths.

## Climbing image and OCINEB

The climbing image (CI-NEB) modification activates after initial convergence:
the highest-energy image climbs along the band to the exact saddle point.

OCINEB (Off-Path Climbing Image NEB) {cite:t}`neb-goswamiEnhancedClimbingImage2026`
extends CI-NEB with Min-Mode Following (MMF) and hessian eigenmode alignment:
once the climbing image stabilizes, a dimer search refines the saddle point.
This is controlled by `NEBOcinebController` which monitors CI stability and
triggers a `MinModeSaddleSearch` when the angular tolerance is met.

## Eigenmode methods

The eigenmode estimation (for dimer/Lanczos) uses `std::variant` instead of
an abstract base class:

```cpp
using EigenmodeStrategy = std::variant<DimerStrategy, ImprovedDimerStrategy,
                                        LanczosStrategy, GPRDimerStrategy>;
```

This eliminates virtual dispatch overhead and enables value semantics. The
variant is constructed by `buildEigenmodeStrategy()` based on configuration.

## Parallel force evaluation

When built with `-Deon_parallel_neb=true` (requires TBB), NEB evaluates image
forces in parallel using `std::for_each` with `std::execution::par`. The
`Potential::isThreadSafe()` virtual method gates parallelism: Python-based
potentials (ASE, CatLearn) return `false` and fall back to serial evaluation.
