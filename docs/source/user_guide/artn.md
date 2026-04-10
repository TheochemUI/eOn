---
myst:
  html_meta:
    "description": "Guide to the ARTn (Activation-Relaxation Technique nouveau) saddle search method in eOn."
    "keywords": "eOn ARTn, pARTn, activation relaxation, saddle point search, min-mode following"
---

# ARTn (Activation-Relaxation Technique nouveau)

```{versionadded} 2.13
```

The Activation-Relaxation Technique nouveau (ARTn)
{cite:t}`artn-barkemaEventBasedRelaxationAtomistic1996,artn-maasMousseau2005` is
a saddle point search method that uses an explicit push phase followed by
Lanczos eigenvalue computation and perpendicular relaxation. It is implemented
via the [pARTn](https://gitlab.com/mammasmias/artn-plugin) Fortran library.

## How ARTn differs from Dimer/Lanczos

Unlike the [dimer](project:dimer.md) and [Lanczos](project:lanczos.md)
methods, which follow the lowest eigenmode from a given starting point, ARTn
actively *pushes* the system away from a local minimum before searching for the
saddle:

1. **Push phase**: displace atoms away from the minimum by `push_step_size`
2. **Lanczos phase**: compute the lowest eigenvalue via internal Lanczos
3. **Perpendicular relaxation**: FIRE relaxation orthogonal to the push direction
4. **Convergence check**: iterate until force threshold is met or saddle found

This push-based exploration can find saddles that gradient-following methods
miss, making ARTn complementary to Dimer/Lanczos rather than a replacement.

## Force call costs

Each ARTn step (one call to `artn_step`) requires exactly 1 external force
evaluation. However, ARTn typically needs more total steps than Lanczos because
of the push and perpendicular relaxation overhead.

Illustrative benchmark on LJ38 TS optimization (100 near-saddle structures,
OptBench):

| Method | Avg force calls | Median | Failed/100 |
|---------|-----------------|--------|------------|
| Lanczos | 238 | 178 | 0 |
| ARTn | 426 | 269 | 10 |
| Dimer | 523 | 359 | 0 |

On this benchmark ARTn is faster than Dimer when it converges but less robust
on non-periodic clusters (10% failure rate from the push phase causing "box
explosion" on systems without periodic boundary conditions). Treat these
numbers as workload-specific rather than universal performance guarantees; for
periodic bulk and surface systems, which are ARTn's design target, the failure
rate is expected to be lower.

## Usage

ARTn requires the `artn-plugin` Fortran subproject source to be available under
`subprojects/` (for example via `meson subprojects download artn-plugin`). When
present, eOn builds it via CMake at configure time as a workaround for Meson's
current Fortran submodule scanner limitations.

ARTn can be used in two ways:

**As a standalone method** (``method = artn``): ARTn starts directly from
``pos.con`` and handles its own push, internal Lanczos eigenmode estimation,
and perpendicular relaxation (FIRE). ``displacement.con`` is ignored in this
mode, and an optional initial mode from ``direction.dat`` biases the push
direction when present.

```{code-block} ini
[Saddle Search]
method = artn

[ARTn]
push_step_size = 0.3
force_threshold = 0.05
max_iterations = 500
```

**As a min-mode drop-in** (``min_mode_method = artn``): eOn's displacement and
epicenter logic seeds the initial structure and mode, then ARTn takes over from
the displaced configuration. This is useful when eOn's displacement strategy
(e.g. least-coordinated weighting) should control where the search starts.

```{code-block} ini
[Saddle Search]
method = min_mode
min_mode_method = artn

[ARTn]
push_step_size = 0.3
force_threshold = 0.05
max_iterations = 500
```

### Parameters

- **`push_step_size`**: Step size (Angstrom) for the initial push away from the
  minimum. Maps to pARTn's ``push_step_size``. Larger values explore further
  but risk instability on clusters.
- **`force_threshold`**: Force convergence criterion (eV/Angstrom). Maps to
  pARTn's ``forc_thr``.
- **`max_iterations`**: Maximum number of ARTn steps.
- **`ninit`**: Number of initial push steps before Lanczos eigenmode estimation.
  Default 0 skips the push and goes directly to Lanczos. Larger values let ARTn
  explore further from the minimum. Maps to pARTn's ``ninit``.
- **`nperp_limitation`**: Comma-separated integers controlling the maximum number
  of perpendicular relaxation steps per Lanczos cycle. ``default`` uses pARTn
  defaults (tuned for exploration from a minimum). ``-1`` for unlimited
  perp-relax (smooth potentials in refinement). ``20,30`` limits steps
  (recommended for corrugated/ML potentials). Maps to pARTn's
  ``nperp_limitation``.
- **`lanczos_min_size`**: Minimum Lanczos iterations before convergence check.
  Default ``-1`` uses pARTn default (3). Set to ``1`` for refinement near a
  known saddle. Maps to pARTn's ``lanczos_min_size``.
- **`nsmooth`**: Number of smooth interpolation steps between push direction and
  eigenvector. ``-1`` uses pARTn default. ``0`` disables smoothing. Maps to
  pARTn's ``nsmooth``.

## Build requirements

ARTn requires a Fortran compiler and LAPACK. Download the subproject first, then
let eOn build it automatically at configure time:

```{code-block} shell
meson subprojects download artn-plugin
meson setup builddir -Dwith_artn=true
```

Or with a prebuilt library:

```{code-block} shell
meson setup builddir -Dwith_artn=true \
  -Dartn_libdir=/path/to/lib -Dartn_includedir=/path/to/include
```

## Comparison with kart

The [kart](https://gitlab.com/groupe_mousseau/kart) kinetic ART code also
integrates pARTn, but differently: kart calls the Fortran API directly and
embeds LAMMPS force computation inside the ARTn loop. eOn uses only the C API
(`artn.h`) and treats ARTn as a black-box step function, keeping the potential
abstraction layer intact. This means eOn's ARTn works with any potential
backend (LJ, EAM, LAMMPS, XTB, metatomic, etc.) without modification.

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: ARTN_
keyprefix: artn-
---
```
