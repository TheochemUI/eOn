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

```{versionchanged} 2.15
ARTn is now loaded at runtime via `dlopen` / `LoadLibrary`, mirroring
the LAMMPS pattern. No `-Dwith_artn` build flag is needed and no
configure-time CMake build of `artn-plugin` is required. A single
eOn binary picks up ARTn iff `libartn.so` is on the library search
path. The pre-2.15 `meson subprojects download artn-plugin` workflow
is no longer required.
```

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
  Default ``-1`` keeps pARTn's own default. Set ``0`` to skip the push and go
  directly to Lanczos (appropriate when eOn supplies the displacement); larger
  values let ARTn explore further from the minimum before engaging Lanczos.
  Maps to pARTn's ``ninit``.
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
- **`filin`**: Path to an ``artn.in`` style input file read by pARTn during
  setup. Empty (the default) keeps pARTn's post ``artn_create`` sentinel
  ``BBBB``, i.e. no file is read and all parameters come from this ``[ARTn]``
  section. Any non-empty value is required to exist; eOn checks before setup
  and aborts with a clear error when it does not. Maps to pARTn's ``filin``.

## Runtime requirements

ARTn needs `libartn.so` (or `libartn.dylib` / `libartn.dll`) on the
library search path at run time. The shim and the SaddleSearch driver
are unconditionally compiled into eonclient; only the dynamic library
is supplied externally.

```{code-block} shell
# Build libartn from the upstream Fortran subproject (one-time):
git clone https://gitlab.com/mammasmias/artn-plugin
cd artn-plugin
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DWITH_LAMMPS=OFF -DWITH_QE=OFF -DWITH_SIESTA=OFF
cmake --build build -j

# Make eonclient see it:
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH      # Linux
export DYLD_LIBRARY_PATH=$PWD/build:$DYLD_LIBRARY_PATH  # macOS
```

If the library is missing, the eonclient banner still says
``ARTn: enabled (dlopen at runtime)`` but `method = artn` (or
`min_mode_method = artn`) raises a runtime error naming the entry
point and pointing at `LD_LIBRARY_PATH`.

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
