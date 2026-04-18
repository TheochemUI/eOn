---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.13.0: ARTn and IRA Fortran subprojects, OCINEB min-mode refinement, NEB strategy decomposition, Highway SIMD, and parallel image evaluation."
    "keywords": "eOn release notes, ARTn, pARTn, IRA, OCINEB, NEB strategy, Highway SIMD, std::thread"
---

# Release notes

## [v2.13.0] - 2026-XX-XX

See the [CHANGELOG](project:../changelog.md) for the full fragment-by-fragment
list; this page picks out the user-facing highlights.

### Saddle search

#### ARTn via pARTn

Added ARTn as a saddle-search method through the `artn-plugin` Fortran
library, exposed in two complementary modalities:

- **Standalone** (`method = artn`): pARTn drives the full push +
  eigenmode + perpendicular-relaxation cycle starting from the minimum.
  `displacement.con` is ignored; an optional initial mode from
  `direction.dat` biases the push.
- **Drop-in min-mode** (`min_mode_method = artn`): eOn's displacement
  and epicenter logic seeds the initial structure and mode, then pARTn
  takes over from the displaced configuration.

Configurable through the `[ARTn]` INI section with `push_step_size`,
`force_threshold`, `max_iterations`, and the pARTn tuning knobs
`ninit`, `nperp_limitation`, `lanczos_min_size`, `nsmooth`, and
`nnewchance`. Integer knobs with a `-1` default respect pARTn's own
upstream default; explicit values are forwarded via `set_param`.

Requires `-Dwith_artn=true` at build time. The subproject tracks
`mammasmias/artn-plugin` directly.

See the [ARTn user guide](project:../../user_guide/artn.md) for
details.

#### OCINEB hybrid (recommended)

OCINEB (Off-Path Climbing Image NEB) runs Min-Mode Following refinement
on the climbing image with Hessian eigenmode alignment, activated via
`ci_mmf = true`. For automated saddle-point refinement in production
workflows this is the recommended path. See Goswami, Gunde, Jónsson,
*Enhanced Climbing Image Nudged Elastic Band Method with Hessian
Eigenmode Alignment*, 2026,
[arXiv:2601.12630](https://arxiv.org/abs/2601.12630)
({cite:t}`neb-goswamiEnhancedClimbingImage2026`).

### NEB

#### Strategy-pattern decomposition

`NudgedElasticBand` is decomposed into modular components: tangent,
projection, spring force, OCINEB controller, spline extrema, initial
paths, and objective function. Each component is selectable via
config and covered by its own test suite.

#### Parallel image evaluation

NEB image forces now evaluate concurrently through `std::thread` (not
`std::jthread` -- Apple Clang's libc++ does not yet ship `jthread`).
Opt in via `parallel = true` in `[Main]` when the potential is
thread-safe or supports per-image instances. Achieves ~2.5x speedup
on a 5-image NEB with a Morse potential on the SVN reference setup.
Replaces the TBB-based `std::execution::par` prototype.

#### Initial paths

IDPP (Image Dependent Pair Potential) path initialization is
available via `neb.initial_path.method = idpp`. A collective IDPP
variant optimizes all movable atoms simultaneously. Optional
oversampling + cubic-spline decimation restores IDPP smoothness
after down-sampling.

#### Onsager-Machlup action paths

Minimum-action paths via the Onsager-Machlup action functional are
available with `onsager_machlup = true`. Spring stiffness can be
optimized automatically (`om_optimize_k = true`).

### Dimer

#### Parallel improved dimer

The two dimer replicas in `ImprovedDimer` evaluate forces
concurrently through `std::thread` when the potential is
thread-safe or supports per-image instances. Opt in via
`parallel = true`.

### Structure comparison

#### IRA (Iterative Rotations and Assignments)

`IRACompare` exposes `match()` (CShDA + SVD alignment),
`matchPBC()` (periodic boundary conditions; note this runs
`cshda_pbc` only, assignment without rotation/SVD), and
`findSymmetry()` (SOFI point-group detection). Requires
`-Dwith_ira=true` at build time.

### Infrastructure

#### Highway SIMD subproject

Highway is added as an optional cmake wrap subproject. When
available, `-DWITH_HIGHWAY` is set so future potential kernels can
opt in. Hand-written SIMD kernels for Morse/LJ/EAM pair loops are
staged for a follow-up release; today the main benefit is
compile-time availability of the subproject for downstream
development.

#### Eigen Fortran layout helpers

New `AtomMatrixF` type alias in `Eigen.h` plus
`from_fortran_layout_vector` covers the zero-copy
RowMajor-`(nat, 3)`-to-ColumnMajor-`(3, nat)` interop that both
ARTn and IRA rely on.

#### Upstream artn-plugin

The `artn-plugin` meson subproject now tracks
`mammasmias/artn-plugin` directly. Downstream builds no longer need
to pin the HaoZeke personal fork. The pARTn `get_error` C wrapper
is used for richer diagnostic messages; older libartn builds that
predate the wrapper fall back to `has_error` via `get_data`.

#### Test infrastructure

- `tests/test_one_pt.py` fixtures gain ARTn and IRA coverage.
- Catch2 test binaries `test_artn`, `test_ira` added.
- JSON serialization for Parameters (`load_json`, `to_json`) via
  `nlohmann/json`.

### Build requirements

- Fortran compiler and LAPACK for ARTn/IRA.
- pARTn `>= c51217b4` on `mammasmias/artn-plugin:DEVEL` (includes
  the `meson dependency('openblas')` CI fix; auto-cloned via
  `subprojects/artn-plugin.wrap`).
- Rust `>= 1.88` and `cbindgen` for the existing `readcon-core`
  subproject (unchanged from v2.12.0).

### Deprecations and breaking changes

None at the eOn public INI/Python surface. Internally,
`MinModeSaddleSearch` moves from a direct member of
`SaddleSearchJob` to a `SaddleSearchMethod` base-class pointer so
ARTn can share the polymorphism slot.
