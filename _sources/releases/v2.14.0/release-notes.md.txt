---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.14.0: structured readcon-core movie metadata, ARTn filin, batched parallel NEB forces, shared-instance threading fix, visualization tutorial."
    "keywords": "eOn release notes, readcon-core, trajectory metadata, ARTn filin, batched NEB forces, MyST-NB tutorial"
---

# Release notes

## [v2.14.0] - 2026-XX-XX

See the [CHANGELOG](project:../changelog.md) for the fragment-by-fragment
list; this page picks out the user-facing highlights.

### Trajectory output

#### Structured per-frame metadata in `.con` movies

Minimization, saddle-search, and NEB movies now embed structured
per-iteration metadata directly in each `.con` frame via
[readcon-core](https://github.com/lode-org/readcon-core) (bumped to
`0.8.0`). Energies, force norms, iteration indices, and
method-specific scalars ride with the frame they describe, so
downstream consumers can parse a single file instead of zipping the
movie against an external `.dat` table.

The legacy sidecar `.dat` files are kept behind
`write_deprecated_outs = true` for one release so pipelines that still
depend on them have a migration window. New code should read the
frame metadata from the `.con` itself. See the
[visualization tutorial](project:../../tutorials/visualization.md)
for a MyST-NB executable walkthrough.

Contributing commits:

- `feat(traj): structured per-iteration .dat output for movies`
- `fix(traj): embed movie metadata in .con frames`
- `fix(io): harden readcon movie outputs`
- `fix(io): keep multi-frame .con writers alive across appends`
- `test(io): cover compatibility outputs and NEB metadata`

#### Visualization tutorial

A new MyST-NB executable tutorial at
`docs/source/tutorials/visualization.md` walks through generating a
trajectory, plotting the energy landscape, and rendering NEB profiles
with `rgpycrumbs`, pinned to a validated plotting stack. Uses an HCN
NEB demo keyed to the `eon_orchestrator` configuration.

### Saddle search

#### ARTn `filin` input file

The `[ARTn]` section gains a `filin` key that maps to pARTn's
namelist filename. Empty (the default) preserves the post
`artn_create` `NAN_STR` (`"BBBB"`) sentinel, so pARTn reads no file
and every parameter comes from the eOn INI. Setting `filin` to a
path wires that file into pARTn setup; eOn now checks the file
exists up front and aborts with a clear message before `setup_artn`
runs, instead of letting the failure hide inside pARTn's generic
`ERR_FILE` code.

See the [ARTn user guide](project:../../user_guide/artn.md).

### Parallelism

#### Batched force evaluation for NEB

NEB parallel image evaluation switches from one-thread-per-image to
batched potential calls, cutting scheduling overhead for
medium-sized bands with cheap potentials. The batched path reuses
the per-image thread-safety contract established in v2.13.0 and does
not change the opt-in (`parallel = true` in `[Main]`).

Contributing commit:

- `feat: batched force evaluation and parallel NEB fix`

#### Shared-instance threading fix for legacy Fortran potentials

Legacy Fortran-backed potentials (EAM_Al, FeHe, Lenosky_Si, Morse_Pt,
SW_Si, Tersoff_Si, EMT, XTB when built without the reentrant flag)
are now gated out of the shared-instance parallel path. The prior
code assumed every `Potential` could be called from multiple threads
against a single instance, which broke potentials whose Fortran
backend carries module-level state. These potentials now fall back
to the serial code path under `parallel = true` rather than
corrupting internal buffers.

Contributing commit:

- `fix(parallel): gate shared-instance threading for legacy Fortran potentials`

### Parameters

#### JSON round-trip for `[Debug]` options

`ParametersJSON::to_json` always emitted a `Debug` section but
`from_json` had no matching handler, so `write_movies`,
`write_movies_interval`, and `write_deprecated_outs` were silently
dropped on any JSON deserialization. The round-trip is now
symmetric and the new `test_params_json` case enforces it.

Contributing commit:

- `fix(params): round-trip debug options through JSON`

### Build

#### cbindgen provisioning

`cbindgen` is now installed into the pixi env's `bin/` and wired
into the `setupeon` dependency graph, so `meson setup` on the
`readcon-core` subproject no longer fails with "Program 'cbindgen'
not found" on fresh clones.

Contributing commit:

- `build(pixi): install cbindgen into env bin and wire setupeon dep`

#### Metatomic runtime RPATH

`eonclient` now encodes the metatomic runtime library paths at link
time, so downstream tutorials and examples running against a pixi
dev env no longer need to hand-massage `LD_LIBRARY_PATH`.

Contributing commit:

- `fix(build): encode metatomic runtime paths in eonclient`

### Deprecations and breaking changes

- `write_deprecated_outs` keeps the legacy `.dat` sidecar tables for
  minimization and saddle-search for one release. The sidecars will
  be removed in v2.15.0; migrate to reading per-frame metadata from
  the `.con` movie output.

### Build requirements

- `readcon-core >= 0.8.0` (auto-resolved by the meson wrap).
- `cbindgen` on `PATH` at configure time (now provisioned by the
  pixi env for dev builds).
- All other requirements unchanged from v2.13.0.
