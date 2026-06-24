---
myst:
  html_meta:
    "description": "Consolidated release notes for eOn v2.15.0 (next cut after v2.14.0): potentials, process search, NEB/MPI, and maintainer release process."
    "keywords": "eOn release notes, v2.15.0, runtime potentials, MPI, LAMMPS, release process, PyPI"
---

# Release notes

## [v2.15.0] - pending cut (post-v2.14.0 on `main`)

See the [CHANGELOG](project:../changelog.md) for the fragment-by-fragment list
after `cog bump` consumes `docs/newsfragments/`; this page consolidates
**user-facing and maintainer** highlights from commits since `v2.14.0` plus the
release-process PR, so the next cut has one place to link from the GitHub
Release and feedstock PR.

**Channel note:** `v2.14.0` may exist as tag/CHANGELOG without finished GH
tarball asset, PyPI, and/or feedstock. This cut is intended to **supersede** that
incomplete multi-channel state (see {doc}`/devdocs/release` §5). Optionally still
attach `eon-v2.14.0.tar.xz` to the old GH release for reproducibility; do not
re-publish a conflicting `2.14.0` wheel if you never owned that semver on PyPI.

### Potentials

#### Runtime-loaded Fortran potentials

Legacy Fortran potential modules can be loaded at runtime via `dlopen` instead of
only at link time. Configure a search path with `[Potential] potentials_path` so
site-specific or experimental pots can ship as shared objects without rebuilding
the full client.

Contributing commits / PRs:

- `feat(potentials): load Fortran potentials at runtime via dlopen` (#342 area)
- `feat(config): add [Potential] potentials_path config key`
- `fix(pot): RTLD_DEEPBIND prevents symbol hijack of Fortran helpers at runtime`

### Process search

#### Fixed atoms preserved after `displacement.con`

Process search restores fixed-atom rows after loading a displacement geometry, so
constrained atoms do not drift when the displacement file omits or resets them.

- `fix(process-search): restore fixed-atom rows after loading displacement.con` (#341)

### NEB and MPI

#### Per-image LAMMPS on private communicators

Parallel NEB no longer shares one LAMMPS world across images in a way that
corrupts neighbor/comm state; each image gets an isolated communicator/instance
where required.

- `fix(neb): isolate per-image LAMMPS instances on private MPI communicators` (#340)

#### MPI C API build path

The MPI integration path is revived on the C API (`with_mpi`), retiring reliance
on obsolete C++ MPI bindings that blocked modern toolchains.

- `feat(mpi): port MPI C++ bindings to the C API and revive the with_mpi build` (#339)

### Windows / conda-forge

win-64 packages enable **in-tree Fortran** (`compiler('fortran')` / m2w64 gfortran,
`-Dwith_fortran=true`, static objects, `/STACK:16M` via meson) per
[conda-forge/eon-feedstock#15](https://github.com/conda-forge/eon-feedstock/issues/15)
and [windows-compat-sci-cpp](https://rgoswami.me/posts/windows-compat-sci-cpp/).
CuH2 is enabled on win with the rest of in-tree Fortran (m2w64 gfortran).

### Maintainer / distribution

#### Release process, Nickel GHA, PyPI name, feedstock checklist

Documented and automated the full cut path: `cog` / towncrier lockstep,
`scripts/release_assert.py`, staged tag workflow (tarball + GH release + PyPI),
`ci/gha/*.ncl` generation (rgpot-style), incomplete-release recovery, conda-forge
feedstock steps, and explicit Doxygen deferral (Sphinx remains the docs gate).

**PyPI:** the project name `eon` on [pypi.org](https://pypi.org/project/eon/) is
already taken by **EoN** (*Epidemics on Networks*, unrelated). This release
targets a **distinct** distribution name **`eon-akmc`** (see `ci/gha/pypi.ncl` /
`pyproject.toml` `project.name` alignment in the process PR)—import package
layout stays under `eon/` in-tree; the **index** name is what trusted publishing
and `pip install` use.

Contributing fragment: `docs/newsfragments/343.dev.md` (release-process PR).

### Housekeeping (non-user-facing, optional mention)

- `chore: cleanup cruft`, `chore(lfs): stop tracking test artifacts`
- `docs(site): load antics tracker` (docs site only)

### How to cut (maintainer)

1. Pre-gates in {doc}`/devdocs/release` §1 (CI green, fragments, cookbook if
   blocking, incomplete-`2.14.0` decision recorded).
2. `cog bump --version 2.15.0` (or `--auto` if commit log warrants) on `main`.
3. `git push origin main --follow-tags` → `release.yml` (from `ci/gha/release.ncl`).
4. Feedstock PR with `eon-v2.15.0.tar.xz` sha256 (template in process PR scratch /
   release.md §4).
