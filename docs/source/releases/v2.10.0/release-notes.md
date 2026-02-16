---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.10.0, covering displacement documentation, ASV benchmarks, and Eigen copy elimination."
    "keywords": "eOn release notes, displacement scripts, ASV benchmarks, performance"
---

# Release notes

## [v2.10.0] - 2026-02-15

This release improves documentation, developer infrastructure, and runtime
performance. There are no breaking API changes.

### Documentation

#### Displacement Strategies Guide

The {doc}`../../user_guide/saddle_search` page now includes a conceptual
explanation of how epicenter selection works during saddle searches, covering
the weight-based probabilistic strategies (`displace_random_weight`,
`displace_listed_atom_weight`, etc.), the `displace_radius` /
`displace_magnitude` mechanism, and the `displace_all_listed` flag.

#### Displacement Scripts Tutorial

A new {doc}`../../tutorials/displacement_scripts` tutorial walks through two
worked examples of targeted displacement:

- **Vacancy diffusion in Cu** using OVITO's Polyhedral Template Matching (PTM)
  to identify non-bulk atoms near a vacancy (`ptmdisp.py`).
- **Adsorbate on a catalyst surface** using ASE to select atoms by element or
  z-coordinate and expand the selection by a distance cutoff
  (`adsorbate_region.py`).

The tutorial also covers the script interface contract (positional `.con` path
argument, comma-separated indices on stdout, per-state caching), the static
`displace_atom_list` alternative, and client-side `listed_atoms` displacement.

#### Schema Descriptions

The `[Saddle Search]` configuration fields `displace_atom_kmc_state_script`,
`displace_all_listed`, `displace_atom_list`, and `client_displace_type` now
have detailed descriptions that appear in the auto-generated documentation
tables.

### GPR-Dimer

- Exposed `gprd_linalg_backend` build option for selecting the linear algebra
  backend used by gpr_optim: `eigen` (default), `cusolver` (NVIDIA GPU),
  `kokkos` (portable CPU/GPU), or `stdpar` (C++17 parallel STL with TBB).
- Updated pinned gpr_optim to latest commit with new backends, subproject build
  fix, and performance improvements.

### Performance

- Eliminated unnecessary Eigen matrix copies in `Matter`, `Potential`, and
  `HelperFunctions` hot paths by using const references and move semantics.

### Developer

- Added ASV benchmark CI workflow with `asv-spyglass` for automated PR
  performance comparison.
- Expanded ASV benchmark suite with point evaluation, LJ minimization, and NEB
  workloads.
- Added macOS arm64 to the metatomic CI matrix using Homebrew gfortran.
- Reduced build times by linking to xTB by default.
- Build cleanup for Windows compatibility.
- Refactored `MetatomicPotential` variant resolution to use upstream
  `metatomic_torch::pick_output`.

### Bug Fixes

- Fixed Windows `STATUS_STACK_OVERFLOW` (`0xC00000FD`) crash in the EAM Al
  potential.  The Fortran subroutine `gagafeDblexp` allocates ~3.6 MB of local
  arrays (`phi`/`phivirst` dimensioned at `MAXPRS=200000`), exceeding the 1 MB
  Windows default stack.  The linker now requests a 16 MB stack on Windows.
- Fixed silent client failure on Windows when stdout is redirected by the Python
  server, by switching spdlog from `stdout_color_sink_mt` (which uses
  `WriteConsole`) to the plain `stdout_sink_mt`.
- Use Goswami & Jonsson 2025 for removing rotations through projections.
