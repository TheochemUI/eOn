---
myst:
  html_meta:
    "description": "Consolidated release notes for eOn v2.16.0: in-process rgpot potential (NWChem/CPMD), Davidson min-mode solver, readcon 0.13 ConFileIO, Metatomic controls, PbcConvention."
    "keywords": "eOn release notes, v2.16.0, rgpot, RgpotPot, Davidson, readcon, Metatomic"
---

# Release notes

## [v2.16.0] - 2026-07-03

See the [CHANGELOG](project:../changelog.md) for the fragment-by-fragment list
after `cog bump` consumes `docs/newsfragments/`; this page consolidates
**user-facing and maintainer** highlights from commits since `v2.15.0`, so the
cut has one place to link from the GitHub Release and feedstock PR.

### Potentials

#### In-process rgpot potential (NWChem / CPMD)

Potential type **`RGPOT`** (build with `-Dwith_rgpot=true`) runs
[rgpot](https://github.com/OmniPotentRPC/rgpot)'s `NWChemPot` / `CPMDPot`
frontends inside `eonclient` and `dlopen`s the `libnwchemc` / `libcpmdc`
engine libraries — no sockets, no potserv. Configuration lives in the
`[RgpotPot]` INI section and the matching Python `RgpotPot` schema model
(backend, basis/theory/scf_type for NWChem, functional/cutoff_ry for CPMD,
charge, multiplicity, engine paths, input_block). Engine-path environment
overrides are backend-scoped (`NWCHEMC_LIBRARY` vs `CPMDC_LIBRARY`). Ships
with a CPMD BLYP example (`examples/rgpot_cpmd_blyp/`), nwchemc + cpmdc smoke
tests, a dedicated CI job, and user-guide pages
({doc}`/user_guide/rgpot_pot`, {doc}`/user_guide/rgpot_integration`).

#### Metatomic output, determinism, and rotation controls

- Explicit `energy_output` / `energy_uncertainty_output` keys for models with
  non-default output names (#215), non-conservative forces via
  `non_conservative` and `variant_force` / `force_output` (#296).
- `deterministic` / `deterministic_strict` knobs for reproducible PyTorch
  evaluation.
- Random SO(3) rotation per evaluation (`random_rotation`, #287) and
  approximate O(3) symmetrization by averaging over `n_symmetry_rotations`
  (#292).
- Builds against metatomic-torch >=0.1.15 / metatensor-torch >=0.10 with
  `sample_kind`-based per-atom outputs.

#### Morse hot path

The Morse pair force loop inlines the potential, uses inverse-r scaling, and
accumulates energy in a register (~60% of instruction count in `Morse::force`
on NEB Morse workloads before the rework).

#### Molecular QM potentials and PBC

NWChem socket, ASE-ORCA, and ASE-NWChem auto-disable PBC on `Matter` attach
and hard-fail if PBC is re-enabled (#188). ASE NWChem takes `mpi_launcher`
(`mpirun` or `srun`) for Slurm-friendly execution (#193).

### Saddle search and min-mode

- **Davidson** minimum-mode solver (`min_mode_method = davidson`): Ritz
  subspace with finite-difference Hessian-vector products and residual
  correction, as an alternative to dimer rotation minimization and Lanczos.
- Min-mode saddle search no longer reports `STATUS_GOOD` while the climb
  objective is unconverged (#20).
- Saddle/process search synthesize a missing `displacement.con` from
  `pos.con` plus `displace_magnitude` times a unit `direction.dat` mode (#79).

### I/O and Matter

- `ConFileIO` targets readcon-core **0.13**: `IoStatus` enum returns, bulk
  zero-copy builder maps, one-pass `writeNebPath`, and force+energy co-loading
  trusted without recomputation.
- `PbcConvention` (Legacy vs MinimumImage) for position wrapping (#176);
  `getForcesFree` is const and `setPositionsFree` applies PBC wrapping (#171).

### NEB and AKMC

- NEB with `initializer = file` loads endpoints from the path frames and no
  longer requires separate `reactant.con` / `product.con` (#278).
- MCAMC `energy_level` superbasin scheme tracks a true run-wide minimum energy
  (#212).

### Maintainer notes

- The `subprojects/rgpot.wrap` dependency is pinned to a released rgpot tag
  (previously floating `head`).
- New CI workflow `ci_rgpot.yml` builds `-Dwith_rgpot=true` and runs the rgpot
  test suite against rgpot's fake engines.
