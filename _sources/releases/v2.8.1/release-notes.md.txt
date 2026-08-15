---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.8.1, covering NEB endpoint minimization, metatomic device selection, and bug fixes."
    "keywords": "eOn v2.8.1 release notes, new features, bug fixes"
---

# Release notes

## [v2.8.1] - 2025-11-03

This release primarily includes bug fixes for Nudged Elastic Band (NEB)
calculations and updates major dependencies, including PyTorch.

### ✨ New Features

* **NEB Endpoint Minimization:** Added a new parameter,
  `minimize_endpoints_for_ipath` (`nebMinimEPIpath`), to allow for endpoint
  minimization even when a full initial path is provided via `initial_path_in`.
* **Metatomic Potential:** Updated to use the newer
  `metatomic_torch::pick_device` API for device selection.

### 🐛 Bug Fixes

* **NEB:** Fixed an issue where `neb.dat` could be incorrectly overwritten. The
  final `neb.dat` file is now generated separately from intermediate iteration
  files (e.g., `neb_000.dat`).
* **NEB:** Prevented zero projected forces from being reported in
  `NudgedElasticBandJob` results when the true force is non-zero.
* **NEB:** The search for the highest energy image now correctly iterates over
  intermediate images only.
* **Matching:** Fixed a bug in the `helper_functions::identical` by correcting a
  loop boundary condition to avoid iterating past the number of atoms.
