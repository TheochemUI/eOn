---
myst:
  html_meta:
    "description": "Detailed release notes for EON v2.8.2, covering Metatomic uncertainty features, build system updates, and dependency changes."
    "keywords": "EON v2.8.2 release notes, new features, enhancements, bug fixes"
---

# Release notes

## [v2.8.2] - 2025-12-01

This release focuses on integrating uncertainty quantification into the
Metatomic potential interface and modernizing the build and dependency
infrastructure.

### ✨ New Features

* **Metatomic Uncertainty:** Added support for checking per-atom energy
    uncertainty in Metatomic potentials.
    * Introduced the `uncertainty_threshold` parameter (default: `0.1`,
        corresponding to \~100 meV/atom).
    * If the model supports the `energy_uncertainty` output (e.g., PET-MAD),
        EON will now warn if atoms exceed this threshold.
    * The **variance** is now automatically populated using the mean of these
        per-atom uncertainties if available.
* **Pixi Build System:** Added new environments (`rel`, `eon`) and tasks to
    `pixi.toml` to streamline the release process and simplify local builds.

### 🚀 Enhancements

* **Device Selection:** The Metatomic potential now utilizes the
    `metatomic_torch::pick_device` API for more robust device selection (CPU/GPU)
    based on model capabilities.
* **Model Loading:** Updated the internal loading mechanism to use
    `metatensor_torch::Module` for better compatibility with recent Metatensor
    versions.

### 📦 Dependencies

* **PyTorch:** Updated to `v2.8`.
* **Metatomic/Metatensor:** Updated dependencies to `metatomic-torch >= 0.1.7`
  and `metatensor-torch >= 0.8.2` to support the new uncertainty APIs.
* **Pydantic:** Pinned specific versions to resolve linting and compatibility
  issues in the build environment.
