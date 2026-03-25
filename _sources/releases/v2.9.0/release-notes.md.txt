---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.9.0, introducing OCI-NEB, Sequential IDPP initialization, and Metatomic variants."
    "keywords": "eOn release notes, OCI-NEB, RONEB, IDPP, Metatomic, Onsager-Machlup"
---

# Release notes

## [v2.9.0] - 2026-01-27

This release introduces significant advancements in Nudged Elastic Band (NEB)
path initialization and saddle point refinement. The update features the novel
Off-path Climbing Image NEB (OCI-NEB) method, a comprehensive suite of Image
Dependent Pair Potential (IDPP) initializers, and enhanced support for
multi-headed Metatomic potentials.

### ✨ Major Features

#### Off-path Climbing Image NEB (OCI-NEB)

Standard Climbing Image NEB constrains the highest energy image to move along
the elastic band tangent. This release implements OCI-NEB (previously referred
to as RONEB), a hybrid method that integrates Min-Mode Following (MMF) directly
into the climbing image phase.

- **Mechanism:** When the climbing image force drops below a specified threshold
  (`trigger_force` or relative `trigger_factor`), the system switches to the
  dimer search. This allows the image to break orthogonality with the band and
  follow the true lowest eigenmode toward the saddle, even if the elastic band
  curvature deviates from the Minimum Energy Path (MEP).
- **Stability:** The implementation includes robust fallback strategies. If the
  dimer mode and the path tangent diverge significantly (controlled by
  `angle_tol`), the system penalizes the trust radius or reverts to standard
  CI-NEB.
- **Curvature Recovery:** The dimer method now caches and restores the
  configuration with the most negative curvature if the rotation fails to
  converge, ensuring partial results contribute to the optimization.

The modalities of this method are described in the [accompanying publication](https://arxiv.org/abs/2601.12630).

#### Advanced Path Initialization (IDPP & S-IDPP)

Linear interpolation often results in high-energy atomic overlaps. This version
provides a suite of interpolation strategies based on the Image Dependent Pair
Potential (IDPP).

- **Collective IDPP:** Solves the IDPP objective function for all images
  simultaneously using the global optimizer (e.g., L-BFGS).
- **Sequential IDPP (S-IDPP):** Grows the path sequentially from the reactant
  and product inward. This method proves superior for complex reaction
  coordinates by resolving clashes at the frontiers before interpolating the
  center.
- **ZBL Repulsion:** An optional Ziegler-Biersack-Littmark (ZBL) repulsive
  potential can now wrap the IDPP objective (`sidpp_zbl`). This prevents atomic
  fusion during the initialization of dense paths.
- **Oversampling:** The system can now generate an initial path with a higher
  density of images (e.g., 3x), relax them via IDPP, and decimate the path back
  to the target image count using cubic Hermite splines.

### 🚀 Enhancements

- **Metatomic Variants:** Added support for multi-headed machine learning
  potentials. Users may now specify `variant_base`, `variant_energy`, or
  `variant_energy_uncertainty` in the configuration to target specific output
  heads (e.g., `energy/pbe0` or `energy_uncertainty/ensemble`).
- **Onsager-Machlup Action:** Implemented spring dynamics based on the
  Onsager-Machlup action. This allows for adaptive spring constants
  (`om_optimize_k`) that scale with the local path curvature and force,
  improving resolution in curved regions.
- **Peak Analysis:** The NEB job now automatically detects local maxima along
  the spline. It writes these configurations (`peakXX_pos.con`) and estimates
  their reaction modes (`peakXX_mode.dat`), facilitating subsequent dimer
  searches.
- **Timing Reports:** The client now reports breakdown of Real, User, and System
  time in `results.dat` and standard logs.
- **Relative triggers:** The NEB now uses relative thresholds along with
  absolute ones for starting the climbing image.


### 🔧 Configuration Changes

New parameters are available in the `[Nudged Elastic Band]` block:

- `initializer`: Options include `linear`, `idpp`, `idpp_collective`, `sidpp`,
  and `sidpp_zbl`.
- `ci_mmf`: Boolean to enable OCI-NEB.
- `onsager_machlup`: Boolean to enable OM-based spring dynamics.
- `setup_mmf_peaks`: Boolean to toggle the writing of peak estimates.

Along with updates for the `[Metatomic]` variants.
