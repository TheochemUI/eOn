---
myst:
  html_meta:
    "description": "Guide to the Nudged Elastic Band (NEB) method in eOn for finding saddle points and minimum energy paths between known reactants and products."
    "keywords": "eOn NEB, Nudged Elastic Band, minimum energy path, saddle point, climbing image"
---

# Nudged Elastic Band

The nudged elastic band (NEB) is a method for finding saddle points and minimum
energy paths between known reactants and products. The method works by
optimizing a number of intermediate images along the reaction path. Each image
finds the lowest energy possible while maintaining equal spacing to neighboring
images. This constrained optimization is done by adding spring forces along the
band between images and by projecting out the component of the force due to the
potential perpendicular to the band.

Details may be found in {cite:t}`neb-jonssonNudgedElasticBand1998`,
{cite:t}`neb-sheppardPathsWhichNudged2011`, and
{cite:t}`neb-asgeirssonExploringPotentialEnergy2018`.

In order to run a nudged elastic band calculation, set **job** to
*nudged_elastic_band* in the **[Main]** section. Details of the optimizer can be
set as per the <project:optimizer.md> document.

```{tip}
For a step-by-step walkthrough contrasting ASE's NEB with eOn's advanced
features (energy-weighted springs, dimer refinement), see the
[atomistic-cookbook tutorial](https://atomistic-cookbook.org/examples/eon-pet-neb/eon-pet-neb.html)
on oxadiazole formation with a PET-MAD metatomic potential.
```

## Variants

- Classic nudged elastic band of {cite:t}`neb-millsQuantumThermalEffects1994` and {cite:t}`neb-schenterReversibleWorkBased1994`.
- Improved tangent method of {cite:t}`neb-henkelmanImprovedTangentEstimate2000`.
- Climbing image NEB of {cite:t}`neb-henkelmanClimbingImageNudged2000`.
- Doubly nudged method of {cite:t}`neb-trygubenkoDoublyNudgedElastic2004`.

```{versionadded} 2.0
- The energy weighted varying springs method of {cite:t}`neb-asgeirssonNudgedElasticBand2021`.
```

```{note}
`eOn`, like many other codes after {cite:t}`neb-sheppardOptimizationMethodsFinding2008` uses one optimizer instance for moving the whole band of images.
```

```{versionadded} 2.8
Via the surrogate potential interface, a native C++ implementation of the Gaussian Process accelerated NEB first described in {cite:t}`neb-koistinenNudgedElasticBand2017` and {cite:t}`neb-koistinenNudgedElasticBand2019`.
```

```{versionadded} 2.12
- Onsager-Machlup action-based NEB for minimum action paths.
- OCINEB (Off-Path Climbing Image NEB) {cite:t}`neb-goswamiEnhancedClimbingImage2026`: hybrid CI-NEB + Min-Mode Following with hessian eigenmode alignment for automated saddle point refinement.
- Parallel image force evaluation (requires TBB, `-Dwith_parallel_neb=true`).
- IDPP (Image Dependent Pair Potential) path initialization.
- Modular strategy pattern for tangent, projection, and spring force components.
```

### Onsager-Machlup NEB

The Onsager-Machlup variant replaces the standard spring force with an
action-based spring that adapts per-image based on the local force magnitude.
Enable with `onsager_machlup = true` in the NEB section.

### OCINEB (hybrid dimer refinement)

OCINEB {cite:t}`neb-goswamiEnhancedClimbingImage2026` activates a Min-Mode
Following (dimer) search on the climbing image after it stabilizes, using
hessian eigenmode alignment to refine the saddle point to higher accuracy
without additional NEB iterations. Enable with `ci_mmf = true`.

### Parallel evaluation

When compiled with TBB support (`-Dwith_parallel_neb=true`), image forces are
evaluated in parallel. Python-based potentials automatically fall back to serial
evaluation.

## Configuration

The NEB section can be specified in `config.ini`:

```{code-block} ini
[Nudged Elastic Band]
images = 7
converged_force = 0.01
climbing_image_method = true
```

Or programmatically via [rgpycrumbs](https://rgpycrumbs.rgoswami.me):

```python
from rgpycrumbs.eon.helpers import write_eon_config

config = {
    "Main": {"job": "nudged_elastic_band"},
    "Nudged Elastic Band": {
        "images": 7,
        "converged_force": 0.01,
        "climbing_image_method": True,
    },
}
write_eon_config(config, Path("config.ini"))
```

See {doc}`/tutorials/dict_config` for the full programmatic workflow.

```{code-block} ini
[Nudged Elastic Band]
```




```{eval-rst}
.. autopydantic_model:: eon.schema.NudgedElasticBandConfig
```

## Outputs

NEB writes the usual `results.dat`, `neb.dat`, and the final band `neb.con`.

With `write_movies = true` in `[Debug]`, eOn also writes per-iteration
`neb_path_*.con` movie files and `neb_maximage.con`. These `.con` outputs now
embed structured frame metadata via `readcon-core`, including fields such as
`energy`, `frame_index`, `neb_bead`, optional `neb_band`,
`reaction_coordinate`, `relative_energy`, and `parallel_force`.

The existing `neb.dat` and `neb_*.dat` outputs are still written and remain the
primary compatibility path for current plotting tools.

## Refinement

```{versionadded} 2.0
```

Far from the minimum energy path, second order optimizers like those using the
LBFGS may not be optimal. In these situations, to traverse uninteresting
sections of the potential energy surface rapidly, it is best to use an
accelerating optimizer like QuickMin to begin with and transition to LBFGS
later. To facilitate this, the `[Refine]` section has been introduced.

```{eval-rst}
.. autopydantic_model:: eon.schema.RefineConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: NEB_
keyprefix: neb-
---
```
