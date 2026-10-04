---
myst:
  html_meta:
    "description": "Guide to the Nudged Elastic Band (NEB) method in eOn for finding saddle points and minimum energy paths between known reactants and products."
    "keywords": "eOn NEB, Nudged Elastic Band, minimum energy path, saddle point, climbing image"
---

# Nudged elastic band

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

To run a nudged elastic band calculation, set **job** to
*nudged_elastic_band* in the **[Main]** section. Details of the optimizer can be
set as per the <project:optimizer.md> document.

```{tip}
**Python API.** Prefer {doc}`pyeonclient` for in-process NEB: build a
`list[Matter]` from ASE images, then `NudgedElasticBand(path, params, pot).compute()`
— same shape as ASE's `NEB` + `LBFGS`, without a workdir.

For a full walkthrough (ASE NEB vs eOn energy-weighted springs + OCI dimer), see
the [atomistic-cookbook PET-MAD example](https://atomistic-cookbook.org/examples/eon-pet-neb/eon-pet-neb.html).

For a **built-in** Morse Pt NEB with current `rgpycrumbs eon plt-neb` 1D/2D
figures (full history, 1:1 reaction-valley panel, structure strip), see
{doc}`/tutorials/systems/morse_pt_neb`.
```

## Variants

- Classic nudged elastic band of {cite:t}`neb-millsQuantumThermalEffects1994` and {cite:t}`neb-schenterReversibleWorkBased1994`.
- Improved tangent method of {cite:t}`neb-henkelmanImprovedTangentEstimate2000`.
- Climbing image NEB of {cite:t}`neb-henkelmanClimbingImageNudged2000`.
- Doubly nudged method of {cite:t}`neb-trygubenkoDoublyNudgedElastic2004`.
- Solid-state band of {cite:t}`neb-sheppardGeneralizedSolidstateNudged2012`,
  enabled with `solid_state = true`. Interior images relax the cell with
  the atoms. The cell stays lower triangular, and the tangent uses one
  Jacobian for atomic displacements and cell strain
  ([doi:10.1063/1.3684549](https://doi.org/10.1063/1.3684549)).

```{versionadded} 2.0
- The energy weighted varying springs method of {cite:t}`neb-asgeirssonNudgedElasticBand2021`.
```

## Dimer seeds from every band peak

After the band is written, eOn walks every spline maximum. Each interior
maximum more than `mmf_peak_tolerance` above
the reactant is written as `peakNN_pos.con` plus `peakNN_mode.dat` (the
interpolated tangent). That is the native gen-dimer seed: run a
`saddle_search` from each pair. `setup_mmf_peaks` defaults to true; set it
false to skip the files.

```ini
[Nudged Elastic Band]
setup_mmf_peaks = true
mmf_peak_tolerance = 0.05
```

Isomer endpoints with scrambled atom order can be aligned before
interpolation with `match_endpoints = true` (default off). The `ira`
method rigid-rotates and permutes the reactant onto the product via
IRACompare. Hungarian assignment is not in-tree; `match_method =
hungarian` currently uses IRA.
`match_endpoints` and `match_method` are client keys.
Neither key is in `eon/config.yaml`, and the server exits on the unknown option.

```{note}
`eOn`, like many other codes after {cite:t}`neb-sheppardOptimizationMethodsFinding2008` uses one optimizer instance for moving the whole band of images.
```

```{versionadded} 2.8
Via the surrogate potential interface, a native C++ implementation of the Gaussian Process accelerated NEB first described in {cite:t}`neb-koistinenNudgedElasticBand2017` and {cite:t}`neb-koistinenNudgedElasticBand2019`.
```

```{versionadded} 2.12
- Onsager-Machlup action-based NEB for minimum action paths.
- OCINEB (Off-Path Climbing Image NEB) {cite:t}`neb-goswamiEnhancedClimbingImage2026`: hybrid CI-NEB + Min-Mode Following with hessian eigenmode alignment for automated saddle point refinement.
- Parallel image force evaluation (`[Main] parallel` on the client). In 2.12 it required TBB and `-Dwith_parallel_neb=true`; see the note under Parallel evaluation.
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

The Frontiers article
([doi:10.3389/fchem.2026.1807063](https://doi.org/10.3389/fchem.2026.1807063))
restores the climbing image only on positive curvature (Algorithm 1:
curvature $> 0$). On alignment failure it keeps the most-negative-curvature
point and raises the trigger with the linear penalty. That is
`ci_mmf_restore_unhelpful = false`, the default.

`ci_mmf_restore_unhelpful = true` also restores after an alignment reject
or a force increase **on the image the dimer moved**. After a downhill
walk `maxEnergyImage` can hop; the band CI force is then a neighbor and
is not the restore score. That is an opt-in deviation from the published
protocol. Do not turn it on by default.

`climbing_image_converged_only` (default true) compares the climbing image
to `converged_force`. That is not a band-wide certificate: a SIDPP path
that lands two images on the same point can put the climber on a
stationary artifact while the rest of the band is still at 2 eV/Å.
`climbing_image_band_slack` (default 10) refuses that report. The job is
not converged while any image exceeds slack times the tolerance. SIDPP
itself throws if adjacent images collapse below \(10^{-6}\) Å.
That key is not in `eon/config.yaml`.
The client reads it.
A file the server loads that sets the key exits on the unknown option.

```{code-block} ini
[Nudged Elastic Band]
climbing_image_method = true
climbing_image_converged_only = true
```

The client sets the climb flag only when climbing is active and the highest
interior image is strictly above the higher endpoint.
A fixed-cell band compares potential energy.
A solid-state band compares enthalpy.
A monotonic band keeps the spring on every interior image.
With `climbing_image_method` on, either force test activates climbing.
One test requires the band force to be under the initial force times `ci_after_rel`.
The other requires it to be under `ci_after`.
`ci_after` is a server option.
`ci_after_rel` is not in `eon/config.yaml`, so the server rejects a file that sets it.

### Zoom-NEB

Zoom-NEB packs the images onto a window around the climbing image once
that image index has stayed fixed for `zoom_ci_stability` iterations and
the band force is below `zoom_after`. A `zoom_after` of 0 uses ten times
`converged_force`. The window endpoints become the fixed ends of the
band, and the images are placed at equal arc length inside it.

`zoom_mode = auto` keeps the contiguous images whose energy is above
`zoom_alpha` of the barrier, measured from the lower endpoint. If that
test leaves only the climbing image, the window falls back to
`zoom_offset` images on each side. `zoom_mode = manual` uses the offset
directly. `zoom_interpolation` is `cubic` (default) or `linear`.

Zoom does not turn off climbing-image NEB or OCINEB. The dimer step is
skipped on the redistribution iteration so it runs on the packed band.
`zoom_max_iterations` limits the steps after the pack; 0 keeps
`max_iterations`.

```ini
[Nudged Elastic Band]
zoom_neb = true
zoom_mode = auto
zoom_alpha = 0.5
zoom_offset = 1
zoom_after = 0.0
zoom_interpolation = cubic
zoom_ci_stability = 5
ci_mmf = true
```

### Solid-state band

`solid_state = true` gives each interior image its own cell. The reactant
and product cells stay fixed. A proper rotation puts every cell into lower
triangular form before the band moves: the first lattice vector lies on x
and the second lies in the xy plane. That removes the three rotations of
the cell. The initial path for `initializer = linear` is a fractional
interpolation of those oriented endpoints. `initializer = file` keeps the
supplied frames after the same rotation.

The spring and the perpendicular force are evaluated in the joint metric of
{cite:t}`neb-sheppardGeneralizedSolidstateNudged2012`. `solid_state_weight`
multiplies the cell block of the Jacobian. The default, 1, gives a cell
strain the same weight as an atomic displacement.
`solid_state_pressure` is a hydrostatic pressure in eV/Angstrom^3 added to
that metric. Positive pressure favors a smaller cell. Zero leaves the
potential-energy surface.

The cell force is the stress tensor. A potential that implements the Cauchy
stress, with the sign `sigma = (1/V) dE/dε` for the right strain
`h <- h (I+ε)` at fixed fractional coordinates, is used directly: the band
reads the stress the potential returned with each image's force call, so an
iteration costs one call per moved image. Any other potential is
differentiated on the six lower strain components with a central difference
of step 1e-5. That costs 12 extra energy calls per image and iteration; they
go to the potential as one batch over the band, which calculator groups or a
batched model evaluate together.

`ci_mmf`, `zoom_neb`, `onsager_machlup`, `neb_doubly_nudged`, and
`neb_elastic_band` are refused. The climbing image itself is the joint-space
reflection of the force. Peak files still carry the atomic block of the tangent.

```ini
[Nudged Elastic Band]
solid_state = true
solid_state_weight = 1.0
solid_state_pressure = 0.0
```

### Parallel evaluation

```{versionchanged} 3.5.0
The image pool is plain `std::thread`, so GCC and Clang builds need no TBB.
`-Dwith_parallel_neb` is deprecated and has no effect.
```

`[Main] parallel` is a client key and defaults to true.
`eon/config.yaml` has no `parallel` key in `[Main]`.
Setting that key makes the server exit before the client runs.
With `parallel = true`, a band that reaches the projection step projects
each image's tangent, spring, and projected force on one pool. Dirty-image
force calls use that pool too. The pool starts once and stays for the
process. On Linux the pool sizes itself to the affinity mask, so a Slurm
or taskset limit sets the thread count. Other systems take the size from
`std::thread::hardware_concurrency()`. The caller is one of those workers,
so a 20-image band on 8 allowed cores keeps 8 threads busy. The pool
reports an error in one image after the other images finish. The pool
needs no extra library.
nvc++ builds with `-Dstdpar=cpu` / `-Dstdpar=gpu` use
`std::execution::par` for those force calls instead; see {doc}`stdpar`. Python-based potentials fall back to serial
evaluation unless they report thread-safe shared instances or per-image
copies. Morse and other host potentials stay on the CPU.

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
later. The `[Refine]` section does that switch.

```{eval-rst}
.. autopydantic_model:: eon.schema.RefineConfig
   :no-index:
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
