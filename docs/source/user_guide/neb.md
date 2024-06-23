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

## Variants

- Classic nudged elastic band of {cite:t}`neb-millsQuantumThermalEffects1994` and {cite:t}`neb-schenterReversibleWorkBased1994`.
- Improved tangent method of {cite:t}`neb-henkelmanImprovedTangentEstimate2000`.
- Climbing image NEB of {cite:t}`neb-henkelmanClimbingImageNudged2000`.
- Doubly nudged method of {cite:t}`neb-trygubenkoDoublyNudgedElastic2004`.

```{versionadded} 2.0
- The energy weighted varying springs method of {cite:t}`neb-asgeirssonNudgedElasticBand2021`.
```

```{note}
`eON`, like many other codes after {cite:t}`neb-sheppardOptimizationMethodsFinding2008` uses one optimizer instance for moving the whole band of images.
```

```{versionadded} 3.0_TBA
Via the surrogate potential interface, a native C++ implementation of the Gaussian Process accelerated NEB first described in {cite:t}`neb-koistinenNudgedElasticBand2017` and {cite:t}`neb-koistinenNudgedElasticBand2019`.
```

## Configuration

```{code-block} ini
[Nudged Elastic Band]
```


```{versionchanged} 2.1_TBA
In TOML, this will be `[NEB]`
```


```{eval-rst}
.. autopydantic_model:: eon.schema.NudgedElasticBandConfig
```

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
