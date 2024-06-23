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

In order to run a nudged elastic band calculation, set **job** to *nudged_elastic_band* in the
**[Main]** section. Details of the optimizer can be set in the **[optimizer]** section.

## Variants

```{versionadded} 2.0
The energy weighted varying springs method of {cite:t}`neb-asgeirssonNudgedElasticBand2021` has been added.
```

- Classic nudged elastic band of {cite:t}`neb-millsQuantumThermalEffects1994` and {cite:t}`neb-schenterReversibleWorkBased1994`.
- Improved tangent method of {cite:t}`neb-henkelmanImprovedTangentEstimate2000`.
- Climbing image NEB of {cite:t}`neb-henkelmanClimbingImageNudged2000`.
- Doubly nudged method of {cite:t}`neb-trygubenkoDoublyNudgedElastic2004`.

## Configuration

```{eval-rst}
.. autopydantic_model:: eon.schema.NudgedElasticBandConfig
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
