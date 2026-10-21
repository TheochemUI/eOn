---
myst:
  html_meta:
    "description": "Guide to the Basin Hopping global optimization method in eOn for exploring potential energy surfaces."
    "keywords": "eOn Basin Hopping, global optimization, Monte Carlo, potential energy surface"
---

# Basin hopping

Basin hopping is a Monte Carlo method in which the energy of each configuration
is taken to be the energy of a local minimum
{cite:p}`bh-walesGlobalOptimizationBasinHopping1997`.

At each basin hopping step the client will print out:
- the current energy (`current`),
- the trial energy (`trial`),
- the lowest energy found (`global min`),
- the number of force calls needed to minimize the structure (`fc`),
- the acceptance ratio (`ar`),
- and the current max displacement (`md`).

## Acceptance

The client minimizes each trial before the test. `de` is the quenched energy
difference between that minimum and the current minimum. During `steps`, an
uphill hop is accepted with probability `exp(-de/(kB*temperature))` when
`temperature` and `kB` are positive. `kB` is in eV/K, and `temperature` is in
kelvin. A temperature or `kB` that is not positive rejects the uphill hop.
Other hops are accepted. `quenching_steps` reject an uphill hop and skip
the factor.

## Notes

- `eOn` defaults to letting displacements occur from minimized structures as per
  the method of {cite:t}`bh-whiteInvestigationTwoApproaches1998`, which is
  controlled by {any}`siginificant_structure
  <eon.schema.BasinHoppingConfig.significant_structure>`
- The occasional jumping variant of {cite:t}`bh-iwamatsuBasinHoppingOccasional2004`
  is also implemented, and is controlled by
  {any}`jump_max <eon.schema.BasinHoppingConfig.jump_max>`.

## Stock quench

Leaving the optimizer block at its defaults (conjugate gradients, 1000
iterations, force tolerance 0.01, maximum move 0.2) selects an L-BFGS quench
for this job. That quench uses energy acceptance, a memory of at least 20,
automatic scaling, and damped curvature. An isolated cluster with no fixed
atoms projects out rigid translation and rotation. The convergence measure is
the root-mean-square force, and the trial step is 0.4. Each hop takes a
loose tolerance of 0.01. A hop that is still far above the best minimum
stops early. A hop whose energy is near the best, and the first quench,
then takes a tolerance of 0.001.

The temperature of a single-element cluster is seeded from one fifth of
the quenched well depth per atom, which is a reduced temperature near 1
for a Lennard-Jones cluster. The step length stays near its initial value.
Acceptance above the target at that step lowers the temperature, and
acceptance below the target at a short step raises it, so one temperature
setting can search a surface whose energy unit is not kelvin. A rejected
hop whose gap is still much larger than kT doubles kT. Once the best
energy has not improved for several steps, a single-element cluster
replaces some Cartesian kicks with a point-group projection or a rotation
of the atoms outside the median radius. A free cluster kept in a bounding
cell is still isolated for those moves. A cluster with more than one element
and `swap_probability` left at 0 swaps elements with probability one half.
A rigid water cluster does not swap, so each hydrogen stays on its oxygen.
`jump_steps` left at 0 takes five steps when a jump fires.

A cluster made only of oxygen and hydrogen, two hydrogens for each oxygen, is
moved as rigid bodies when every oxygen can be paired with two hydrogens
inside 1.6 Å. The pair is the one whose bond lengths sit nearest 0.9572 Å and
whose angle sits nearest 104.52 degrees, so a contact shorter than a bond
does not take a hydrogen from another molecule. The whole molecule is
translated by the configured displacement and rotated by 0.9 radians, and
the quench varies the center of mass and the angle-axis rotation. The Monte
Carlo temperature is seeded from the depth of the quenched well and stays
there while hops are still being accepted. It rises only when acceptance
falls below the target. Molecules that would sit on top of
each other are pushed apart before the quench. Any other geometry keeps the
Cartesian displacement.

## Configuration

```{code-block} ini
[Basin Hopping]
```


```{eval-rst}
.. autopydantic_model:: eon.schema.BasinHoppingConfig
```

## References


```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: BH_
keyprefix: bh-
---
```
