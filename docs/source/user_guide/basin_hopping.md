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
