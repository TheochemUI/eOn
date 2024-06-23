# Basin Hopping

Basin hopping is a Monte Carlo method in which the energy of each configuration
is taken to be the energy of a local minimum
{cite:p}`walesGlobalOptimizationBasinHopping1997`.

At each basin hopping step the client will print out the current energy
(current), the trial energy (trial), the lowest energy found (global min), the
number of force calls needed to minimize the structure (fc), the acceptance
ratio (ar), and the current max displacement (md).

## Notes

- `eON` defaults to letting displacements occur from minimized structures as per the method of {cite:t}`whiteInvestigationTwoApproaches1998`.
- The occational jumping variant of {cite:t}`iwamatsuBasinHoppingOccasional2004` is also implemented, and is controlled by {any}`eon.schema.BasinHoppingConfig.jump_max`.

## Configuration

```{eval-rst}
.. autopydantic_model:: eon.schema.BasinHoppingConfig
```

## References

```{bibliography}
```
