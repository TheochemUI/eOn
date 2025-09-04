# Hyperdynamics

The hyperdynamics method uses a bias potential which should be zero at
transition states and positive in minima in order to accelerate the rate of
transitions as noted by {cite:t}`hd-voterHyperdynamicsAcceleratedMolecular1997`.

The hyperdynamics time step {math}`\delta t` can be obtained from the molecular
dynamics simulation time step {math}`\delta t^b` multiplied by a boost factor
{math}`e^{\beta \Delta V}`, where {math}`{\Delta V}` is the bias potential.

There are several possible forms of bias potential. In EON, we have implemented
the bond-boost method of {cite:t}`hd-mironAcceleratedMolecularDynamics2003`
where the bias potential is controlled by the maximal (fractional) change in any
bond length in the system. This is a good bias potential for systems in which
the dynamics is governed by bond breaking and forming events.

You can run a hyperdynamics job by:
- Setting the {any}`bias_potential` option.
- Within a ``parallel_replica`` {any}`job`.

## Configuration

```{code-block} ini
[Hyperdynamics]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.HyperdynamicsConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: HD_
keyprefix: hd-
---
```
