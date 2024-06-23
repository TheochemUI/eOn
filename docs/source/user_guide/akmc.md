# Adaptive Kinetic Monte Carlo

Unlike traditional Kinetic Monte Carlo (KMC) methods, adaptive (a), or
off-lattice KMC {cite:p}`ak-trochetOffLatticeKineticMonte2020` doesn't require a
lattice of predetermined points in configuration space. The method dynamically
builds a rate table based on local atomic movements and energy barriers. This
allows for more accurate modeling of complex systems where atomic positions are
not fixed to a grid. The confidence scheme and thermal accessibility settings in
AKMC ensure that transitions between states are based on well-defined criteria,
enhancing the precision of simulations. 

## Configuration

```{code-block} ini
[akmc]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.AKMCConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: AK_
keyprefix: ak-
---
```
