# Optimizers

The key difference in design between the optimizers and Job / Potential is that,
similar to the Potentials, the inputs are always the same, and thus they are
implemented the in the same fashion, but with a simpler form, since we do not
keep count of Optimizers.

Note that by design, this means that the points of extension are (at runtime):
- Potentials
- Optimizers

So bindings can be used to provide either of these, while the Job types are
implemented as variants and cannot be overridden at runtime.
