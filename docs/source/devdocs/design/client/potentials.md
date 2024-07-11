# Potential design

Potentials in `eON` were designed for inter-operability from the start, with the business end of the implementation handed off to a generic `C-style` call:

```{code-block} cpp
AtomMatrix force(long nAtoms, AtomMatrix positions,
                 VectorXi atomicNrs, double *energy, Matrix3d box);
```

```{versionchanged} 2.x
- To support surrogate potentials, the signature include `double *variance` which is typically `nullptr` for exact potentials.
- A factory function is used for thread-safety instead of the static raw pointer `getPotential`
- `static` usage has been curtailed for thread-safety
- `parameters` are no longer required to consrtuct a potential, just the type of potential
```

