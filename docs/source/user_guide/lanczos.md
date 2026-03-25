---
myst:
  html_meta:
    "description": "Guide to the Lanczos eigenmode method in eOn for saddle point searches."
    "keywords": "eOn Lanczos, lowest eigenmode, saddle point search, min-mode following"
---

# Lanczos

```{versionchanged} 2.12
Force call tracking fixed to use the job's potential instance.
```

The Lanczos method determines the lowest curvature mode of the potential energy
surface using the iterative Lanczos algorithm
{cite:t}`lcz-malekDynamicsLennardJonesClusters2000`. It is an alternative to
the [dimer method](project:dimer.md) for min-mode following saddle searches.

## When to Use Lanczos vs Dimer

Both methods find the same lowest eigenmode with comparable cost. Either can
be used with OCINEB climbing image refinement and other min-mode following
saddle searches. The choice is largely a matter of preference:

- **Lanczos**: 1 gradient evaluation per iteration. Natural for large systems.
- **Improved Dimer**: 2 evaluations per iteration (finite difference). Integrates
  with the [AtomicGPDimer](project:dimer.md) for GP-accelerated searches.

## Usage

The Lanczos method is selected in the saddle search configuration:

```{code-block} ini
[Saddle Search]
min_mode_method = lanczos

[Lanczos]
tolerance = 0.01
max_iterations = 20
```

`max_iterations` controls how many Lanczos iterations are performed per
eigenmode computation. Unlike the dimer's `rotations_max`, each Lanczos
iteration costs exactly 1 gradient evaluation.

## Configuration

```{code-block} ini
[Lanczos]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.LanczosConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: LCZ_
keyprefix: lcz-
---
```
