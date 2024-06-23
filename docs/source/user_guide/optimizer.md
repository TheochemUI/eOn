# Optimizer

```{note}
All the optimizers in `eON` **do not** use a line search to determine steepest descent!!
```

There are several other ways in which the `eON` implementations differ from a
standard optimizer software suite[^1] . Some prominent reasons are:
- They are meant to be used with the <project:neb.md> which does not
  have well defined closed form energy surface
   + Especially in the global band optimization approach {cite:p}`opt-sheppardOptimizationMethodsFinding2008`.
- They are specialized for atomic systems
- The optimizers **only** see the moving atoms, the frozen atoms are omitted
  before being passed to the optimizer

## Configuration

```{code-block} ini
[Optimizer]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.OptimizerConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: OPT_
keyprefix: opt-
---
```

[^1]: e.g. as may be found in `scipy` or `ceres` for instance
