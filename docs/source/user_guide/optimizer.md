---
myst:
  html_meta:
    "description": "Configuration options for the various geometry optimizers in eOn, including LBFGS, QuickMin, FIRE, CG, and SD."
    "keywords": "eOn optimizer, LBFGS, QuickMin, FIRE, CG, SD, geometry optimization"
---

# Optimizer

```{note}
The eOn optimizers skip a line search when taking the steepest-descent step.
```

There are several other ways in which the `eOn` implementations differ from a
standard optimizer software suite[^1] . Some prominent reasons are:
- They are written for the <project:neb.md>, whose energy surface has no
  closed form
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

Each of the optimizer methods have their own settings as well.

### LBFGS

```{code-block} ini
[LBFGS]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.LBFGSConfig
```

### QuickMin

```{code-block} ini
[QuickMin]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.QuickMinConfig
```

### FIRE

```{code-block} ini
[FIRE]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.FIREConfig
```

### CG

```{code-block} ini
[CG]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.CGConfig
```

### SD

```{code-block} ini
[SD]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.SDConfig
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
