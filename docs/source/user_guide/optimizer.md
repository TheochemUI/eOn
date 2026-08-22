---
myst:
  html_meta:
    "description": "Configuration options for the various geometry optimizers in eOn, including LBFGS, QuickMin, FIRE, CG, SD, and xtsci."
    "keywords": "eOn optimizer, LBFGS, QuickMin, FIRE, CG, SD, xtsci, Newton, RFO, geometry optimization"
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

### Xtsci

``opt_method = xtsci`` selects the xtsci-optimize engine. The solver
inside that engine is ``[Xtsci] method``. Build with
``-Dwith_xtsci=true``. The adapter holds an ``xts_solver_t`` session;
each ``step()`` is one outer iteration so L-BFGS pairs and NLCG
conjugacy survive the host loop.

``[Xtsci] qn_step`` and ``[Xtsci] precon`` are the first-class knobs
for an L-BFGS session. They map onto ``xts_solver_set_qn_step`` and
``xts_solver_step_hess``. The pair / Lindh matrix is still built in
eOn. Two-loop plus that matrix is \(H_0 = P^{-1}\);
``qn_step = newton`` is regularized Newton on the same ``P``.
``method = newton`` is a different, memoryless Newton solver.

The native ``[LBFGS]`` optimizer is unchanged. ``lbfgs_step`` and
``lbfgs_precon`` still drive it, and the xtsci adapter falls back to
them when the ``[Xtsci]`` fields are left at their defaults.

```{code-block} ini
[Optimizer]
opt_method = xtsci

[Xtsci]
method = lbfgs
qn_step = newton
precon = pair
```

```{eval-rst}
.. autopydantic_model:: eon.schema.XtsciConfig
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
