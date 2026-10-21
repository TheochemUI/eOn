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

For min-mode saddle searches (for example, `min_mode_method = dimer`), the
dimer mode curvature (`|C_tau|`) is not used as a known first-step Hessian
scale for translation. The translation objective climbs on an effective
force, so the initial L-BFGS scale still follows the usual `[LBFGS]`
controls (`lbfgs_inverse_curvature`, `lbfgs_auto_scale`) and its own
finite-difference probe, with `max_move` as the final cap.

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

When ``eindir-core`` is on the pkg-config path, the same backend borrows
an ``eindir_objective_t`` around the ``ObjectiveFunction``. The adapter
checks the eindir ABI stamp, the gradient feature bit, and the semantic
descriptor (eV, angstrom, analytic forces) before dispatch.
``Optimizer::run`` then calls ``xts_minimize_eindir`` unless the step
needs a host Hessian. ``step()`` stays one ``xts_solver_t`` iteration.
LBFGS, CG, and the other built-in solvers are unchanged.

``[Xtsci] qn_step``, ``precon``, and ``accept`` are the first-class
knobs for an L-BFGS session. They map onto
``xts_solver_set_qn_step``, ``xts_solver_step_hess_fg``, and
``xts_solver_set_accept``. Energy and forces share one
``potential->force()`` call. ``accept = none`` takes the
maxmove-clipped step (one geometry). ``nonmonotone`` is the Grippo
window used by the rematch INI. The pair / Lindh matrix is still
built in eOn. Two-loop plus that matrix is \(H_0 = P^{-1}\);
``qn_step = newton`` is regularized Newton on the same ``P``.
``method = newton`` is a different, memoryless Newton solver.

The native ``[LBFGS]`` optimizer is unchanged. ``lbfgs_step``,
``lbfgs_precon``, and ``lbfgs_accept`` still drive it, and the
xtsci adapter falls back to them when the ``[Xtsci]`` fields are
left at their defaults.

```{code-block} ini
[Optimizer]
opt_method = xtsci

[Xtsci]
method = lbfgs
qn_step = newton
precon = pair
accept = nonmonotone
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

### Rgsaddle band and min-mode

``[Nudged Elastic Band] opt_method = xtsci`` with ``-Dwith_rgsaddle=true`` steps the band
through a rgsaddle session instead of the flattened band objective.
``[Xtsci] method = fire`` is the band stepper. ``method = lbfgs``
selects the L-BFGS session. Matter still owns the images, and the NEB
loop still owns climbing policy, resampling, and the spline.
``-Dwith_rgsaddle=true`` needs rgsaddle 0.2.0 or newer. The surface
callback returns a status. The first band request carries every image.
Later requests carry the interior images, or one image. When the
coordinates already match, that image stays put and its energy and
forces stay with it. The images that move go out as one batch. The
client resends the band only after a reparameterization, or after
another move outside the session.

``min_mode_method = xtsci`` rotates the lowest mode in that session
with a zero translation cap. The min-mode callback reads the centre
from its cached evaluation and runs a displaced probe on a copy. The
min-mode saddle search still climbs with the selected optimizer.
In a build with ``-Dwith_rgsaddle=true``, ``min_mode_method = dimer``
and ``lanczos`` rotate in the same session, ``lanczos`` as a Lanczos
session. Davidson and the GP dimer keep their own solvers.
Against rgsaddle ABI minor 5 or newer one session serves every step of the
search: the client moves it to each new centre with the gradient it
already has, so the centre costs no extra call, and the previous
mode seeds the next rotation.

The session takes its settings from the dimer and Hessian keys:

| rgsaddle setting | eOn key | Conversion |
|---|---|---|
| ``dr`` | ``[Main] finite_difference`` | Angstrom; 0.001 when not positive |
| ``rotation_tol`` | ``[Dimer] torque_min`` | none |
| ``rotation_angle_tol`` | ``[Dimer] converged_angle`` | degrees to radians (ABI minor 5) |
| ``max_rotations`` | ``[Dimer] rotations_max`` | at least 1 |
| ``difference`` | ``[Hessian] fd_scheme`` | ``one_sided`` is forward; ``central`` and ``fourth`` are central (ABI minor 5) |
| ``method`` | ``[Xtsci] method`` | ``lbfgs`` gives L-BFGS; any other gives FIRE, the fast inertial relaxation engine |
| ``krylov_dim`` | none | 12 |

rgsaddle has no fourth-order stencil, so ``fd_scheme = fourth`` gives
central differences here. ``[Dimer] rotation_angle`` is a trial angle
of the host dimer and does not reach the session.
