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

A classic dimer (`[Dimer] improved = false`) reports `|C_tau|` as a known
first-step scale when the mode is negative. L-BFGS keeps `H0 = 1/|C_tau|`
when that step fits in `max_move`, and also when the mode is milder than
30 eV/Å² and the cap would otherwise set the length. A steeper mode whose
uncapped step exceeds `max_move` takes one finite-difference probe of the
effective force. Improved dimer and Lanczos do not report `|C_tau|`. Those
searches climb on the effective force, so the first step is the probe and
`max_move` is only the cap. Later steps still use the `[LBFGS]` secant scale.
A build that rotates the mode in rgsaddle uses the same rule.
`improved = false` reports `|C_tau|`, and the improved dimer and Lanczos
do not.

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

### rgmin backend

eOn's minimizers can run inside rgmin, the Rust minimizer library
(formerly xtsci-optimize), instead of eOn's own C++ code. Two keys select
it:

- ``[Optimizer] opt_backend = rgmin`` keeps ``opt_method`` and hands the
  step to the rgmin solver of the same name.
- ``[Optimizer] opt_method = xtsci`` runs whichever rgmin solver
  ``[Xtsci] method`` names, including those eOn has no native version of
  (``bfgs``, ``newton``, ``rfo``, ``dogleg``, ``fire2``, ...).

| ``opt_method`` | rgmin solver under ``opt_backend = rgmin`` |
|---|---|
| ``lbfgs`` | ``lbfgs`` |
| ``fire`` | ``fire`` |
| ``cg`` | ``polak_ribiere`` |
| ``sd`` | ``steepest`` |
| ``qm`` | ``quickmin`` |

```{code-block} ini
[Optimizer]
opt_method = lbfgs
opt_backend = rgmin
```

Both need a client that links rgmin. ``-Dwith_rgmin=auto`` (the
default) links an rgmin found by pkg-config and changes nothing until an
input selects it; ``-Dwith_rgmin=enabled`` makes a missing rgmin a
configure error. ``-Dwith_xtsci=true`` is the deprecated spelling of
``enabled``. The client needs rgmin 0.4.0 or newer, installed with

```{code-block} bash
cargo cinstall --release --features capi --prefix=$PREFIX
# cinstall regenerates rgmin.h; the C++ client needs the curated header
cp include/rgmin.h $PREFIX/include/
cp -r include/dlpack $PREFIX/include/
```

``eonclient --features`` lists ``rgmin minimizers`` when it is linked.
A client without rgmin rejects ``opt_backend = rgmin`` with an error
naming ``-Dwith_rgmin``.

``opt_backend`` only applies where the optimizer minimizes an energy
whose gradient is the force: relaxations in every job, NEB endpoint
minimization, the per-image IDPP path and BGSD. The NEB band, the
collective IDPP band and the min-mode translation step are not
minimizations of an energy, so they keep eOn's optimizer; the rgsaddle
section below covers the band.

Matter evaluates every energy and force rgmin asks for, so
``total_force_calls`` counts the same calls as with ``opt_backend = eon``.
The host loop applies ``converged_force``, ``max_iterations`` and
``max_move`` (as eOn's largest per-atom move) as before. rgmin reads
``time_step`` for ``fire`` and ``qm``, ``[LBFGS] lbfgs_memory`` for the
memory, and the ``[Xtsci]`` or ``[LBFGS]`` preconditioner and acceptance
keys described below. The other ``[FIRE]``, ``[CG]``, ``[QuickMin]``
and ``[SD]`` keys belong to eOn's own solvers and do not reach rgmin.

### Xtsci

``opt_method = xtsci`` selects the rgmin engine. The solver inside that
engine is ``[Xtsci] method``. The adapter holds an ``rgmin_solver_t``
session; each ``step()`` is one outer iteration so L-BFGS pairs and NLCG
conjugacy survive the host loop.

When ``eindir-core`` is on the pkg-config path, the same backend borrows
an ``eindir_objective_t`` around the ``ObjectiveFunction``. The adapter
checks the eindir ABI stamp, the gradient feature bit, and the semantic
descriptor (eV, angstrom, analytic forces) before dispatch.
``Optimizer::run`` then calls ``rgmin_minimize_eindir`` unless the step
needs a host Hessian. ``step()`` stays one ``rgmin_solver_t`` iteration.

``[Xtsci] qn_step``, ``precon``, and ``accept`` are the first-class
knobs for an L-BFGS session. They map onto
``rgmin_solver_set_qn_step``, ``rgmin_solver_step_hess_fg``, and
``rgmin_solver_set_accept``. Energy and forces share one
``potential->force()`` call. ``accept = none`` takes the
maxmove-clipped step (one geometry). ``nonmonotone`` is the Grippo
window used by the rematch INI. The pair / Lindh matrix is still
built in eOn. Two-loop plus that matrix is \(H_0 = P^{-1}\);
``qn_step = newton`` is regularized Newton on the same ``P``.
``method = newton`` is a different, memoryless Newton solver.

The native ``[LBFGS]`` optimizer is unchanged. ``lbfgs_step``,
``lbfgs_precon``, and ``lbfgs_accept`` still drive it, and the
rgmin adapter falls back to them when the ``[Xtsci]`` fields are
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

### Rgsaddle band and min-mode

rgsaddle, the Rust saddle-search library, steps the NEB band and rotates
the min-mode dimer when the client is built with ``-Dwith_rgsaddle=true``
(rgsaddle 0.2.0 or newer, found by pkg-config). There is no runtime
switch: in such a build the existing keys pick the rgsaddle path.

| Input | Build without rgsaddle | Build with ``-Dwith_rgsaddle=true`` |
|---|---|---|
| ``[Nudged Elastic Band] opt_method = lbfgs``, ``fire`` or ``xtsci``, non-geometric spring | eOn optimizer on the band objective | rgsaddle band session |
| ``[Nudged Elastic Band] opt_method = cg``, ``sd``, ``qm``, or a geometric spring | eOn optimizer | eOn optimizer |
| ``[Saddle Search] min_mode_method = dimer``, ``lanczos`` or ``xtsci`` | eOn dimer or Lanczos (``xtsci`` is an error) | rotation in an rgsaddle session |
| ``min_mode_method = davidson``, ``gprdimer`` | eOn | eOn |

``eonclient --features`` lists ``rgsaddle band and min-mode`` in such a
build. In the band session ``[Xtsci] method = fire`` is the band stepper
and ``method = lbfgs`` selects the L-BFGS session; ``opt_method = lbfgs``
and ``fire`` pick the matching stepper. ``[Optimizer] refine_*`` does not
apply to a band stepped by rgsaddle. Matter still owns the images, and
the NEB loop still owns climbing policy, resampling, and the spline.
The surface callback returns a status. The first band request carries
every image. Later requests carry the interior images, or one image.
When the coordinates already match, that image stays put and its energy
and forces stay with it. The images that move go out as one batch. The
client resends the band only after a reparameterization, or after
another move outside the session. Every image evaluation goes through
Matter, so ``force_calls_neb`` keeps its meaning.

``min_mode_method = xtsci`` rotates the lowest mode in that session
with a zero translation cap. The min-mode callback reads the centre
from its cached evaluation and runs a displaced probe on a copy. The
min-mode saddle search still climbs with the selected optimizer.
``min_mode_method = dimer`` and ``lanczos`` rotate in the same session,
``lanczos`` as a Lanczos session. Davidson and the GP dimer keep their
own solvers. Against rgsaddle ABI minor 5 or newer one session serves
every step of the search: the client moves it to each new centre with
the gradient it already has, so the centre costs no extra call, and the
previous mode seeds the next rotation.

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
