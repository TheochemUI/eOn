---
myst:
  html_meta:
    "description": "Instructions for running a local energy minimization calculation in eOn to find the nearest stable structure."
    "keywords": "eOn minimization, geometry optimization, energy minimization, optimizer"
---

# Minimization

Local minimization relaxes a structure to the nearest potential energy minimum
using one of several optimization algorithms.

To run a minimization, set **job** to *minimization* in the **[Main]** section:

```{code-block} ini
[Main]
job = minimization

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 1000
```

Or programmatically via [rgpycrumbs](https://rgpycrumbs.rgoswami.me):

```python
from rgpycrumbs.eon.helpers import write_eon_config

config = {
    "Main": {"job": "minimization"},
    "Optimizer": {"opt_method": "lbfgs", "converged_force": 0.01},
}
write_eon_config(config, Path("config.ini"))
```

The optimizer settings control convergence. See <project:optimizer.md> for the
full list of available optimizers and their parameters.

```{tip}
Enable `[Debug] write_movies = true` (and optionally `write_deprecated_outs =
true`) so `rgpycrumbs eon plt-min` can build energy profiles, 2D landscapes, and
convergence panels. For two endpoints, use **separate** landscape plots (one
`--job-dir` / `--label` each). See {doc}`/tutorials/systems/lj_minimization` and
{doc}`/tutorials/systems/index`.
```

## Optimizer selection

| Optimizer | Best for | Key parameter |
|---|---|---|
| **LBFGS** (default) | Most minimizations, fast convergence near minima | `lbfgs_memory` (default 20); optional `lbfgs_curvature` / `lbfgs_project_rigid` |
| **FIRE** | Systems far from equilibrium, tolerant of bad initial guesses | `time_step` |
| **CG** | Large systems where LBFGS memory is a concern | `cg_line_search` |
| **QuickMin** | Simple dynamics-based relaxation | `time_step` |
| **SD** | Debugging, guaranteed descent direction | `sd_alpha` |

## Convergence

The minimization converges when `converged_force` is met. The number is
always in eV/A; `convergence_metric` chooses **which** force scalar is
compared. The default is eOn's historical `norm`, not OPTIM's RMS.

- **norm** (default): Euclidean `||F||_2` of the free-atom force vector
  (`sqrt(sum_i F_i^2)`). This is **not** an RMS: it grows with system
  size. Chill / OPTIM `||F|| < 1e-2` is a different scalar.
- **rms**: `||F||_2 / sqrt(3 N_free)`, matching OPTIM MYLBFGS
  (`SQRT(SUM(VNEW**2)/(3*NATOMS))`). Opt-in for an OPTIM-comparable
  stop; it is not the eOn default.
- **max_atom**: maximum `||F||` on any single free atom
- **max_component**: maximum Cartesian component of any free-atom force

L-BFGS memory updates are controlled by `lbfgs_curvature` (default
`reset`, the ASE wipe). `skip` is the Li–Fukushima cautious rule
(keep older pairs when `s·y` is not safely positive). `damped` is
Powell's 1978 replacement `ŷ = θ y + (1-θ) B0 s` with
`θ = 0.8 sᵀB0s / (sᵀB0s − sᵀy)` when `sᵀy < 0.2 sᵀB0s` (Nocedal
and Wright, *Numerical Optimization*, §18.3). Isolated clusters can
set `lbfgs_project_rigid = true` to remove the six rigid-body modes
from the force and the step (those modes are Hessian-null; they
should not enter the limited-memory pairs).

## Refinement

For paths far from the minimum, a two-stage optimization can be faster: start
with QuickMin or FIRE and switch to LBFGS after the forces
drop below a threshold.

```{code-block} ini
[Optimizer]
opt_method = fire

[Refine]
opt_method = lbfgs
threshold = 0.5
```

This runs FIRE until the max force drops below 0.5 eV/A, then switches to LBFGS
for final convergence.

## Output

The minimization writes:
- `min.con`: the minimized structure
- `results.dat`: energy, force calls, convergence status

With `write_movies = true` (in `[Debug]`), `minimization.con` is written as a
concatenated structure movie (one frame per iteration). Each frame stores
structured JSON metadata in the second CON header line via `readcon-core`, including `energy`,
`frame_index`, `step_size`, and `convergence`.

Set `write_deprecated_outs = true` in `[Debug]` to also emit the legacy
`minimization.dat` sidecar during the compatibility window.

## Configuration

```{eval-rst}
.. autopydantic_model:: eon.schema.MinimizationConfig
```
