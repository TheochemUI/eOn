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

## Optimizer selection

| Optimizer | Best for | Key parameter |
|---|---|---|
| **LBFGS** (default) | Most minimizations, fast convergence near minima | `lbfgs_memory` (default 20) |
| **FIRE** | Systems far from equilibrium, robust for bad initial guesses | `time_step` |
| **CG** | Large systems where LBFGS memory is a concern | `cg_line_search` |
| **QuickMin** | Simple dynamics-based relaxation | `time_step` |
| **SD** | Debugging, guaranteed descent direction | `sd_alpha` |

## Convergence

The minimization converges when the force criterion is met. Three metrics are
available via `convergence_metric`:

- **norm** (default): root-mean-square force across all free atoms
- **max_atom**: maximum force on any single atom
- **max_component**: maximum force component (x, y, or z)

## Refinement

For paths far from the minimum, a two-stage optimization can be faster: start
with a robust optimizer (QuickMin or FIRE) and switch to LBFGS after the forces
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
structured JSON metadata on line 2 via `readcon-core`, including `energy`,
`frame_index`, `step_size`, and `convergence`.

Set `write_deprecated_outs = true` in `[Debug]` to also emit the legacy
`minimization.dat` sidecar during the compatibility window.

## Configuration

```{eval-rst}
.. autopydantic_model:: eon.schema.MinimizationConfig
```
