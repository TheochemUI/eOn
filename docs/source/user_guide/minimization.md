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
| **LBFGS** (default) | Most minimizations, fast convergence near minima | `lbfgs_memory` (default 20); optional `lbfgs_secant` / `lbfgs_precon` / `lbfgs_curvature` |
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
`reset`, the ASE wipe). `skip` drops a pair when `s·y < 0.2 sᵀB0s`.
`damped` is Powell's 1978 replacement `ŷ = θ y + (1-θ) B0 s` with
`θ = 0.8 sᵀB0s / (sᵀB0s − sᵀy)` when `sᵀy < 0.2 sᵀB0s` (Nocedal
and Wright, *Numerical Optimization*, §18.3). `cautious` is the
Li–Fukushima (2001) rule: keep the pair only if
`s·y ≥ ε ||s||² ||g||^α`. Isolated clusters can set
`lbfgs_project_rigid = true` to remove the six rigid-body modes
from the force and the step (those modes are Hessian-null; they
should not enter the limited-memory pairs).

Further L-BFGS knobs, all optional and off by default:

- `lbfgs_secant = zhangxu` replaces `y` by the Zhang–Deng–Chen /
  Zhang–Xu modified secant, which folds the two function values into
  the pair (JOTA 1999 / 2001). No extra potential call. On
  Lennard-Jones and Morse pair potentials this pair inflates the
  force-call count; leave it at `standard` there.
- `lbfgs_precon` selects the two-loop `H0` metric.
  `none` is `H0 I`.
  `exp` / `c1` is the Packwood–Kermode–Csanyi pair Laplacian
  (J. Chem. Phys. 2016, doi:10.1063/1.4947024). ASE falls back
  below about 100 atoms; `exp` on LJ38 costs force calls.
  `lindh` is the Lindh stretch model Hessian
  (Chem. Phys. Lett. 1995, doi:10.1016/0009-2614(95)00646-L):
  `k = μ exp(α (r_nn² - r²))` along each pair, `α = lbfgs_precon_A`.
  `pair` is the Mones–Ortner–Csanyi positive-definite pair Hessian
  (Sci. Rep. 2018, doi:10.1038/s41598-018-32105-x) of the running
  potential: analytic `V''` and `V'/r` of LJ (`4ε[(σ/r)¹²-(σ/r)⁶]`,
  `ε=σ=1`) or Morse_Pt (`D=0.7102`, `α=1.6047`, `r0=2.897`), with
  negative pieces dropped. That is the right metric when the
  energy is itself a pair potential (Morse bulk). `pair_abs` uses
  `|V''|` and `|V'/r|` so attractive LJ neighbours still contribute;
  on LJ38 it does not beat an unpreconditioned two-loop.

Newton and RFO are not L-BFGS. They live in
[xtsci-optimize](https://github.com/HaoZeke/xtsci-optimize) as
`Method::Newton`. Build eOn with `-Dwith_xtsci=true` and select them
from the INI:

```{code-block} ini
[Optimizer]
opt_method = xtsci
converged_force = 0.01

[Xtsci]
method = newton

[LBFGS]
lbfgs_precon = pair
```

`method = rfo` is Banerjee RFO. The pair Hessian is still built in eOn
(`lbfgs_precon`) and handed to Rust through `xts_minimize_hess`. The
engine also runs the other xtsci-optimize methods (`lbfgs`, `bfgs`,
`sr1`, `sr2`, `steepest`, `adam`, `pso`, and the NLCG conjugacies).
- `lbfgs_h0 = sy_yy` (Nocedal–Wright 7.20), `ss_sy` (Barzilai–Borwein
  2), or `adaptive` (the smaller of the two when both are positive).
- `lbfgs_accept = none` (default) takes the two-loop step, matching
  ASE and historical eOn. `energy` refuses a rise (OPTIM-style
  shrink; each trial is a potential call). `nonmonotone` is the
  Grippo–Lampariello–Lucidi window of the last five energies
  (SIAM J. Numer. Anal. 1986).
- `lbfgs_extra_updates = p` reuses the newest pair `p` extra times
  in the two-loop recursion (Al-Baali, 2000/2014).

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
