---
myst:
  html_meta:
    "description": "Fixed-energy line-integral nudged elastic band (LI-NEB) in eOn through rgsaddle: optimal tunneling paths and optional harmonic quantum TST prefactors and actions."
    "keywords": "eOn LI-NEB, lineb, instanton, tunneling path, rgsaddle, Ásgeirsson."
---

# Line-integral NEB (LI-NEB)

`job = lineb` finds the optimal tunneling path at a fixed energy {math}`E`
with the line-integral nudged elastic band of Ásgeirsson, Arnaldsson and
Jónsson, [J. Chem. Phys. 148, 102334 (2018)](https://doi.org/10.1063/1.5007180).
The path is the mass-weighted band that minimizes the abbreviated action
{math}`W(E) = \int \sqrt{2(V-E)}\,ds`, with endpoints held on {math}`V = E`.

The search runs in rgsaddle through the lineb C API
(`rgsaddle_lineb_solve`). Build with `-Dwith_rgsaddle=true` against rgsaddle with C ABI 1.46 or compatible later minor.

Energies are in eV, lengths in Å and masses in amu. Defaults for
`hbar` and `boltzmann` match eOn's tunneling constants so the reported
temperature is in kelvin.

## Inputs

| File / key | Role |
|---|---|
| `[Main] job = lineb` | Select the job |
| `[Lineb] initial_path` | Multi-frame seed band (for example `neb.con`), at least three images, interior images above `energy` |
| `[Lineb] energy` | Required tunneling energy {math}`E` in eV |
| `[Lineb] beads` | At least two harmonic quantum TST replicas, or `0` for the path only (`log_rate_over_reactant` is NaN) |
| `[Lineb] molecular` | `false` expands atom masses to coordinates (adsorbate / fixed field); `true` uses molecule alignment and needs every atom free |

Frozen atoms are dropped from the path degrees of freedom when
`molecular = false`. A per-axis freeze (only some of x, y, z fixed) is
refused. Keep atom order, masses, cell, fixed masks and frozen substrate
positions consistent across the seed. Under periodic boundaries the job joins
adjacent images by their minimum-image displacements. Use a nonperiodic system
with all atoms free for molecular alignment.

Set `endpoint_kappa` to a positive value for the fixed turning-surface force
coefficient in Eqs. (19) and (21) of the paper. The coefficient multiplies
the energy error in eV to give a force in mass-weighted action units. Its
default of zero selects the gradient-normalized `endpoint_scale` restraint.
A positive value requires rgsaddle ABI minor 47 or higher.

For the initialization and convergence test used in the paper’s Cu(110)
calculation, supply a full path between the minima and set both
`linear_seed = true` and `endpoint_convergence = true`. The first option
joins the two energy crossings with a straight line. The second tests the
nudged force at each endpoint and the perpendicular action gradient at
each interior image. Each option can be used independently with rgsaddle
ABI minor 48 or higher.

Tolerances `opt_tol` and `perp_tol` are mass-weighted force norms. The
remaining keys (`spring`, `endpoint_scale`, `time_step`, `max_move`,
`turning_width`, `quadrature_pairs`, `max_steps`) match rgsaddle's
`LinebConfig` defaults unless you change them.

## Outputs

| File | Contents |
|---|---|
| `lineb.con` | Optimized band |
| `results.dat` | Action, both half periods and temperatures, optional `log_rate_over_reactant`, call counts |

For positive replica counts, `log_rate_over_reactant` is the natural
logarithm of the rate times the reactant partition function, {math}`\ln(k Z_r)`.
To compare with a rate in inverse seconds, subtract {math}`\ln Z_r` evaluated
at the same temperature, energy zero and rigid-mode convention, then convert
from the eOn time unit (about 10.18 fs). A path-only run writes NaN for this
field. A solver failure appears in the job log with the rgsaddle status. When the solver exhausts `max_steps`, the diagnostic reports both final force norms and the tolerance for each.
