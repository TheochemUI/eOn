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
calculation, supply a path that crosses the target energy on both sides
of the saddle and set both
`linear_seed = true` and `endpoint_convergence = true`. The first option
joins the two energy crossings with a straight line. The second tests the
nudged force at each endpoint and the perpendicular action gradient at
each interior image. Each option can be used independently with rgsaddle
ABI minor 48 or higher.

The eOn job uses the number of frames in `initial_path` as the number of
LI-NEB images. For a ten-image calculation with a densely sampled MEP,
locate its two crossings of the target energy, interpolate ten frames
along the straight line between them, and pass those frames as the seed.
The accuracy of the crossings matters near the saddle. For the frozen
Cu(110) reconstruction, following the intrinsic reaction coordinate
from the saddle in both directions resolved the first pair of crossings
at 0.05 eV below the barrier.

For a frozen Cu substrate with two mobile hydrogen atoms of mass
1.00793 amu, the paper uses a spring of 20 Å⁻² and an endpoint coefficient
of 1000 eV⁻¹ Å⁻¹ for the dimensionless action. With eOn’s default
`hbar`, the corresponding input values are

```ini
[Lineb]
linear_seed = true
endpoint_convergence = true
spring = 1.28291
endpoint_kappa = 64.3993
opt_tol = 0.000643993
perp_tol = 0.000643993
```

Supply `initial_path` and an absolute `energy` from the same potential.
These numbers use the fact that every mobile coordinate has the same mass.
They convert the published controls and do not establish a converged rate.

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
