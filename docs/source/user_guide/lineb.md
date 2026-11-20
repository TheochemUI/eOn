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
| `[Lineb] energy` | Tunneling energy {math}`E` in eV |
| `[Lineb] beads` | Harmonic quantum TST replicas; `0` writes the path only (`log_rate_over_reactant` is NaN) |
| `[Lineb] molecular` | `false` expands atom masses to coordinates (adsorbate / fixed field); `true` uses molecule alignment and needs every atom free |

Frozen atoms are dropped from the path degrees of freedom when
`molecular = false`. A per-axis freeze (only some of x, y, z fixed) is
refused. With `molecular = true`, any fixed atom is refused.

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
field. A solver failure appears in the job log with the rgsaddle status.
