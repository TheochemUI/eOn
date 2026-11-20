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

## Run a Morse barrier example

From the repository root, copy the example to a run directory and invoke
the client built with rgsaddle and the `morse_pt` potential.

```sh
mkdir -p run/lineb-morse
cp examples/lineb-morse/{config.ini,seed.con} run/lineb-morse/
(cd run/lineb-morse && eonclient)
```

Two fixed Pt atoms sit 10 Å apart. A third Pt atom moves between two
equivalent binding sites. The three seed frames bracket the barrier at
the specified energy of −0.55436031 eV. The job writes the optimized
path to `run/lineb-morse/lineb.con`. With the counted solve API, the first frame
gives the action and the cost for the entire band:

```python
import readcon

frame = readcon.read_con("run/lineb-morse/lineb.con")[0]
print(round(frame.metadata["lineb"]["action"], 4))
print(frame.metadata["evaluation_counts"]["potential_calls"])
```

A reference run gave an action of 36.6004 and used 1084 potential
evaluations for the solve and three to write the image energies.
The example requests the path only. Its JSON rate field is null.

## Inputs

| File / key | Role |
|---|---|
| `[Main] job = lineb` | Select the job |
| `[Lineb] initial_path` | Multi-frame seed band (for example `neb.con`), at least three images, interior images above `energy` |
| `[Lineb] energy` | Required tunneling energy {math}`E` in eV |
| `[Lineb] beads` | At least two harmonic quantum TST replicas, or `0` for the path only (JSON rate is null) |
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
A positive value requires the rgsaddle fixed-kappa solve API.

For the initialization and convergence test used in the paper’s Cu(110)
calculation, supply a path that crosses the target energy on both sides
of the saddle and set both
`linear_seed = true` and `endpoint_convergence = true`. The first option
joins the two energy crossings with a straight line. The second tests the
nudged force at each endpoint and the perpendicular action gradient at
each interior image. Each option can be used independently with a library that exposes the rgsaddle
seed options.

The paper also uses velocity-projected velocity Verlet for optimizing the
band. Set `projected_verlet = true` to request that optimizer with
a library that defines `RGSADDLE_LINEB_PROJECTED_VERLET`. The default optimizer is FIRE 2.0.

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
`hbar`, the corresponding inputs for projected velocity Verlet are

```ini
[Lineb]
linear_seed = true
endpoint_convergence = true
projected_verlet = true
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
| `lineb.con` | Optimized band, per-image energy and job-wide force-call counts in the second-line JSON of each frame |
| `lineb-orbit.con` | Equal-time rate replicas, mobile-atom potential gradients, and finite-bead stationarity for a rate run with an orbit-capable backend |
| `results.dat` | Action, both half periods and temperatures, optional `log_rate_over_reactant`, call counts |

Read the counters from any frame with readcon. Every frame of the optimized
band repeats the same job totals and has its own `neb_bead` and image energy.

```python
import readcon

frame = readcon.read_con("lineb.con")[0]
ledger = frame.metadata["evaluation_counts"]
counts = ledger["surface_calls"]
backend = ledger["potential_calls"]
print(counts)
```

With `rgsaddle_lineb_solve_counted`, the six disjoint counts sum to `total`.
Each surface call evaluates one atomic configuration.

| Key | Evaluations counted |
|---|---|
| `seed` | Inspecting the supplied band and locating the initial energy crossings |
| `optimization` | Relaxing the band, including trial steps |
| `path_preparation` | Preparing the converged path and refining its turning points |
| `quadrature` | Integrating the action and building the imaginary-time table |
| `dynamics` | Integrating the independent constrained trajectory |
| `rate` | Evaluating the rate prefactor, zero for a path-only job |

The sibling `potential_calls` object records changes to eOn's potential counter
during `setup`, `solve` and `output`, plus their `total`. Caching can make
these counts differ from the solver's requests. An MEP computed by a separate
job has its own cost. Older rgsaddle libraries supply the combined
`optimizer_reported` and `measure_reported` fields, plus `dynamics`, `rate`
and the observed remainder in `unclassified`.

To measure work inside rgpot, set `result_cache` in `[Lineb]` to a
RocksDB directory and use an rgpot potential with cache support
(`-Drgpot:with_cache=true` for the bundled dependency). A new directory
starts cold. Reusing the directory lets rgpot serve matching geometries.
The readcon metadata field `evaluation_counts.rgpot_results` gives
`computed` and `served` for the job and for each of `setup`, `solve`,
and `output` under `phases`. Without an attached cache this field is null.
For cost comparisons, record whether the database was cold or reused.

For positive replica counts, `log_rate_over_reactant` is the natural
logarithm of the rate times the reactant partition function, {math}`\ln(k Z_r)`.
To compare with a rate in inverse seconds, subtract {math}`\ln Z_r` evaluated
at the same temperature, energy zero and rigid-mode convention, then convert
from the eOn time unit (about 10.18 fs). A path-only run writes NaN for this
field. A solver failure appears in the job log with the rgsaddle status. When the solver exhausts `max_steps`, the diagnostic reports both final force norms and the tolerance for each.

The readcon field `lineb.spectrum` reports the number of retained ring
modes and the number with negative eigenvalues when a backend with
`RGSADDLE_HAS_LINEB_RATE_DIAGNOSTICS` evaluates the rate. It is null for a path-only run or a backend
that does not return the spectrum. A first-order instanton has one negative
mode. Check the force residual around the ring as well as the mode count
and rate when increasing the number of replicas.

For a rate job (`beads > 0`) with `RGSADDLE_HAS_LINEB_ORBIT`, open
`lineb-orbit.con` with readcon. Each frame holds one replica of the
full ring in continuous Cartesian coordinates, including the fixed atoms.
The second-line JSON records the energy gradients of the mobile atoms at
those coordinates, in the order of their CON atom IDs. Returning replicas
reuse the measured outbound coordinates and gradients. The diagnostic
uses the potential evaluations already counted in the rate.

```python
import readcon

ring = readcon.read_con("lineb-orbit.con")
frame = ring[0].metadata
print(frame["lineb_orbit"]["mobile_atom_ids"])
print(frame["lineb_orbit"]["gradients_eV_per_angstrom"])
print(frame["lineb"]["stationarity"])
```

For replica $j$ and mobile coordinate $a$, the discrete residual is
$R_{ja}=\partial_a V(x_j)+m_a(P/\tau)^2(2x_{ja}-x_{j-1,a}-x_{j+1,a})$.
Indices wrap around the ring, $P$ is the number of beads, and $\tau$ is
twice the reported half period. `max_abs_eV_per_angstrom` and
`rms_eV_per_angstrom` summarize all mobile Cartesian components.
`max_bead` and `max_axis` are zero-based. `max_atom_id` identifies
the CON atom. Compare these values as the path images and bead count
increase. A continuum residual requires differentiating the path in
continuous imaginary time.
