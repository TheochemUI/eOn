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

Run `eonclient --features` on the installed client. A build with
phase-resolved counts and rate orbit output reports

```text
  rgsaddle fixed-energy LI-NEB: enabled
  LI-NEB phase-resolved counts: enabled
  LI-NEB rate orbit output: enabled
  LI-NEB measured turning point: enabled
  LI-NEB image Hessians and rate counts: enabled
```

The search runs in rgsaddle. Build eOn with `-Dwith_rgsaddle=true` and
an rgsaddle library that exports `rgsaddle_lineb_solve`. Meson checks
the installed library during configuration. Look for the following result
in the setup log. A cached check also ends in `YES (cached)`.

```text
Checking for function "rgsaddle_lineb_solve" with dependency rgsaddle: YES
```

For phase-resolved counts, the library must also export
`rgsaddle_lineb_solve_counted`. For bead-resolved gradients and residuals,
use a matching header and library with `RGSADDLE_HAS_LINEB_ORBIT` and
`rgsaddle_lineb_solve_counted_orbit`. These checks distinguish builds
that share an ABI minor version but provide different entry points.
The last two entries report support for `rgsaddle_lineb_solve_counted_turning`
and `rgsaddle_lineb_solve_counted_rate`. Check the latter before submitting
a rate job with `interpolated_hessian = true`.

If Meson finds the pkg-config package but cannot link
`rgsaddle_abi_stamp`, inspect the linker error in
`meson-logs/meson-log.txt`. Use matching headers and a library built
for the C++ compiler's target and C library. For example, a system
rgsaddle build may require a newer glibc than the compiler in a
packaged environment can link against.

Energies are in eV, lengths in Å and masses in amu. Defaults for
`hbar` and `boltzmann` match eOn's tunneling constants so the reported
temperature is in kelvin.

## Build the Linux client from source

Start in the eOn checkout containing `examples/lineb-morse`. Open the
default Pixi shell as described in {doc}`../install/index`. In that shell,
install the rgsaddle C library under a private prefix. The revision below
is the one selected by the `LI-NEB` row of
`.github/workflows/ci_rg_backends.yml` in this checkout.

```bash
prefix="$PWD/build-lineb-deps"
mkdir -p "$prefix"
git clone https://github.com/HaoZeke/rgsaddle.git "$prefix/src/rgsaddle"
git -C "$prefix/src/rgsaddle" checkout --detach 9ac3607724ea97f40bb27c5172b2492e7ee782d6
cargo install cargo-c --version 0.10.23 --locked --root "$prefix/tools"
export PATH="$prefix/tools/bin:$PATH"
(cd "$prefix/src/rgsaddle" && cargo cinstall --locked --release --features capi \
  --prefix="$prefix" --libdir="$prefix/lib")
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="$prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
meson setup build-lineb --buildtype=release -Dwith_rgsaddle=true
meson compile -C build-lineb
export PATH="$PWD/build-lineb/client:$PATH"
eonclient --features
```

Keep this shell open for the Morse example so the client can load the
same library that Meson checked.

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

To try the rate diagnostic, set `beads = 8` in the copied `[Lineb]`
section and run the client in the same directory. The reference run
wrote eight frames to `lineb-orbit.con`. It used 35 rate evaluations,
bringing the solver total to 1,119 and the job total to 1,122.

For the image-Hessian approximation in section II.D of the 2018 paper,
set `interpolated_hessian = true` in `[Lineb]`. The solver computes
Cartesian Hessians by central differences at the optimized path images
and interpolates them to the rate beads. The default computes the
curvature at the beads. The choice is independent of the orbit clock.
Set `image_hessian_step` to a positive displacement in Å to control
the finite differences at the images. Zero selects the solver default.

The job records the chosen method in `lineb.hessian`. With the counted
rate interface, `lineb.rate_calls.curvature` and `lineb.rate_calls.bead`
separate curvature evaluations from bead energy/gradient evaluations.
These are subdivisions of `surface_calls.rate`. Include that rate total
once when comparing the cost of a reactive flux.

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

To select the imaginary-time propagation described in section II.C, set
`interpolated_force_clock = true` with a library that defines
`RGSADDLE_LINEB_INTERPOLATED_FORCE_CLOCK`. This uses a cubic interpolation
of the signed tangential image forces to determine the half-period,
temperature and bead positions. The default obtains them by quadrature.
The `lineb.replica_clock` field in the CON frame metadata records the selection.
With a library exporting `rgsaddle_lineb_solve_counted_turning`, the
`lineb.turning` object also gives the zero-based launch image, the
mass-weighted path length and return arc coordinate (√amu Å), the largest
Cartesian component of the displacement from the opposite image (Å), and
the return energy relative to the target and to the opposite image (eV).
These measurements show how closely the propagated orbit reaches the
optimized turning point.

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

For an unconstrained nonlinear molecule, `rgpycrumbs` can read the full
orbit through readcon and normalize the rate. The command below pins the
revision with native eOn orbit support.
Supply the optimized reactant in extended XYZ format with matching species,
atom order and isotope masses. Its spectrum JSON supplies the electronic
energy in `electronic_energy_eV` and the $3N-6$ positive internal frequencies
in `signed_wavenumbers_cm-1`. Use the same potential, energy zero and unit
conversion for the reactant and the instanton. The native frames record the
run's `lineb.hbar`, `lineb.boltzmann` and `lineb.molecular` settings.

```sh
uvx --from git+https://github.com/HaoZeke/rgpycrumbs.git@0c8cbb32a6060b5670af91a93d790f234e54ef9e \
    rgpycrumbs eon molecular-instanton-rate lineb-orbit.con \
    reactant_spectrum.json reactant.extxyz --confirm-same-energy-zero
```

The `log10_rate_per_s` output includes the finite-bead harmonic reactant
partition function and the rotational inertia ratio for a single reaction
path. Apply any symmetry-equivalent path multiplicity separately when
comparing to a total experimental or published rate.

The readcon field `lineb.spectrum` reports the number of retained ring
modes and the number with negative eigenvalues when a backend with
`RGSADDLE_HAS_LINEB_RATE_DIAGNOSTICS` evaluates the rate. It is null for a path-only run or a backend
that does not return the spectrum. A first-order instanton has one negative
mode. Check the force residual around the ring as well as the mode count
and rate when increasing the number of replicas.

For a rate job (`beads > 0`) with `RGSADDLE_HAS_LINEB_ORBIT`, open
`lineb-orbit.con` with readcon. Each frame holds one replica of the
full ring in continuous Cartesian coordinates, including the fixed atoms.
The fields `temperature_K` and `tau_over_beta` give the temperature and
the bead's imaginary time as a fraction of the full period.
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

To check both residuals from the saved gradients, run the following command
in the job directory. It reads the ring through readcon and uses twice the
recorded half-period. It requires no calculator or electronic evaluations.

```sh
uvx --from git+https://github.com/HaoZeke/rgpycrumbs.git@0c8cbb32a6060b5670af91a93d790f234e54ef9e \
    rgpycrumbs eon orbit-residual lineb-orbit.con --compare-discretizations
```

For the eight-bead Morse example, the largest finite-bead residual is
0.0229891 eV/Å, which agrees with the native `lineb.stationarity` field.
Fourier differentiation on the same eight coordinates gives 0.0234458 eV/Å.
