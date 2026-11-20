---
sd_hide_title: true
myst:
  html_meta:
    "description": "Line-integral tunnelling rates in eOn: energy-continuation LI-NEB paths, imaginary-time orbits, and harmonic quantum transition-state theory."
    "keywords": "eOn line-integral tunnelling, LI-NEB, HQTST, quantum transition-state theory, tunnelling rate."
---

# Line-integral tunnelling

`job = line_integral_tunneling` implements the LI-NEB path and harmonic
quantum transition-state theory (HQTST) of Ásgeirsson, Arnaldsson and Jónsson,
[*J. Chem. Phys.* **148**, 102334](https://doi.org/10.1063/1.5007180).
The path optimization was introduced by Einarsdóttir *et al.*,
[*Lecture Notes in Computer Science* **7134**, 45](https://doi.org/10.1007/978-3-642-28145-7_5).

The job starts from a converged minimum-energy path in `neb.con`. At energy
{math}`E=E^\ddagger-\Delta E`, the two MEP segments that bracket {math}`E`
give the turning points and the initial band. Each subsequent, lower energy
uses the preceding optimal tunnelling path and the newly exposed MEP
segments as its seed. Continuation ends when the orbit temperature reaches
the requested `lowest_temperature`. `energy_step` caps each energy decrease;
the job can take a smaller decrease so that the predicted change in
{math}`\ln k` does not exceed `log_rate_tolerance`.

## LI-NEB path

The optimized quantity is the abbreviated JWKB action

```{math}
W(E)=\int_{\rm OTP}\sqrt{2[V(\mathbf q)-E]}\,d s ,
```

in mass-weighted coordinates {math}`\mathbf q_i=\sqrt{m_i}\mathbf r_i`.
The masses therefore remain the real atomic masses. At an interior image the
discrete action gradient is projected perpendicular to the NEB tangent. A
spring acts on mass-weighted path lengths, only along the tangent, and keeps
the images evenly spaced. At each endpoint the contour force likewise acts
on mass-weighted path length: the spring force is projected onto the local
{math}`V=E` contour and the normal component pulls the endpoint onto that
contour.

The default optimizer is velocity-projection velocity Verlet. `fire`,
`lbfgs`, `cg`, `qm`, `sd`, and `xtsci` select the corresponding eOn
optimizer. Convergence requires both the largest interior perpendicular
action gradient and the largest endpoint optimization force to be below
`force_tolerance`.

For a finite system with no fixed atoms, every path evaluation removes
mass-weighted translation and applies mass-weighted Kabsch RMSD alignment to
remove rotation. Periodic structures retain their physical rotational
coordinates.

## Period, replicas, and rate

Each converged OTP is interpolated by cubic polynomials. The imaginary
half-period is evaluated directly by adaptive quadrature,

```{math}
\frac{\tau}{2}=\int_{\rm OTP}
\frac{|d\mathbf q/ds|}{\sqrt{2[V(\mathbf q(s))-E]}}\,ds .
```

A squared-sine change of variable removes the integrable square-root
singularities at both turning points. The quadrature error is estimated by
Gauss--Legendre refinement and replicas are placed by inverting its cumulative
integral. The reported temperature is

```{math}
T=\frac{\hbar}{k_B\tau}.
```

The job places {math}`P` replicas at {math}`t_j=j\tau/P` and evaluates Eq. 6
of the JCP paper. When the potential provides an exact mass-weighted
Hessian-vector product, the image and ring operators use it directly. For a
potential without that capability, eOn assembles each image Hessian with the
shared finite-difference adapter, using `[Main] finite_difference`, reports
that fallback as a warning, and reuses the resulting matrix for ring products.
The rate always uses the scalable cyclic block-ring determinant, with no dense
diagonalization branch or size-dependent method switch: all
{math}`P N-1` nonzero frequencies enter, where {math}`N` is the number of
active mass-weighted coordinates, and the imaginary-time zero mode is removed.
The closed orbit traverses the OTP twice, so the exponent uses
{math}`S(E)=2W(E)+E\tau`. Reactant rigid modes and the corresponding ring
rigid modes are removed together. For a finite molecule with no fixed atoms,
the reported value is the fixed-frame vibrational rate used by eOn's
instanton implementation. Eq. 6's separate rotational partition-function
quotient is not implemented. Fixed or periodic systems do not need that
quotient.

## Configuration

```ini
[Main]
job = line_integral_tunneling

[Line Integral Tunneling]
reactant_filename = reactant.con
saddle_filename = saddle.con
mep_filename = neb.con

images = 10
replicas = 128
energy_offset = 0.01
energy_step = 0.01
lowest_temperature = 50
log_rate_tolerance = 0.5

spring = 20
contour = 1000
force_tolerance = 0.01
max_iterations = 5000
optimizer = velocity_verlet
li_neb_time_step = 0.01
period_tolerance = 0.01
determinant_estimated_error = 0.1
```

`energy_offset` and `energy_step` are in eV. `spring` is in
Å{sup}`-2`, and `contour` is in eV{sup}`-1` Å{sup}`-1`.
`li_neb_time_step` is the velocity-projection optimizer step, while
`period_tolerance` is the absolute tolerance on the imaginary-time period
quadrature. Both use the internal mass-weighted time unit
{math}`\sqrt{\mathrm{amu}}\,\text{Å}/\sqrt{\mathrm{eV}}`.

The job rejects an MEP that does not cross the requested energy on both sides
of its highest image, or that contains more than one disjoint above-energy
barrier. Each input MEP must describe one elementary barrier. The job also
rejects partially constrained atoms, an interior image below {math}`E`, an
unconverged OTP, an orbit whose period uncertainty overlaps the crossover
period, or a ring Hessian that is not a first-order saddle. The crossover
period is obtained from the saddle Hessian's lowest mass-weighted mode. Whole
atoms may be fixed. It does not change methods or adjust user parameters in
response.

## Output

`line_integral_tunneling.dat` contains one row per energy:

- {math}`E`, {math}`W(E)`, {math}`T`, {math}`\tau`, and its estimated
  quadrature-plus-discretization error;
- {math}`\ln k` in the internal time unit and
  {math}`\ln(k/\mathrm{s}^{-1})`, plus {math}`k` in s{sup}`-1`;
- force calls for the first LI-NEB and lower-energy continuation, plus
  Hessian-vector products (or matrix products after the reported
  finite-difference fallback);
- the two convergence measures and ring-Hessian negative-mode count.

`line_integral_path_NNN.con` contains the corresponding OTP. Its first frame
repeats the action, period, period error, temperature, rate, and split
force-call counts as metadata. `results.dat` summarizes the final energy,
crossover temperature, period error, cumulative stage counts, and whether the
finite-difference Hessian fallback was used.

For a direct check at one of the reported temperatures, run
{doc}`instanton` in rate mode with the same saddle and `temperature`. PI-QTST
from that job provides the independent constrained-ring calculation when
`pi_planes` is nonzero.

## Regression checks

The two-dimensional regression uses the canonical four-term surface of
[Müller and Brown](https://doi.org/10.1007/BF00547608) and the JCP paper's
stated modification: three Gaussians shift and deepen the right minimum. The
article and preliminary paper do not tabulate those three terms. The test
therefore labels its surface as a reconstruction rather than claiming
unpublished coefficients. It places a positive Gaussian at the canonical
right basin and two negative Gaussians at the shifted basin read from Fig. 2:

| amplitude | {math}`a` | {math}`b` | {math}`c` | {math}`x_0` | {math}`y_0` |
|---:|---:|---:|---:|---:|---:|
| 100 | -10 | 0 | -10 | 0.623 | 0.028 |
| -100 | -10 | 0 | -10 | 1.0 | 0.0 |
| -50 | -5 | 0 | -5 | 1.2 | 0.0 |

Here each term is
{math}`A\exp[a(x-x_0)^2+b(x-x_0)(y-y_0)+c(y-y_0)^2]`; 600 surface-energy
units equal one eV in the test. A steepest-descent MEP generated on this
fixed surface is committed with the test. Below the intermediate minimum,
LI-NEB reproduces the corner cutting shown in Fig. 2. Digitizing the red
energy profile in Fig. 3 at 0.1-displacement intervals gives
{math}`W=0.81\pm0.04` for the figure's unit-mass surface; the regression gives
{math}`W=0.7911`. The uncertainty combines the 0.01 eV ordinate reading and
0.05 displacement reading from the printed tick spacing.

When the built-in CuH2 embedded-atom potential is enabled, a second regression
uses 216 fixed Cu atoms in a six-layer Cu(110) slab and two mobile H atoms. A
separately converged eOn NEB supplies its MEP. The test uses 10 LI images,
the paper defaults, and a 150-replica independently optimized eOn rate
instanton. That calculation is a useful implementation comparison, but it is
not the paper's reference calculation:

```{include} ../fig/generated/line_integral_cuh2_rates.md
```

The two lowest-temperature points include explicit image-count, replica-count,
and Hessian-interpolation checks. The {math}`\ln k` values are in the internal
time unit:

```{include} ../fig/generated/line_integral_cuh2_convergence.md
```

The tables, and every number quoted about them in this section, are generated
by {file}`scripts/line_integral_cuh2_tables.py` from the CSV files in
{file}`docs/source/fig/data/`. To refresh them, run
`meson test -C build test_cuh2_line_integral --print-errorlogs`; it writes
`cuh2_line_integral_*.csv` into `eon_cuh2_line_integral` under the system
temporary directory, and those files replace the committed ones.
The regression context sizes its bead-product workers from the available
hardware, with one cloned CuH2 potential per worker. The production job and
the documented cost measurements use one.

The 2018 article does not contain a Cu rate table or rate-against-temperature
figure from which absolute published rates can be recovered. Section II E
only states that the natural logarithm of its LI flux was within 2% of its
150-replica reference from 260 K to 50 K, except at 50 K. Consequently, that
unit-dependent statement is not used as a tolerance on the rate ratio above.
The rate plot is instead Fig. 4(b) of the preliminary 2012 paper. It compares
HQTST with RAW-QTST for a different model containing 432 active coordinates,
and its caption says that the LI path uses 10 images at high energy and 40 at
low energy. Values digitized from that HQTST curve are:

| {math}`T` (K) | {math}`\log_{10}(k/\mathrm{s}^{-1})` |
|---:|---:|
| 255.3 | 9.3 |
| 234.9 | 7.4 |
| 214.3 | 5.3 |
| 193.1 | 3.1 |
| 174.3 | 0.9 |
| 156.7 | -1.4 |
| 138.9 | -3.8 |
| 120.5 | -6.4 |
| 101.0 | -9.2 |
| 78.6 | -11.9 |

The values were read at the original eOn comparison temperatures;
their uncertainty is about ±0.5 in {math}`\log_{10}k`, rising to ±0.6 on the
steep high-temperature segment. Because the fixed-Cu fixture has only six
active coordinates, no paper-derived acceptance tolerance is applied to this
non-equivalent comparison.

```{include} ../fig/generated/line_integral_cuh2_calls.md
```
