---
myst:
  html_meta:
    "description": "Tunnelling splittings between two minima in eOn, and the thermal rate through a saddle below the crossover: WKB along a NEB band and the ring-polymer instanton job."
    "keywords": "eOn instanton, tunnelling splitting, two-level system, WKB, ring polymer, instanton rate"
---

# Tunnelling splittings

A two-level system (TLS) in a glass is a pair of adjacent minima that the
structure tunnels between at about one kelvin. Its tunnelling splitting
{math}`\Delta_0` and asymmetry {math}`\Delta` set the TLS energy
{math}`E = \sqrt{\Delta^2 + \Delta_0^2}`. eOn estimates {math}`\Delta_0` in two
ways:

| | NEB band (WKB) | `job = instanton` |
|---|---|---|
| Path | the minimum energy path | the path of least imaginary-time action |
| Dimensions | one, along the mass-weighted band | every free degree of freedom |
| Cost | the NEB itself | one batch of forces per iteration over the beads, plus a Hessian per bead |
| Output | first frame of `neb.con` | `instanton.con` and `results.dat` |

With `mode = rate` the same job estimates the thermal rate through a saddle
instead of a splitting. That path is below.

Energies are in eV, lengths in Å and masses in amu throughout, so
mass-weighted lengths are in amu^0.5 Å.

## WKB along a NEB band

Every NEB job writes `reaction_coordinate_mw`, the mass-weighted arc length, on
each frame of `neb.con`. The first frame also carries the band's WKB estimate:

| Key | Meaning |
|---|---|
| `hbar_omega_reactant`, `hbar_omega_product` | {math}`\hbar\omega` of each well along the band, from a fit of {math}`a s^2 + b s^3` to the images within half the barrier |
| `tunnel_action` | {math}`S = \hbar^{-1} \int \sqrt{2 (V(s) - E)}\, ds` over the forbidden region |
| `tunnel_splitting` | {math}`\Delta_0 = (\hbar\omega / \pi) e^{-S}`, with {math}`\omega` the geometric mean of the wells |
| `tls_energy` | {math}`\sqrt{\Delta^2 + \Delta_0^2}` |
| `tunnel_deep_wells` | 1 when both barriers exceed {math}`\hbar\omega`; below that, WKB is the wrong tool |

The profile between images is a monotone cubic, so it cannot dip below the
data. The level {math}`E` is the higher of the two harmonic ground states. A
structure without masses, or a band whose end is flat, leaves these keys out.

WKB along the band is exact in one dimension up to its semiclassical error.
When the path curves, the tunnelling cuts the corner, and the transverse
zero-point energy changes along the way. In the two-dimensional test valley
below, both effects together put the band estimate a factor of 2.8 below the
exact splitting.

## The instanton job

The instanton is a path of {math}`P + 1` beads in imaginary time
{math}`\beta\hbar`, with its ends fixed at the two minima. It minimises the
discretised Euclidean action

```{math}
S = \sum_j \frac{|q_{j+1} - q_j|^2}{2\,\delta\tau} + \delta\tau \sum_j V(q_j),
\qquad \delta\tau = \beta\hbar / P,
```

in mass-weighted coordinates {math}`q`. The splitting comes from the ratio of
the off-diagonal to the diagonal imaginary-time propagator, both taken in the
same steepest-descent approximation:

```{math}
\Delta_0 = 2\hbar \sqrt{\frac{S_0}{2\pi\hbar\,\delta\tau}}
\sqrt{\frac{\det J_\text{well}}{\det' J}}\; e^{-(S - S_\text{well})/\hbar}.
```

Here {math}`J` is the Hessian of {math}`S` over the interior beads. The
prime leaves out its zero mode, the kink's position in imaginary time.
{math}`S_0 = \int |\dot q|^2 d\tau`, and {math}`J_\text{well}` is the same
Hessian with every bead at a minimum.

```{code-block} ini
[Main]
job = instanton

[Potential]
potential = rgpot

[Instanton]
reactant_filename = reactant.con
product_filename = product.con
; start from a converged band instead of the straight line
initial_path = neb.con
beads = 256
beta_hbar_omega = 30
force_tolerance = 1e-3
; one finite-difference Hessian every 4 beads, linear in between
hessian_stride = 4
```

`beta_hbar_omega` sets the imaginary time in units of {math}`1/\omega` of
the stiffer minimum along the line between them. It must be long enough for
the kink to relax into both wells. `results.dat` reports
`instanton_mode_separation`, the second eigenvalue of {math}`J` over its zero
mode. Values above {math}`10^3` mean the kink is isolated; small values mean
{math}`\beta\hbar` is too short. With no atom fixed, the product is aligned
to the reactant first: its mass-weighted mean displacement is removed, and
for a cluster its best rotation as well.

Each iteration evaluates every interior bead in one call. Under
`[RgpotPot] ranks_per_image`, that call spreads the beads over the CPMD
calculator groups the same way a NEB spreads its images.

`instanton.con` holds one frame per bead, with `imaginary_time_fs`. Its first
frame and `results.dat` carry:

| Key | Meaning |
|---|---|
| `tunnel_splitting_instanton` | {math}`\Delta_0`, eV |
| `instanton_action` | {math}`(S - S_\text{well})/\hbar` |
| `tls_energy_instanton` | {math}`\sqrt{\Delta^2 + \Delta_0^2}`, eV |
| `tunnel_asymmetry` | {math}`V(\text{product}) - V(\text{reactant})`, eV |
| `instanton_temperature_K` | {math}`1/(k_B \beta)` for the imaginary time used |
| `instanton_mode_separation` | how well the kink's translation separates from the other modes |
| `instanton_symmetric` | 1 when {math}`\beta|\Delta| < 0.1` |
| `instanton_beta_asymmetry` | {math}`\beta|\Delta|` |

The propagator ratio measures {math}`\Delta_0` when the two wells lie within
a small fraction of {math}`k_B T` of each other. `instanton_symmetric = 0`
flags a pair outside that window. The job still writes the path and the
action, but no `tunnel_splitting_instanton`, and reports success: the flag
says why. For such a pair, set `mode = rate`, give the saddle and a
temperature below the crossover, and read the rate in the next sections.

## Which path object

NEB images and ring-polymer beads are different objects. An image is a point
on a path in configuration space between two minima. Its springs are
fictitious, and only the force perpendicular to the path is kept. A bead is
one imaginary-time slice of a single quantum system. Its springs are
physical, with stiffness fixed by the temperature and the number of beads.

| Object | Points | Springs | What it returns |
|---|---|---|---|
| `mode = splitting` | open string between two minima, optionally started from a band | Euclidean action, no tangent projection | tunnelling splitting when the wells are close in energy |
| `mode = rate` | closed ring through one saddle | same action, stiffness set by {math}`T` and {math}`N` | thermal rate below the crossover temperature |
| Centroid PMF | one ring per image, centroid held on the image | sampled, not minimised | quantum free-energy barrier along the path |
| Harmonic centroid string | one ring per image, optimised | local harmonic quantum correction | a free-energy estimate only as good as that harmonic well |

The first two rows are `job = instanton`. The centroid potential of mean
force is a constrained path-integral molecular dynamics sample, one
thermostatted ring per image. That sample is not an optimisation, and it
does not belong in this job. A string of harmonically corrected rings is a
different calculation again, and this page does not implement it.

## The rate below the crossover

`mode = rate` reads the reactant and `saddle_filename` (default `saddle.con`).
The instanton is a closed ring, a first-order saddle of the ring-polymer
potential, with one negative eigenvalue and one zero eigenvalue that cycles
the beads. Richardson and Althorpe, J. Chem. Phys. 131, 214106 (2009), give

```{math}
k Z_r = \frac{1}{\beta_N \hbar}
\sqrt{\frac{B_N}{2\pi \beta_N \hbar^2}}
\prod_k' \frac{1}{\beta_N \hbar |\omega_k|}
\exp(-\beta_N U_N).
```

{math}`B_N` is the sum of squared steps around the ring, {math}`\omega_k^2`
are the eigenvalues of the mass-weighted ring Hessian, and the prime leaves
out the cyclic zero mode and the rigid translations and rotations.
{math}`Z_r` is the harmonic ring-polymer partition function of the reactant.
{math}`U_N` is the ring-polymer potential, the bead potentials plus the springs.

The crossover temperature is {math}`T_c = \hbar \omega_b / (2\pi k_B)`, with
{math}`\omega_b` the imaginary frequency at the saddle. At or above {math}`T_c`
the ring collapses onto the saddle. The rate there is a parabolic barrier
correction, and classical transition-state theory is only the one-bead limit
of that correction. The job does not evaluate it: it stops and says so.

```{code-block} ini
[Main]
job = instanton

[Instanton]
mode = rate
reactant_filename = reactant.con
saddle_filename = saddle.con
temperature = 5
beads = 256
hessian_stride = 1
```

`temperature` is in kelvin. It must be positive and below {math}`T_c`. The
default 0 means the temperature was not set, and `mode = rate` then refuses
to run. `beads` defaults to 256, the same default as the splitting. A ring
of 32 beads at 5 K does not resolve a stiff bond: the path integral starts
to converge once the bead count exceeds {math}`\beta \hbar \omega` of the
stiffest mode. `hessian_stride` of 1 takes a Hessian on every bead. A larger
stride keeps a Hessian on every stride-th bead and interpolates linearly
between those anchors. That interpolation is an approximation: the rate
formula uses the Hessian of every bead.

With no atom fixed, the three translations are omitted on both sides. A
rotation is omitted when it is a zero mode of the reactant Hessian, which a
free cluster has and a crystal does not. A cluster in a large periodic cell
is told apart by that Hessian, not by the periodic flag. The springs along
those directions stay, so they cancel between the instanton and the reactant.

`results.dat` then carries:

| Key | Meaning |
|---|---|
| `rate_instanton` | {math}`k` in s^{-1} |
| `rate_instanton_log` | {math}`\ln(k\,/\,\mathrm{s}^{-1})` |
| `rate_htst_log` | classical harmonic transition-state theory, the same logarithm |
| `instanton_crossover_K` | {math}`T_c`, K |
| `instanton_negative_modes` | negative eigenvalues of the ring Hessian; a first-order saddle has 1 |
| `instanton_zero_mode` | the eigenvalue left out |

Habershon, Manolopoulos, Markland and Miller, Annu. Rev. Phys. Chem. 64, 387
(2013), doi:10.1146/annurev-physchem-040412-110122, expect the sampled
ring-polymer rate to lie within about a factor of two of the exact quantum
rate between {math}`T_c` and {math}`T_c/2`. That bound compares the sampled
rate with the exact rate. It is not a comparison of this instanton with a
free-energy profile.

The cubic metastable well, {math}`V = \omega_0^2 q^2/2 - g q^3/3`, is the
check. Deep below the crossover its rate approaches the zero-temperature
decay of Caldeira and Leggett, Ann. Phys. 149, 374 (1983).

## Checks

The Catch2 cases `Instanton splitting in a curved valley matches the exact
gap` and `The instanton cuts the corner the minimum energy path takes` use
{math}`V = V_0 (x^2 - 1)^2 + \tfrac{K}{2} (y - C (1 - x^2))^2` at unit mass
with {math}`K = 4` eV/Å². The exact gap comes from a fourth-order
finite-difference Hamiltonian, converged to {math}`10^{-5}`.

| {math}`V_0` / eV | {math}`C` / Å | exact {math}`\Delta_0` / eV | instanton / exact | WKB on the valley floor / exact |
|---|---|---|---|---|
| 0.12 | 0.35 | 2.520e-5 | 1.17 | 0.35 |
| 0.30 | 0.35 | 9.818e-8 | 1.08 | |
| 0.12 | 0 | 2.040e-5 | 1.12 | 1.02 |
| 0.30 | 0 | 1.195e-7 | 1.07 | |

The instanton's error falls as the barrier deepens, the regime glass TLS sit
in. The same cases tie the C++ path and splitting to an independent
implementation of the discretisation to {math}`2 \times 10^{-3}`.
