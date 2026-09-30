---
myst:
  html_meta:
    "description": "Tunnelling splittings between two minima in eOn: WKB along a NEB band and the ring-polymer instanton job."
    "keywords": "eOn instanton, tunnelling splitting, two-level system, WKB, ring polymer"
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

The propagator ratio measures {math}`\Delta_0` when the two wells lie within
a small fraction of {math}`k_B T` of each other. `instanton_symmetric = 0`
flags a pair outside that window. The action and the path still hold there,
but the splitting prefactor does not.

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
