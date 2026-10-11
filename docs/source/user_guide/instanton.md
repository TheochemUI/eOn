---
myst:
  html_meta:
    "description": "Tunnelling splittings between two minima in eOn, and the thermal rate through a saddle: WKB along a NEB band, the ring-polymer instanton below the crossover, and the parabolic barrier factor above it."
    "keywords": "eOn instanton, tunnelling splitting, two-level system, WKB, ring polymer, instanton rate."
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

Every nudged elastic band (NEB) job writes `reaction_coordinate_mw`, the mass-weighted arc length, on
each frame of `neb.con`. The first frame also carries the band's Wentzel-Kramers-Brillouin (WKB) estimate.
Every band variant eOn runs carries them, the off-path climbing image band
{cite:p}`inst-goswamiEnhancedClimbingImage2026` among them. How eOn finds the
minima and saddles these jobs start from is set out in
{cite:t}`inst-goswamiEfficientExplorationChemical2025`.
The keys are:

| Key | Meaning |
|---|---|
| `hbar_omega_reactant`, `hbar_omega_product` | Wells: {math}`\hbar\omega` of each well along the band, from a fit of {math}`a s^2 + b s^3` to the images within half the barrier |
| `tunnel_action` | Action: {math}`S = \hbar^{-1} \int \sqrt{2 (V(s) - E)}\, ds` over the forbidden region |
| `tunnel_splitting` | Estimate: {math}`\Delta_0 = (\hbar\omega / \pi) e^{-S}`, with {math}`\omega` the geometric mean of the wells |
| `tunnel_asymmetry_zpe` | Asymmetry: {math}`\Delta = V_p - V_r + (\hbar\omega_p - \hbar\omega_r)/2`, the difference of the two local ground states |
| `tls_energy` | Energy: {math}`\sqrt{\Delta^2 + \Delta_0^2}` |
| `tunnel_deep_wells` | Flag: 1 when both barriers exceed {math}`\hbar\omega`; below that, WKB is the wrong tool |

The profile between images is a monotone cubic, so it cannot dip below the
data. The level {math}`E` is the higher of the two harmonic ground states. A
structure without masses, or a band whose end is flat, leaves these keys out.

{math}`\Delta` is the diagonal term of the two-level Hamiltonian
{cite:p}`inst-andersonAnomalousLowtemperatureThermal1972`, the eq. S9 of
{cite:t}`inst-khomenkoDepletionTwoLevelSystems2020`. A tilt changes the
curvature of each well as well as its depth: on a quartic double well tilted
by 1 meV, the difference of the minima alone puts the energy 6 percent above
the exact gap, and {math}`\Delta` within 0.6 percent.

## The band solved in one dimension

The same frame carries the exact levels of the one-dimensional problem
along the band, {math}`-\tfrac{\hbar^2}{2} \partial_s^2 + V(s)` in the
mass-weighted coordinate, the protocol of the glass two-level-system
studies {cite:p}`inst-khomenkoDepletionTwoLevelSystems2020,inst-mocanuMicroscopicObservationTwolevel2023`.
The keys are:

| Key | Meaning |
|---|---|
| `tls_energy_dvr` | Energy: {math}`E_2 - E_1` of the two lowest levels, eV |
| `tunnel_asymmetry_dvr` | Asymmetry: {math}`E \cos 2\theta`, positive when the product well lies higher, eV |
| `tunnel_splitting_dvr` | Matrix element: {math}`\Delta_0 = E \sin 2\theta`, eV |
| `tunnel_wall_sampled` | Flag: 1 when the outer walls were sampled, below |
| `tunnel_double_well` | Flag: 1 when the images rise to one maximum and fall, a double well and no more |

The kinetic energy is the fourth-order central difference on a uniform grid
of 20 points per oscillator length {math}`\hbar / \sqrt{\hbar\omega}` of the
stiffer well, at most 1500 points, reaching seven such lengths past each
minimum. On a quartic double well, symmetric or tilted, the gap lands
within {math}`10^{-5}` of a converged grid. {math}`\theta` rotates the two
lowest states into the pair most localised on either side of the barrier
top. That rotation cancels the tail each state leaves across the barrier,
so on a tilted quartic {math}`\Delta_0` stays within {math}`10^{-3}` of the
gap of the level wells with the asymmetry thousands of times larger, where
the ground state's weight on each side alone would overshoot fourfold.

A band runs from minimum to minimum and never sees the outer walls. Past
each end the profile continues by the fit of that well, its cubic or its
parabola, whichever is stiffer; on a quartic double well that keeps the gap
within 6 percent. `[Nudged Elastic Band] tunnel_wall_points` above 0 samples
the real walls instead, the linear extrapolation of the reaction path of
{cite:t}`inst-khomenkoDepletionTwoLevelSystems2020`: that many force calls
past each minimum along the band's end segment continued outward, out to
seven oscillator lengths, once the band has converged. The fit stays as the
model and the residual of the samples, a cubic spline flat at the minimum,
is added to it. Eight samples bring the quartic within 1.5 percent of the
exact gap, the rest being the band's own images.

```{code-block} ini
[Nudged Elastic Band]
images = 13
tunnel_wall_points = 8
```

Either solution stays one-dimensional; WKB along the band is exact in one
dimension up to its semiclassical error. When the path curves, the
tunnelling cuts the corner, and the transverse zero-point energy changes
along the way. In the two-dimensional test valley below, both effects
together put the WKB band estimate a factor of 2.8 below the exact
splitting. For both, run the instanton job on the pair.

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
same steepest-descent approximation.

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
springs = trotter
beta_hbar_omega = 30
force_tolerance = 1e-3
; one finite-difference Hessian every 4 beads, linear in between
hessian_stride = 4
```

`springs = eco` is refused. The instanton is a Trotter discretisation of
the ring.

`beta_hbar_omega` sets the imaginary time in units of {math}`1/\omega` of
the stiffer minimum along the line between them. It must be long enough for
the kink to relax into both wells. `results.dat` reports
`instanton_mode_separation`, the second eigenvalue of {math}`J` over its zero
mode. Values above {math}`10^3` mean the kink is isolated; small values mean
{math}`\beta\hbar` is too short. With no atom fixed, the product is aligned
to the reactant first: its mass-weighted mean displacement is removed, and
for a cluster its best rotation as well.

With that reactant rotation removed, each iteration evaluates every interior bead of the kink in one call. Under
`[RgpotPot] ranks_per_image`, that call spreads the beads over the Car-Parrinello molecular dynamics (CPMD)
calculator groups the same way a NEB spreads its images.

`instanton.con` writes one frame per bead, with `imaginary_time_fs`. The keys are:

| Key | Meaning |
|---|---|
| `tunnel_splitting_instanton` | Splitting: {math}`\Delta_0`, eV |
| `instanton_action` | Action: {math}`(S - S_\text{well})/\hbar` |
| `tls_energy_instanton` | Energy: {math}`\sqrt{\Delta^2 + \Delta_0^2}` with {math}`\Delta` the zero-point asymmetry below, eV |
| `tunnel_asymmetry` | Asymmetry: {math}`V(\text{product}) - V(\text{reactant})`, eV |
| `tunnel_asymmetry_zpe` | Asymmetry of the local ground states: {math}`V_p - V_r + \tfrac{\hbar}{2}\left(\sum_k \omega_k^{p} - \sum_k \omega_k^{r}\right)` over every vibration of the two minimum Hessians, eV |
| `instanton_temperature_K` | Temperature: {math}`1/(k_B \beta)` for the imaginary time used |
| `instanton_mode_separation` | Separation: how well the kink's translation separates from the other modes |
| `instanton_symmetric` | Symmetry: 1 when {math}`\beta|V_p - V_r| < 0.1` |
| `instanton_beta_asymmetry` | Magnitude: {math}`\beta|V_p - V_r|` |
| `instanton_symmetrized` | Levelled: 1 when the instanton ran on the levelled surface below |

Two minima of one energy can still hold different zero-point energies,
when an isotope or the cage around a defect makes one well stiffer than the
other. Then {math}`\Delta_0` is the tunnelling matrix element taken with the
geometric mean of the two well kernels, and the levels split by
{math}`\sqrt{\Delta^2 + \Delta_0^2}` with the harmonic zero-point
difference in {math}`\Delta`
{cite:p}`inst-jahrInstantonTheoryTunneling2020`. On a valley whose
transverse stiffness differs by 10 percent between the wells, that
difference splits the levels 140 times further than tunnelling does, and
{math}`\sqrt{\Delta^2 + \Delta_0^2}` lands within 6 percent of the exact
gap.

Most two-level systems in a glass also have minima of different energy,
with the asymmetry comparable to the splitting. On such a surface no kink
joins the two minima at rest, and the propagator ratio measures
{math}`\Delta_0` only while the wells lie within a small fraction of
{math}`k_B T` of each other. The job therefore levels the wells first. It
subtracts {math}`(V_p - V_r)\,\sigma(\xi)`, with {math}`\xi` the projection
of {math}`q` on the line from the reactant to the product (0 at one, 1 at
the other) and {math}`\sigma = 6\xi^5 - 15\xi^4 + 10\xi^3`. Both minima
stay stationary with their Hessians, since {math}`\sigma`, {math}`\sigma'`
and {math}`\sigma''` vanish at both ends, and the instanton of the levelled
surface gives the tunnelling matrix element. The difference of the minima
goes back into {math}`\Delta`, so moving the bias out of the path and into
the asymmetry, as
{cite:t}`inst-erakovicVibrationalTunnelingSpectra2022` do for the coupling
of nondegenerate well states. That holds while {math}`|V_p - V_r|` is small
against the barrier, the regime of a two-level system. On a quartic double
well tilted until {math}`\beta |V_p - V_r| = 0.47`, the matrix element stays
within 1 percent of the level wells' and {math}`\sqrt{\Delta^2 + \Delta_0^2}`
within 2 percent of the exact gap. `instanton_action` is then the action on
the levelled surface, and the frames of `instanton.con` keep the energies of
the real one.

`symmetrize = false` keeps the bare surface. `instanton_symmetric = 0` then
flags a pair outside the window. The job still writes the path and the
action, but no `tunnel_splitting_instanton`, and reports success: the flag
says why. For such a pair, set `mode = rate`, give the saddle and a
temperature below the crossover, and read the rate in the next sections.

## Which path object

NEB images and ring-polymer beads are different objects. An image is a point
on a path in configuration space between two minima. Its springs are
fictitious. The parallel force is removed. A bead is
one imaginary-time slice of a single quantum system. Its springs are
physical, with stiffness fixed by the temperature and the number of beads.

The columns are:

| Object | Points | Springs | What it returns |
|---|---|---|---|
| `mode = splitting` | open string between two minima, started from a band when one is present | Euclidean action, no tangent projection | tunnelling splitting when the wells are close in energy |
| `mode = rate` | closed ring through one saddle | same action, stiffness set by {math}`T` and {math}`N` | thermal rate: the ring below the crossover, the parabolic factor above it |
| Centroid potential of mean force (PMF) | one ring per image, centroid held on the image | sampled, not minimised | quantum free-energy barrier along the path |
| Harmonic centroid string | a ring at each image, optimised | local harmonic quantum correction | a free-energy estimate as good as that harmonic well |

The first two rows are `job = instanton`. The centroid potential of mean
force is a constrained path-integral molecular dynamics sample, one
thermostatted ring per image. That sample is not an optimisation, and it
does not belong in this job. A string of harmonically corrected rings is a
different calculation again, and this page does not implement it.

## The rate below the crossover

`mode = rate` reads the reactant and `saddle_filename` (default `saddle.con`).
Any eOn saddle search can supply it, the Gaussian process accelerated dimer
{cite:p}`inst-goswamiEfficientImplementationGaussian2025,inst-goswamiAdaptivePruningIncreased2025b`
among them.
The instanton is a closed ring, a first-order saddle of the ring-polymer
potential, with one negative eigenvalue and one zero eigenvalue that cycles
the beads. {cite:t}`inst-richardsonRingpolymerMolecularDynamics2009` give

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
the ring collapses onto the saddle. Above {math}`T_c` the rate is the
parabolic barrier factor

```{math}
\kappa = \frac{\pi T_c / T}{\sin(\pi T_c / T)}
```

times the quantum harmonic transition-state theory rate from the same
Hessians. The job writes the parabolic barrier factor. It does not search
a ring at that temperature.
{math}`\kappa` tends to 1 at high temperature, which is the one-bead limit.
The factor diverges as {math}`T` approaches {math}`T_c` from above.

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

`temperature` is in kelvin. It must be positive. Below {math}`T_c` the job
optimises the ring. Above {math}`T_c` the job writes the parabolic barrier
factor and does not optimise a ring. At {math}`T_c` the factor diverges
and that temperature records no rate. The optimizer is not what keeps the path unused. Each
bead still needs a force. The rate uses the fluctuation prefactor from the
bead Hessians, and above the crossover the ring collapses onto the saddle.
The default 0 means the temperature was not set, and
`mode = rate` then refuses to run. `beads` defaults to 256, the same default as the splitting. A ring
of 32 beads at 5 K does not resolve a stiff bond: the path integral starts
to converge once the bead count exceeds {math}`\beta \hbar \omega` of the
stiffest mode. `hessian_stride` of 1 takes a Hessian on every bead. A larger
stride keeps a Hessian on every stride-th bead and interpolates linearly
between those anchors. That interpolation is an approximation: the rate
formula uses the Hessian of every bead.

`half_ring` defaults to enabled. On an even bead count the potential is
evaluated from one turning point to the other and copied onto the mirror.
The reaction coordinate on that half stays monotone, so the search cannot
settle on an out-and-back bounce.
A half ring that has converged, or that has stopped with the gradient
under the force tolerance while the negative-mode count is not 1, is
checked for an unstable mode that is odd under the mirror when the
interior beads still match their mirrors. Two copies of the instanton
on one ring are that mode. The search steps along it and continues on
the whole ring.
An odd count keeps every bead. `energy_shift` (default 0, in eV) is
subtracted from every bead potential and from the reactant and saddle
energies in the rate. An empty `discretization` keeps every imaginary-time
step at {math}`\beta_N \hbar`. A comma-separated list of one positive weight
per bead makes the ring an adaptive grid
{cite:p}`inst-rommelAdaptiveIntegrationGrids2011`: the weights are scaled to
a mean of one, the link from bead {math}`j` to bead {math}`j+1` lasts
{math}`w_j \beta_N \hbar`, so its spring is {math}`c / w_j`, and bead
{math}`j` carries {math}`(w_{j-1} + w_j)/2` of its potential, the
trapezoidal rule over unequal steps. The rate keeps the same grid: the
ring Hessian, the speed {math}`B_N = \sum_j |q_{j+1} - q_j|^2 / w_j^2` of
the time shift, and the harmonic reactant ring of {math}`Z_r`. Below about
a fifth of {math}`T_c` the beads of a uniform ring crowd the turning
points; longer steps there and shorter ones where the ring crosses the
barrier spread them along the path, and on the Eckart barrier at
{math}`0.2\,T_c` steps {math}`1 + 0.3 \cos(4\pi (j + 1/2)/N)` cut the
error of 48 beads fourfold. Weights keep the whole ring (no `half_ring`),
need the Newton search, and take no friction bath, whose frequencies
belong to the uniform ring.

`friction` adds the electronic-friction bath of
{cite:t}`inst-litmanDissipativeTunnelingRates2022` to the ring,
{math}`\sum_{l \ne 0} \tfrac{\omega_l}{2} |G_l|^2` over the ring's normal
modes with {math}`\omega_l = 2\omega_P |\sin(\pi l / N)|`,
{math}`\omega_P = 1/(\beta_N \hbar)`, and {math}`G` the normal modes of the
line integral of {math}`\sqrt{\eta}` along the ring (their Eqs. 20 and 35).
{math}`\eta` is a friction per unit mass in inverse time units of
{math}`\sqrt{\mathrm{amu}\,\mathrm{Å}^2/\mathrm{eV}}` (10.18 fs): one value
on every bead with `friction = implicit` and `friction_eta`, or one per
bead with `friction = explicit` and `friction_eta_beads`, whose links are
closed around the ring so that no bead is its start. The bath couples every
bead to every other, so the ring is searched with the dimer on its
gradients, and the rate takes the bath's curvature from the dense ring
Hessian and adds {math}`\eta\,\omega_k` to every {math}`k > 0` mode of
{math}`Z_r` (their Eq. 36). A bead-wise bath gives the reactant no friction
of its own, so with `friction = explicit` the ring is written but no rate.

The search is an index-1 Newton step on the ring Hessian while the active
coordinate count is within the Newton limit. Past that limit the ring is one
structure and the dimer follows its unstable mode, with the beads in one force
batch. Without an
`initial_path`, the rate job traces a steepest-descent path out of the
saddle along both signs of its unstable mode and seeds the ring from it by
the period condition below. That path is written to `instanton_sd_path.con`,
one frame per point, with the energy and the arc length. A later rate job
reads it as `initial_path` and does not repeat those force calls. On LJ13
that path costs 586 gradient calls. Cooling a cosine
ring from 0.85 of the crossover is the fallback when no path can be
built; it finds the ring only where the ring grows continuously out of
the saddle as the temperature drops. Where it does not, the search walks
to a neighbouring saddle, so a converged ring that fails the channel test
described with the `results.dat` keys below is refused and no rate is
written. The step climbs the mode that overlaps
the last climb and turns every other negative curvature downhill. The
rigid motions of the whole ring, rebuilt from the current beads, are left
out of the step; the imaginary-time cycle keeps its component, on its own
curvature, and is lifted to a spring-sized curvature only while that
curvature is near zero. A
converged gradient is classified with finite-difference bead Hessians
before it counts: the Bofill blocks learn curvature only along the steps
taken, so they can carry negative curvatures the surface does not have and
miss ones it does. A second negative curvature that survives is a
higher-index stationary ring, and the search steps down that mode by a
trust-sized displacement. The ring Hessian is block cyclic
tridiagonal in the beads, and every solve, determinant and inertia count
goes through a block LU of the open chain plus a low-rank Woodbury
correction for the closure, the cycle, the rigid modes and each
eigenvector-following flip: {math}`O(N f^3)` for {math}`f` degrees of
freedom, and the {math}`Nf \times Nf` matrix is never formed, so the same
step serves a seven-atom cluster and a 254-atom cell. The lowest ring
modes come from Lanczos on matrix-vector products. The bead curvature
blocks start from the saddle Hessian (`initial_hessians = saddle`) and
follow accepted moves with a Bofill update, rebuilt from finite
differences when their count of negative curvatures grows and up to three
times when the trust radius reaches its floor. Below the crossover the
copied blocks already hold a negative curvature in the ring's slowest
mode, so the first entry rebuilds them, and `saddle` costs the same force
calls as `finite_difference`. Near-zero eigenvectors of the saddle Hessian are held at a
spring-sized curvature in the Newton step, so a rigid displacement does
not singularize the chain. `initial_hessians = finite_difference` takes
{math}`2f` gradient calls per
bead first. The rate uses the bead Hessians chosen by `hessian_final`, not
that update. On the one-dimensional Eckart barrier the search converges
in 4 to 7 steps from either seed.

`temperatures` is a comma-separated list in kelvin. The search starts at
the highest and each ring starts the next, colder one. An empty list uses
`temperature`. `bead_ladder` (default off) starts a ring of at least 16
beads at a quarter of that count and doubles. `hessian_final` is
`recomputed`: the rate takes a Hessian on every stride-th bead.

On a rate calculation, `initial_path` is a band over the barrier. The ring
starts on the closed orbit of that band whose period is {math}`\beta \hbar`,
and the same band carries a one-dimensional WKB rate.
`rate_instanton.dat` has one row per temperature, with columns `T_K`,
`T_c_K`, `beads`, `converged`, `iterations`, `U_N_eV`, `negative_modes`,
`ln_k_per_s`, `k_per_s`, `ln_k_htst_per_s`, `barrier_effective_eV`,
`ln_k_wkb_path_per_s`, `ln_k_parabolic_per_s` and `parabolic_factor`.
Above {math}`T_c` the instanton columns are empty and the last two hold
the parabolic rate. Below {math}`T_c` those two are empty. At {math}`T_c`
the parabolic columns are empty. A run at more than one temperature also writes
`instanton_<T>K.con`. The coldest temperature is written to `instanton.con`.

With no atom fixed, both the reactant and the instanton omit the three
translations. A
rotation is omitted when it is a zero mode of the reactant Hessian, which a
free cluster has and a crystal does not: when its curvature r^T H r / r^T r
is at most 0.1 of the softest vibration, the lowest eigenvalue of the
Hessian once the translations and rotations are projected out. Both sides
scale with the Hessian and neither grows with the atom count. A cluster in a large periodic cell
is told apart by that Hessian, not by the periodic flag. The springs along
those directions stay, so they cancel between the instanton and the reactant.
On the ring the omitted rotations are those of the beads themselves about
the ring's centre of mass: a rotation moves each bead by a different
amount, so the reactant's rotation copied to every bead is not a zero mode
of the ring. The bead Hessians keep their rotational curvature, which
balances the springs on a bead that is not a minimum, and lose only their
translations. The omitted rotations come back as the ratio of the
classical rotational partition functions: the square root of the product
of the ring's principal moments of inertia, every bead's atoms about the
ring's centre of mass, over the same product for the reactant repeated on
every bead. That is the moment of inertia of the ring polymer averaged
over its beads, as the ring-polymer instanton rate takes it for a
molecule that turns freely. `instanton_rotation_ratio_log` reports its
logarithm. The classical comparison `rate_htst_log` and the parabolic rate
above {math}`T_c` take the saddle's ratio, so all three rates count the
rotations alike. A linear structure keeps the two moments perpendicular to
its axis.

`results.dat` reports the rate. The keys are:

| Key | Meaning |
|---|---|
| `rate_instanton` | {math}`k` in s^{-1} |
| `rate_instanton_log` | {math}`\ln(k\,/\,\mathrm{s}^{-1})` |
| `rate_htst_log` | classical harmonic transition-state theory, the same logarithm |
| `parabolic_factor` | {math}`\kappa`, above {math}`T_c` |
| `rate_parabolic` | {math}`\kappa` times the quantum harmonic TST rate, s^{-1} |
| `rate_parabolic_log` | {math}`\ln(k\,/\,\mathrm{s}^{-1})` of that rate |
| `instanton_crossover_K` | {math}`T_c`, K |
| `instanton_negative_modes` | negative eigenvalues of the ring Hessian; a first-order saddle has 1 |
| `instanton_zero_mode` | the eigenvalue left out |
| `instanton_rotation_ratio_log` | ln of the ring's rotational partition function over the reactant's; 0 when no rotation is free |
| `instanton_s_min`, `instanton_s_max` | the turning points along the saddle's unstable mode, amu^0.5 Angstrom from the saddle |
| `instanton_chord_overlap` | cosine of the angle between the chord joining the turning points and the unstable mode |
| `instanton_collapsed` | 1 when the search stopped because the beads fell onto one point (B_N below 1e-4 of the starting ring's) |
| `instanton_crossing_offset` | largest distance from the saddle, across the mode, at which the ring crosses the dividing plane |

A converged ring is given a rate only when it belongs to the seeded saddle:
it straddles the plane through the saddle normal to the unstable mode,
crosses that plane no farther from the saddle than its own span
`instanton_s_max - instanton_s_min`, and its chord lies within 60 degrees
of the mode. A ring that slid to a neighbouring saddle can still straddle
the plane and have one negative mode, and then its rate would belong to a
different reaction.

{cite:t}`inst-habershonRingpolymerMolecularDynamics2013` expect the sampled
ring-polymer rate to lie within about a factor of two of the exact quantum
rate between {math}`T_c` and {math}`T_c/2`. That bound compares the sampled
rate with the exact rate. It is not a comparison of this instanton with a
free-energy profile.

The cubic metastable well, {math}`V = \omega_0^2 q^2/2 - g q^3/3`, is the
check. Deep below the crossover its rate approaches the zero-temperature
decay of {cite:t}`inst-caldeiraQuantumTunnellingDissipative1983`.

## Path-integral quantum TST on planes

With `pi_planes` above 0, `mode = rate` follows the instanton with
path-integral quantum transition-state theory (PI-QTST,
{cite:t}`inst-vothRigorousFormulationQuantum1989`; review in
{cite:t}`inst-vothFeynmanPathIntegral1993`). A ring polymer is sampled with
its centroid held on each of a set of parallel planes. The centroid
potential of mean force along the plane coordinate gives a free-energy
barrier and a rate at every temperature in `temperature` or `temperatures`.

```ini
[Instanton]
mode = rate
saddle_filename = saddle.con
temperatures = 150, 105
pi_planes = 21
pi_beads = 32
pi_equilibration_steps = 500
pi_sampling_steps = 4000
pi_time_step = 0.5
pi_thermostat = pile
pi_direction = mode
pi_reactant_extent = 0.5
```

The coordinate. With {math}`q` the mass-weighted displacement from the
reactant over the free atoms and {math}`\hat n` a unit vector in those
coordinates, the plane coordinate is {math}`s = \hat n \cdot q`, in
amu^0.5 Å. The reactant sits at {math}`s = 0` and the saddle at
{math}`s^* = \hat n \cdot q_\mathrm{saddle}`. `pi_direction = mode` (the
default) takes {math}`\hat n` from the unstable eigenvector of the saddle's
mass-weighted Hessian, oriented toward the saddle, so the last plane is the
dividing surface normal to the barrier mode. When that mode makes more than
60 degrees with the reactant-saddle line, the job warns and uses
`pi_direction = line`, the straight mass-weighted line from the reactant to
the saddle. (A saddle with no negative eigenvalue has no crossover
temperature, and the job stops before the planes.) One normal serves every
plane, so {math}`s` is a linear coordinate and its mean force integrates to
its free energy with no metric correction. The planes are fixed in the reactant's frame. Translations
of a structure with no atom fixed lie within every plane, because the
unstable mode of the projected Hessian has no rigid component, but a
rotation of a free cluster changes {math}`s`; fix an atom, or use a cell,
when the sampling is long enough for the cluster to turn. The job warns
when no atom is fixed in an aperiodic cell. `pi_planes` planes are spaced
uniformly from {math}`s_0 = -\,\mathtt{pi\_reactant\_extent}\; s^*`, behind
the reactant, to {math}`s^*`.

The sampling. On each plane the ring's centroid starts where
`initial_path` crosses the plane, or on the reactant-saddle line when no
band is given. Projecting the centroid position and momentum holds it on
the plane. The ring is carried from plane to plane: its centroid moves to
the next plane and its internal modes keep their thermalised state.
The free-ring springs are propagated exactly, so `pi_time_step` is
limited by the physical vibrations, as for classical dynamics.
`pi_equilibration_steps` are discarded, then `pi_sampling_steps` steps of
`pi_time_step` fs record {math}`\hat a \cdot f_c`, the centroid force along
the plane's Cartesian normal {math}`a = M^{1/2} \hat n`. The bead forces of
a step are one batch, so a calculator group carries the beads. `pi_thermostat = pile` puts PILE on the internal
modes and a Langevin thermostat of time `pi_pile_tau` fs on the centroid
within the plane, and `pi_pile_scale` (default 1) scales the critical
damping of the internal modes. Below the crossover the lowest ring modes
are soft at the barrier and overdamped at critical damping; on the Eckart
check 0.5 cuts the mean-force error at 0.7 {math}`T_c` by a third. `piglet` reads a normal-mode GLE from `pi_gle_file`, in
the same format and with the same meaning as `[Dynamics] path_gle_file`.
`pi_seed` seeds the noise. PIGLET tunes the internal modes so that bead
averages converge with fewer beads; the centroid free energy assumes the
ring polymer's own Boltzmann distribution, so treat a `piglet` barrier as
an approximation and check it against `pile`. The transmission factor's
parents always take `pile`.

The free energy. The mean force on the plane at {math}`s` is

```{math}
F'(s) = -\frac{\left\langle a \cdot f_c \right\rangle_s}{|a|^2} ,
\qquad a = M^{1/2} \hat n ,
```

in eV per amu^0.5 Å. Any vector {math}`u` with {math}`a \cdot u = 1`,
{math}`M^{-1/2} \hat n` among them, gives the same mean, because the force
within the plane averages to zero; the normal component is the estimator
the job records. {math}`F(s)` is the trapezoid integral of {math}`F'` from
{math}`s_0`. The production run is cut into ten equal blocks. The standard
error of the block means is the error of {math}`F'(s)`. The errors of
{math}`F` and of the rate follow by linear propagation, with the planes
taken as independent. The quantum free-energy barrier is
{math}`\Delta F = F(s^*) - \min_{s < s^*} F(s)`, the minimum over the
planes before {math}`s^*`; the job warns when it is not positive, since
{math}`s^*` is then no barrier on the centroid free energy. The log and `results.dat`
give it beside the classical barrier
{math}`V(\mathrm{saddle}) - V(\mathrm{reactant})` and the instanton's
effective barrier.

The rate. With {math}`s` a unit-mass coordinate,

```{math}
k_\mathrm{PI\text{-}QTST} = \frac{1}{2}\sqrt{\frac{2}{\pi\beta}}\;
\frac{e^{-\beta F(s^*)}}{\int_{s_0}^{s^*} e^{-\beta F(s)}\,ds},
```

with {math}`\beta = 1/k_B T` in eV^{-1}. The prefactor is
{math}`\tfrac12\langle|\dot s|\rangle`, half the mean speed of a free unit
mass, in amu^0.5 Å per time unit of
{math}`\sqrt{\mathrm{amu}\,\mathrm{Å}^2/\mathrm{eV}}` (10.18 fs). The
integral, again the trapezoid rule over the planes, is the reactant's
centroid density along {math}`s`. The job reports {math}`k` in s^{-1}. The
integral stops at the first plane, so that plane must sit several
{math}`k_B T` above the reactant minimum of {math}`F`; the job warns below
5 {math}`k_B T`. It also warns when {math}`F'(s^*)` differs from zero by
more than three standard errors: the maximum of the centroid free energy
then lies off the classical saddle's plane.

The formula is classical TST on the centroid free-energy surface. With one
bead, {math}`F` is the classical free energy and the rate is classical TST
along {math}`s`. Near and above the crossover it carries the tunnelling
correction of a symmetric barrier. For strongly asymmetric barriers well
below the crossover it fails {cite:p}`inst-vothFeynmanPathIntegral1993`,
and the instanton gives the rate there.

Output. `rate_piqtst.dat` has one row per temperature and plane, with
columns `T_K`, `s_amu05A`, `dF_ds_eV_per_amu05A`, `dF_ds_error`, `F_eV`,
`F_error_eV` and `spread_max_A`. `piqtst_planes.con` holds one frame per
plane at the coldest temperature, at the production-averaged centroid. Its
readcon `spreads` section holds each atom's root-mean-square bead
displacement from the centroid along x, y and z in Å. The scalars carry
{math}`s`, {math}`F'` and {math}`F` with their errors, the temperature and
the bead count. A run at more than one temperature also writes
`piqtst_planes_<T>K.con`. `results.dat` gains, for the coldest temperature:

| Key | Meaning |
|---|---|
| `barrier_piqtst` | {math}`\Delta F`, eV |
| `barrier_piqtst_error` | its standard error, eV |
| `rate_piqtst` | {math}`k_\mathrm{PI\text{-}QTST}`, s^{-1} |
| `rate_piqtst_log` | {math}`\ln(k\,/\,\mathrm{s}^{-1})` |
| `rate_piqtst_log_error` | standard error of that logarithm |
| `piqtst_s_star` | {math}`s^*`, amu^0.5 Å |
| `piqtst_dF_ds_star`, `piqtst_dF_ds_star_error` | {math}`F'(s^*)` and its error |
| `piqtst_first_plane_kT` | {math}`\beta (F(s_0) - \min F)` |
| `piqtst_temperature_K`, `piqtst_planes`, `piqtst_beads` | the run |

Cost: `pi_planes` times (`pi_equilibration_steps` plus `pi_sampling_steps`
plus 1) force batches of `pi_beads` beads per temperature: one per step,
and one when the ring moves to a plane.

### Recrossing and the RPMD rate

PI-QTST counts every centroid that reaches {math}`s^*` moving forward as
reactive. Some of those trajectories turn back. With
`pi_recrossing_parents` above 0 the job measures that fraction, the
Bennett-Chandler transmission factor of ring-polymer molecular dynamics
(RPMD, {cite:t}`inst-craigChemicalReactionRates2005`; the two-step scheme
of {cite:t}`inst-suleimanovRPMDrateBimolecularChemical2013`), and reports

```{math}
k_\mathrm{RPMD} = \kappa\, k_\mathrm{PI\text{-}QTST} .
```

```ini
[Instanton]
pi_recrossing_parents = 100
pi_recrossing_children = 20
pi_recrossing_time = 100
pi_recrossing_spacing = 50
```

The parents. A ring with its centroid held on the top plane {math}`s^*`,
the dividing surface of the scan, is thermostatted with PILE for
`pi_equilibration_steps`, as the scan's `pile` is: Bennett-Chandler needs
the ring polymer's own constrained distribution, which PIGLET does not
sample. Every `pi_recrossing_spacing` steps after that its
beads are one parent, `pi_recrossing_parents` in all.

The children. Each parent launches `pi_recrossing_children` momentum draws.
A draw takes every ring normal mode from the Maxwell-Boltzmann distribution
at {math}`\beta_P = \beta / P` with the plane constraint removed, and runs
twice, with {math}`p` and with {math}`-p`, which halves the variance. Each
child runs `pi_recrossing_time` fs of thermostat-free RPMD in steps of
`pi_time_step`: velocity Verlet on the physical forces with the free ring
propagated exactly in normal modes. Along the way the job records the
centroid coordinate {math}`s(t)`. The centroid velocity along the coordinate is
{math}`\dot s = \hat n \cdot M^{1/2} v_c`, and

```{math}
\kappa(t) = \frac{\langle \dot s(0)\, h(s(t) - s^*) \rangle}
                  {\langle \dot s(0)\, h(\dot s(0)) \rangle} ,
```

with {math}`h` the step function. At {math}`t = 0` the side is that of
{math}`\dot s(0)`, the limit {math}`t \to 0^+`, so {math}`\kappa(0) = 1`.
{math}`\kappa(t)` falls as children recross and levels off once they have
committed to a side. The reported {math}`\kappa` is the mean of
{math}`\kappa(t)` over the last quarter of `pi_recrossing_time`. Its
standard error is the jackknife over parents. Lengthen
`pi_recrossing_time` when `kappa_piqtst.dat` has not levelled off by the last
quarter.

What {math}`\kappa` means. {math}`\kappa` lies between 0 and 1. A value
of 1 means no trajectory that crosses {math}`s^*` forward returns, and
PI-QTST is the RPMD rate. A value below 1 means the plane is not the
dynamical bottleneck: either it is tilted from the barrier's own dividing
surface, or a bath coupling turns trajectories back. One bead gives the
classical transmission through the plane. On a harmonic saddle whose plane
normal makes an angle {math}`\theta` with the unstable mode it is
{math}`\sqrt{\cos^2\theta - \sin^2\theta\,\omega_b^2/\omega_\perp^2}`.
The RPMD rate is independent of the choice of {math}`s^*` in the
long-time limit. The PI-QTST rate is not.

Output. `kappa_piqtst.dat` has the coldest temperature's curve, columns
`t_fs` and `kappa`. A run at more than one temperature also writes
`kappa_piqtst_<T>K.dat`. `rate_piqtst.dat` gains the columns `kappa`,
`kappa_error`, `ln_k_rpmd_s` and `ln_k_rpmd_s_error`, the same on every
plane of a temperature. `results.dat` gains:

| Key | Meaning |
|---|---|
| `piqtst_kappa`, `piqtst_kappa_error` | {math}`\kappa` and its standard error |
| `rate_rpmd` | {math}`k_\mathrm{RPMD}`, s^{-1} |
| `rate_rpmd_log` | {math}`\ln(k_\mathrm{RPMD}\,/\,\mathrm{s}^{-1})` |
| `rate_rpmd_log_error` | its standard error, both errors in quadrature |

| Option | Default | Meaning |
|---|---|---|
| `pi_recrossing_parents` | 0 | parent configurations; 0 is off, otherwise at least 2 |
| `pi_recrossing_children` | 20 | momentum draws per parent, each run forward and reversed |
| `pi_recrossing_time` | 100 | fs per child, at least four `pi_time_step` |
| `pi_recrossing_spacing` | 50 | thermostatted steps between parents |

Cost per temperature: `pi_equilibration_steps` plus
`pi_recrossing_parents` times `pi_recrossing_spacing` constrained steps,
one force batch each, and `2 * pi_recrossing_parents * pi_recrossing_children`
child runs of `pi_recrossing_time / pi_time_step + 1` force batches each,
the first at the parent's beads.

## Centroid and spread

Both modes also write `instanton_centroid.con` (and
`instanton_centroid_<T>K.con` per temperature): one frame at the mean of
the beads, with the readcon `spreads` section holding each atom's
root-mean-square displacement from that mean along x, y and z in Å. This is
the delocalised configuration as a centroid plus a per-atom spread, the same
representation a path-integral trajectory collapses to, so the atoms that
tunnel and how far they spread read off one frame. Below the crossover the
density is bimodal along the reaction path and the spread there is a width,
not a Gaussian; the beads in `instanton.con` keep the full path. The frame
carries `spread_max`, the largest entry, beside the temperature.

## Checks

These checks use two Catch2 cases. The curved-valley case is `Instanton splitting in a curved valley matches the exact gap`. The corner case is `The instanton cuts the corner the minimum energy path takes`.
The potential is {math}`V = V_0 (x^2 - 1)^2 + \tfrac{K}{2} (y - C (1 - x^2))^2` at unit mass
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

For the rate, `The Eckart rate instanton matches the exact flux to its
semiclassical error` compares {math}`k Z_r` through the symmetric Eckart
barrier {math}`V_0 / \cosh^2(x/a)` ({math}`V_0 = 0.425` eV, {math}`a =
0.734` amu^0.5 Å, {math}`T_c = 150` K) with the exact flux
{math}`(2\pi\hbar)^{-1}\int P(E) e^{-\beta E} dE` from Eckart's transmission
probability, at {math}`T = 0.5\,T_c` and {math}`0.35\,T_c`:

| beads | instanton / exact |
|---|---|
| 64 | 0.94 to 0.96 |
| 128 | 0.93 to 0.94 |
| {math}`N \to \infty` (1/N² extrapolation) | 0.928 |

In one dimension the instanton is the steepest-descent evaluation of the
WKB thermal integral, so its limit shares the uniform WKB error; the Kemble
integral along the path gives the same 0.928. `The ring spectrum from the
block chain matches the dense Hessian` ties the chain's determinant and
inertia to a dense eigendecomposition, `A rigid mode leaves the
instanton rate unchanged` checks the rigid-mode bookkeeping through the
search, and `The rate lifts the rotations of a free diatomic's ring`
checks the rotations against the dense spectrum of a stretching bond.

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: INST_
keyprefix: inst-
---
```
