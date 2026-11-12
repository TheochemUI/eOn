---
myst:
  html_meta:
    "description": "Variable-cell minimization and the solid-state dimer: min-mode saddle search on atoms and cell strain."
    "keywords": "eOn solid-state dimer, variable cell, Cauchy stress, elastic dipole, Jacobian"
---

# Solid-state dimer

A solid-state dimer is a min-mode saddle search on the atoms and the cell
together. The min-mode translation is the saddle search in
{cite:t}`ssd-goswamiEfficientImplementationGaussian2025` and
{cite:t}`ssd-goswamiBayesianHierarchicalModels2025a`, using the dimer estimate
of {cite:t}`ssd-henkelmanDimerMethodFinding1999` and the comparison in
{cite:t}`ssd-olsenComparisonMethodsFinding2004`. The cell block is the
generalized solid-state coordinate of
{cite:t}`ssd-sheppardGeneralizedSolidstateNudged2012`, the same Jacobian the
solid-state band uses. AKMC, and the clients that drive the same saddle
search (amsel, pyKMC), then find a saddle whose cell has relaxed.

The cell stays frozen unless `[Solid State] enabled` is true. A mask of all
zeros also leaves the cell frozen. Either way the fixed-cell dimer and the
fixed-cell minimizer are unchanged.

## Coordinate

The cell is kept lower triangular: the first lattice vector lies on x and
the second lies in the xy plane, which removes the three rotations of the
cell. The six remaining strains are xx, yx, yy, zx, zy, zz. Relative to the
cell at the start of the search, h0,

```{math}
h = h_0 (I + \varepsilon)
```

with {math}`\varepsilon` lower triangular. Atomic coordinates in the chart
are the Cartesian positions in that reference cell, so a pure cell step
keeps fractional coordinates fixed. The scaled strain is

```{math}
q = J \varepsilon, \qquad J = V_0^{1/3} N^{1/6} w
```

{math}`V_0` is the reference volume, {math}`N` the number of atoms, and
{math}`w` is `weight`. This is
{math}`(V_0/N)^{1/3} N^{1/2} w`, the Jacobian of Sheppard, Xiao and
Henkelman. {math}`J` is fixed for the search. The generalized force on
{math}`q` is the derivative of the energy with respect to that scaled
strain, from the Cauchy stress the potential reports,

```{math}
\sigma = \frac{1}{V} \frac{\partial E}{\partial \varepsilon}
```

for the right strain {math}`h \leftarrow h(I+\varepsilon)` at fixed
fractional coordinates. An external stress, including a hydrostatic
pressure, is added before that derivative. A potential that does not report
the stress is refused.

ASE's `FrechetCellFilter` ({cite:t}`ssd-larsenAtomicSimulationEnvironment2017`)
uses the matrix logarithm of the deformation gradient and the Fréchet
derivative of the exponential, and scales the cell block by the atom count.
The dimer keeps Sheppard's Jacobian so a solid-state saddle and a
solid-state band are in the same metric. The mask and the external stress
follow the filter's interface: any of the six strains can be frozen, and
the external stress is a full tensor, not only a pressure.

`min_mode_method` may be `dimer` or `lanczos`. The dimer rotation is the
classical finite-difference rotation in the joint space. Lanczos builds the
same joint Hessian-vector product. Other min-mode methods are refused, as
are an LBFGS pair preconditioner and rigid-body projection. This is not a
change to the solid-state band: dynamic NEB and image freezing are not part
of this search.

`index_check` is on by default for this search, including a mask that
freezes every strain. A historical dimer with `enabled` false does not run
it. When the force tolerance is met, the search builds the joint Hessian
by central differences in the chart `J = V^{1/3} N^{1/6}` and projects out
rigid translations. The point is kept only when exactly one eigenvalue lies
below -1e-4. A higher-index point is left along its softest negative mode,
and the search continues. If nothing is accepted, the result is no barrier
rather than that point. The potential calls spent on the Hessian and on
that walk are counted apart from the dimer steps. `flat` is 0.005 by
default. While the curvature {math}`\lambda` is negative, the largest
force component must also fall to {math}`0.005|\lambda|`, as well as under
the usual tolerance. Set `flat` to 0 to leave the usual tolerance alone,
or `index_check` to false to accept the minimum-mode test by itself.

## Minimizer

The same chart is the objective of a local minimization when `enabled` is
true. Conjugate gradients, LBFGS, FIRE, steepest descent and quickmin all
see the joint gradient. `converged_force` is compared with that gradient.
The default `norm` is the Euclidean norm of the atomic forces and the cell
forces together.

## Elastic dipole

When the cell is free, `results.dat` carries the Cauchy stress, the cell
volume, and the elastic dipole {math}`P = V \sigma` at the minimum and at
the saddle. A minimization writes `stress_xx` ... `stress_zz`,
`cell_volume`, and `elastic_dipole_xx` ... `elastic_dipole_zz`. A saddle
search prefixes those keys with `saddle_` and `reactant_`. A process search
writes `reactant_`, `saddle_` and `product_`. Callers that want the dipole
can read `elastic_dipole_*` or multiply `stress_*` by `cell_volume`.

## Copper interstitial

The shear benchmark uses the analytic Mishin EAM1 of Table I in
Phys. Rev. B 63, 224106. An fcc crystal at 3.615 Angstrom has energy
-3.5402183104868774 eV per atom. The setfl file is that table rounded
to -3.54 eV per atom, so `potential = eam_cu` does not reproduce these
barriers. The defect is a <100> dumbbell on one site of a 3x3x3 cell,
109 atoms. The shear is the xy map `h = L * [[1, 0, 0], [gamma, 1, 0],
[0, 0, 1]]`, and the positions move with it, so `x' = x + gamma y`.
The fixed-cell barrier falls as that shear grows. A variable-cell
search sets `stress` to the negative of the stress at the minimum, so
the enthalpy is stationary there. Both searches use `index_check` and
`flat = 0.005`. The fixed-cell barriers are 0.4738130973207717,
0.2435460895545134, 0.14829201811522807, 0.067862165555311549 and
0.0077937374149428251 eV. The variable-cell enthalpy barriers at 0, 2, 3
and 4% are 0.34700693547728179, 0.16255571725275786, 0.088051942122888249
and 0.027924583692424676 eV. Each of those points has one projected
eigenvalue below -1e-4 and is the shear-aligned <110> switch. At 5%
variable the reactant is a minimum (lowest projected eigenvalue
0.0028631813667541446) and the search stops on an index-1 point whose
enthalpy is 8.1226496035924356e-6 eV above that reactant and which is
not the switch, so that row is not a barrier. Raw force calls and
`n_neg` are in `benchmarks/solid_dimer_cu.csv`.

`max_strain_step` limits how far any strain component may move in one
update. `max_strain` is the largest strain relative to the cell at the
start of that run. Past it the update stops with an error, so a saddle
search cannot shear the cell without a bound. The imposed shear belongs
in the starting cell. Freeze that component in `mask` when it must stay
at the applied value. Stiff cell directions are scaled to an atomic
curvature on top of the Jacobian, so the dimer rotation does not spend
its force calls there. The scale does not change the stationary points.

```{code-block} ini
[Main]
job = process_search

[Potential]
potential = eam_cu

[EAM]
file = Cu01.eam.alloy

[Solid State]
enabled = true
weight = 1.0
pressure = 0.0
stress = 0 0.001 0 0.001 0 0 0 0 0
mask = 1 1 1 1 1 1
max_strain_step = 0.05
max_strain = 0.20
index_check = true
flat = 0.005

[Saddle Search]
min_mode_method = dimer

[Optimizer]
opt_method = cg
```

`stress` is row-major. The example adds a yx shear of 0.001 eV/Angstrom^3.
`mask` lists xx, yx, yy, zx, zy, zz.

## Configuration

```{eval-rst}
.. autopydantic_model:: eon.schema.SolidStateConfig
```

```{eval-rst}
.. autopydantic_model:: eon.schema.EamConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: SSD_
keyprefix: ssd-
---
```
