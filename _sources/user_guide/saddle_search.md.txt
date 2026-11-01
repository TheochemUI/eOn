---
myst:
  html_meta:
    "description": "Configuration for initiating saddle searches in eOn by displacing atoms from a minimum energy state."
    "keywords": "eOn saddle search, transition state, displacement, potential energy surface"
---

# Saddle search

A saddle search is initiated by making a local displacement of atoms from their
position at the minimum of the current state. This displacement can be done
using the different strategies indicated by the {any}`client_displace_type`
option, and the following parameters. If the user knows something about the
local environment where reactions are likely to take place in the system, this
information can be used to make saddle searches more efficient by getting them
started in the right part of configuration space.

## Displacement strategies

### Epicenters and Weight-Based Selection

Each saddle search begins by choosing an **epicenter** — a single atom around
which the initial displacement is constructed. The epicenter is selected
probabilistically from several strategies, each controlled by a relative weight:

| Weight parameter                       | Strategy                                                       |
| -------------------------------------- | -------------------------------------------------------------- |
| `displace_random_weight`               | Pick any atom at random                                        |
| `displace_least_coordinated_weight`    | Pick an atom with the lowest coordination number               |
| `displace_not_FCC_HCP_weight`          | Pick an atom whose local environment is neither FCC nor HCP    |
| `displace_under_coordinated_weight`    | Pick an atom with coordination ≤ `displace_max_coordination`   |
| `displace_listed_atom_weight`          | Pick an atom from `displace_atom_list`                         |
| `displace_listed_type_weight`          | Pick an atom whose type is in `displace_type_list`             |

The weights do not need to sum to 1; they are normalised internally. Setting a
weight to 0 disables that strategy. For example, setting only
`displace_listed_atom_weight = 1.0` and all other weights to 0 starts
every saddle search from an atom in the explicit list.

### Radius and Magnitude

Once the epicenter is chosen, all atoms within `displace_radius` of the
epicenter are displaced. Each displaced atom receives a random perturbation
drawn from a Gaussian distribution with standard deviation
`displace_magnitude` (in Ångströms) independently in each Cartesian direction.

### Displacing all listed atoms

When `displace_all_listed` is **true** and a listed-atom strategy is selected,
*every* atom in `displace_atom_list` (or `displace_type_list`) is displaced,
not a single random pick. Atoms within `displace_radius` of *any*
displaced atom are also included. Set `displace_radius` to 0 to restrict the
displacement strictly to the listed atoms.

### Dynamic atom lists via scripts

For systems where the relevant atoms change from state to state (e.g. a
migrating vacancy), a static list is insufficient. The
`displace_atom_kmc_state_script` option lets you specify a Python script that
is executed once per new AKMC state to determine the atom list dynamically. See
the {doc}`../tutorials/displacement_scripts` tutorial for worked examples
covering vacancy diffusion and adsorbate-on-surface scenarios.

### Client-Side Displacement

The `client_displace_type` option selects how the **client** (C++ code) picks
the epicenter when the displacement is performed client-side rather than by the
server:

- `random` — uniform random atom
- `last_atom` — the last atom in the configuration
- `least_coordinated` — the atom with the fewest neighbours
- `not_fcc_hcp_coordinated` — an atom whose local structure is neither FCC nor HCP
- `listed_atoms` — pick via `listedAtomEpiCenter` from `displace_atom_list`
  (CON file-order, remapped; lone `-1` is every free atom) and apply a
  `displace_radius` / `displace_magnitude` kick. No server displacement
  file needed. `random`, `last_atom`, and `least_coordinated` use the
  same client helper.
- `load` — read a displacement vector from a file written by the server

## Basin hopping

`method = basin_hopping` quenches the hop, builds a band, and starts a
dimer on the highest interior image. The dimer direction is the
minimum-image `(next - prev) / 2`. A bead that has crossed a periodic
face takes the short step. The acceptance test is on the
[basin hopping](project:basin_hopping.md) page.

A classic dimer (`[Dimer] improved = false`) translated with L-BFGS reports
the mode curvature `|C_tau|` once that mode is negative. The first translation
step keeps `H0 = 1/|C_tau|` when the step fits in `max_move`, and also when
the mode is milder than 30 eV/Å² and the cap would otherwise choose the
length. A steeper mode whose uncapped step exceeds `max_move` is measured
with one finite-difference probe of the effective force. Improved dimer and
Lanczos do not report `|C_tau|`. Their first step is that probe, and
`max_move` is only the cap. The adatom AKMC templates use the classic dimer.
The LJ13 and Pt heptamer force-call counts follow the probe. A build with
`-Dwith_rgsaddle` rotates every dimer and Lanczos search in rgsaddle. The
scale still follows `improved = false`.

## Configuration

```{code-block} ini
[Saddle Search]
```




```{eval-rst}
.. autopydantic_model:: eon.schema.SaddleSearchConfig
```

## Gaussian process dimer

`min_mode_method = gprdimer` runs the Gaussian process dimer from a
`gpr_optim` checkout in `subprojects/gpr_optim`, for example
`rsync -a ../gpr_optim/ subprojects/gpr_optim/`. Under `-Dwith_gprd=auto`
(the default) a Linux build links that checkout when it is there and leaves
the dimer off when it is not. `-Dwith_gprd=enabled` stops configuration
instead. `-Dwith_gprd=disabled` omits the method. The install page has the same `rsync` next to the
default `meson setup`.

## Output

The saddle search writes:
- `saddle.con`: the saddle point structure
- `mode.dat`: eigenvector at the saddle (Nx3 whitespace-separated)
- `results.dat`: energy, force calls, convergence status, eigenvalue

With `write_movies = true` (in `[Debug]`), `climb.con` is written as a
concatenated structure movie (one frame per iteration). Each frame stores
structured JSON metadata in the second CON header line via `readcon-core`, including `energy`,
`frame_index`, `step_size`, `delta_e`, `convergence`, `eigenvalue`, `torque`,
`angle`, and `rotations`.

Set `write_deprecated_outs = true` in `[Debug]` to also emit the legacy
`climb.dat` sidecar during the compatibility window.

## Refining a converged saddle

This section covers partitioned rational function optimization (P-RFO) refinement.
It takes a saddle from a dimer or Lanczos search in eOn and lowers its max-atom force further.
The searches it continues are those of {cite:t}`ss-goswamiEOnV3602026` and {cite:t}`ss-goswamiEfficientImplementationGaussian2025`, and for a band, the climbing image of {cite:t}`ss-goswamiEnhancedClimbingImage2026`.
The step is the P-RFO step of {cite:t}`ss-bakerAlgorithmLocationTransition1986` in the trust region of the Sella optimizer of {cite:t}`ss-hermesSellaOpenSourceAutomationFriendly2022`, as rgsaddle implements it.
The refinement is off by default.

### How to refine a dimer saddle

This recipe takes a dimer search that converges at 1e-3 eV/Å and refines its saddle to 1e-5 eV/Å.

1. Build the client against rgsaddle with `meson setup bbdir -Dwith_rgsaddle=true`.
   The refinement needs rgsaddle C application binary interface (ABI) 1.15 or newer, which the pinned rgsaddle 0.2.2 provides.
   A client built without it stops with an error once the search converges.
   On a band or another large free set the refinement is faster with rgsaddle 0.3.2, as the measurements below show, but eOn's continuous integration (CI) tests only the pinned rgsaddle.
2. Add the refinement keys to `[Saddle Search]`, next to the search's own settings.

   ```{code-block} ini
   [Saddle Search]
   min_mode_method = dimer
   converged_force = 0.001
   prfo_refine = true
   prfo_refine_force = 1e-5
   ```

3. Run the job as usual.
   After the search converges, the log reports the refinement in one line.
   For `benchmarks/data/one_pt_saddle_search` at the pinned rgsaddle it reads as follows.

   ```{code-block} text
   [P-RFO] 7 steps, 7 force calls, max atom force 2.434e-04 -> 1.093e-08 eV/A
   ```

4. Read the outputs.
   `saddle.con` is the refined geometry.
   `mode.dat` and `final_eigenvalue` in `results.dat` come from one more min-mode estimate at the refined point.
   `force_calls_saddle` and `total_force_calls` include the refinement and that estimate.

If the log line ends with `kept the starting point`, the refinement did not reach a smaller force and the search's own saddle is written unchanged.
A nudged elastic band (NEB) job reads the same keys.
It refines its climbing image, as {ref}`neb-prfo-refine` describes.

### Reference

The keys are `prfo_refine`, `prfo_refine_force` and `prfo_refine_max_iterations`, listed in the configuration model above.
The refinement follows these rules.

- It runs only after a search ends with status good and a negative eigenvalue.
  The Gaussian process dimer (`gprdimer`) is not refined.
  Process searches in adaptive kinetic Monte Carlo (AKMC) use the same search, so the keys apply there too.
- The search hands over its geometry, its lowest mode and that mode's curvature.
  The session starts from an identity Hessian model with that one mode set to that curvature.
  It updates the model with the TS-BFGS update (a Broyden-Fletcher-Goldfarb-Shanno variant for transition states).
- Only free atoms move.
  A fixed coordinate of a partly fixed atom keeps its starting value.
  There is no projection of rigid-body motion, because a periodic cell or a slab with fixed atoms has none to remove.
- The refinement stops when the max-atom force is at most `prfo_refine_force`, or after `prfo_refine_max_iterations` steps.
  At the pinned rgsaddle each step costs one force call.
  rgsaddle 0.3.2 also spends force calls on Hessian-vector products when it refreshes its curvature estimate, so some steps cost three.
- If the last point has a larger max-atom force than the starting saddle, the starting saddle is kept.
- Each step diagonalizes dense matrices with one row per free coordinate.
  At the pinned rgsaddle all of them use Jacobi rotations, so the time per step grows quickly on larger free sets.
  rgsaddle 0.3.2 moves the curvature estimate to a faster dense solver above 128 free coordinates.
  Its proposed P-RFO eigensolve series does the same for the P-RFO step above 32.
  The TS-BFGS update still diagonalizes the whole model with a cyclic Jacobi solver after every force call, and on a band of hundreds of free atoms this sets the time per step.
- The log can show the floating-point exception notice `FPE (continuing, masking further): overflow` once during the refinement.
  With rgsaddle 0.3.2 the faster dense solver raises it.
  That solver forms a scaled sum of squares that can overflow, then discards the sum.
  At the pinned rgsaddle it can appear late in a long refinement, from a rotation angle in the Jacobi solver that is too large to square, where the masked result is a zero rotation.
  The client masks the trap and continues.
  With the eigensolve series, runs with the trap turned off wrote the same outputs.

### When refinement pays off

Refinement starts from a converged saddle, so it is not a cold start.
It reuses the mode of the search.
Its benefit depends on how cheaply the search itself would continue.
The tables compare one run of each case on one machine, one thread per run.
"Search" is the repository input at 1e-3 eV/Å, "refined" adds `prfo_refine_force = 1e-5`, and "tightened" sets the search's own tolerance to 1e-5 eV/Å.
Each cell gives force calls, the max-atom force in eV/Å of the written geometry, evaluated again with a `point` job, and the wall time of the whole client run.

The first table gives the search, refined and tightened runs with rgsaddle 0.3.2 and its proposed P-RFO eigensolve series.
That series builds on the main branch at 276a984 and is not in a release yet.

| Case | Free atoms | Search | Refined | Tightened |
| --- | --- | --- | --- | --- |
| Dimer, `benchmarks/data/one_pt_saddle_search` | 1 | 39, 2.4e-4, 0.06 s | 43, 1.1e-7, 0.06 s | 54, 8.9e-6, 0.06 s |
| Dimer, Pt heptamer test system | 7 | 130, 2.9e-4, 0.11 s | 163, 2.4e-6, 0.16 s | 136, 4.0e-6, 0.11 s |
| Dimer, 13-atom Lennard-Jones (LJ13) test cluster | 13 | 158, 3.5e-4, 0.02 s | 201, 9.9e-6, 0.27 s | 176, 4.3e-6, 0.02 s |
| NEB, `benchmarks/data/neb_morse_pt` | 175 | 277, 1.7e-4, 3.1 s | 465, 7.4e-6, 1207 s | 427, 1.8e-6, 3.2 s |
| NEB, `examples/neb-al` | 401 | 240, 2.2e-4, 0.87 s | not finished after 3600 s | 373, 3.4e-6, 0.90 s |

With the pinned rgsaddle 0.2.2:

| Case | Free atoms | Search | Refined | Tightened |
| --- | --- | --- | --- | --- |
| Dimer, `benchmarks/data/one_pt_saddle_search` | 1 | 39, 2.4e-4, 0.07 s | 47, 1.1e-8, 0.06 s | 54, 8.9e-6, 0.06 s |
| Dimer, Pt heptamer test system | 7 | 180, 4.5e-4, 0.17 s | 213, 1.6e-6, 0.21 s | 192, 3.9e-6, 0.21 s |
| Dimer, 13-atom Lennard-Jones (LJ13) test cluster | 13 | 836, 3.4e-4, 0.02 s | 891, 6.2e-6, 0.61 s | 854, 2.0e-6, 0.02 s |
| NEB, `benchmarks/data/neb_morse_pt` | 175 | 277, 1.7e-4, 3.1 s | not finished; a run capped at 3 steps took 41 s | 427, 1.8e-6, 3.0 s |
| NEB, `examples/neb-al` | 401 | 240, 2.2e-4, 0.86 s | not finished; a run capped at 3 steps took 1583 s | 373, 3.4e-6, 0.97 s |

With one free atom, refinement reached 1.1e-7 eV/Å in 43 calls, against 54 calls for the tightened dimer.
With 7 and 13 free atoms it reached the target, but the tightened dimer needed fewer calls, 136 against 163 and 176 against 201.
On the 175-atom band the refinement converged in 181 steps, but it took 465 calls and 1207 s, against 427 calls and 3.2 s for the tightened band.
On the 401-atom band it did not finish within 3600 s, and a run capped at 3 steps took 1140 s.
rgsaddle 0.3.2 without the series took the same steps and calls on the four refinements that finished, and 1295 s on the 175-atom band.
The dimer searches also run through rgsaddle, so their calls depend on its version, and at 0.3.2 the heptamer and LJ13 searches need fewer calls than at the pin.
At the pinned rgsaddle a run capped at 3 steps took 41 s on the 175-atom band and 1583 s on the 401-atom band.
The uncapped 175-atom refinement at the pin had not finished after 8074 s of wall time.
Use the refinement on small free sets, after a converged search, and compare it with tightening the search itself.

## See also

- <project:dimer.md> for the eigenmode method used during the saddle search
- <project:process_search.md> for automated saddle search + minimization
- <project:neb.md> for finding paths when both endpoints are known
- {doc}`/devdocs/design/client/neb-modularization` for the eigenmode variant
  implementation (dimer, improved dimer, Lanczos, GPR dimer)

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: SS_
keyprefix: ss-
---
```
