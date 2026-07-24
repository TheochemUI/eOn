---
myst:
  html_meta:
    "description": "pyeonclient Python API: Matter, min-mode (Dimer/Lanczos/Davidson), MinModeSaddleSearch, NEB, Hessian, Prefactor."
    "keywords": "pyeonclient, Dimer, saddle search, NEB, Hessian, Matter, ASE"
---

# pyeonclient (Python client API)

```{versionadded} 2.16
In-process client: first-class algorithm objects on :class:`~pyeonclient.Matter`,
not a working directory and not the `eonclient` binary alone.
```

[pyeonclient](https://pypi.org/project/pyeonclient/) exposes the eOn **C++
client algorithms** to Python with nanobind. Geometries are
{py:class}`~pyeonclient.Matter`. Each method is a live class you construct and
call — the same engines the binary uses.

| Algorithm | Python class / function | Binary analogue |
|-----------|-------------------------|-----------------|
| Local min | `Matter.relax()` | `job = minimization` |
| Min-mode | `Dimer`, `ImprovedDimer`, `Lanczos`, `Davidson` | dimer / lanczos / davidson |
| Saddle search | `MinModeSaddleSearch` | `job = saddle_search` |
| NEB | `NudgedElasticBand` | `job = nudged_elastic_band` |
| Hessian | `Hessian` | `job = hessian` |
| Prefactor | `get_prefactors` | `job = prefactor` |
| Dynamics | `run_dynamics` | `job = dynamics` |
| Monte Carlo | `run_monte_carlo` | `job = monte_carlo` |
| Basin hopping | `run_basin_hopping` | `job = basin_hopping` |
| Process search | `process_search` | `job = process_search` |
| Structure compare | `structures_equal` / `structure_distance` | `job = structure_comparison` |
| TAD | `run_tad` → `TADJob` | `job = tad` |
| Parallel replica | `run_parallel_replica` → `ParallelReplicaJob` | `job = parallel_replica` |
| Safe hyperdynamics | `run_safe_hyperdynamics` → `SafeHyperJob` | `job = safe_hyperdynamics` |
| Replica exchange | `run_replica_exchange` → `ReplicaExchangeJob` | `job = replica_exchange` |
| GP surrogate NEB | `run_gp_surrogate_neb` → `GPSurrogateJob` (gated) | `WITH_GP_SURROGATE` |
| Opaque batch | `make_job` + `Job.run` | any `JobType` via factory |

`make_job(params)` still constructs every `JobType` the C++ factory supports
(process search, dynamics, basin hopping, TAD, …) for workdir-oriented runs.
Prefer the first-class classes when you want Matter in / Matter out.

## Install

```{code-block} bash
pip install pyeonclient
pip install 'pyeonclient[ase]'
pip install 'rgpot>=2.5.2'   # multi-ABI engines for Metatomic via RGPOT
```

Typed Dimer/NEB specs come from the shared **[eon-schema](project:../devdocs/eon-schema.md)**
package (`pip install 'pyeonclient[models]'` → `eon-schema>=0.2`). Job-config
pydantic models used in docs live there too (`eon.schema` is a re-export for
eon-akmc).



## Metatomic force backends

```{versionadded} 2.17
```

Named force factories live in ``pyeonclient.backends`` (``make_backend`` /
``list_backends``). Fat Metatomic, RGPOT ``libmetatomic_engine.so``, and ASE
``MetatomicCalculator`` all produce a pyeonclient ``Potential`` for the same
``.pt`` model. See [Metatomic backends](project:rgpot_metatomic.md) for the
API, packaging notes, and a PET-MAD single-point benchmark plot.


## Mutation policy (`inplace`)

Algorithms that change geometry **do not mutate caller-owned Matter by default**.
Pass ``inplace=True`` to write into the input object.

```{code-block} python
out, ok = matter.relax()                 # copy-then-relax; matter unchanged
out, ok = matter.relax(inplace=True)     # mutates matter

saddle, status = pyec.min_mode_saddle_search(
    matter, mode, E0, params, pot)       # matter unchanged
saddle, status = pyec.min_mode_saddle_search(
    matter, mode, E0, params, pot, inplace=True)
```

Same for ``run_dynamics``, ``run_monte_carlo``, ``run_basin_hopping``,
``process_search``, ``run_tad``, ``run_parallel_replica``,
``run_safe_hyperdynamics``, and ``run_replica_exchange``. There are **no**
silent mutators: every algorithm takes ``inplace=False`` by default.

```{code-block} python
out = pyec.TAD(matter, params, pot).run(matter, params, pot)                 # TADJob; matter unchanged
out = pyec.run_parallel_replica(matter, params, pot)    # ParallelReplicaJob
out = pyec.run_safe_hyperdynamics(matter, params, pot)
out = pyec.run_replica_exchange(matter, params, pot)
out = pyec.TAD(matter, params, pot).run(matter, params, pot, inplace=True)   # updates matter
```

## Shared setup

Matter-first (preferred for new code). One ``Parameters`` is shared with
``make_backend`` so pot type and job knobs stay consistent:

```{code-block} python
import numpy as np
from pyeonclient import Matter, Parameters
from pyeonclient.backends import make_backend

params = Parameters()
params.opt_converged_force = 0.01
params.opt_max_iterations = 1000
pot = make_backend("lj", params=params)   # or "metatomic", "ase_metatomic", …

matter = Matter(pot, params)
matter.resize(n)
matter.positions = ...                    # (n, 3) float64
matter.cell = np.eye(3) * L
matter.atomic_numbers = ...
matter.masses = ...
# or: from pyeonclient import from_ase; matter = from_ase(atoms, pot, params)
```

INI / PotType path (eonclient workdir parity) still works::

    params = load_parameters("config.ini")
    pot = make_potential(params.potential, params)

## Minimization

```{code-block} python
out, ok = matter.relax()                  # like LBFGS(atoms).run(fmax=...)
print(ok, matter.potential_energy, matter.max_force)
```

## Min-mode (Dimer)

```{code-block} python
params.dimer_rotations_max = 10
params.dimer_opt_method = "cg"            # cg | lbfgs | sd

direction = np.random.default_rng(0).normal(size=matter.positions.shape)
# Default method is "improved"
d = pyec.Dimer(matter, params, pot)
# d = pyec.Dimer(matter, params, pot, method="classic")
# d = pyec.Dimer(matter, params, pot, method="improved", accelerant="gp")  # WITH_GPRD
d.compute(matter, direction)
print(d.eigenvalue, d.eigenvector)
```

Low-level engines remain: ``ClassicDimer``, ``ImprovedDimer``, ``Lanczos``,
``Davidson``. Prefer ``Dimer(method=..., accelerant=...)``.

Parameters: `dimer_improved`, `dimer_rotation_angle`, `dimer_converged_angle`,
`dimer_max_iterations`, `dimer_opt_method`, `dimer_rotations_max` /
`_min`, `dimer_torque_max` / `_min`, `dimer_remove_rotation`.

## Single-ended saddle search

```{code-block} python
params.saddle_minmode_method = "dimer"    # dimer | lanczos | davidson
params.saddle_max_iterations = 1000
params.saddle_max_energy = 20.0

# Displace from the minimum, then search
mode = dimer.eigenvector                  # or a random (n, 3) guess
E_react = matter.potential_energy
ss = pyec.MinModeSaddleSearch(matter, mode, E_react, params, pot)
status = ss.run()                         # mutates matter → saddle geometry
print(pyec.saddle_status_message(status))
print(ss.eigenvalue, ss.iteration, ss.force_calls)
saddle = matter                           # now at the saddle
```

`SaddleStatus` enumerates the same codes as the C++ client (`GOOD`,
`BAD_MAX_ITERATIONS`, `DIMER_LOST_MODE`, …).

Parameters: `saddle_method`, `saddle_minmode_method`, `saddle_max_iterations`,
`saddle_max_energy`, `saddle_displace_magnitude`, `saddle_displace_type`,
`saddle_converged_force`, `saddle_remove_rotation`.

## NEB

Use **eOn-native** path initializers (IDPP / SIDPP / linear). Do not pull
ASE ``interpolate("idpp")`` into a pyeonclient band — the same engines live
in ``helpers::neb_paths`` and are already what the endpoint NEB constructor
uses.

Typed knobs live on :class:`~pyeonclient.NebSpec`. Compose the same steps as a
workdir job without a package alias:

```{code-block} python
from pyeonclient import (
    NEB,
    NebSpec,
    Parameters,
    PathInit,
    append_timing,
    from_ase,
    pot_registry_total_force_calls,
    steady_clock_now,
    write_neb_results,
    write_potcall_summary,
)
from pyeonclient.backends import make_backend

params = Parameters()
spec = NebSpec(
    n_images=10,
    path_init=PathInit.idpp,
    energy_weighted=True,
    ci_mmf=True,
    max_iterations=1000,
    force_tolerance=0.01,
)
spec.apply_to_parameters(params)
pot = make_backend("metatomic", model_path="pet-mad.pt", params=params)
initial = from_ase(reactant, pot, params)
final = from_ase(product, pot, params)

t0 = steady_clock_now()
f0 = pot_registry_total_force_calls()
neb = NEB(initial, final, params, pot, spec=spec)
status = neb.compute()
if status.name == "GOOD":
    neb.find_extrema()
write_neb_results(neb, params, pot_registry_total_force_calls() - f0)
write_potcall_summary("_potcalls.json")
append_timing("results.dat", t0)
```

Explicit path object still available::

    path = neb_idpp_path(initial, final, n_intermediate=10, parameters=params)
    band = NudgedElasticBand(path, params, pot)
    band.compute()

Helpers: ``neb_linear_path``, ``neb_idpp_path``, ``neb_idpp_collective_path``,
``neb_sidpp_path``, ``neb_initial_path``. See {doc}`neb` for energy-weighted
springs and OCINEB (``neb_ci_mmf``, …).

GP-surrogate NEB is still the same **NEB** idea with a GP accelerant on the
band engine (``WITH_GP_SURROGATE``), not a separate product noun — and **not**
prefactor / Parallel Replica.

## Hessian and prefactors

```{code-block} python
atoms = pyec.all_free_atoms(matter)         # int64 indices
hess = pyec.Hessian(params, matter)
H = hess.get_hessian(matter, atoms)       # (3n, 3n)
freqs = hess.get_freqs(matter, atoms)

pref1, pref2 = pyec.get_prefactors(params, min1, saddle, min2)
```

Parameters: `hessian_phva_atoms`, `hessian_zero_freq_value`, `prefactor_rate`
(`htst` / `qqhtst`), `prefactor_filter_scheme`, `prefactor_within_radius`, …

### Free/fixed versus PHVA active (Lanczos / Davidson)

**free/fixed** is the optimizer mask. **`phva_atoms`** (or `compute(..., atoms=)`)
is the PHVA mobile set for dense Hessian and matrix-free min-mode. Default
`All` keeps historical behavior (Krylov over all free atoms). An explicit list
bounds the Krylov dimension to \(3 N_\mathrm{active}\) without rewriting free
flags:

```{code-block} python
mobile = pyec.resolve_mobile_atoms(matter, "12,15,18")  # ∩ free
lanczos = pyec.Lanczos(matter, params, pot)
lanczos.compute(matter, direction, atoms=mobile)
# or params.lanczos_phva_atoms = "12,15,18" then compute(matter, direction)
```

Same polarity for `params.davidson_phva_atoms` / `Davidson.compute(..., atoms=)`.

## Job factory (full JobType surface)

Every type in `JobType` that the C++ `makeJob` factory implements is available
as a workdir-oriented job:

```{code-block} python
params.job = pyec.JobType.Saddle_Search     # or Process_Search, Dynamics, …
job = pyec.make_job(params)                 # reads cwd files like eonclient
files = job.run()
pyec.write_potcall_summary()
```

Use this for process search, dynamics, basin hopping, TAD, parallel replica,
Monte Carlo, GP surrogate, structure comparison, etc., until those gain
dedicated Matter-first classes. The algorithms above are already first-class.

## ASE bridge

Geometry (always array-copy, fast)::

```{code-block} python
matter = pyec.from_ase(atoms, pot, params)
atoms = pyec.to_ase(matter)
```

**Seamless calculator** — wrap a live ASE Calculator as the Matter PEF
(no file script, no ``-Dwith_ase``)::

```{code-block} python
from ase.calculators.emt import EMT

atoms.calc = EMT()
matter = pyec.from_ase(atoms)                 # uses atoms.calc
# or
pot = pyec.potential_from_ase(atoms.calc)
matter = pyec.from_ase(atoms, pot, params)
```

Requires `pip install 'pyeonclient[ase]'` and ASE installed at runtime.
Not thread-safe across images (one calculator object).

## Related

- {doc}`dimer`, {doc}`neb`, {doc}`hessian`, {doc}`prefactor`, {doc}`saddle_search`
- {doc}`/devdocs/pyeonclient-pypi` — wheels
- [Cookbook: eOn + PET-MAD NEB](https://atomistic-cookbook.org/examples/eon-pet-neb/eon-pet-neb.html)
