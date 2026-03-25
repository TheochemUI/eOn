---
myst:
  html_meta:
    "description": "Tutorial on using Python dictionaries to generate eOn configuration files programmatically, as an alternative to hand-written config.ini files."
    "keywords": "eOn dictionary config, Python configuration, write_eon_config, rgpycrumbs, programmatic setup"
---

# Dictionary-Style Configuration

eOn reads its parameters from a `config.ini` file in INI format. While writing
these files by hand works well for one-off runs, programmatic workflows benefit
from generating the configuration from Python dictionaries. The
{func}`~rgpycrumbs.eon.helpers.write_eon_config` function in
[rgpycrumbs](https://github.com/HaoZeke/rgpycrumbs) handles this conversion.

All configuration options are documented in {mod}`eon.schema`, which is the
single source of truth for parameter names, types, defaults, and allowed values.
The user guide pages (e.g. {doc}`/user_guide/main`, {doc}`/user_guide/neb`)
render these Pydantic models directly.

## Installation

`rgpycrumbs` is available on PyPI:

```bash
pip install rgpycrumbs
# or
uv add rgpycrumbs
```

## Basic Usage

A dictionary maps INI section names to their key-value pairs. The section names
match the headers in `config.ini` (and the headings in the user guide). For
example, an AKMC run on a Pt heptamer:

```python
from pathlib import Path
from rgpycrumbs.eon.helpers import write_eon_config

settings = {
    "Main": {
        "job": "akmc",
        "temperature": 300,
    },
    "Potential": {
        "potential": "morse_pt",
    },
    "Optimizer": {
        "converged_force": 0.001,
        "max_iterations": 1000,
    },
    "AKMC": {
        "confidence": 0.95,
    },
    "Process Search": {
        "minimize_first": True,
    },
    "Communicator": {
        "type": "local",
        "number_of_CPUs": 2,
        "num_jobs": 2,
    },
    "Saddle Search": {
        "displace_least_coordinated_weight": 1.0,
        "displace_radius": 3.3,
        "displace_magnitude": 0.1,
        "min_mode_method": "dimer",
        "max_energy": 10.0,
    },
}

write_eon_config(Path("."), settings)
```

This writes a `config.ini` in the current directory. You can then run
`python -m eon.server` (or `eonclient` directly) as usual.

```{tip}
Pass a directory path and `write_eon_config` creates `config.ini` inside it.
Pass a full file path to control the output name.
```

## Why Dictionaries?

Compared to editing INI files by hand, the dictionary approach provides:

- **Reproducibility**: the script *is* the configuration, checked into version
  control alongside the structure files.
- **Parameterization**: loop over temperatures, spring constants, or image
  counts without duplicating INI files.
- **Validation reference**: option names and types are defined in
  {mod}`eon.schema`. Typos that would silently fall back to defaults in an INI
  file become obvious when compared against the schema documentation.
- **Notebook integration**: generate and run eOn configurations within Jupyter
  notebooks. The
  [atomistic cookbook](https://atomistic-cookbook.org/examples/eon-pet-neb/eon-pet-neb.html)
  has a worked NEB example using this approach.

## Examples by Job Type

Each example directory under `examples/` now contains both a `config.ini` and
a Python equivalent (`run_*.py`). The scripts are self-contained and can be
executed with:

```bash
cd examples/akmc-pt
python run_akmc_pt.py      # generates config.ini
python -m eon.server       # runs the simulation
```

### AKMC

:::::{tab-set}

::::{tab-item} Dictionary (Python)
```{literalinclude} ../../../examples/akmc-pt/run_akmc_pt.py
:language: python
```
::::

::::{tab-item} INI
```{literalinclude} ../../../examples/akmc-pt/config.ini
:language: ini
```
::::

:::::

### NEB

:::::{tab-set}

::::{tab-item} Dictionary (Python)
```{literalinclude} ../../../examples/neb-al/run_neb_al.py
:language: python
```
::::

::::{tab-item} INI
```{literalinclude} ../../../examples/neb-al/config.ini
:language: ini
```
::::

:::::

For advanced NEB options (climbing image, energy-weighted springs, IDPP
initialization, off-path CI with MMF), see
`examples/neb-al/run_neb_advanced.py`.

### Basin Hopping

:::::{tab-set}

::::{tab-item} Dictionary (Python)
```{literalinclude} ../../../examples/basin-hopping/run_basin_hopping.py
:language: python
```
::::

::::{tab-item} INI
```{literalinclude} ../../../examples/basin-hopping/config.ini
:language: ini
```
::::

:::::

### Parallel Replica Dynamics

:::::{tab-set}

::::{tab-item} Dictionary (Python)
```{literalinclude} ../../../examples/parallel-replica/run_parrep.py
:language: python
```
::::

::::{tab-item} INI
```{literalinclude} ../../../examples/parallel-replica/config.ini
:language: ini
```
::::

:::::

### AKMC with Displacement Script

The Cu vacancy example shows how dictionary config works alongside a
displacement script (`ptmdisp.py`). The script path is just another string
parameter:

```python
"Saddle Search": {
    "displace_atom_kmc_state_script": "ptmdisp.py",
    "displace_all_listed": True,
    # ...
},
```

See `examples/akmc-cu-vacancy/run_akmc_cu.py` for the full configuration and
{doc}`displacement_scripts` for details on writing displacement scripts.

## Parameter Sweeps

The dictionary approach makes parameter sweeps straightforward:

```python
from pathlib import Path
from rgpycrumbs.eon.helpers import write_eon_config

base = {
    "Main": {"job": "nudged_elastic_band"},
    "Potential": {"potential": "eam_al"},
    "Optimizer": {
        "opt_method": "lbfgs",
        "max_move": 0.1,
        "converged_force": 0.001,
        "max_iterations": 1000,
    },
}

for n_images in [5, 7, 11, 15]:
    run_dir = Path(f"neb_{n_images}img")
    run_dir.mkdir(exist_ok=True)
    settings = {
        **base,
        "Nudged Elastic Band": {
            "images": n_images,
            "spring": 5.0,
        },
    }
    write_eon_config(run_dir, settings)
```

## Schema Reference

The authoritative documentation for every configuration option lives in the
Pydantic models in {mod}`eon.schema`. The user guide pages render these models
automatically:

- {doc}`/user_guide/main` -- general simulation parameters
- {doc}`/user_guide/akmc` -- adaptive kinetic Monte Carlo
- {doc}`/user_guide/neb` -- nudged elastic band
- {doc}`/user_guide/saddle_search` -- saddle search methods
- {doc}`/user_guide/optimizer` -- optimization algorithms
- {doc}`/user_guide/potential` -- interatomic potentials
- {doc}`/user_guide/dynamics` -- molecular dynamics
- {doc}`/user_guide/parallel_replica` -- parallel replica dynamics
- {doc}`/user_guide/basin_hopping` -- basin hopping
- {doc}`/user_guide/communicator` -- job communicators
