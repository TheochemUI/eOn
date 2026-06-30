---
myst:
  html_meta:
    "description": "Learn how to set up and use the Metatomic interface in eOn for machine learning potentials with metatensor and PyTorch."
    "keywords": "eOn metatomic, machine learning potential, metatensor, PyTorch, MLIP"
---

# Metatomic Interface

```{versionadded} 2.0
```

```{admonition} conda-forge availability
:class: tip
Included in the `conda-forge` package. No additional build flags required.
```

The Metatomic interface allows eOn to use machine learning potentials developed
with the `metatensor` and `pytorch` libraries.

## Setup and Compilation

The most robust way to handle the dependencies is to use `-Dwith_metatomic=True
-Dpip_metatomic=True -Dtorch_version=2.10` inside a `pixi` environment with the
`metatomic` feature (for example `pixi run -e dev-mta` / `setupeon` + `mkeon`).
That path expects **metatomic-torch >= 0.1.15** and **metatensor-torch >= 0.10**
(pip layout `metatensor_torch/torch-<ver>/` for the shared library).

Some notes about this implementation:

- Unlike the `python` potentials, using `pip_metatomic` does not share a `guard`
  across the client
- These models are typically stored in `.pt` files.
- Prefer `device = cuda` (or `auto`) on a GPU host; CPU remains the default when
  unset in the INI reader legacy path.

For more details including building from source, refer to the [upstream documentation](https://docs.metatensor.org/latest/index.html).

## Basic Configuration

```{code-block} ini
[Potential]
potential = Metatomic

[Metatomic]
model_path = lennard-jones.pt
device = cpu
deterministic = true
```

The `model_path` is the path to the PyTorch model file (`.pt`). As with other
external files, it is highly recommended to provide the full absolute path to
the model. Configuration keys are defined in `eon/schema.py` (`Metatomic`
Pydantic model) and mirrored in `eon/config.yaml` for the YAML front end;
`config.ini` / `[Metatomic]` INI keys are read in `ParametersINI.cpp`.

## Usage example

### Lennard Jones Baseline

Here is a complete workflow, from generating a simple Lennard-Jones test model
to running a calculation with eonclient.

The metatomic-lj-test package provides a simple way to create a sample model
file. Run this Python script:

```{code-block} python
import metatomic_lj_test

model = metatomic_lj_test.lennard_jones_model(
    atomic_type=1,
    cutoff=3.4,
    sigma=1.5,
    epsilon=23.0,
    length_unit="Angstrom",
    energy_unit="eV",
    with_extension=False,
)

model.save("lennard-jones.pt", collect_extensions="extensions/")

print("Saved model to lennard-jones.pt")

```

This will create the `lennard-jones.pt` file in your current directory.

<!-- This needs to be re-exported every time the Torch version changes. -->


This can be loaded checked with a configuration file for eOn:

```{code-block} ini
[Main]
job = point
temperature = 300
random_seed = 706253457

[Potential]
potential = metatomic

[Metatomic]
model_path = lennard-jones.pt
```

With coordinates from the `lj13.con` file:

```{code-block} bash
cp docs/lj13.con pos.con
eonclient
# Energy:         98374.877530582715
# Max atom force: 109591.908376179999
```

Validated with the `metatomic` wrapper.

```{code-block} python
import numpy as np
import ase.io as aseio
from metatomic.torch import load_atomistic_model
from metatomic_ase import MetatomicCalculator
atomistic_model = load_atomistic_model("lennard-jones.pt")
mta_calculator = MetatomicCalculator(atomistic_model)
atoms = aseio.read("pos.con")
atoms.calc = mta_calculator
# NOTE(rg): Normally needs a mask to remove fixed atoms, lj13 has no fixed atoms
print(atoms.get_potential_energy())
print(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
```

Which yields the expected result.

```{code-block} python
In [2]: atoms.get_potential_energy()
Out[2]: np.float64(98374.87753058573)
In [6]: np.max(np.linalg.norm(atoms.get_forces(), axis=1))
Out[6]: np.float64(109591.90837618345)
```


### PET-MAD

To use `pet-mad` we can use `metatrain` to grab the model.

```{code-block} sh
mtt export https://huggingface.co/lab-cosmo/pet-mad/resolve/v1.1.0/models/pet-mad-v1.1.0.ckpt
```

## Variance

```{versionadded} 2.2
```

To enable per-atom energy uncertainty checks, set the `uncertainty_threshold`
parameter to a positive value. This parameter serves simultaneously as the
activation switch and the warning threshold.

This functionality supports models which expose an `energy_uncertainty` [output
key](https://docs.metatensor.org/metatomic/latest/outputs/energy.html#energy-uncertainty).
If the configuration specifies a `variant_base` or `variant_energy_uncertainty`,
eOn will automatically target the corresponding variant key (e.g.,
`energy_uncertainty/ensemble`).

````{margin}
```{note}

* These are **not** force uncertainities and **may not** be correlated to accuracy.

```
````

## Variants

```{versionadded} 2.9.0
```

Metatomic models frequently act as multi-headed neural networks, capable of
predicting properties corresponding to different levels of theory (e.g.,
`energy/pbe0` versus `energy/r2scan`) or auxiliary outputs within a single model
file. The eOn interface permits precise selection of these output heads through
the `[Metatomic]` configuration block.

The implementation follows the [upstream Metatomic
specification](https://docs.metatensor.org/metatomic/latest/outputs/variants.html#output-variants)
where variants appear as suffixes to the base quantity (e.g.,
`energy/<variant>`), via `metatomic_torch::pick_output`.

### Variant configuration

Four keys control variant selection:

- **`variant_base`**: Default variant suffix for outputs that do not set an
  override. For example, setting this to `pbe` targets `energy/pbe` and, when
  uncertainty is enabled, `energy_uncertainty/pbe`.
- **`variant_energy`**: Overrides the energy selection specifically. Setting
  this takes precedence over `variant_base`. Use the special value `off` to
  revert to the default `energy` output even if `variant_base` remains active
  for other quantities.
- **`variant_energy_uncertainty`**: Overrides the uncertainty selection
  specifically, functioning identically to `variant_energy`.
- **`variant_force`**: Override for the non-conservative force head (see
  below). Defaults to the energy variant when empty.

### Explicit output keys

```{versionadded} 2.15
```

When a model uses non-standard output names (for example a custom energy head),
set the literal keys instead of relying on `pick_output`:

- **`energy_output`**: Literal energy key in the model capabilities.
- **`energy_uncertainty_output`**: Literal uncertainty key (used when
  `uncertainty_threshold` is positive).
- **`force_output`**: Literal non-conservative force key (used when
  `non_conservative = true`). Missing explicit keys raise at potential
  construction.

### Example

To select a specific functional (e.g., PBE0) for both energy and uncertainty:

```{code-block} ini
[Metatomic]
model_path = universal-potential.pt
variant_base = pbe0

```

To use a baseline potential for energy but a specific variance head for
uncertainty quantification:

```{code-block} ini
[Metatomic]
model_path = active-learning.pt
# Uses 'energy' (default)
variant_energy = off
# Uses 'energy_uncertainty/ensemble'
variant_energy_uncertainty = ensemble
```

## Non-conservative forces

```{versionadded} 2.16
```

By default, forces are the negative gradient of the requested energy output
(conservative). Set **`non_conservative = true`** to read forces from the
model's `non_conservative_force` (or variant) output instead, matching
`metatomic_ase.MetatomicCalculator(non_conservative=...)`.

- **`non_conservative`**: `false` (default) uses autograd on energy; `true`
  requests the NC force output.
- **`variant_force`**: Variant suffix for `non_conservative_force`
  (defaults to the energy variant).
- **`force_output`**: Literal key override for the NC force output.

```{code-block} ini
[Metatomic]
model_path = pet-mad.pt
non_conservative = true
variant_force =   # empty -> same as energy variant / variant_base
```

Use non-conservative forces for speed only when the model documents them; they
are not guaranteed to be energy-conserving in MD or long NEB runs.

## Random rotations and symmetrization

```{versionadded} 2.16
```

Broken rotational symmetry in some ML potentials can be mitigated by evaluating
structures in a randomly rotated frame and rotating forces back (see Langer,
Pozdnyakov, Ceriotti, *Mach. Learn.: Sci. Technol.* **5**, 4LT01 (2024)).

- **`random_rotation`**: If `true`, each `force()` call draws a uniform SO(3)
  rotation, evaluates the model on rotated positions (and cell), and maps
  forces back to the lab frame.
- **`n_symmetry_rotations`**: If greater than zero, average energy and forces
  over that many independent random rotations (approximate O(3)
  symmetrization). When set, it supersedes single-shot `random_rotation`.

```{code-block} ini
[Metatomic]
model_path = pet-mad.pt
random_rotation = true
# or, for averaging:
# n_symmetry_rotations = 8
```

These options increase cost proportionally to the number of rotations and are
orthogonal to non-conservative forces (NC forces are rotated the same way).

## Determinism (PyTorch)

```{versionadded} 2.16
```

Parallel NEB and multi-image workloads need bit-stable forces across model
instances. eOn applies the same mitigations discussed in
[PyTorch deterministic regression](https://rgoswami.me/snippets/pytorch-deterministic-regression/)
when determinism is enabled:

- Disable TorchScript profiling-guided optimization (`getProfilingMode() = false`).
- Enable `at::globalContext().setDeterministicAlgorithms(...)`.
- Disable cuDNN benchmarking.

Keys:

- **`deterministic`**: `true` (default) applies the above. Set `false` only for
  exploratory runs where speed matters more than reproducibility across
  images or processes.
- **`deterministic_strict`**: When `true`, nondeterministic CUDA ops fail
  instead of warning. Requires `CUBLAS_WORKSPACE_CONFIG` in the environment
  (for example `:4096:8`). If that variable is set, strict mode is also
  inferred even when `deterministic_strict` is left `false`.

```{code-block} ini
[Metatomic]
model_path = pet-mad.pt
deterministic = true
deterministic_strict = false
```

```{code-block} bash
export CUBLAS_WORKSPACE_CONFIG=:4096:8
# then set deterministic_strict = true for hard failures on nondeterministic ops
```

## Device selection

Set **`device`** to `cpu`, `cuda`, `mps`, or leave empty / `auto` to let
`metatomic_torch::pick_device` choose from the model's
`supported_devices`. GPU builds must link a CUDA-enabled libtorch; the
cookbook and conda packages typically ship a CPU libtorch with optional CUDA
at runtime when the host provides it.

```{code-block} ini
[Metatomic]
model_path = /abs/path/pet-mad.pt
device = cuda
deterministic = true
```
