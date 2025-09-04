# Metatomic Interface

```{versionadded} 2.0
```

The Metatomic interface allows EON to use machine learning potentials developed
with the `metatensor` and `pytorch` libraries.

## Setup and Compilation

The most robust way to handle the dependencies is to use `-Dwith_metatomic=True
-Dpip_metatomic=True -Dtorch_version=2.7` inside the `pixi s -e meta-dev`
environment.

Some notes about this implemetnation:

- Unlike the `python` potentials, using `pip_metatomic` does not share a `guard`
  across the client
- These models are typically stored in `.pt` files.

For more details including building from source, refer to the [upstream documentation](https://docs.metatensor.org/latest/index.html).

## Configuration

```{code-block} ini
[Potential]
potential = metatomic

[Metatomic]
model_path = lennard-jones.pt
```

The model_path is the path to the PyTorch model file (`.pt`). As with other
external files, it is highly recommended to provide the full absolute path to
the model.

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

This will create the lennard-jones.pt file in your current directory.


This can be loaded checked with a configuration file for EON:

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
from metatomic.torch.ase_calculator import MetatomicCalculator
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
