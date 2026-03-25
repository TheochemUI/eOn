---
myst:
  html_meta:
    "description": "Guide to using the Atomic Simulation Environment (ASE) interface in eOn, allowing the use of any ASE calculator as a potential."
    "keywords": "eOn ASE, Atomic Simulation Environment, ASE calculator, Python potential"
---

# ASE Interface

```{admonition} conda-forge availability
:class: warning
**Not** included in the `conda-forge` package. Requires building from source
with `-Dwith_ase=True`. For a simpler alternative that works with any eOn
install, see the [external potential](project:ext_pot.md) guide.
```

The ASE interface embeds a Python interpreter inside the eOn client, letting you
use **any** [ASE calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html)
as a potential without per-call process overhead.

## Building from source

### Prerequisites

A conda (or mamba/pixi) environment with at minimum:

```{code-block} bash
conda create -n eon python numpy ase pybind11 eigen spdlog fmt meson ninja pkg-config compilers
conda activate eon
```

Install any additional calculator packages you need (e.g. `mace-torch`,
`chgnet`, etc.) into the same environment.

### Build

```{code-block} bash
git clone https://github.com/TheochemUI/eOn.git
cd eOn
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype=release -Dwith_ase=True
meson install -C bbdir
```

```{note}
The `-Dwith_python=True` flag is no longer needed; `-Dwith_ase=True` implies it.
The Python interpreter is started lazily (only when the ASE potential is
actually used), so builds with this flag carry zero overhead for non-ASE runs.
```

## Configuration

```{code-block} ini
[Potential]
potential = ase_pot
ext_pot_path = /full/path/to/the/ase_script.py
```

The `ext_pot_path` value is the path to a Python script that defines the
calculator (see next section). When running compound jobs (aKMC, process search,
etc.), **always use an absolute path** because eOn copies configuration into
per-job scratch directories.

## Writing the calculator script

eOn imports this script at startup and calls two functions:

- `ase_calc()` -- returns an initialized ASE calculator object.
- `_calculate(R, atomicNrs, box, calc)` -- evaluates energy and forces.

Only `ase_calc()` should be edited. The `_calculate` function is a fixed
bridge between eOn and ASE; do not modify it.

### Template

```{code-block} python
:caption: my_calculator.py

from ase import Atoms
from ase.calculators.lj import LennardJones

def ase_calc():
    # --- customize this section ---
    calc = LennardJones(epsilon=0.0103, sigma=3.40, rc=10.0, ro=0.0, smooth=True)
    return calc

#=======================================================================
# DO NOT EDIT below this line
def _calculate(R, atomicNrs, box, calc):
    system = Atoms(symbols=atomicNrs, positions=R, pbc=True, cell=box)
    system.calc = calc
    forces = system.get_forces()
    energy = system.get_potential_energy()
    return energy, forces
#=======================================================================


if __name__ == '__main__':
    # Quick sanity check: python3 my_calculator.py
    from ase.io import read
    atoms = read("pos.con")       # or any structure file
    pos = atoms.positions
    atomicNrs = atoms.get_atomic_numbers()
    cell = atoms.cell
    calc = ase_calc()
    e, f = _calculate(pos, atomicNrs, cell, calc)
    print(f"E = {e}")
    print("F =")
    print(f)
```

### MACE example

```{code-block} python
def ase_calc():
    from mace.calculators import MACECalculator
    return MACECalculator(model_paths="/absolute/path/to/model.pt", device="cuda")
```

## Alternatives

If you installed eOn from conda-forge or prefer not to build from source, the
[external potential](project:ext_pot.md) interface provides the same calculator
flexibility through file-based communication. It works out of the box with any
eOn install at the cost of per-call process overhead.
