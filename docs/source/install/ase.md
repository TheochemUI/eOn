# ASE Interface

In order to build eON client with the ability to use ASE calculators as a
potential, the appropriate option must be setup when using `meson` which is
`-Dwith_python=True -Dwith_ase=True`. 

Some notes about this implementation:

- It is registered as an external potential
- The `conda` environment must have `ase` present
- Like other `python` potentials, it shares the `guard` across the client

## Configuration

```{code-block} ini
[Potential]
potential = ase
ext_pot_path = /full/path/to/the/ase_script.py
```

The last line (`ext_pot_path`) is required, and it may be the relative or the
absolute path to the Python script containing the desired ASE calculator (see
next section).

When running EON (e.g. AKMC jobs), which calls `eonclient` in the back, **it is
highly recommended to write the full absolute path to the script**.

### Scripting `ASE`

EON client imports this Python script to get the energy and forces.
The script should look like this (e.g. using ASE's Lennard Jones calculator)::

```{code-block} python
from ase import Atoms
from ase.calculators.lj import LennardJones

def ase_calc():
    # MODIFY THIS SECTION
    calc = LennardJones(epsilon=0.0103, sigma=3.40, rc=10.0, ro=0.0, smooth=True)
    return calc

#=======================================================================
# DO NOT EDIT
def _calculate(R, atomicNrs, box, calc):
    system = Atoms(symbols=atomicNrs, positions=R, pbc=True, cell=box)
    system.calc = calc
    forces = system.get_forces()
    energy = system.get_potential_energy()
    return energy, forces
#=======================================================================


if __name__ == '__main__':
    # This is just to verify that the script runs properly with `python3 <this_script>`
    from ase.io import read
    atoms = read("...")  # structure file name
    pos = atoms.positions  # array of positions (n x 3)
    atomicNrs = atoms.get_atomic_numbers()  # array of atomic numbers (n)
    cell = atoms.cell  # cell array (3 x 3)
    calc = ase_calc()
    e, f = _calculate(pos, atomicNrs, cell, calc)
    print(f"E = {e}")
    print("F =")
    print(f)
```

The user should customize the content of the function `ase_calc()` to use the
desired ASE calculator, **but not that of `_calculate()` nor the function
names**.

Often, the ASE calculator may require an external potential file. When running
EON (e.g. AKMC jobs) and an external file name is used in this script, **it is
highly recommended to write the full absolute path to the external file**.

Before using it, make sure that this script contains no errors on the Python
side.
