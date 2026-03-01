---
myst:
  html_meta:
    "description": "Guide to using the external potential (ext_pot) interface in eOn for wrapping arbitrary calculators via file-based communication."
    "keywords": "eOn external potential, ext_pot, MLIP, DeePMD, wrapper script, file interface"
---

# External Potential

```{admonition} conda-forge availability
:class: tip
Included in the `conda-forge` package. Always compiled in; no build flags required.
```

The external potential (`ext_pot`) interface allows eOn to use **any** energy and
force calculator by communicating through files and a system call. Because it
requires no compile-time dependencies, `ext_pot` is always available in every eOn
build, including the `conda-forge` package.

```{tip}
If you installed eOn from `conda-forge` and want to use an MLIP (DeePMD, MACE,
etc.) that is not available through the [Metatomic](project:metatomic_pot.md)
interface, `ext_pot` is the recommended path. The
[LAMMPS](project:lammps_pot.md) and [ASE](project:ase_pot.md) potentials require
additional compile-time flags that are **not** enabled in the `conda-forge`
build.
```

## Protocol

When the client needs an energy/force evaluation it:

1. Writes the file `from_eon_to_extpot` in the working directory.
2. Calls the executable (or script) specified by `ext_pot_path` via `system()`.
3. Reads the file `from_extpot_to_eon` that the script must create.

### Input file (`from_eon_to_extpot`)

The first three lines are the 3x3 box matrix (one row per line, tab-separated).
Each subsequent line contains one atom: `atomic_number  x  y  z`
(tab-separated, double precision).

```
box_00  box_01  box_02
box_10  box_11  box_12
box_20  box_21  box_22
6   1.234   5.678   9.012
6   ...
```

### Output file (`from_extpot_to_eon`)

The first line is the total energy (scalar).
Each subsequent line contains the force on one atom: `fx  fy  fz`
(space-separated). The atom order must match the input.

```
-42.12345
0.123  -0.456  0.789
...
```

## Configuration

```{code-block} ini
[Potential]
potential = ext_pot
ext_pot_path = /full/path/to/your_wrapper
```

```{important}
When running compound jobs (aKMC, process search, etc.) eOn copies
configuration into per-job scratch directories. Use **absolute paths** for
`ext_pot_path` and for any model files referenced inside your wrapper.
```

## Examples

### DeePMD (PyTorch) wrapper

The following Python script wraps a DeePMD-kit v3 PyTorch model and can be used
directly as the `ext_pot_path` target. Save it as, e.g.,
`/home/user/scripts/deepmd_extpot` and make it executable (`chmod +x`).

```{code-block} python
:caption: deepmd_extpot

#!/usr/bin/env python
"""eOn ext_pot wrapper for DeePMD-kit (PyTorch backend)."""
import numpy as np
from deepmd.infer import DeepPot

# --- user settings ---
MODEL = "/absolute/path/to/model.pth"
# ----------------------

dp = DeepPot(MODEL)

# Read input
lines = open("from_eon_to_extpot").readlines()
box = np.array([[float(v) for v in l.split()] for l in lines[:3]])
atoms = []
coords = []
for l in lines[3:]:
    tok = l.split()
    atoms.append(int(tok[0]))
    coords.append([float(tok[1]), float(tok[2]), float(tok[3])])
atoms = np.array(atoms)
coords = np.array(coords)

# Evaluate
energy, forces, _ = dp.eval(coords.reshape(1, -1, 3),
                             box.reshape(1, 9),
                             atoms)

# Write output
with open("from_extpot_to_eon", "w") as f:
    f.write(f"{energy[0]:.15f}\n")
    for fx, fy, fz in forces[0]:
        f.write(f"{fx:.15f} {fy:.15f} {fz:.15f}\n")
```

### Generic ASE calculator wrapper

Any ASE-compatible calculator can be wrapped in a similar fashion:

```{code-block} python
:caption: ase_extpot

#!/usr/bin/env python
"""eOn ext_pot wrapper using an ASE calculator."""
import numpy as np
from ase import Atoms

# --- user settings ---
from mace_mp import MACECalculator
calc = MACECalculator(model="/path/to/model.pt", device="cuda")
# ----------------------

lines = open("from_eon_to_extpot").readlines()
box = np.array([[float(v) for v in l.split()] for l in lines[:3]])
numbers = []
positions = []
for l in lines[3:]:
    tok = l.split()
    numbers.append(int(tok[0]))
    positions.append([float(tok[1]), float(tok[2]), float(tok[3])])

system = Atoms(numbers=numbers, positions=positions, cell=box, pbc=True)
system.calc = calc

energy = system.get_potential_energy()
forces = system.get_forces()

with open("from_extpot_to_eon", "w") as f:
    f.write(f"{energy:.15f}\n")
    for fx, fy, fz in forces:
        f.write(f"{fx:.15f} {fy:.15f} {fz:.15f}\n")
```

## Performance considerations

Each force call spawns a new process and re-initializes the calculator. For
GPU-accelerated MLIPs this overhead is typically small compared to the inference
cost, but for very fast potentials the startup cost can dominate. In that
scenario consider:

- The built-in [ASE interface](project:ase_pot.md) (requires building from
  source with `-Dwith_python=True -Dwith_ase=True`), which embeds the
  interpreter and avoids per-call process overhead.
- The [Metatomic](project:metatomic_pot.md) interface for models that support
  the metatensor format (compiled into the `conda-forge` build).
- The [serve mode](project:serve_mode.md) to expose any eOn potential over RPC
  to external callers.
