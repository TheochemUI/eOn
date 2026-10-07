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
interface, `ext_pot` or [LAMMPS](project:lammps_pot.md) is the recommended path.
The [ASE](project:ase_pot.md) potentials require additional compile-time flags
that are **not** enabled in the `conda-forge` build.
```

## CPMD method file

A CPMD script behind `ext_pot` should read one job TOML, the same file the
[rgpot page](project:rgpot_pot.md) compiles to a CPMDParams message. The
file route renders that message. The script does not keep a second cutoff.

## Protocol

When the client needs an energy/force evaluation it:

1. Creates a private exchange directory named `extpot_<pid>_<n>` inside the
   directory eOn runs in, where `<pid>` is the client's process id and `<n>`
   counts the potentials it built.
2. Writes the file `from_eon_to_extpot` in that directory.
3. Calls the executable (or script) specified by `ext_pot_path` via
   `system()`, with the exchange directory as the command's working
   directory.
4. Reads the file `from_extpot_to_eon` that the script must create there.

The two filenames do not change, so a wrapper that opens them by relative
name, as every example below does, needs no edit. The exchange directory is
what keeps two clients started in one directory from reading each other's
numbers. The directory lives as long as the potential: each call reuses it,
so a wrapper can keep state there between calls (the CPMD wrapper below
keeps its wavefunction). When the potential is destroyed the client
removes the two exchange files and then the directory, if nothing else is
left in it. A failed call leaves everything in place for inspection.

## One directory per image

`[Main] parallel` defaults to true. A dimer, and the two product
minimizations in a process search, then build a second potential when the
potential asks for one instance per image. `ext_pot` does. Each instance
takes the next `extpot_<pid>_<n>`, so the two calls do not share
`from_eon_to_extpot`. An error on the second thread is stored, the first
thread is joined, and the stored error is raised. `parallel = false` runs
both images on the one instance, one after the other.

```{important}
`ext_pot_path` is resolved against the directory eOn runs in when it names an
existing file, so `./ext_pot` and an absolute path both work. A command line
with arguments (`python wrapper.py`) or a name looked up on `PATH` is passed
to the shell unchanged and therefore has to be written so it resolves from
the exchange directory: give the script an absolute path. On Windows,
`cmd.exe` does not honor a shebang: a resolved script whose first line names
python, or whose name ends in `.py`, is invoked as `python <script>`. A
wrapper that needs to reach files in the directory eOn runs in finds that
directory in the environment variable `EON_EXTPOT_RUN_DIR`.
```

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

### CPMD on several MPI ranks

`examples/akmc-cpmd-slurm/potfiles/cpmd_extpot.py` answers each call with a
`cpmd.x` wavefunction optimisation on several MPI ranks. It reads its
settings from the environment:

| Variable | Meaning | Default |
| --- | --- | --- |
| `CPMD_BIN` | path to `cpmd.x` | required |
| `CPMD_LAUNCH` | MPI launcher put in front of `cpmd.x` | `mpirun -np 4` |
| `PP_LIBRARY_PATH` | pseudopotential directory, passed to `cpmd.x` | required |
| `CPMD_PP` | `Z:file[:options]` per element, comma separated, e.g. `14:Si_MT_BLYP.psp:LMAX=D,7:N_MT_BLYP.psp:LMAX=P` | required |
| `CPMD_CUTOFF` | plane-wave cutoff in Ry | `30` |
| `CPMD_FUNCTIONAL` | `&DFT FUNCTIONAL` | `BLYP` |
| `CPMD_CONVERGENCE` | `CONVERGENCE ORBITALS` | `1e-6` |
| `CPMD_MAXITER` | `MAXITER` | `300` |
| `CPMD_RESTART_POOL` | shared directory of converged `RESTART.1` files that seed a job's first call | unset |

A single point, a minimization or a nudged elastic band on one machine
needs nothing else than the script as `ext_pot_path`:

```{code-block} bash
export CPMD_BIN=/path/to/cpmd.x PP_LIBRARY_PATH=/path/to/pseudopotentials/
export CPMD_LAUNCH="mpirun -np 4 --bind-to none" CPMD_CUTOFF=70
export CPMD_PP="14:Si_MT_BLYP.psp:LMAX=D,7:N_MT_BLYP.psp:LMAX=P"
cat > config.ini <<'EOF'
[Main]
job = point

[Potential]
potential = ext_pot
ext_pot_path = /absolute/path/to/eOn/examples/akmc-cpmd-slurm/potfiles/cpmd_extpot.py
EOF
eonclient
```

Inside a Slurm job, ask for as many tasks as `CPMD_LAUNCH` starts ranks
(`--ntasks=4`), or OpenMPI finds one slot and refuses `-np 4`.
`--bind-to none` keeps it from refusing to bind when the allocation is
hardware threads rather than whole cores.

`examples/akmc-cpmd-slurm` runs adaptive kinetic Monte Carlo on a silicon
vacancy with the same script. The `cluster` communicator submits each
process search as a Slurm job (see [Communicator](project:communicator.md)),
and inside the job `ext_pot` calls `cpmd_extpot.py`. The client stays one
process; the ranks belong to CPMD alone.

What the script hands CPMD and reads back:

- The deck says `ANGSTROM`, so CPMD converts eOn's positions with its own
  Bohr (0.529177210859 A). CPMD answers in Hartree and Hartree/Bohr, and the
  script converts with CODATA 2018 (27.211386245988 eV, 0.529177210903 A),
  as the in-process [RGPOT](project:rgpot_pot.md) route does. The two Bohr
  values differ by 8.3e-11 relative; `client/validation/cpmd_units.py`
  derives the chain.
- The table under `GRADIENTS (-FORCES)` holds forces, not gradients:
  `wrgeo` prints the ionic force `fion`. The first table, printed before the
  wavefunction optimisation, is all zeros; the last one is the result.
- That table carries three significant digits. `GEOMETRY` holds the same
  forces in Ha/Bohr to twelve decimals in columns 4 to 6, so the script
  reads forces from `GEOMETRY` and checks them against the table.
- The energy is the last `TOTAL ENERGY` of `cpmd.out`, which CPMD prints
  with eight decimals: the file route's energy is good to 5e-9 Ha
  (1.4e-7 eV). The in-process route returns it to full precision.
- CPMD lists atoms species by species, in the order each element first
  appears. The script maps them back to eOn's order.
- eOn's cell goes to `CELL VECTORS` unchanged, one lattice vector per row.
  CPMD refuses `SYMMETRY` beside it, so the deck omits `SYMMETRY`. A
  system without a cell reaches the script as a zero box, which it
  refuses: CPMD needs a periodic cell.

The script refuses a result whose wavefunction optimisation printed
`BUT NO CONVERGENCE`, after one retry from an atomic guess with
`PCG MINIMIZE`. It keeps `RESTART.1` in the exchange directory, so every call
after the first in a job restarts from the previous wavefunction.

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
