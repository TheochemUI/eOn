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

## Start with a complete working example

This example needs Python and the eOn client only. It evaluates a harmonic
potential, so it needs no electronic-structure executable, model download or
pseudopotentials. The commands below target Linux, including Ubuntu under WSL.
Run calculations inside your site's compute allocation, not on its login node.

Install [Pixi](https://pixi.sh), then create a separate workspace:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"
mkdir eon-external-demo
cd eon-external-demo
pixi init --channel conda-forge
pixi add eon python
pixi run eonclient --version
```

`pixi add` installs the client and Python and records their resolved packages
in `pixi.lock`. Keep `pixi.toml` and `pixi.lock` together; a recipient can
recreate the environment with `pixi install --locked`. No compiler or MPI
installation is needed for this file-interface example.

Save the following executable as `harmonic_extpot.py`:

```python
#!/usr/bin/env python3
"""Harmonic energy in eV and its negative derivative in eV/Angstrom."""
import math
from pathlib import Path

rows = [line.split() for line in
        Path("from_eon_to_extpot").read_text().splitlines() if line.strip()]
if len(rows) < 4 or any(len(row) != 3 for row in rows[:3]):
    raise ValueError("Expected three box rows followed by at least one atom")
positions = []
for row in rows[3:]:
    if len(row) != 4:
        raise ValueError("Expected atomic_number x y z")
    int(row[0])
    xyz = [float(value) for value in row[1:]]
    if not all(math.isfinite(value) for value in xyz):
        raise ValueError("Non-finite coordinate")
    positions.append(xyz)
k = 1.0  # eV/Angstrom squared
energy = 0.5 * k * sum(x*x for xyz in positions for x in xyz)
forces = [[-k*x for x in xyz] for xyz in positions]
with Path("from_extpot_to_eon").open("w") as output:
    output.write(f"{energy:.15g}\n")
    for force in forces:
        output.write(" ".join(f"{value:.15g}" for value in force) + "\n")
```

Create `pos.con`, with one silicon atom at (1, 0, 0) Angstrom in a 20 Angstrom
cube. The two blank lines between the cell angles and component count are
part of this structure format:

```text
One atom harmonic example
0.0
20.0 20.0 20.0
90.0 90.0 90.0


1
1
28.0855
Si
Coordinates of Component 1
1.0 0.0 0.0 0 0
```

Make the wrapper executable, create the configuration and evaluate the point:

```bash
chmod +x harmonic_extpot.py
cat > config.ini <<EOF
[Main]
job = point

[Potential]
potential = ext_pot
ext_pot_path = $PWD/harmonic_extpot.py
EOF
pixi run eonclient > eon.log 2>&1
cat results.dat
```

The `potential_energy` result is **0.5 eV**. The wrapper returns one force,
**(-1, 0, 0) eV/Angstrom**. This potential is a protocol demonstration:
its origin is fixed and it does not model periodic interatomic interactions.
The exact result-file formatting depends on the installed eOn release.

For example, the point check prints these result entries:

```text
GOOD termination_reason_text
EXT_POT potential_type
1 total_force_calls
5.000000000000e-01 potential_energy
1.000000000000e+00 Max_Force
```

To check the wrapper independently, run it from a separate directory
containing this `from_eon_to_extpot` file:

```text
20 0 0
0 20 0
0 0 20
14 1 0 0
```

Its `from_extpot_to_eon` must contain:

```text
0.5
-1 -0 -0
```

Signed zero is equivalent to zero. A successful wrapper-only check tests
the file protocol; the client run also tests the structure, configuration
and executable lookup.

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

Positions and box vectors are in Angstrom. The first three lines are the
3x3 box matrix (one lattice vector per row, tab-separated).
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

The first line is the total energy in eV (scalar).
Each subsequent line contains the force on one atom in eV/Angstrom: `fx  fy  fz`
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

### Use an existing HPC CPMD installation

Keep the eOn client in its Pixi environment and use the site's existing CPMD
executable for the expensive force evaluations. You need all of:

1. The absolute path to a working `cpmd.x` (or site-named CPMD executable).
2. Its compatible MPI runtime and the modules used to run it.
3. The licensed pseudopotentials for every element in `pos.con`, with the
   correct functional and angular-momentum choices.
4. A periodic `pos.con`, the wrapper script and `config.ini`.

For example, clone the public eOn sources using Pixi's Git. This supplies
the CPMD wrapper and a small silicon-vacancy structure; it does not build
or install another eOn client:

```bash
pixi add git
pixi run git clone --branch develop https://github.com/TheochemUI/eOn.git eon-source
pixi run git -C eon-source rev-parse HEAD
cp eon-source/examples/akmc-cpmd-slurm/potfiles/cpmd_extpot.py ./cpmd_extpot.py
cp eon-source/examples/akmc-cpmd-slurm/pos.con ./pos.con
chmod +x cpmd_extpot.py
```

Record the printed source commit alongside the Pixi lockfile. The silicon
example contains seven atoms and requires a silicon pseudopotential; it is
separate from the harmonic demonstration.

An executable built with Intel MPI must use its matching Intel MPI launcher.
Do not send it through Open MPI installed for another calculator. A launcher
script can restore the site's compiler/MPI modules in the CPMD child process
without replacing the eOn environment. Save this as `launch-cpmd-mpi.sh` and
replace the module initialization and module names with those from the
working CPMD batch script:

```bash
#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module purge
unset LD_LIBRARY_PATH LD_PRELOAD LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH
module load intel/oneAPI-2022 intelmpi/oneAPI-2022
exec mpirun "$@"
```

The module initialization path is site-specific. Add further modules only
when that CPMD binary needs them. Verify this launcher on the allocated node;
the `mpirun` it selects must belong to the MPI used to link CPMD.
The library-path cleanup belongs in this child launcher before loading its
modules; applying it to the parent process can break the eOn installation.

From the workspace root, set the actual executable and pseudopotential
locations. These are examples of account-relative paths, not files supplied
by Pixi:

```bash
chmod +x launch-cpmd-mpi.sh
export CPMD_BIN="$HOME/myCPMD/bin/cpmd4.3-oneAPI-2022-AMD.x"
export PP_LIBRARY_PATH="$HOME/pseudo/"
export CPMD_LAUNCH="$PWD/launch-cpmd-mpi.sh -np 4"
export CPMD_PP="14:Si_MT_BLYP.psp:LMAX=D LOC=D"
export CPMD_CUTOFF=70
export CPMD_FUNCTIONAL=BLYP
export CPMD_CONVERGENCE=1e-6
export OMP_NUM_THREADS=1
test -x "$CPMD_BIN"
test -f "${PP_LIBRARY_PATH}Si_MT_BLYP.psp"
cat > config.ini <<EOF
[Main]
job = point

[Potential]
potential = ext_pot
ext_pot_path = $PWD/cpmd_extpot.py
EOF
pixi run eonclient > eon-cpmd.log 2>&1
cat results.dat
```

Run this command inside a compute allocation with four MPI tasks available
for CPMD. The eOn client is one process; its wrapper launches those four
CPMD processes. For Slurm, request a node with `--ntasks=4` and enter the
node using the site's supported interactive-job procedure. A site that
requires `srun` instead of `mpirun` needs a launcher adapted to its MPI/PMI
configuration. Do not start four copies of `eonclient` for this file route.

For the Intel oneAPI 2022 setup, the source tree also supplies a complete
Slurm job at `examples/akmc-cpmd-slurm/run_point.sbatch`. From the workspace
root, copy it and edit its `#SBATCH` partition/account and the child
launcher's module names to match the site. Submit through the Pixi
environment so the installed client and Python are inherited by the job:

```bash
cp eon-source/examples/akmc-cpmd-slurm/run_point.sbatch ./run_point.sbatch
# Edit the scheduler directives and site module names in run_point.sbatch.
export EON_REPO="$PWD/eon-source"
pixi run sbatch ./run_point.sbatch "$PWD/pos.con"
```

Keep `CPMD_BIN`, `PP_LIBRARY_PATH` and the method variables exported as
above. The script's own defaults are the public Si7 vacancy demonstration
at 20 Ry with `LMAX=P`; an exported override changes those defaults.
Slurm prints `Submitted batch job <jobid>`. The calculation's record is
`eon-cpmd-point-<jobid>/`, including the copied geometry, configuration,
wrapper and `eon.log`. The scheduler output prints the converted finite
energy when the point calculation succeeds. The script leaves that
directory in the submission filesystem for inspection.

The result contains a finite `potential_energy`. Its numerical value depends
on the supplied pseudopotential, cutoff, geometry and convergence settings;
the relaxed structure's recorded example energy is not a universal check
at another cutoff. Inspect the exchange directory's `cpmd.out`, `input`,
`GEOMETRY` and `RESTART.1` if the call fails. Keep these with the result when
reporting an error.

The supplied wrapper builds a neutral, single-point BLYP deck by default.
Its environment controls listed below do **not** provide a complete method
description for a charged, spin-polarized, dispersion-corrected or constrained
system. Such a calculation requires a wrapper/deck that explicitly retains
those physical settings. A dynamics or restart deck cannot be substituted
without checking its task, charge, cell, functional, dispersion, species,
occupations and SCF convergence. The in-process
[RGPOT route](project:rgpot_pot.md) uses a separate serialized CPMD method
message; it does not read this wrapper's environment as that message.

### Units and force sign

The external interface expects the **total energy in eV**, positions and box
vectors in **Angstrom**, and forces in **eV/Angstrom**. For CPMD:

```text
energy_eV       = energy_Hartree * 27.211386245988
force_eV_per_A  = force_Hartree_per_Bohr * 51.4220674763
force           = -gradient
```

The last line applies when the source quantity is an energy gradient.
CPMD's `GEOMETRY` columns 4–6 already contain the ionic **force** `fion`:
convert the units without negating them. Returning a positive gradient as
a force drives minimization uphill. The header `GRADIENTS (-FORCES)` in
the printed CPMD table does not change what `fion` contains.

### CPMD on several MPI ranks

`examples/akmc-cpmd-slurm/potfiles/cpmd_extpot.py` answers each call with a
`cpmd.x` wavefunction optimisation on the MPI ranks that `CPMD_LAUNCH`
starts. It reads its
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
| `CPMDC_RESTART` | path to the `cpmdc-restart` tool; a pooled `RESTART.1` seeds a call only when `cpmdc-restart info` reports the same atom count, species counts and cutoff | unset, no check |

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
hardware threads rather than physical cores.

`examples/akmc-cpmd-slurm` runs adaptive kinetic Monte Carlo on a silicon
vacancy with the same script. The `cluster` communicator submits each
process search as a Slurm job (see [Communicator](project:communicator.md)),
and inside the job `ext_pot` calls `cpmd_extpot.py`. The client stays one
process; the ranks belong to CPMD alone.

What the script hands CPMD and reads back:

- The deck says `ANGSTROM`. CPMD then converts eOn's positions with its
  own Bohr, 0.529177210859 A. CPMD answers in Hartree and Hartree/Bohr. The
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

## Troubleshooting

| Symptom | Check |
| --- | --- |
| Wrapper cannot be executed or exit status 126/127 | Make it executable; check its shebang/interpreter under `pixi run`; use an absolute `ext_pot_path`. |
| No `from_extpot_to_eon` | Run the wrapper manually from its exchange directory and inspect the calculator log and exit status. Write the response only after the calculation succeeds. |
| Too few force rows | Return exactly one force triple per input atom, in the same order; restore species-grouped calculator output to the input order. |
| `cpmd.x` cannot load a shared library or MPI fails | Use the compiler/MPI runtime matching that executable inside the CPMD launcher, not the Pixi MPI for a different build. |
| SCF does not converge | Inspect `cpmd.out`, the rendered input and restart compatibility. Do not pass an unconverged energy to eOn. |
| Minimization moves uphill | Check force units and sign; forces are the negative energy derivative. |
| Several concurrent image calculations exhaust the allocation | Count calculator ranks per concurrent image, rather than counting the single eOn process alone. |

An external wrapper must return a nonzero status for a failed calculation;
it must not invent an energy or reuse a stale response. Return finite values
and preserve diagnostics. A point job validates setup before a minimization,
NEB or aKMC run makes many expensive calls.

## Performance considerations

Each force call spawns a new process and re-initializes the calculator. For
calculators with substantial initialization or model-loading costs, startup
can dominate the force evaluation. Measure it for the actual calculator. In that
scenario consider:

- The built-in [ASE interface](project:ase_pot.md) (requires building from
  source with `-Dwith_python=True -Dwith_ase=True`), which embeds the
  interpreter and avoids per-call process overhead.
- The [Metatomic](project:metatomic_pot.md) interface for models that support
  the metatensor format (compiled into the `conda-forge` build).
- The [serve mode](project:serve_mode.md) to expose any eOn potential over RPC
  to external callers.
