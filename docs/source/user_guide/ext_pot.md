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

## Run CPMD through eOn

Use the runnable example in `examples/ext_pot_cpmd` on the **main** branch.
It evaluates a geometry with your existing CPMD executable and
pseudopotentials. Pixi supplies eOn, Python and uv; no cpmdc library or eOn
source compilation is needed for this external route.

The example contains:

| File | Purpose |
| --- | --- |
| `pixi.toml` | Installs eOn, Python and uv; defines `run` and `neb` tasks. |
| `.env.example` | Executable, pseudopotential, MPI and physical-method settings; copy it to `.env`. |
| `config.ini` | A point job with `ext_pot_path = ./ext_pot`. |
| `ext_pot` | Loads `.env` and invokes the CPMD file-interface wrapper. |
| `cpmd_extpot.py` | Renders a single-point deck, runs CPMD and converts its results. |
| `cpmd_launch.sh` | Starts CPMD with the compatible MPI runtime. |
| `run.sh` | Converts an input if needed and preserves a separate point-job record. |

### Install the example environment

Use a Linux terminal, including Ubuntu under WSL. Install
[Pixi](https://pixi.sh), clone the public repository, then install the
example's environment:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"
git clone --branch main https://github.com/TheochemUI/eOn.git
cd eOn/examples/ext_pot_cpmd
pixi install
pixi run eonclient --version
cp .env.example .env
```

Git is needed for the initial clone. Pixi installs the scientific
environment without administrator rights. Keep its generated `pixi.lock`
with `pixi.toml` when sharing a resolved environment; another account can
recreate it using `pixi install --locked`. On an HPC system, run calculations
inside a compute allocation according to the site policy. A workstation
does not require Slurm.

### Set the existing CPMD paths

Edit `.env`. For an existing Intel oneAPI CPMD installation, the essential
settings look like:

```bash
EON_CLIENT=eonclient
CPMD_BIN=$HOME/myCPMD/bin/cpmd4.3-oneAPI-2022-AMD.x
PP_LIBRARY_PATH=$HOME/pseudo
CPMD_PP='14:Si_MT_BLYP.psp:LMAX=D LOC=D,7:N_MT_BLYP.psp:LMAX=P LOC=P'
CPMD_MPI_LAUNCH='mpirun -np 4'
CPMD_MODULES='intel/oneAPI-2022 intelmpi/oneAPI-2022'
CPMD_FUNCTIONAL=BLYP
CPMD_CUTOFF=70
CPMD_CONVERGENCE=1e-6
```

Replace the paths, module names, pseudopotential mapping and ranks with
the actual CPMD installation. The executable and licensed pseudopotentials
are supplied by your group, not downloaded by Pixi. Every element in the
structure needs a mapping in `CPMD_PP`. Its third field is CPMD's
angular-momentum/options line; the wrapper uses Kleinman–Bylander
pseudopotentials.

Leave `CPMD_MODULES=''` on a workstation where the compatible runtime is
already available. A binary linked with Intel MPI needs the matching
Intel MPI launcher; another environment's Open MPI is not a substitute.
When modules are requested, only the CPMD child clears inherited library
overrides and loads them. The parent eOn client keeps its Pixi libraries.

`.env` is sourced by Bash: values set there take precedence over the
environment inherited by the launcher. Keep account paths in `.env`,
which is ignored by Git.

### Evaluate a geometry

From `examples/ext_pot_cpmd`, use the same command for XYZ or CON input:

```bash
pixi run run /absolute/path/input.xyz
pixi run run /absolute/path/pos.con
```

Each invocation prints `Record: .../runs/point-...`, evaluates one geometry,
and prints a finite energy in eV when successful. The record retains
`pos.con`, `config.ini`, `eon.log` and `results.dat`. CPMD input, output and
`GEOMETRY` remain in its external-potential exchange directory. A failed
client or a non-finite energy causes a nonzero exit; inspect the preserved
logs rather than using a partial result.

For the supplied CON example, the static configuration can also be run
directly from this directory after editing `.env`:

```bash
pixi run eonclient > eon.log 2>&1
cat results.dat
```

This direct command uses the existing `pos.con` and `config.ini` in place;
`pixi run run` creates a separate record and handles an arbitrary input.
For a file route using four CPMD ranks, reserve four ranks for each
concurrent calculator. The eOn client itself stays one process.

### XYZ frames, units and the periodic cell

The XYZ launcher uses `uv run --script xyz_to_con.py` with the pinned
`readcon-chemfiles==0.17.1` dependency to read XYZ and write CON. The first
conversion needs network access unless the wheel is cached. Prepare that
cache with an actual input on a host where dependency installation is
permitted:

```bash
pixi run uv run --script xyz_to_con.py /absolute/path/input.xyz /tmp/input.con
```

Set `XYZ_FRAME=-1` for the last frame, or a zero-based nonnegative frame
index. `XYZ_UNITS=angstrom` is the default; `XYZ_UNITS=bohr` converts
coordinates and a metadata cell to Angstrom. An explicit
`CELL_ANGSTROM='9.984 19.968 14.976'` is always in Angstrom. Edit these
values in `.env` before launching.

A positive cell in extended XYZ metadata is preserved unless
`CELL_ANGSTROM` overrides it. This wrapper requires an orthorhombic box.
For plain XYZ with no cell, conversion prints a **LOUD WARNING** and
infers a box from the coordinate bounds with `VACUUM_ANGSTROM=10` on
each side. It translates the coordinates into that artificial box;
**this does not recover the original bulk cell**. Supply the actual
periodic cell for a material calculation. `pos.source.json` records
the frame, input units, cell source, translation, atom IDs and the
mapping between CON and original XYZ order.

### Run a nudged elastic band

Use two endpoints with the same atoms, atom correspondence and periodic
cell, relaxed on the same physical Hamiltonian:

```bash
pixi run neb /absolute/path/reactant.xyz /absolute/path/product.xyz
```

The `neb` task uses `run_neb.sh` and `config-neb.ini`; XYZ endpoints use
the same conversion and cell settings. Confirm endpoint identities and
the real cell before launching. Reserve resources for the number of
concurrent images times the CPMD rank count, and inspect the band's
termination and force convergence before treating its peak as a barrier.

## CPMD method requirements

The distributed settings use a neutral BLYP single-point calculation.
A charged, spin-polarized, dispersion-corrected or constrained system
needs a deck/wrapper that explicitly retains those settings. Replacing
its geometry alone does not reproduce that Hamiltonian. Check charge,
occupations, functional, dispersion, pseudopotentials, cell, cutoff and
SCF convergence against the actual source calculation.

The [RGPOT route](project:rgpot_pot.md) loads cpmdc in process and uses a
separate CPMD method message. This external executable route reads its
own settings from `.env`; the two configurations are distinct inputs.

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

## Calculator units

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

The runnable point and NEB commands above use this wrapper. The older
`examples/akmc-cpmd-slurm` campaign uses the same protocol for process
searches submitted through the cluster communicator (see
[Communicator](project:communicator.md)). Its scheduler setup is separate
from the workstation launcher in `examples/ext_pot_cpmd`.

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

### Empty CPMD output or a site that requires an allocation

The CPMD example launcher writes `CPMD launcher: initializing`, the module
loading phase and the selected executable to `cpmd.out`. If these lines
are absent, check that the checkout contains the current launcher and
that the inspected file belongs to the failed call. eOn names exchange
directories `extpot_<process-id>_<call-id>`; use the directory printed in
the error rather than guessing a hidden directory name. If initialization
appears but execution does not, inspect the module-loading error. If
execution appears, inspect the MPI/runtime diagnostics that follow it.
An empty file alone does not identify an SCF failure.

A workstation can run the example directly. On a cluster whose MPI launch
requires a compute allocation, submit the entire eOn point job from
`examples/ext_pot_cpmd`, with the site's partition and account:

```bash
sbatch --partition=YOUR_PARTITION --account=YOUR_ACCOUNT submit.sbatch
```

The script requests four MPI ranks and runs `run.sh`; its job output is
`eon-cpmd-<job-id>.out`, and `run.sh` prints the record directory. Edit the
rank request and `CPMD_MPI_LAUNCH` together if the calculation needs a
different rank count. A batch script that runs a generated `cpmd.inp`
directly checks CPMD startup, but does not resume an eOn process that has
already exited. Run eOn inside the allocation to complete the external
potential exchange.

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
