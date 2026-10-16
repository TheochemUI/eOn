# CPMD through the external potential interface

Use an existing `cpmd.x` executable to evaluate an XYZ or CON geometry through
eOn. This example uses the external potential interface and needs no cpmdc
library. The Pixi environment supplies eOn, Python and uv; `cpmd.x` and its
pseudopotentials come from your CPMD installation.

Install Pixi, then create the example environment:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"
git clone --branch main https://github.com/TheochemUI/eOn.git
cd eOn/examples/ext_pot_cpmd
pixi install
cp .env.example .env
```

Edit `.env` with your executable and potential paths: `CPMD_BIN` names
`cpmd.x`, `PP_LIBRARY_PATH` names its pseudopotential directory, and
`CPMD_PP` maps every element's atomic number to its pseudopotential filename
and options. For example, a Si/N calculation needs both Si and N entries.
`CPMD_MPI_LAUNCH` specifies the launcher and ranks appropriate to that CPMD
binary. Leave `CPMD_MODULES` empty on a workstation; a site using Intel
modules can list them there. Only the CPMD child loads these modules and
clears inherited library overrides; eOn keeps the Pixi environment.

Run a geometry on a workstation:

```bash
pixi run run /absolute/path/input.xyz
```

A CON input uses the same command:

```bash
pixi run run /absolute/path/pos.con
```

CON inputs are read and written through the pinned readcon converter before
evaluation. It preserves the cell, coordinates, atom IDs, masses and fixed
flags, and removes stored energies and forces so the point uses your configured
CPMD executable and pseudopotentials. The source file is preserved.
`CON_FRAME=0` selects the first CON frame; set `CON_FRAME=-1` in `.env` to
select the last. `pos.source.json` records the selected frame and file hashes.

For a nudged elastic band between two endpoints:

```bash
pixi run neb /absolute/path/reactant.xyz /absolute/path/product.xyz
```

`run.sh` has no Slurm dependency. On an HPC system run the command inside a
compute allocation according to the site's policy. If a site supplies a
submission script, that script can call the same `pixi run run` command.

## XYZ frames, units and cells

The launcher calls `uv run --script xyz_to_con.py`, which installs the pinned
`readcon-chemfiles==0.17.1` dependency and uses readcon for both XYZ reading
and CON writing. The first conversion needs network access to fetch that
wheel unless uv already has it cached. To prepare the conversion before a
compute node loses network access, run it with your actual input on a host
where installing dependencies is permitted:

```bash
pixi run uv run --script xyz_to_con.py /absolute/path/input.xyz /tmp/input.con
```

For a CON source, prime the same readcon dependency with:

```bash
pixi run uv run --script con_to_input.py /absolute/path/source.con /tmp/point-input.con
```

`XYZ_FRAME=-1` selects the last frame; nonnegative values select a zero-based
frame. `XYZ_UNITS=angstrom` is the default. `XYZ_UNITS=bohr` converts both
coordinates and a cell read from XYZ metadata to Angstrom. An explicit
`CELL_ANGSTROM` is always in Angstrom. To select the first frame of a
Bohr-coordinate input with a known periodic cell, edit these lines in `.env`:

```bash
XYZ_FRAME=0
XYZ_UNITS=bohr
CELL_ANGSTROM='9.984 19.968 14.976'
```

Then run `pixi run run /absolute/path/input.xyz`. The launcher sources
`.env`, so its values take precedence over one-command environment
assignments. A positive cell supplied by XYZ metadata is retained
unless `CELL_ANGSTROM` overrides it. This generic CPMD wrapper requires an
orthorhombic cell and rejects non-right cell angles.

For a plain XYZ without cell metadata, the converter prints a **LOUD
WARNING** and infers a box from the coordinate bounds, adding
`VACUUM_ANGSTROM=10` on each side. It translates the coordinates into this
box. The warning includes its lengths and states that this artificial box
is **not the original bulk cell**. Supply the real `CELL_ANGSTROM` to
reproduce a periodic material. Changing an unknown box changes the physical
calculation.

Each call creates a fresh `runs/point-*` directory. `pos.source.json` records
the selected frame, units, cell source, translation, atom IDs and mapping
between CON order and the original XYZ order. `pos.con`, `config.ini`,
`eon.log` and `results.dat` stay with that calculation. CPMD input, output
and `GEOMETRY` stay in its external-potential exchange directory. A failed
client or non-finite energy produces a nonzero exit. `.env`, downloaded
environments and run outputs are ignored by Git.

## Hamiltonian and units

The distributed settings use neutral BLYP at 70 Ry, with Si `LMAX=D LOC=D`
and N `LMAX=P LOC=P`. Keep the pseudopotential mapping consistent with the
species in your input. These settings do not reproduce the charged BLYP/D2
amorphous-cell Hamiltonian merely by substituting its XYZ. Use that system's
validated charge, dispersion and cell settings for that calculation.

The wrapper converts total energy from Hartree to eV using
`27.211386245988`, and forces from Hartree/bohr to eV/Angstrom using
`51.4220674763`. CPMD `GEOMETRY` contains forces, so it applies no sign
reversal. The coordinates passed by eOn are in Angstrom.
