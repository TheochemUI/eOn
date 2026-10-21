# Structure input formats

eOn reads `.con` and `.convel` files through readcon. Builds that link
[chemfiles](https://chemfiles.org) also read `.xyz`, `.extxyz`, `.pdb`,
`.cif`, `.gro`, LAMMPS data, `POSCAR` and the other chemfiles formats. Output
stays `.con`, so restarts and the server keep working unchanged.

## Install

The pixi environment carries `chemfiles-lib`. The default
`-Dwith_chemfiles=auto` finds it, and builds the bundled CMake wrap when the
environment lacks it.

```bash
pixi shell
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype=release \
  --force-fallback-for=nlohmann_json,hwy -Drgpot:with_mpi=enabled
meson compile -C bbdir
meson test -C bbdir --suite eon
```

Pass `-Dwith_chemfiles=enabled` to stop configure when chemfiles is missing,
or `-Dwith_chemfiles=disabled` for a `.con`-only client. A client without
chemfiles refuses a foreign file with a message saying how to enable it.

## Usage

The extension picks the reader. `.con`, `.convel` and their `.gz`/`.zst`
forms go to readcon. Known chemfiles extensions (`.xyz`, `.extxyz`, `.pdb`,
`.cif`, `.gro`, `.data`, `.lmp`, `.mol2`, `.sdf`, `.vasp`, plus `POSCAR` and
`CONTCAR` names) go to chemfiles.

Jobs name their inputs `pos.con` or `reactant.con`. Set
`input_format` and eOn reads `pos.xyz` in place of `pos.con` whenever
`pos.xyz` exists:

```ini
[Main]
job = minimization
input_format = xyz
```

The job writes `min.con` as before.

### What gets mapped

- Positions, elements, masses, and velocities when the file has them.
- Fixed atoms from `.extxyz` per-atom properties: `fixed` (true means
  frozen) or ASE's `move_mask` (true means free). A three-column property
  freezes each direction on its own.

### Cell

1. A cell in the file wins, such as the `.extxyz` `Lattice=` or a `.pdb`
   `CRYST1` record. The system then runs periodic.
2. With no cell in the file, `[Main] cell` supplies one: three numbers
   `a b c` for an orthorhombic box, or nine for the lattice vectors as rows.
   `[Main] pbc` (default `true`) sets periodicity.
3. With neither, eOn logs a warning that names the file and these keys, then
   builds a non-periodic orthorhombic box around the atoms with
   `[Main] default_vacuum` (default 10 Å) on every side and shifts the atoms
   into it.

```ini
[Main]
input_format = xyz
cell = 20 20 20
pbc = true
default_vacuum = 10.0
```

The four keys are `[Main]` parameters on the client `Parameters` object.
They round-trip through JSON and appear in the config schema and in
{any}`eon.schema.MainConfig`. A config load, and a pyeonclient assignment,
publishes them process-wide: a structure read does not carry a Parameters
object. `con` or an empty `input_format` keeps readcon.
