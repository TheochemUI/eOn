---
myst:
  html_meta:
    "description": "Using LAMMPS potentials with eOn via runtime dynamic loading."
    "keywords": "eOn LAMMPS, LAMMPS potential, molecular dynamics, dlopen"
---

# LAMMPS Potential

```{versionchanged} 2.12
LAMMPS is now loaded at runtime via `dlopen`/`LoadLibrary`. No compile-time
flag is needed. A single eOn binary can use LAMMPS potentials if `liblammps`
is installed, without requiring LAMMPS at build time.
```

```{admonition} conda-forge availability
:class: tip
The `conda-forge` eOn package works with LAMMPS out of the box. Simply install
LAMMPS in the same environment: `conda install lammps`. eOn will find and load
`liblammps` automatically at runtime.
```

## How It Works

eOn uses a singleton `LammpsLoader` that searches for the LAMMPS shared library
at the first potential evaluation. The search order is:

- **Linux**: `liblammps.so`, `liblammps.so.0`
- **macOS**: `liblammps.dylib`, `liblammps.0.dylib`
- **Windows**: `lammps.dll`, `liblammps.dll`

If the library is found, seven LAMMPS C API functions are loaded via `dlsym`.
If not found, eOn prints a clear error message with installation instructions.

The cross-platform dynamic loading is handled by `DynLib.h`, a lightweight
header-only abstraction over `dlopen`/`LoadLibrary`.

## Setup

Install LAMMPS in the same conda/pixi environment as eOn:

```{code-block} bash
conda install -c conda-forge lammps
```

Or build LAMMPS from source and install the shared library:

```{code-block} bash
cmake -D CMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
      -D BUILD_MPI=no \
      -D BUILD_SHARED_LIBS=yes \
      ../cmake
make -j4
make install
```

No eOn rebuild is needed. The `potential = lammps` configuration option works
immediately once `liblammps` is on the library search path.

## Usage

Set the potential to `lammps` in the configuration file and place a LAMMPS
input file named `in.lammps` in the `potfiles` directory. This file specifies
which LAMMPS potential to use. Example for the Morse potential:

```ini
pair_style morse 9.5 #morse potential with 9.5 Angstrom cutoff
pair_coeff * * 0.7102 1.6047 2.797 #specify parameters
pair_modify shift yes #shift the potential to be zero at the cutoff
```

## Troubleshooting

If eOn reports "LAMMPS library not found", ensure that:

1. `liblammps.so` (or `.dylib`/`.dll`) is on the library search path
   (`LD_LIBRARY_PATH` on Linux, `DYLD_LIBRARY_PATH` on macOS)
2. The LAMMPS library was built as a shared library (`BUILD_SHARED_LIBS=yes`)
3. The LAMMPS version is compatible (tested with LAMMPS 2Aug2023 and later)

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: LMP_
keyprefix: lmp-
---
```
