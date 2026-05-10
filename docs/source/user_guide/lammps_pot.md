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

Set the potential to `lammps` in the configuration file and provide a LAMMPS
input file named `in.lammps` plus any potential parameter files (e.g.
`Cu01.eam.alloy`). The `in.lammps` file specifies which LAMMPS potential
to use. Example for the Morse potential:

```ini
pair_style morse 9.5 #morse potential with 9.5 Angstrom cutoff
pair_coeff * * 0.7102 1.6047 2.797 #specify parameters
pair_modify shift yes #shift the potential to be zero at the cutoff
```

### File placement

Two ways to feed LAMMPS its run-time files. Pick one:

#### 1. Bundle (recommended) -- `lammps_bundle = path/to/run.eonlpb`

```{versionadded} 2.15
```

Pack `in.lammps` and every file it references (pair_coeff data,
custom `pair_style` `.so` plugins, KIM tables, `read_data` inputs,
shell helpers, ...) into one `.eonlpb` blob. Point eOn at the bundle
and the eonclient CWD becomes irrelevant for LAMMPS file lookups.

```{code-block} bash
# Pack the directory holding in.lammps + all referenced files
python -m eon.lammps_bundle pack potfiles/ run.eonlpb

# Inspect what got packed
python -m eon.lammps_bundle list run.eonlpb
```

```{code-block} ini
[Potential]
potential = lammps
lammps_bundle = run.eonlpb
```

`LAMMPSPot` extracts the bundle into a per-instance scratch dir
under `$TMPDIR` and issues `shell cd <scratch>` to liblammps before
sourcing `in.lammps`, so every relative reference inside `in.lammps`
resolves there. The scratch dir is removed when the potential is
destroyed. Multiple `eonclient` instances on the same host get their
own private scratch dirs (PID + 64-bit random suffix).

#### 2. Legacy CWD mode (no bundle)

If `lammps_bundle` is unset, `in.lammps` and every file it references
must live in the directory where `eonclient` runs:

- **Multi-job drivers** (`job = akmc`, `job = parallel_replica`,
  `job = basin_hopping`, ...): keep them in `potfiles/` next to your
  `config.ini`. The Python driver copies the contents into each
  per-job scratch directory (`jobs/scratch/<wuid>/`) before launching
  `eonclient`.

- **Direct single-job `eonclient`** (`job = minimization`,
  `job = nudged_elastic_band`, `job = single_point`,
  `job = process_search`): place `in.lammps` and every referenced file
  in the SAME directory as `config.ini`. eOn now refuses to construct
  a `LAMMPSPot` if `in.lammps` is missing from CWD and surfaces the
  CWD path in the error message; `lammps_has_error` is also checked
  after `lammps_file` so LAMMPS-side syntax / pair_coeff / pair_style
  errors no longer hide.

```{admonition} Why a bundle?
:class: tip
Pre-2.15 eOn coupled the eonclient CWD to LAMMPS's CWD: any file
liblammps reads (pair_coeff, read_data, custom plugin `.so`, ...) had
to be in eonclient's CWD or absolute. Job-scratch cleanup, tmpfs
eviction, or a wrapper script chdir-ing somewhere unexpected would
silently break the next force call. The bundle decouples them; one
blob travels with the run.
```

## Troubleshooting

If eOn reports "LAMMPS library not found", ensure that:

1. `liblammps.so` (or `.dylib`/`.dll`) is on the library search path
   (`LD_LIBRARY_PATH` on Linux, `DYLD_LIBRARY_PATH` on macOS).
   For pixi/conda envs, exporting
   `LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH` is usually enough.
2. The LAMMPS library was built as a shared library (`BUILD_SHARED_LIBS=yes`)
3. The LAMMPS version is compatible (tested with LAMMPS 2Aug2023 and later;
   conda-forge `lammps=2024.08.29` works with eOn 2.14)

If `eonclient` reports `LAMMPSPot: in.lammps not found in eonclient CWD`,
either point `lammps_bundle` at a `.eonlpb` blob (recommended) or place
`in.lammps` and every file it references next to `config.ini` (legacy
CWD mode). If liblammps surfaces a syntax / pair_coeff / pair_style
error after sourcing `in.lammps`, eOn now reports the LAMMPS-side
message verbatim instead of segfaulting on the next force call.

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: LMP_
keyprefix: lmp-
---
```
