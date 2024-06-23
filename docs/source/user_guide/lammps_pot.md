# LAMMPS Potential

```{versionchanged} 2.0
Instructions now reflect using the `cmake` and `meson` build system.
```

## Setup

In order to use LAMMPS, reviewed recently in
{cite:t}`lmp-thompsonLAMMPSFlexibleSimulation2022`, as a potential you must first
build the serial library version of LAMMPS. This can be done by following the
instructions in the [lammps
documentation](https://docs.lammps.org/Build_basics.html).

For most cases executing these commands from the LAMMPS `src` folder should
work:

```{code-block} bash
cmake -D CMAKE_INSTALL_PREFIX=$CONDA_PREFIX -D BUILD_MPI=no -D BUILD_SHARED_LIBS=yes ../cmake
make
make install
```

If the same environment is used as `eon`, as described in the [installation
instructions](project:../install/index.md) then integration is via:

```{code-block} bash
meson reconfigure bbdir -Dwith_lammps=True
meson install -C bbdir
```

```{note}
With this setup, there is no need to explicitly let `eON` know about other libraries, since they are all installed together in the environment.

For backwards compatibility, we also support copying in `liblammps.so` into `potentials/LAMMPS/liblammps.so`.
```

## Usage

After setting the potential to `lammps` in the configuration file you need to
place a LAMMPS input file in the `potfiles` directory in your simulation. This
file should be named `in.lammps` and it needs to specify what potential LAMMPS
should use.  Here is an example `in.lammps` file that uses the morse potential::

```ini
pair_style morse 9.5 #morse potential with 9.5 Angstrom cutoff
pair_coeff * * 0.7102 1.6047 2.797 #specify parameters
pair_modify shift yes #shift the potential to be zero at the cutoff
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: LMP_
keyprefix: lmp-
---
```
