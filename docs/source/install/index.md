# Installation

eOn is divided up into two separate programs: a server and a client. The client
does most of the computation (e.g. saddle searches, minimizations, and molecular
dynamics) while the server creates the input for the client and processes its
results.

## Getting started

The simplest way to hit the ground running is with the `conda` package:

```{code-block} bash
micromamba install -c "https://prefix.dev/channels/rg-forge" eon
```

At this point any of the many examples should be good to go.

The server is accessed through `python -m eon.server`, and the `eonclient`
binary is automatically made available in the activated environment..

```{versionchanged} 2.0
While reading older documentation, calls to `eon` must now be `python -m
eon.server`.
```

# Obtaining sources

```{versionadded} 2.0
`eOn` is now developed and distributed primarily via GitHub.
```

Once git is present[^1]:

```{code-block} bash
git clone https://github.com/TheochemUI/eOn.git
cd eOn
```

````{margin}
```{note}

* The [GitHub CLI tool](https://cli.github.com/) makes authentication much easier.
* Micromamba ([installation](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)) is a faster conda-helper
```
````

## Building from source

We provide a `conda` environment and `pixi` setup, with dependencies handled by `conda-lock`.

```{code-block} bash
micromamba create -n eongit -f conda-lock.yml
micromamba activate eongit
# Or
pixi shell
```

This leads to the most robust installation approach:

```{code-block} bash
# conda-compilers may try to install to $CONDA_PREFIX/lib/x86_64-linux-gnu without --libdir
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype=release
meson install -C bbdir
```

Some additional performance can be gained with `ccache` and `mold`, which can be
passed with `--native-file`:

- With `ccache` installed, add `--native-file nativeFiles/ccache_gnu.ini`
- With `mold` installed, add `--native-file nativeFiles/mold.ini`

### Optional packages

The full listing of options is found in the `meson_options.txt` file. These can
all be turned on and off at the command line. As an example see the [LAMMPS
integration instructions](project:../user_guide/lammps_pot.md).

# Additional topics

This page lists generally applicable installation instructions, for specific
systems, follow the sub-parts of this document.

```{toctree}
:maxdepth: 2
:caption: Special topics

windows
lammps
metatomic
ase
svn
```

# Licenses

`eOn` is released under the [BSD 3-Clause
License](https://opensource.org/license/BSD-3-Clause).

## Vendored
Some libraries[^2] are distributed along with `eOn`, namely:

- `mcamc` which contains `libqd` :: BSD-3-Clause license
```{versionadded} 2.0
- `cxxopts` :: MIT License
- `magic_enum` :: MIT License
- `catch2` :: Boost Software License, Version 1.0
- `ApprovalTests.cpp` :: Apache 2.0 License
```
```{deprecated} 2.0
- Eigen 2.x :: Mozilla Public License
```


<!-- pipx run pdm run sphinx-build -b html docs/source docs/build/html -->

[^1]: Installation instructions [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
[^2]: All with compatible licenses
