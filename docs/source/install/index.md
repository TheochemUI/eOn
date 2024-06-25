# Installation

Eon is divided up into two separate programs: a server and a client. The client
does most of the computation (e.g. saddle searches, minimizations, and molecular
dynamics) while the server creates the input for the client and processes its
results.

# Obtaining sources

```{versionadded} 2.0
`eON` is now developed and distributed primarily via GitHub.
```

Assuming access has been granted to the Github repository:

```{code-block} bash
git clone https://github.com/TheochemUI/EONgit.git
cd EONgit
```

````{margin}
```{note}

* The [GitHub CLI tool](https://cli.github.com/) makes authentication much easier.
* Micromamba ([installation](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)) is a faster conda-helper
```
````

## Setup

We provide a `conda` environment.

```{code-block} bash
micromamba create -f environment.yml
micromamba activate eongit
```

This leads to the most robust installation approach:

```{code-block} bash
meson setup bbdir --prefix=$CONDA_PREFIX
meson install -C bbdir
```

Occasionally, the library is installed to the wrong path[^1], in which case the
output of the above command should be checked for the location and then that
should be added to `$PATH`.

```{code-block} bash
meson install -C bbdir | grep libeon
Installing client/libeoncbase.so to /micromamba/envs/eongit/lib
Installing client/libeonclib.so to /micromamba/envs/eongit/lib
```

At which point we should add the library path:

```{code-block} bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/micromamba/envs/eongit/lib
```

The server is accessed through `python -m eon.server`, and the `eonclient`
binary is automatically made available in the activated environment..

```{versionchanged} 2.0
While reading older documentation, calls to `eon` must now be `python -m
eon.server`. 
```

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
svn
```

# Licenses

`eON` is released under the [BSD 3-Clause
License](https://opensource.org/license/BSD-3-Clause).

## Vendored
Some libraries[^2] are distributed along with `eON`, namely:

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

[^1]: Due to the [post-fix labeling](https://conda-forge.org/docs/user/faq/#why-dont-the-cc-compilers-automatically-know-how-to-find-libraries-installed-by-conda) of `conda-compilers`
[^2]: All with compatible licenses
