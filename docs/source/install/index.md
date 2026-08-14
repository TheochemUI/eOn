---
myst:
  html_meta:
    "description": "Installation guide for the eOn software package, using pixi/conda or source builds."
    "keywords": "install eOn, build eOn, conda, meson, compilation"
---

# Installation

eOn is divided up into two separate programs: a server and a client. The client
does most of the computation (e.g. saddle searches, minimizations, and molecular
dynamics) while the server creates the input for the client and processes its
results.

## Getting started

The simplest way to hit the ground running is with the `conda` package:

```{code-block} bash
# best with pixi
pixi init
pixi add eon
# or with conda/micromamba
micromamba install -c conda-forge eon
```

At this point any of the many examples should be good to go.

The conda package is a maximalist build with the following potentials and
features enabled:

- [Metatomic](project:../user_guide/metatomic_pot.md) (machine-learned potentials via libtorch)
- [xTB](https://xtb-docs.readthedocs.io/) (semi-empirical tight-binding)
- [rgpot integration](project:../user_guide/rgpot_integration.md) (direct dlopen vs serve vs potserv client)
- [RgpotPot / RGPOT](project:../user_guide/rgpot_pot.md) (`-Dwith_rgpot`: in-process NWChemPot/CPMDPot)
- [Serve mode](project:../user_guide/serve_mode.md) (`-Dwith_serve`: eOn as rgpot-compatible RPC *server*)

The server is accessed through `python -m eon.server`, and the `eonclient`
binary is automatically made available in the activated environment.

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
* [Pixi](https://pixi.sh/) is now recommended
```
````

## Building from source

We provide a `conda` environment and `pixi` setup, with dependencies handled by `conda-lock`.

```{code-block} bash
pixi shell
# or
pixi s -e dev-lite
```

Other environments can be found by inspecting the `pixi.toml` file.

This leads to the most robust installation approach:

```{code-block} bash
# conda-compilers may try to install to
# $CONDA_PREFIX/lib/x86_64-linux-gnu
# without --libdir
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype=release
meson install -C bbdir
```

If that `meson setup` stops on glibc errors raised from inside `<cmath>`, add
`--force-fallback-for=nlohmann_json` and read the rolling distro section below
for what causes it.

Some additional performance can be gained with `ccache` and `mold`, which can be
passed with `--native-file`:

- With `ccache` installed, add `--native-file nativeFiles/ccache_gnu.ini`
- With `mold` installed, add `--native-file nativeFiles/mold.ini`

### Troubleshooting: rolling distros (Arch, Fedora)

On rolling-release distributions with newer system packages, the conda-forge
compiler sysroot conflicts with system headers. `meson setup` stops during the
compiler checks with errors such as `__iseqsigf128 was not declared` or
`__fpclassify has not been declared`, raised from inside `<cmath>`. Nothing in
the output names the package responsible, so the failure reads as a broken
toolchain.

By default meson prefers an installed `nlohmann_json` module (pkg-config or
CMake, from EasyBuild and distro packages) and falls back to the wrap,
the same pattern as readcon. The system CMake config for `nlohmann_json`
exports `-I/usr/include`, which mixes glibc headers into the conda sysroot.
Force the wrap to keep that include path out of the build:

```{code-block} bash
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib \
  --force-fallback-for=nlohmann_json
```

`client/meson.build` carries the same note beside the `nlohmann_json`
dependency lookup.

### Optional packages

The full listing of options is found in the `meson_options.txt` file. These can
all be turned on and off at the command line. As an example see the [LAMMPS
integration instructions](project:../user_guide/lammps_pot.md).

For optional wrapped dependencies such as ARTn and IRA, make sure the
subproject sources are present before configuring:

```{code-block} bash
meson subprojects download artn-plugin ira
```

# Licenses

`eOn` is released under the [BSD 3-Clause
License](https://opensource.org/license/BSD-3-Clause).

## Vendored

Some libraries[^2] are distributed along with `eOn`, namely:

- `mcamc` which contains `libqd` :: BSD-3-Clause license

```{versionadded} 2.0
- `magic_enum` :: MIT License
- `catch2` :: Boost Software License, Version 1.0
- `ApprovalTests.cpp` :: Apache 2.0 License
```

```{deprecated} 2.0
- Eigen 2.x :: Mozilla Public License
```

[^1]: Installation instructions [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
[^2]: All with compatible licenses
