---
myst:
  html_meta:
    "description": "Installation guide for the eOn software package, using pixi/conda or source builds."
    "keywords": "install eOn, build eOn, conda, meson, compilation"
---

# Installation

```{toctree}
:hidden:

eessi
lammps
```

eOn is divided up into two separate programs: a server and a client. The client
does most of the computation (e.g. saddle searches, minimizations, and molecular
dynamics) while the server creates the input for the client and processes that
output.

## Getting started

The shortest start is the `conda` package:

```{code-block} bash
# best with pixi
pixi init
pixi add eon
# or with conda/micromamba
micromamba install -c conda-forge eon
```

Those commands install the conda-forge release package. They do not check out `develop`. A source build of `develop` is the section below. The examples in this book run against the release package. `pixi.toml` on `develop` records version 3.6.0. The conda-forge package does not contain the CPMD engine. It is built without `-Drgpot:with_mpi=enabled`, so it does not ship calculator groups or in-process CPMD.

The conda package is a maximalist build with the following potentials and
features enabled:

- [Metatomic](project:../user_guide/metatomic_pot.md) (machine-learned potentials via libtorch)
- [xTB](https://xtb-docs.readthedocs.io/) (semi-empirical tight-binding)
- [rgpot integration](project:../user_guide/rgpot_integration.md) (direct dlopen vs serve vs potserv client)
- [RgpotPot / RGPOT](project:../user_guide/rgpot_pot.md) (in-process NWChemPot/CPMDPot on a source build; the conda-forge package does not enable calculator groups)
- [Serve mode](project:../user_guide/serve_mode.md) (`-Dwith_serve`: eOn as an rgpot-compatible remote procedure call server)

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

Once git is present[^1], clone `develop`. A clone with no branch name follows the repository default.

```{code-block} bash
git clone -b develop https://github.com/TheochemUI/eOn.git
cd eOn
```

```{note}
* Authentication is easier with the [command line tool](https://cli.github.com/).
* [Pixi](https://pixi.sh/) is now recommended
```

## Building from source

Pixi installs the dependencies recorded in `pixi.lock`. `pixi shell` opens the default environment. `dev-lite` is the lighter shell: the `lint` and `develop` features from `pixi.toml`.

```{code-block} bash
pixi shell
# or
pixi shell -e dev-lite
```

Other environments can be found by inspecting the `pixi.toml` file.

This is the installation path that fails least often:

```{code-block} bash
# conda-compilers may try to install to
# $CONDA_PREFIX/lib/x86_64-linux-gnu
# without --libdir
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype=release \
  --force-fallback-for=nlohmann_json,hwy \
  -Drgpot:with_mpi=enabled
meson compile -C bbdir
meson test -C bbdir --suite eon
meson install -C bbdir
```

`-Drgpot:with_mpi=enabled` turns on calculator groups. eOn's own `-Dwith_mpi` is the AKMC client, not those groups. In-process CPMD also needs `libcpmdc` on the loader path. The client tests skip that engine when the library is absent.

`min_mode_method = gprdimer` links the GP dimer from a checkout of
[gpr_optim](https://github.com/TheochemUI/gpr_optim) in
`subprojects/gpr_optim`. `-Dwith_gprd=auto` is the default. On Linux
that checkout is linked when it is present. `-Dwith_gprd=enabled`
stops configuration when the checkout is absent.
A build directory that stored `with_gprd` as `true` or `false` rejects
`meson setup --reconfigure` (`Option "with_gprd" value auto is not boolean`).
Run `python scripts/migrate_with_gprd_option.py <builddir>` once: `true`
becomes `enabled` and `false` becomes `disabled`.
Copy a sibling checkout into place with:

```{code-block} bash
rsync -a ../gpr_optim/ subprojects/gpr_optim/
```

A build configured with `-Dwith_gprd=disabled` has no GP dimer. Asking
for `gprdimer` then stops, and the message names `-Dwith_gprd=enabled`
and that `rsync`.

The setup line already passes `--force-fallback-for=nlohmann_json,hwy`. The rolling distro section below says why a host `nlohmann_json` breaks the conda compiler, and why the Highway name is `hwy` rather than `libhwy`. The test line should finish with a fail count of 0.

Some additional performance can be gained with `ccache` and `mold`, which can be
passed with `--native-file`:

- With `ccache` installed, add `--native-file nativeFiles/ccache_gnu.ini`
- With `mold` installed, add `--native-file nativeFiles/mold.ini`

### Troubleshooting on rolling distros

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
  --force-fallback-for=nlohmann_json,hwy
```

`hwy` is the dependency name `highway.wrap` provides. `--force-fallback-for=libhwy`
does not select that wrap: pkg-config can find `libhwy` and the cmake
subproject is never configured. `client/meson.build` looks up `hwy`.

### Troubleshooting: a global cargo linker setting

readcon-core is built by cargo, which reads `~/.cargo/config.toml` in addition
to the environment. A `rustflags` entry there applies to every crate eOn builds.
That includes a conda or pixi environment. The compiler there is conda-forge's.

One failure is a linker override. A system GCC that can resolve `mold` by name accepts `-C link-arg=-fuse-ld=mold`. conda-forge's GCC 13 rejects `-fuse-ld=/usr/bin/mold`:

```text
error: linking with `x86_64-conda-linux-gnu-cc` failed: exit status: 1
  = note: x86_64-conda-linux-gnu-cc: error: unrecognized command-line option
          '-fuse-ld=/usr/bin/mold'
error: could not compile `zmij` (build script)
FAILED: subprojects/readcon-core/libreadcon_core.a
```

`RUSTFLAGS` overrides `build.rustflags` from the config file, so clearing it for
the build leaves the global setting alone:

```{code-block} bash
RUSTFLAGS= meson compile -C bbdir
```

The same applies to `scripts/pyeonclient_build_wheel.sh`, which builds
readcon-core through cargo as part of the wheel.

### Optional packages

The full listing of options is found in the `meson_options.txt` file. These can
all be turned on and off at the command line. As an example see the [LAMMPS
integration instructions](project:../user_guide/lammps_pot.md).

For an optional wrapped dependency, download the subproject sources before configuring:

```{code-block} bash
meson subprojects download artn-plugin ira
```

# European environment for scientific software installations (EESSI)

A `develop` checkout on this 2026.06 release uses EESSI-extend, the `eOn-devel` bundle, CapnProto 1.4.0, and `foss/2026.1`. The run-path list and the calculator groups are on the [build page](eessi.md).

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
