---
myst:
  html_meta:
    "description": "Build an eOn develop checkout on EESSI 2026.06 with EESSI-extend, CapnProto 1.4.0, the buildenv run-path wrappers, and calculator groups."
    "keywords": "eOn, EESSI, EasyBuild, Meson, develop, CapnProto"
---

# Build develop on EESSI

This page builds a `develop` checkout on the European Environment for Scientific Software Installations (EESSI), release 2026.06. The check is `meson test --suite eon`. The fail count is 0.

EESSI 2026.06 is the tree under `/cvmfs/software.eessi.io/versions/2026.06`. A clone with no branch name follows the repository default. This build uses `develop`.

```{code-block} bash
git clone -b develop https://github.com/TheochemUI/eOn.git
cd eOn
```

The first `meson setup` downloads the wraps in `subprojects/`. That step needs a network.

## Install the modules

`eOn-devel/2026.06-GCCcore-15.2.0` is one EasyBuild bundle. It loads Meson, Ninja, CMake, pkgconf, Eigen, Python, and Rust. It depends on `foss/2026.1`, which loads OpenMPI, FFTW, and OpenBLAS. The bundle builds Cap'n Proto 1.4.0 into the same prefix. After `module load eOn-devel`, `pkg-config` finds `capnp-rpc` and `mpicc` is on `PATH`. `CC`, `CXX`, and `FC` stay `gcc`, `g++`, and `gfortran`.

Set the install prefix before the module load. If `EESSI_USER_INSTALL` is unset, the prefix is `$HOME/eessi`. If you set it, the directory must already exist. The load stops when that path is missing.

```{code-block} bash
source /cvmfs/software.eessi.io/versions/2026.06/init/bash
mkdir -p "${EESSI_USER_INSTALL:-$HOME/eessi}"
module load EESSI-extend
```

Some sites set `LMOD_EXACT_MATCH=yes`; Snellius is one. A short name then fails with "These module(s) or extension(s) exist but cannot be loaded". Use full names, such as `EESSI-extend/2026.06-easybuild`, or run `export LMOD_EXACT_MATCH=no` before the loads on this page.

A user prefix sets `EASYBUILD_UMASK` to `077`. Lmod, running as root, then reports the new module as unknown, because the module directory is not searchable by others. Pass `--umask=022` on every `eb` command.

A shell with nounset (`set -u`) aborts during this load. `EESSI-extend` reads `LD_LIBRARY_PATH`, and the stack leaves that variable unset. Load the module with nounset off.

```{code-block} bash
eb --umask=022 eessi/eOn-devel-2026.06-GCCcore-15.2.0.eb
```

The direct in-process rgpot arm builds rgpot's remote procedure call stack. That stack needs the `capnp` program and the `capnp-rpc` pkg-config file. EESSI 2026.06 has no separate CapnProto module. The bundle builds Cap'n Proto 1.4.0 from `https://capnproto.org/capnproto-c++-1.4.0.tar.gz` and installs it into the `eOn-devel` prefix.

## Configure and build

Open a new shell. Load the bundle, then `buildenv/default-foss-2026.1`. The bundle pulls `foss/2026.1` and puts Cap'n Proto on `PATH`. `buildenv` comes last, for the run paths below.

```{code-block} bash
source /cvmfs/software.eessi.io/versions/2026.06/init/bash
module load EESSI-extend
module load eOn-devel/2026.06-GCCcore-15.2.0
module load buildenv/default-foss-2026.1
hash -r
command -v gcc
command -v mpicc
command -v capnp
pkg-config --exists capnp-rpc
```

`capnp --version` should print `Cap'n Proto version 1.4.0`. `pkg-config --exists capnp-rpc` should return 0.

`foss/2026.1` provides `mpicc`, `mpicxx`, and `mpif90`. Export those three as `CC`, `CXX`, and `FC`. Meson then stores the Message Passing Interface (MPI) compiler wrappers.

`command -v gcc` should print a path ending in `buildenv/default-foss-2026.1/bin/rpath_wrappers/gcc_wrapper/gcc`.

Loaded modules do not set `LD_LIBRARY_PATH`. Leave it unset. `buildenv` puts EasyBuild's run-path wrappers for `gcc`, `g++`, `gfortran` and `ld` first on `PATH`. At link time they add a `-Wl,-rpath` for every `-L` directory and every `LIBRARY_PATH` entry, and every loaded module sets `LIBRARY_PATH`. So `eonclient` finds Cap'n Proto, OpenMPI and FFTW at run time, before and after `meson install`. `mpicc` calls `gcc` from `PATH`, so it goes through the same wrapper. Keep the `LDFLAGS` that `buildenv` sets; replacing it drops the toolchain's `-L` directories.

```{code-block} bash
export CC=mpicc CXX=mpicxx FC=mpif90
```

`readelf -d build-eessi/client/eonclient | grep RPATH` lists the module library directories after the build.

Point `CARGO_HOME` at an empty directory outside the checkout. Cargo then skips a host cargo configuration in its home. Cargo also reads `.cargo/config.toml` from the checkout and from every parent directory. A tree under `$HOME` still reads `$HOME/.cargo/config.toml`. Move the checkout in that case.

Clear host linker flags and compiler wrappers before configure.

```{code-block} bash
unset RUSTC_WRAPPER CARGO_BUILD_RUSTC_WRAPPER RUSTC_WORKSPACE_WRAPPER
unset RUSTFLAGS CARGO_ENCODED_RUSTFLAGS
export CARGO_HOME="$PWD/../eon-cargo-home"
mkdir -p "$CARGO_HOME"
```

A git `insteadOf` rule that rewrites `https://github.com/` to SSH also rewrites Meson wrap downloads. Drop that rule for the configure, or the wrap fetch fails. An empty global git configuration drops it for one shell:

```{code-block} bash
touch "$PWD/../eon-empty-gitconfig"
export GIT_CONFIG_GLOBAL="$PWD/../eon-empty-gitconfig"
```

That configure has two Message Passing Interface (MPI) switches. `with_mpi` in `meson_options.txt` is the eOn client switch, and this build leaves it off. Calculator groups are the rgpot subproject switch, `-Drgpot:with_mpi=enabled`. `ranks_per_image` greater than 0 throws when that switch is off.

```{code-block} bash
meson setup build-eessi \
  --prefix "$PWD/../eon-prefix" \
  --libdir lib \
  --buildtype release \
  -Dwith_tests=true \
  -Dpython.install_env=prefix \
  -Drgpot:with_mpi=enabled
meson compile -C build-eessi
```

## Run the suite

Give the ranks a writable directory, and put that full path in `TMPDIR`. Keep `OMP_NUM_THREADS` at 1 so OpenMP does not multiply the ranks. When the host has fewer free cores than ranks, set `OMPI_MCA_rmaps_base_oversubscribe` to 1.

Those ranks, on a host whose user id is 0, need two further OpenMPI variables. OpenMPI refuses to start without them. A rootless namespace is that case. `PMIX_MCA_psec=none` skips the PMIx security handshake, which that namespace does not provide. A normal account can leave the two `OMPI_ALLOW_RUN_AS_ROOT` variables unset.

```{code-block} bash
mkdir -p "$PWD/../eon-tmp"
export TMPDIR="$PWD/../eon-tmp"
export OMP_NUM_THREADS=1
export OMPI_MCA_rmaps_base_oversubscribe=1
export PMIX_MCA_psec=none
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
meson test -C build-eessi --suite eon
```

Those ranks, on 2026-10-02, for `develop` commit `6f0199ee4`, gave 57 tests: 54 passed, 3 skipped and 0 failed, in 84 s on 8 cores. The three skips are `test_rgpot_mpi_fault`, `test_rgpot_mpi_abort` and `test_rgpot_mpi_params`. Each exits 77 when the CPMD engine library `libcpmdc` is not installed. Later commits on `develop` add tests, so a newer checkout can report a higher pass count.

On 2026-11-11 the `buildenv` route above ran on an AMD Rome node of Snellius with 16 cores, for `develop` commit `0bb29da9e`. The suite ran 80 tests: 74 passed, 5 skipped, and `test_bench_parity` hit the 30 s default timeout. Alone, that test takes 22 s to 27 s on the same node. A parallel run on a loaded node pushes it past 30 s, so it now has a 120 s timeout. With that timeout the same run gave 75 passed, 5 skipped and 0 failed, in 5 min 26 s. A normal account leaves the two `OMPI_ALLOW_RUN_AS_ROOT` variables unset.
