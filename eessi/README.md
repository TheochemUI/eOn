# European environment for scientific software installations (EESSI)

`eOn-devel-2026.06-GCCcore-15.2.0.eb` is an EasyBuild bundle. The module loads
Meson, Ninja, CMake, pkgconf, Eigen, Python, and Rust from EESSI 2026.06.
It depends on `foss/2026.1`, which loads OpenMPI, FFTW, and OpenBLAS.
The same `eb` command builds Cap'n Proto 1.4.0 into the bundle prefix, so
`pkg-config` finds `capnp-rpc` after `module load eOn-devel`. quill, inih,
nlohmann_json, Highway, rgpot, and readcon-core come from the wraps in
`subprojects/`, so the first `meson setup` needs network.

`docs/source/install/eessi.md` is the annotated sequence, including the test count.
The commands here match that page.

Install the bundle once. Create the prefix first. An unset
`EESSI_USER_INSTALL` means `$HOME/eessi`. A set value must already be a
directory. A user prefix sets `EASYBUILD_UMASK` to `077`, so pass
`--umask=022`. Without that flag, Lmod running as root reports the new
module as unknown.

```bash
source /cvmfs/software.eessi.io/versions/2026.06/init/bash
mkdir -p "${EESSI_USER_INSTALL:-$HOME/eessi}"
module load EESSI-extend
eb --umask=022 eessi/eOn-devel-2026.06-GCCcore-15.2.0.eb
```

In each new shell, from the `develop` checkout:

```bash
source /cvmfs/software.eessi.io/versions/2026.06/init/bash
module load EESSI-extend
module load eOn-devel/2026.06-GCCcore-15.2.0
hash -r
command -v mpicc
pkg-config --exists capnp-rpc
unset RUSTC_WRAPPER CARGO_BUILD_RUSTC_WRAPPER RUSTC_WORKSPACE_WRAPPER
unset RUSTFLAGS CARGO_ENCODED_RUSTFLAGS
export CARGO_HOME="$PWD/../eon-cargo-home"
mkdir -p "$CARGO_HOME"
EON_LDFLAGS=""
for v in $(env | grep -o '^EBROOT[A-Z0-9_]*'); do
  for l in lib lib64; do
    if [ -d "${!v}/$l" ]; then
      EON_LDFLAGS="$EON_LDFLAGS -Wl,-rpath,${!v}/$l"
    fi
  done
done
export LDFLAGS="$EON_LDFLAGS"
export CC=mpicc CXX=mpicxx FC=mpif90
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc
meson setup build-eessi \
  --prefix "$PWD/../eon-prefix" \
  --libdir lib \
  --buildtype release \
  -Dwith_tests=true \
  -Dpython.install_env=prefix \
  -Drgpot:with_mpi=enabled
meson compile -C build-eessi
mkdir -p "$PWD/../eon-tmp"
export TMPDIR="$PWD/../eon-tmp"
export OMP_NUM_THREADS=1
export OMPI_MCA_rmaps_base_oversubscribe=1
export PMIX_MCA_psec=none
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
meson test -C build-eessi --suite eon
```

The module also sets:

- `CC`, `CXX` and `FC` to the gcc compilers, so Meson does not wrap them in a
  `ccache` or `sccache` it finds on `PATH`;
- `CARGO_HTTP_CAINFO` to the EESSI `CURL_CA_BUNDLE` when that is set. cargo
  ignores `CURL_CA_BUNDLE`. On Red Hat Enterprise Linux (RHEL) the default CA
  path does not exist, so readcon-core's crates fail to download without it;
- `LDFLAGS` to `-Wl,-rpath,$EESSI_SOFTWARE_PATH/software/GCCcore/15.2.0/lib64`
  when `EESSI_SOFTWARE_PATH` is set.

This module's `LDFLAGS` is the GCCcore `lib64` run path above. EESSI 2026.06
leaves `LD_LIBRARY_PATH` unset, and it stays unset after
`module load FFTW/3.3.10-GCC-15.2.0`. That run path selects the GCCcore
`libstdc++`. The compatibility copy stops at `GLIBCXX_3.4.33` (gcc 14). gcc
15.2 needs `GLIBCXX_3.4.34`.

A program that needs `libfftw3.so.3`, linked with that gcc and with `RUNPATH`
set to the GCCcore `lib64` alone, exits 127. The loader reports
`libfftw3.so.3: cannot open shared object file`. The same program with
`software/FFTW/3.3.10-GCC-15.2.0/lib` added to `RUNPATH` exits 0 and prints
`fftw fftw-3.3.10-sse2-avx-avx2-avx2_128`. Each further module library needs
its own `-Wl,-rpath`.

Host `ldd` resolves that binary's `libfftw3.so.3` to the distro copy and
reports `GLIBC_2.34` missing, because it uses the host loader. The exit status
above is from the binary whose interpreter is the EESSI
`ld-linux-x86-64.so.2`.

`buildenv/default-foss-2026.1` is part of this release. It depends on
`foss/2026.1`, adds `rpath_wrappers` entries for `gcc`, `g++`, `gfortran`,
`ld` and `ld.bfd` to `PATH`, and sets `LDFLAGS` to `-L` directories for
ScaLAPACK, the fastest Fourier transform in the west (FFTW), FlexiBLAS,
OpenMPI and GCCcore. Those `-L` flags are the link search path. `RUNPATH` is
the `-Wl,-rpath` list the loader reads. Loading this module after `buildenv`
replaces `LDFLAGS` with the single GCCcore entry. The exit codes above come
from passing `-Wl,-rpath` on the link line.

Two host settings the module cannot fix:

- cargo reads `.cargo/config.toml` from every parent of the build directory
  as well as from `$CARGO_HOME`. A linker or `rustc-wrapper` set in
  `~/.cargo/config.toml` then applies to readcon-core's build for any checkout
  under `$HOME`, whatever `CARGO_HOME` says. Build outside `$HOME`, and set
  `CARGO_HOME` to an empty directory.
- A git `url.<base>.insteadOf` rule that rewrites `https://github.com/` to SSH
  applies to Meson's wrap downloads as well. `GIT_CONFIG_GLOBAL` pointed at an
  empty file drops it for one shell.

## In-process CPMD point job

The bundle load puts the schema compiler on `PATH`. `capnp` answers
`command -v capnp`, and `pkg-config --exists capnp-rpc` returns 0.

`scripts/ci/opencpmd_point.sh` builds `libcpmdc.so` against OpenCPMD commit
`062582b7cfd832d36f88f504cd08e4ead42eb404` and cpmdc commit
`8439c25cf0bd5f4caefd6ed07bebf142e2b181f5`. `KIT` is a checkout of
`TheochemUI/sige_repro`. The script copies `structures/si3n4_isomer1.con`
from that tree and runs the isomer 1 point job. The energy must sit within
1e-6 eV of -1396.269526. The pseudopotential files stay in `KIT`. They are
not in this repository.

`scripts/build-client.sh` checks `eonclient` and `libcpmdc.so` after the
build. Their run path must not contain a `/home/` directory. `--sanitize`
rewrites a run path that does. The stack directory is an argument. It is
not taken from a home scratch tree.

`EONCLIENT` is `build-eessi/client/eonclient` from the configure above.

```bash
module load eOn-devel/2026.06-GCCcore-15.2.0
hash -r
command -v capnp
pkg-config --exists capnp-rpc
export OPENCPMD_COMMIT=062582b7cfd832d36f88f504cd08e4ead42eb404
export CPMDC_COMMIT=8439c25cf0bd5f4caefd6ed07bebf142e2b181f5
export STACK="$PWD/../eon-opencpmd-stack"
export KIT=/path/to/sige_repro
export EONCLIENT="$PWD/build-eessi/client/eonclient"
scripts/ci/opencpmd_point.sh
```
