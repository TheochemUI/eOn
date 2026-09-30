# European environment for scientific software installations (EESSI)

`eOn-devel-2026.06-GCCcore-15.2.0.eb` is an EasyBuild bundle. The module loads
Meson, Ninja, CMake, pkgconf, Eigen, Python, and Rust from EESSI 2026.06.
CapnProto 1.4.0 and `foss/2026.1` are separate loads. quill, inih,
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

The same prefix then receives CapnProto 1.4.0. The direct rgpot arm needs
the `capnp` program and `capnp-rpc`. EasyBuild fetches that module from
easybuild-easyconfigs pull request 26480. The filename keeps the other
recipes in that pull request unused.

```bash
eb --umask=022 --from-pr 26480 CapnProto-1.4.0-GCCcore-15.2.0.eb
```

A full 40-character commit id works in place of the pull request number.
This id is the pull request head. `--from-commit` stays on it.

```bash
eb --umask=022 --from-commit 66fa89934f0476cd4f9ff14154ee4c87ae5c5d82 CapnProto-1.4.0-GCCcore-15.2.0.eb
```

If that commit is rejected because toolchain GCCcore 15.2.0 is unsupported
and the supported list is empty, export the variable that message prints,
then rerun `eb`:

```bash
export EESSI_SITE_TOP_LEVEL_TOOLCHAINS_2026_06='[{"name": "GCCcore", "version": "15.2.0"}]'
```

In each new shell, from the `develop` checkout:

```bash
source /cvmfs/software.eessi.io/versions/2026.06/init/bash
module load EESSI-extend
module load foss/2026.1 eOn-devel/2026.06-GCCcore-15.2.0 CapnProto/1.4.0-GCCcore-15.2.0
hash -r
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

- `CC`, `CXX` and `FC` to the GCC compilers, so Meson does not wrap them in a
  `ccache` or `sccache` it finds on `PATH`;
- `CARGO_HTTP_CAINFO` follows EESSI's `CURL_CA_BUNDLE` when that is set.
  cargo ignores `CURL_CA_BUNDLE`. The default CA path on Red Hat Enterprise
  Linux (RHEL) does not exist, so readcon-core's crates fail to download
  without this variable.

Two host settings the module cannot fix:

- cargo reads `.cargo/config.toml` from every parent of the build directory.
  A linker or `rustc-wrapper` set in `~/.cargo/config.toml` then applies to
  readcon-core's build. Build outside `$HOME`, or set `CARGO_HOME` to an empty
  directory.
- A git `url.<base>.insteadOf` rule that rewrites `https://github.com/` to SSH
  applies to Meson's wrap downloads as well.
