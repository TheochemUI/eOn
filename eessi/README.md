# Developing eOn on the European environment for scientific software installations (EESSI)

`eOn-devel-2026.06-GCCcore-15.2.0.eb` is an EasyBuild bundle: one module that
loads everything a Meson build of this checkout needs from the European
environment for scientific software installations (EESSI) 2026.06
(Meson, Ninja, CMake, pkgconf, Eigen, Python, Rust). quill, inih,
nlohmann_json, Highway, rgpot and readcon-core come from the wraps in
`subprojects/`, so the first `meson setup` needs network.

Install it once into your EESSI-extend prefix:

```bash
source /cvmfs/software.eessi.io/versions/2026.06/init/bash
module load EESSI-extend
eb eessi/eOn-devel-2026.06-GCCcore-15.2.0.eb
```

Then, in each new shell:

```bash
source /cvmfs/software.eessi.io/versions/2026.06/init/bash
module load EESSI-extend eOn-devel/2026.06-GCCcore-15.2.0
meson setup bbdir -Dwith_tests=true
meson compile -C bbdir
meson test -C bbdir --suite eon
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

- cargo reads `.cargo/config.toml` from every parent of the build directory.
  A linker or `rustc-wrapper` set in `~/.cargo/config.toml` then applies to
  readcon-core's build. Build outside `$HOME`, or set `CARGO_HOME` to an empty
  directory.
- A git `url.<base>.insteadOf` rule that rewrites `https://github.com/` to SSH
  applies to Meson's wrap downloads as well.
