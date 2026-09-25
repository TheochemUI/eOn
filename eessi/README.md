# Developing eOn on EESSI

`eOn-devel-2026.06-GCCcore-15.2.0.eb` is an EasyBuild bundle: one module that
loads everything a Meson build of this checkout needs from EESSI 2026.06
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

- `CC`, `CXX` and `FC` to the GCC compilers, so Meson does not wrap them in a
  `ccache` or `sccache` it finds on `PATH`;
- `CARGO_HTTP_CAINFO` to EESSI's `CURL_CA_BUNDLE` when that is set. cargo
  ignores `CURL_CA_BUNDLE`, and on RHEL-family systems its default CA path does
  not exist, so readcon-core's crates fail to download without it.

Two host settings the module cannot fix:

- cargo reads `.cargo/config.toml` from every parent of the build directory.
  A linker or `rustc-wrapper` set in `~/.cargo/config.toml` then applies to
  readcon-core's build. Build outside `$HOME`, or set `CARGO_HOME` to an empty
  directory.
- A git `url.<base>.insteadOf` rule that rewrites `https://github.com/` to SSH
  applies to Meson's wrap downloads as well.
