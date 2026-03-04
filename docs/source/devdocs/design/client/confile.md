---
myst:
  html_meta:
    "description": "How eOn reads and writes .con files using the readcon-core Rust library."
    "keywords": "eOn, readcon-core, .con file, Rust FFI, con2matter, matter2con"
---

# .con file I/O

eOn uses [readcon-core](https://github.com/lode-org/readcon-core), a Rust
library with C/C++ FFI and Python bindings, for all `.con` file reading and
writing.  The library is included as a Meson subproject via
`subprojects/readcon-core.wrap` for the C++ client, and as a PyPI dependency
(`readcon`) for the Python server.

## C++ client

### Reading (`con2matter`)

`Matter::con2matter(std::string)` calls `readcon::read_first_frame()` (mmap
based) and passes the result to `Matter::con2matter(const readcon::ConFrame &)`
which populates positions, cell, masses, atomic numbers and fixed-atom flags.

### Writing (`matter2con`)

`Matter::matter2con(std::string)` builds a `readcon::ConFrameBuilder` from the
current Matter state, calls `build()`, and writes the resulting frame via
`readcon::ConFrameWriter` with precision 17 (matching the previous `%22.17f`
format for lossless floating-point roundtrip).

### Velocity files (`convel2matter` / `matter2convel`)

The same readcon-core reader and builder handle `.convel` files, which extend
`.con` with a velocity section per atom type.

## Python server

`eon/fileio.py` uses the `readcon` Python package:

- `loadcon()` / `loadcons()` call `readcon.read_con()` or
  `readcon.read_con_string()` and convert frames to eOn `Atoms` objects.
- `savecon()` converts `Atoms` to `readcon.ConFrame` and writes with
  `readcon.write_con()` or `readcon.write_con_string()` (for StringIO).

## Build requirements

readcon-core requires **Rust >= 1.88** and **cbindgen >= 0.29**.  The Rust
toolchain is provided by `conda-forge` (`rust` package); cbindgen is installed
via `cargo install cbindgen` if not already present (see the `ensure_cbindgen`
pixi task).

```{versionadded} 2.11
Replaced the hand-written C `sscanf`/`fgets` parser with readcon-core for
reading.
```

```{versionchanged} 2.12
Replaced `FILE *`/`fprintf` writing with readcon-core `ConFrameBuilder` +
`ConFrameWriter` for both C++ and Python.  The Python server now uses the
`readcon` PyPI package instead of custom parsing in `fileio.py`.
```
