---
myst:
  html_meta:
    "description": "How eOn reads and writes .con files using the readcon-core Rust library."
    "keywords": "eOn, readcon-core, .con file, Rust FFI, con2matter, matter2con"
---

# .con file I/O

eOn uses [readcon-core](https://github.com/lode-org/readcon-core), a Rust
library with C/C++ FFI, for reading `.con` files.  The library is included as a
Meson subproject via `subprojects/readcon-core.wrap`.

## Reading (`con2matter`)

`Matter::con2matter(std::string)` opens the file through
`readcon::ConFrameIterator` and extracts the first frame.
`Matter::con2matter(const readcon::ConFrame &)` then populates positions, cell,
masses, atomic numbers and fixed-atom flags from the frame's cached data.

## Writing (`matter2con`)

Writing is still handled by `Matter::matter2con` using `FILE *` I/O, because
readcon-core's writer can only serialise opaque `RKRConFrame` handles that
originated from a previous read.  Once readcon-core gains a frame-builder API,
the write path can be migrated as well.

## Build requirements

readcon-core requires **Rust >= 1.88** and **cbindgen >= 0.29**.  The Rust
toolchain is provided by `conda-forge` (`rust` package); cbindgen is installed
via `cargo install cbindgen` if not already present (see the `ensure_cbindgen`
pixi task).

```{versionadded} 2.11
Replaced the hand-written C `sscanf`/`fgets` parser with readcon-core.
```
