# Con File I/O via readcon-core

Starting from v2.12, eOn uses the [readcon-core](https://github.com/lode-org/readcon-core)
Rust library for all `.con` and `.convel` file I/O in both the C++ client and the
Python server.

## C++ Client

### Reading

```cpp
auto frame = readcon::read_first_frame(filename);
matter.con2matter(frame);
```

`read_first_frame` uses memory-mapped I/O for efficient parsing. The returned
`ConFrame` is then used to populate the `Matter` object via the
`con2matter(const readcon::ConFrame&)` overload.

### Writing

```cpp
readcon::ConFrameBuilder builder(lengths, angles, prebox_header, postbox_header);
for (long i = 0; i < nAtoms; i++) {
    builder.add_atom(symbol, x, y, z, is_fixed, atom_id, mass);
}
auto frame = builder.build();
readcon::ConFrameWriter writer(filename, precision);
writer.extend({frame});
```

Position precision is 17 digits (matching the previous `%22.17f` format).
Velocity files use 6-digit precision.

### Multi-frame (append)

For NEB path output (multiple images in one `.con` file), frames are appended:

```cpp
for (long i = 0; i <= neb->numImages + 1; i++) {
    neb->path[i]->matter2con(nebFilename, /*append=*/i > 0);
}
```

The `append=true` mode reads existing frames, appends the new one, and rewrites.

## Python Server

The `eon/fileio.py` module uses the `readcon` PyPI package:

```python
import readcon

frames = readcon.read_con(filename)      # read all frames
frame = readcon.read_first_frame(filename)  # single frame
readcon.write_con(filename, [frame])     # write
```

## Build Requirements

- **Rust** >= 1.88 (provided by pixi via conda-forge)
- **cbindgen** >= 0.29 (installed via `cargo install cbindgen`)
- **Meson** >= 1.8.0

The `readcon-core` library is fetched as a Meson subproject via
`subprojects/readcon-core.wrap` and linked statically.
