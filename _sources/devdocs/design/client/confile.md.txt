# Con File I/O via readcon-core

Starting from v2.12, eOn uses the [readcon-core](https://github.com/lode-org/readcon-core)
Rust library for `.con` and `.convel` file I/O in the C++ client, with the
Python side using the companion `readcon` package.

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

Movie-like outputs can also carry structured frame metadata:

```cpp
eonc::io::ConFrameMetadata metadata;
metadata.frame_index = 7;
metadata.energy = matter.getPotentialEnergy();
metadata.scalars.push_back({"step_size", step_size});
metadata.scalars.push_back({"convergence", convergence});
matter.matter2con("minimization.con", /*append=*/true, &metadata);
```

This maps onto `readcon-core` metadata fields such as `energy`,
`frame_index`, `neb_bead`, `neb_band`, arbitrary scalar metadata, string
metadata, and an escape hatch for raw JSON.

### Multi-frame (append)

For multi-frame output, eOn currently uses correctness-first append semantics:

```cpp
if (!path[i]->matter2con(filename, /*append=*/i > 0, &metadata)) {
    throw std::runtime_error("Failed to append frame");
}
```

`append=true` reads the existing frames, appends the new frame in memory, and
rewrites the file. If the target file exists but cannot be parsed, the append
fails instead of silently truncating the history.

NEB path writers use a small helper wrapper so that per-image metadata stays
consistent across `neb.con`, `neb_path_*.con`, and related outputs:

```cpp
eonc::neb::writePathCon(path, tangent, eigenmode_solvers, numImages,
                        estimateEigenvalues, "neb.con");
```

That helper currently stores per-frame fields such as `neb_bead`, optional
`neb_band`, `reaction_coordinate`, `relative_energy`, and `parallel_force`.

## Python Server

The `eon/fileio.py` module uses the `readcon` PyPI package:

```python
import readcon

frames = readcon.read_con(filename)      # read all frames
frame = readcon.read_first_frame(filename)  # single frame
readcon.write_con(filename, [frame])     # write
```

At the time of writing, the Python orchestration layer is still primarily a
consumer of `.con` data; the richer frame metadata is currently produced on the
C++ side for trajectory outputs and downstream visualization tooling.

## Build Requirements

- **Rust** >= 1.88 (provided by pixi via conda-forge)
- **cbindgen** >= 0.29 (installed via `cargo install cbindgen`)
- **Meson** >= 1.8.0

The `readcon-core` library is fetched as a Meson subproject via
`subprojects/readcon-core.wrap` and linked statically.
