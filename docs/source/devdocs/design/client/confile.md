# CON file I/O via readcon-core

Starting from v2.12, eOn uses the [readcon-core](https://github.com/lode-org/readcon-core)
Rust library for `.con` and `.convel` file I/O in the C++ client, while
Python uses the companion `readcon` package.

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
A `.convel` file uses binary64 round-trip text, so a kinetic energy reads
back to a relative 1e-12. `matter2convel` takes a decimal count for the
older fixed layout. Six decimal places lose about 5.8e-7 relative kinetic
energy.

Velocities are stored in `(amu*angstrom^2/eV)^0.5` (10.1805055 fs). A frame
that carries them is written as `con_spec_minor` 1 with that spelling in
`units.time`, plus length angstrom, energy eV, and mass amu. eOn reads
spec versions 1, 2, and 3.0 in that internal unit even when the label says
`fs`. From minor 1, `units.time` is applied. Files from eOn 3.6.0 and
earlier (doi:10.5281/zenodo.23166734) used the natural unit under an `fs`
label and have no minor. readcon 0.17 leaves those values unscaled.
readcon-db `export_h5md` takes `velocity_time_unit` set to
`(amu*angstrom^2/eV)^0.5` when a later file mislabels the same numbers
as `fs`.

Movie-like outputs can also carry structured frame metadata:

```cpp
eonc::io::ConFrameMetadata metadata;
metadata.frame_index = 7;
metadata.energy = matter.getPotentialEnergy();
metadata.scalars.push_back({"step_size", step_size});
metadata.scalars.push_back({"convergence", convergence});
matter.matter2con("minimization.con", /*append=*/true, &metadata);
```

The metadata maps onto `readcon-core` fields such as `energy`,
`frame_index`, `neb_bead`, `neb_band`, arbitrary scalar metadata, string
metadata, and an escape hatch for raw JSON.

### Multi-frame (append)

For multi-frame output, `append=true` extends the file in place:

```cpp
if (!path[i]->matter2con(filename, /*append=*/i > 0, &metadata)) {
    throw std::runtime_error("Failed to append frame");
}
```

The new frame is serialized on its own and its bytes are concatenated onto the
target, so an N-frame trajectory costs N frame writes rather than N(N+1)/2.
`ConFrameWriter` emits self-contained frames with no file-level preamble, which
makes the concatenated result identical to the same frames written in one call.
Serialization goes through a scratch file next to the target because
readcon-core exposes neither a writer over a caller-owned stream nor a flush on
an open writer; the target itself is opened in append mode and flushed before
the call returns, so every intermediate state on disk is a complete multi-frame
`.con` that a reader or a `tail` can parse.

If the target file exists but cannot be parsed, the append fails rather than
truncating the history. That check runs when eOn did not write the file, or
when its size or modification time moved since eOn wrote it; frames this
process wrote and nobody touched since are extended without a re-parse.
`eonc::io::resetConAppendState()` drops that bookkeeping for a file replaced
with same-sized content inside one filesystem timestamp tick.

Gzip and zstd targets keep the read-all-and-rewrite path. A compressed member
cannot be extended in place, and readcon-core reads a single gzip member, so
appended members would be invisible on read-back.

NEB path writers use a small helper wrapper so that per-image metadata stays
consistent across `neb.con`, `neb_path_*.con`, and related outputs:

```cpp
eonc::neb::writePathCon(path, tangent, eigenmode_solvers, numImages,
                        estimateEigenvalues, "neb.con");
```

That helper currently stores per-frame fields such as `neb_bead`, optional
`neb_band`, `reaction_coordinate`, `relative_energy`, and `parallel_force`.

## Python server

The `eon/fileio.py` module uses the `readcon` PyPI package:

```python
import readcon

frames = readcon.read_con(filename)      # read all frames
frame = readcon.read_first_frame(filename)  # single frame
readcon.write_con(filename, [frame])     # write
```

At the time of writing, the Python orchestration layer is still primarily a
consumer of `.con` data; the richer frame metadata is currently produced in
the C++ client for trajectory outputs and downstream visualization tooling.

## Build requirements

- **Rust** >= 1.88 (provided by pixi via conda-forge)
- **cbindgen** >= 0.29 (installed via `cargo install cbindgen`)
- **Meson** >= 1.8.0

The `readcon-core` library is fetched as a Meson subproject via
`subprojects/readcon-core.wrap` and linked statically. The wrap pins
v0.17.0, and `client/meson.build` asks for `>=0.17.0` on both resolution
paths. That release writes `con_spec_minor` and keeps the natural time
spelling `(amu*angstrom^2/eV)^0.5`. `eon/fileio.py` needs the companion
`readcon` package at `>=0.14.5`.

## Updating the subproject

Meson uses an existing `subprojects/readcon-core/` directory unchanged and does
not re-fetch it when the wrap revision changes. After pulling a commit that
moves the wrap, reset the checkout:

```bash
meson subprojects update --reset readcon-core
```

Skipping the reset leaves the build configuring against whatever tag sits on
disk. The version constraint on the `subproject()` call catches that at
configure time and names `readcon-core`, so a version complaint about the
subproject points at a stale checkout rather than at anything in
`client/ConFileIO.cpp`.
