# eOn 3.2.0

An audit release. Of the 331 commits since `v3.1.0`, 170 are fixes to the
client and server IO paths, the potential backends, and the Python bindings.
Three changes need a decision from anyone upgrading; the rest are repairs.

See monorepo `CHANGELOG.md` section 3.2.0 for the full fragment list.

## Compatibility: eOn writes CON spec 3

The C++ wrap and the Python `readcon` pin both move to **0.14.5**, and the
0.14 writer emits `con_spec_version = 3`. A reader linked against readcon
0.13.x rejects those files with `UnsupportedSpecVersion`.

What this means in practice:

- A `.con` written by eOn 3.2.0 is not readable by an eOn 3.1.x client, nor by
  a downstream tool pinned to `readcon>=0.13,<0.14`.
- Pin downstream consumers to `readcon>=0.14.5` in the same change that moves
  them to eOn 3.2.0. `chemparseplot` and anything reading trajectory output
  through it are the ones to check.
- Reading older spec-2 files is unaffected; the change is on the write side.

## Constraints are stored per Cartesian axis

`Matter` holds the CON column-4 bitmask as an `N x 3` matrix rather than one
flag per atom, and `getFree()` derives the optimizer mask per axis instead of
broadcasting. A slab or reaction coordinate fixed along one axis keeps the two
degrees of freedom it declared free.

This depends on readcon 0.14.5: the 0.13 decoder read column-4 value 1 as
fully fixed, so an x-only constraint could not survive a write and read
(readcon-core #24, fixed in #25). eOn 3.2.0 round-trips every mask value.

`getFixed(atom)` keeps its whole-atom meaning and now answers 1 only when all
three axes are fixed. Code branching on it for a partially constrained atom
sees a different answer than it did in 3.1.x.

## Movie writers no longer grow quadratically

Appending a frame to a `.con` serializes that frame and concatenates it
instead of reading and rewriting the whole file. An N-frame trajectory costs N
frame writes rather than N(N+1)/2. `dynamics.con`, `movie.con` and the basin
hopping trial movies are the visible beneficiaries. Output bytes are
unchanged. Gzip and zstd targets keep the rewrite path, because a compressed
member cannot be extended in place.

## External-program potentials check their IO

`AMS`, `AMS_IO`, `ExtPot` and `VASP` now check every file open, read and shell
command. A missing, truncated or stale result file used to leave the force and
energy arrays holding whatever they held before, and those values went into
the optimizer unremarked.

Two behavioural corrections in the same area:

- `ExtPot` and `VASP` no longer declare themselves safe to call from several
  threads on one instance. Both exchange structures through files at fixed
  names, so parallel NEB images were overwriting each other's input and
  reading each other's forces. `ExtPot` asks for an instance per image and
  stays parallel; `VASP` drives one process through one set of files, so its
  image evaluation is now sequential.
- The `VASP` constructor no longer deletes the working directory's contents.
  It removed `WAVECAR`, `CHGCAR` and `TMPCAR` among others, which is precisely
  what an `ISTART`/`ICHARG` restart reads, and it never removed the `STOPCAR`
  that aborts the next run.

## Batched force evaluation

`RgpotAdapter` reports `PotCaps::batched` through `supportsBatchEvaluation()`
and forwards a batch to rgpot's `forceBatchImpl`, so a kernel that evaluates
many systems in one pass gets the chance to. Kernels that do not override it
take a per-system loop, so every existing potential answers a batch call
correctly. Requires rgpot carrying `ForceBatch`.

## Tests run on Unix CI

`tests/` runs after the client install on Linux and macOS. It had not run in
CI before, which is why several of the fixes in this release are for defects
that existed for some time: `test_one_pt` died at collection on a shell
builtin, and `min_mode_method = gprdimer` in a build without `-Dwith_gprd`
silently fell through to ImprovedDimer and reported a nonnegative-mode abort
rather than refusing.

## Upgrade checklist

1. Move `readcon` to `>=0.14.5` wherever a downstream tool reads eOn output.
2. Re-check any code branching on `Matter::getFixed(atom)` for partially
   constrained atoms.
3. If you drive `ExtPot` from parallel NEB, note that each image now runs in
   `extpot_<pid>_<n>`; a wrapper reaching for the client's working directory
   should read `EON_EXTPOT_RUN_DIR`.
