---
myst:
  html_meta:
    "description": "How eOn handles parallel force evaluation with thread-safe and per-image potential instances."
    "keywords": "eOn, parallel, threading, NEB, potential, MetatomicPotential, XTB"
---

# Parallel force evaluation

eOn supports parallel force evaluation in NEB, Dimer/ImprovedDimer, and
ProcessSearchJob. The threading model uses `std::thread` with per-image
potential ownership.

## Thread-safe interface

`Potential` and `Matter` are not internally synchronized except where a
mutex is documented (Metatomic inference lock). The contract is:

- After `makePotential`, the pot is used from one thread unless
  `isThreadSafe()` is true.
- `Matter` is not shared across threads. Parallel NEB uses one Matter
  (or one pot instance) per image.
- `layoutFlags()` says whether the pot is in-process, needs cwd, or is a
  subprocess. Do not share a `NeedsWorkingDirectory` pot across threads
  that would race on cwd.

## Threading model

Callers should use the nonmember queries in `eon/PotCapabilities.h`
(`potIsThreadSafe`, `potAllowsSharedInstance`) rather than naming the
virtuals at every NEB/dimer site. The virtuals stay the override point
on each backend; the free functions are a C++20 concept-constrained
adapter so a later trait cut does not touch those TUs. The header is
header-only (no extra object file).

Two virtual methods on `Potential` control the behavior:

### `isThreadSafe()`

Returns whether the *same* potential instance can be called from multiple
threads concurrently. Most empirical potentials (LJ, Morse, SW, etc.)
return `true` (default). Python-based potentials (ASE, CatLearn) return
`false` because CPython has the GIL.

When `true`, NEB spawns one thread per image and all threads call
`force()` on the shared potential instance.

### `needsPerImageInstance()`

Returns whether NEB should create a *separate* `Potential` instance per
image via `makePotential()`. Needed for potentials where:

- The same instance cannot be called concurrently (internal state, caches)
- But separate instances CAN run in parallel (each has its own state)

Examples:
- **MetatomicPotential**: PyTorch model has internal caches. Same instance
  needs a mutex; separate instances run independently. Returns
  `needsPerImageInstance() = true`.
- **XTB**: `RgpotAdapter` over `rgpot::XTBPot`. The kernel's caps are
  per-instance, so `isThreadSafe()` is false and `needsPerImageInstance()`
  is true. Each image owns its own xtb handles.

When `needsPerImageInstance()` is `true`, NEB gives each interior image
its own instance at construction time, whether or not `[Main] parallel`
is set; the endpoints keep the job's instance. A kernel that keeps a
neighbour table then follows one image instead of rebuilding for every
image of every step: the 601-atom EAM Al band runs 1.68 times faster on
one core with the same iterations. `[Main] parallel` only decides whether
the images run on the image pool, and the threaded evaluation proceeds
lock-free.

Every per-image instance comes from `eonc::helpers::makeImagePotential`,
which clones the job's potential or builds a new one from the parameters
inside a `PotentialConstructionScope` on the job's registry. The force
calls of every image are then part of `total_force_calls`.

## Decision table

| `isThreadSafe()` | `needsPerImageInstance()` | Behavior | Examples |
|:-:|:-:|:--|:--|
| `true` | `false` | Shared instance, parallel threads | LJ, Morse, LJCluster, EMT |
| `false` | `true` | Per-image instances, parallel threads | ASE, CatLearn |
| `true` | `true` | Per-image instances, parallel threads | MetatomicPotential (mutex fallback) |
| `false` | `false` | Sequential evaluation | XTB (`restart.f90` units), SW, EDIP, Lenosky, Tersoff, EAM-Al, FeHe, CuH2, TIP4P-H |

### Where the answers come from

Potentials that reach eOn through `RgpotAdapter` take both flags from
`caps().reentrancy` on the rgpot kernel, so the answer lives next to the
physics.
`SharedInstance` maps to `isThreadSafe()`, `PerInstance` to
`needsPerImageInstance()`, and `ProcessSerial` to neither, which is the
last row above.

The Fortran kernels all declare `ProcessSerial` today. That matches how
they were treated before the port -- they were on the old hard-coded
thread-safety blacklist -- and it is deliberately the conservative
starting point. The kernels
carry no mutable module state any more, so promoting one is a per-kernel
exercise in checking the kernel and adding a test.

The parallel check in NEB is:
```cpp
bool canParallel = pot->isThreadSafe() || perImagePotentials_;
if (numImages > 1 && params.main_options.parallel && canParallel) { ... }
```

The rgsaddle band (`opt_method` `lbfgs`, `fire` or `xtsci` in a
`-Dwith_rgsaddle=true` build) asks the client for the images it carries
through one surface callback. The callback hands them to
`eonc::evaluateTogether`, which with `[Main] parallel` runs them on the
image pool when no two of them share an instance that is unsafe to share.
A potential that batches still takes them in one `forceBatch`. Each image
computes exactly what it computes in a serial run, so energies, forces
and iteration counts do not depend on the number of threads.

## Affected code paths

| Component | Parallel Units | Per-Image Potential |
|:--|:--|:--|
| NEB (host optimizers and the rgsaddle band) | N images | Each interior `path[i]` Matter |
| Dimer | center + forward | `matterDimer` |
| ImprovedDimer | x0 + x1 | `x1` |
| ProcessSearchJob | min1 + min2 | `min2` |

## Performance

With the Morse empirical potential (337-atom Pt, 5 NEB images), parallel
force evaluation gives a **2.3x speedup** over SVN sequential.

With PET-MAD-S ML potential (14-atom Claisen, 10 NEB images):
- Mutex-serialized (shared instance): 192 seconds
- Per-image instances (lock-free): 69 seconds (**2.8x speedup**)

## Adding a New Potential

If your potential has internal state that prevents concurrent calls on the
same instance but supports independent instances:

1. Override `isThreadSafe()` to return `true` (with internal mutex as
   fallback) or `false`
2. Override `needsPerImageInstance()` to return `true`
3. The constructor (called by `makePotential()`) must create an
   independent instance (no shared static state)

The `[Main] parallel = true` config option (default) enables threading.
Set `parallel = false` to force sequential evaluation.
