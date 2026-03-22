---
myst:
  html_meta:
    "description": "How eOn handles parallel force evaluation with thread-safe and per-image potential instances."
    "keywords": "eOn, parallel, threading, NEB, potential, MetatomicPotential, XTB"
---

# Parallel Force Evaluation

eOn supports parallel force evaluation in NEB, Dimer/ImprovedDimer, and
ProcessSearchJob. The threading model uses `std::thread` with per-image
potential ownership.

## Threading Model

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
image via `makePotential()`. This is needed for potentials where:

- The same instance cannot be called concurrently (internal state, caches)
- But separate instances CAN run in parallel (each has its own state)

Examples:
- **MetatomicPotential**: PyTorch model has internal caches. Same instance
  needs a mutex; separate instances run truly in parallel. Returns
  `needsPerImageInstance() = true`.
- **XTBPot**: Fortran library has per-instance state (`xtb_TEnvironment`,
  `xtb_TCalculator`). Same instance is not thread-safe; separate instances
  are independent. Returns `needsPerImageInstance() = true`.

When `needsPerImageInstance()` is `true`, NEB creates N+2 potential
instances (one per image) at construction time. The parallel force
evaluation then proceeds lock-free.

## Decision Table

| `isThreadSafe()` | `needsPerImageInstance()` | Behavior | Examples |
|:-:|:-:|:--|:--|
| `true` | `false` | Shared instance, parallel threads | LJ, Morse, SW, EMT |
| `false` | `true` | Per-image instances, parallel threads | XTB, ASE, metatomic |
| `true` | `true` | Per-image instances, parallel threads | MetatomicPotential (mutex fallback) |
| `false` | `false` | Sequential evaluation | (none currently) |

The parallel check in NEB is:
```cpp
bool canParallel = pot->isThreadSafe() || perImagePotentials_;
if (numImages > 1 && params.main_options.parallel && canParallel) { ... }
```

## Affected Code Paths

| Component | Parallel Units | Per-Image Potential |
|:--|:--|:--|
| NEB | N images | Each `path[i]` Matter |
| Dimer | center + forward | `matterDimer` |
| ImprovedDimer | x0 + x1 | `x1` |
| ProcessSearchJob | min1 + min2 | `min2` |

## Performance

With the Morse empirical potential (337-atom Pt, 5 NEB images), parallel
force evaluation gives a **2.3x speedup** over SVN sequential.

With PET-MAD-S ML potential (14-atom Claisen, 10 NEB images):
- Mutex-serialized (shared instance): 192 seconds
- Per-image instances (true parallel): 69 seconds (**2.8x speedup**)

## Adding a New Potential

If your potential has internal state that prevents concurrent calls on the
same instance but supports independent instances:

1. Override `isThreadSafe()` to return `true` (with internal mutex as
   fallback) or `false`
2. Override `needsPerImageInstance()` to return `true`
3. Ensure the constructor (called by `makePotential()`) creates an
   independent instance (no shared static state)

The `[Main] parallel = true` config option (default) enables threading.
Set `parallel = false` to force sequential evaluation.
