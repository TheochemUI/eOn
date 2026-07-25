---
myst:
  html_meta:
    "description": "Neighbor lists in eOn: vesin as the single geometry backend; performance via ASV eonclient benches."
    "keywords": "eOn, vesin, neighbor list, ASV, Morse, LJ, QSC"
---

# Neighbor lists (vesin)

eOn uses **[vesin](https://luthaf.fr/vesin/)** as the neighbor-list engine
wherever the code owns pair finding itself. Do not reintroduce pot-local
Verlet tables or O(N²) MIC loops for new work.

## Layers

| Layer | API | Backend |
|-------|-----|---------|
| Python server geometry | `eon.geometry.neighbors.neighbor_list` | `vesin.NeighborList` |
| Process-atom shells / displace | `eon.atoms` / `get_process_atoms` | same |
| Classical C++ pair pots (LJ, Morse, LJCluster, QSC) | `eonc::PairListCache` (Verlet-skin cache over vesin) | `vesin_neighbors` (C) |
| Metatomic | pot-local call into vesin C | same `libvesin` / vendored TU |
| External engines (LAMMPS, VASP, ASE, …) | engine-owned | not eOn's NL |
| Legacy Fortran pots (SW, FeHe, …) | still pot-local | migrate via [vesin Fortran](https://luthaf.fr/vesin/) (`module vesin`) |

## Python

```python
from eon.geometry import neighbor_list
nl = neighbor_list(structure, cutoff=4.0)
```

`brute=True` is API-compat only; the algorithm is always vesin.

## C++

Pair potentials go through the Verlet-skin cache, not raw `VesinNeighbors`:

```cpp
#include "eon/VesinNeighbors.h"
eonc::CachedPairList::Options opt;
opt.cutoff = 5.0; // skin defaults to 1.0 Angstrom
eonc::PairListCache::global().ensureVisit(
    R, nAtoms, box9, opt,
    [&](int32_t i, int32_t j, double dx, double dy, double dz, double r2) {
      // d = r_i - r_j (minimum image applied), r2 = |d|^2 <= cutoff^2
    });
```

The candidate list is built at `cutoff + skin`; while every atom stays
within `skin/2` of its build position the cached pairs are re-used and the
evaluation derives exact vectors from the *current* positions, filtered at
the true cutoff — results match a fresh build bit-for-bit in pair content.
On a slot miss in the MIC regime (orthorhombic box, cutoff within half the
smallest periodic width) the kernel inlines into the build's single
brute-force scan (`vesin_visit.hpp`); everything else builds through
vesin's cell list and evaluates with the stored cell shifts. The pool
(`PairListCache::global()`) matches slots by geometry proximity, so several
NEB images sharing one pot instance each keep a live list on any
thread-to-image assignment — including NEB's per-iteration worker threads,
which would destroy any thread_local cache. `ensureVisit` returns the slot
for further passes over the same list (`forEach`), e.g. QSC's force pass.

`forEach` hands out the historical eOn convention **`d = r_i − r_j`** (the
negated vesin vector). Raw `eonc::VesinNeighbors` (vesin convention
`r_ij = r_j − r_i + S·H`) remains for consumers that need vesin's own
buffers, e.g. Metatomic. `eoncbase` always links vesin (pkg-config / system
/ vendored).

**Do not** call `vesin_free` (or destroy a stack-local list) before every
`compute`. Upstream documents that the same `VesinNeighborList` should be
re-used across calls so allocations are recycled; free only when finished.

## How we measure (ASV, not one-off scripts)

Institutional timing is the **existing ASV suite** under `benchmarks/`:

| ASV class | Fixture | Hits vesin-backed pot? |
|-----------|---------|------------------------|
| `TimePointMorsePt` | `data/point_morse_pt` | Morse force pairs |
| `TimeMinimizationLJCluster` | `data/min_lj_cluster` | LJCluster force pairs (many force evals) |
| `TimeSaddleSearchMorseDimer` | `data/one_pt_saddle_search` | Morse + saddle path |
| `TimeNEBMorsePt` | `data/neb_morse_pt` | Morse band forces |

Mechanism (already in CI):

1. **PR** — `Benchmark PR` workflow builds `eonclient` on **base SHA** and
   **head SHA** with the **PR’s** `benchmarks/` tree, runs
   `asv run --set-commit-hash … --quick`, then asv-perch comments the
   main→PR comparison.
2. **main** — `ASV dashboard` restores orphan branch `asv-results`, runs
   history, publishes with asv-tachyon →
   [eondocs.org/bench](https://eondocs.org/bench/) /
   [bench.eondocs.org](https://bench.eondocs.org/).

Local (same as `benchmarks/README.md`):

```bash
pixi run meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype release
pixi run meson install -C bbdir
pixi run bash -c 'pip install asv asv-tachyon'
pixi run asv machine --yes
pixi run asv run -E "existing:$(which python)" --quick
```

To compare a vesin-NL branch to `main` without waiting for GHA: build+install
both commits’ clients (or sequential install) and use the same ASV class names;
the PR workflow is the authoritative pair-wise run.

## Not yet unified

Fortran Stillinger–Weber / FeHe / Aluminum / EDIP still build neighbors
inside `.f`/`.f90`. Upstream vesin provides `fortran/src/vesin.f90`; wire under
`with_fortran` and port pot generators next.
