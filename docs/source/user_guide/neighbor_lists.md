---
myst:
  html_meta:
    "description": "Neighbor lists in eOn: vesin as the single geometry backend; performance via ASV eonclient benches."
    "keywords": "eOn, vesin, neighbor list, ASV, Morse, LJ, rgpot"
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
| Classical C++ pair pots (LJ, Morse, LJCluster, ZBL) | `rgpot::nlist::PairListCache` (Verlet-skin cache over vesin, in the [rgpot](https://github.com/OmniPotentRPC/rgpot) repo) | `vesin_neighbors` (C) |
| Metatomic | pot-local call into vesin C | same `libvesin` / vendored TU |
| External engines (LAMMPS, VASP, ASE, …) | engine-owned | not eOn's NL |
| Fortran pots (SW, EDIP, Lenosky, EAM-Al, CuH2, FeHe, Tersoff, TIP4P-H) | `rgpot_neighbors`, a CSR full-list table (in rgpot) | `vesin_fortran`, vendored in rgpot |

## Python

```python
from eon.geometry import neighbor_list
nl = neighbor_list(structure, cutoff=4.0)
```

`brute=True` is API-compat only; the algorithm is always vesin.

## C++

The classical pair pots (LJ, LJCluster, Morse, ZBL) live in **rgpot** and
reach vesin through `rgpot::nlist::PairListCache`, the Verlet-skin cache in
that repo. eOn consumes them through `RgpotAdapter`; the cache is not part
of eOn's own surface. See the rgpot docs for `Options`, `evaluate`, and
`ensureVisit`.

The candidate list is built at `cutoff + skin`; while every atom stays
within `skin/2` of its build position the cached pairs are re-used and the
evaluation derives exact vectors from the *current* positions, filtered at
the true cutoff — results match a fresh build in pair content. The pool
matches slots by geometry proximity, so several NEB images sharing one pot
instance each keep a live list on any thread-to-image assignment —
including NEB's per-iteration worker threads, which would destroy any
thread_local cache.

`eonc::VesinNeighbors` (vesin convention `r_ij = r_j − r_i + S·H`) remains
in eOn as the RAII wrapper for consumers that need vesin's own buffers:

```cpp
#include "eon/VesinNeighbors.h"
eonc::VesinNeighbors nl;
eonc::VesinNeighbors::Options opt;
opt.cutoff = 5.0;
nl.compute(R, nAtoms, box9, opt);
```

`eoncbase` always links vesin. It resolves in this order: an installed
vesin (pkg-config, cmake, or plain library), else the copy the rgpot
subproject already builds, else eOn's own vendored translation unit. The
middle case matters when rgpot comes from the wrap: building both copies
would put two definitions of every vesin symbol in one link.

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

## Fortran

eOn compiles no Fortran. The kernels and the vendored vesin Fortran
interface both live in rgpot, where every one of them reaches vesin the
same way: `rgpot_neighbors` builds a CSR full neighbour list carrying pair
vectors and distances, with the C API's Verlet `skin` option, and the
kernels read rows out of it. The per-pot CSR helper modules that used to
sit beside the fixed-form kernels here (`vesin_al`, `vesin_cuh2`,
`vesin_fehe`) are gone with them.

Every Fortran pot remains orthorhombic-only (they take the box diagonal).
For the porting rules and the neighbour-table cell convention, see
`docs/orgmode/reference/fortran_potentials.org` in rgpot.
