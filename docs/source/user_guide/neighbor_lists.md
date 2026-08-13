---
myst:
  html_meta:
    "description": "Neighbor lists in eOn: vesin as the geometry backend; performance via the ASV eonclient suite."
    "keywords": "eOn, vesin, neighbor list, ASV, Morse, LJ, rgpot"
---

# Neighbor lists (vesin)

eOn uses [vesin](https://luthaf.fr/vesin/) for pair finding whenever the code
owns that step. New work should not add pot-local Verlet tables or O(N²) MIC
loops.

## Layers

| Layer | API | Backend |
|-------|-----|---------|
| Python server geometry | `eon.geometry.neighbors.neighbor_list` | `vesin.NeighborList` |
| Process-atom shells / displace | `eon.atoms` / `get_process_atoms` | same |
| Classical C++ pair pots (LJ, Morse, LJCluster, ZBL) | `rgpot::nlist::PairListCache` (Verlet-skin cache over vesin, in [rgpot](https://github.com/OmniPotentRPC/rgpot)) | `vesin_neighbors` (C) |
| Metatomic | pot-local call into vesin C | same `libvesin` / vendored TU |
| External engines (LAMMPS, VASP, ASE, …) | engine-owned | not eOn's NL |
| Fortran pots (SW, EDIP, Lenosky, EAM-Al, CuH2, FeHe, Tersoff, TIP4P-H) | `rgpot_neighbors`, CSR full-list table in rgpot | `vesin_fortran`, vendored in rgpot |

## Python

```python
from eon.geometry import neighbor_list
nl = neighbor_list(structure, cutoff=4.0)
```

`brute=True` keeps API compatibility only. The algorithm is always vesin.

## C++

Classical pair pots (LJ, LJCluster, Morse, ZBL) live in **rgpot** and reach
vesin through `rgpot::nlist::PairListCache`. eOn consumes them through
`RgpotAdapter`; the cache is not part of eOn's public surface. See the rgpot
docs for `Options`, `evaluate`, and `ensureVisit`.

The candidate list is built at `cutoff + skin`. While every atom stays within
`skin/2` of its build position, the cache reuses pairs and derives exact
vectors from current positions, filtered at the true cutoff. Pair content
matches a fresh build. The pool matches slots by geometry proximity, so
several NEB images that share one pot instance each keep a live list under
any thread-to-image assignment, including NEB's per-iteration worker threads.
A `thread_local` cache would not survive that assignment.

`eonc::VesinNeighbors` (vesin convention `r_ij = r_j − r_i + S·H`) is the RAII
wrapper for consumers that need vesin's own buffers:

```cpp
#include "eon/VesinNeighbors.h"
eonc::VesinNeighbors nl;
eonc::VesinNeighbors::Options opt;
opt.cutoff = 5.0;
nl.compute(R, nAtoms, box9, opt);
```

`eoncbase` always links vesin. Resolution order: an installed vesin
(pkg-config, cmake, or plain library); else the copy the rgpot subproject
already builds; else eOn's own vendored translation unit. When rgpot comes
from the wrap, reusing its vesin avoids two definitions of every vesin symbol
in one link.

Do not call `vesin_free` (or destroy a stack-local list) before every
`compute`. Upstream documents reusing the same `VesinNeighborList` across
calls so allocations recycle; free only when finished.

## ASV benchmarks

Wall-clock timing for client changes goes through the ASV suite under
`benchmarks/`, not ad-hoc scripts.

| ASV class | Fixture | vesin pot path |
|-----------|---------|----------------|
| `TimePointMorsePt` | `data/point_morse_pt` | Morse force pairs |
| `TimeMinimizationLJCluster` | `data/min_lj_cluster` | LJCluster force pairs (many force evals) |
| `TimeSaddleSearchMorseDimer` | `data/one_pt_saddle_search` | Morse + saddle path |
| `TimeNEBMorsePt` | `data/neb_morse_pt` | Morse band forces |

### CI

1. PR: the `Benchmark PR` workflow builds `eonclient` at the base SHA and the
   head SHA with the PR's `benchmarks/` tree, runs
   `asv run --set-commit-hash … --quick`, then asv-perch comments the
   `main`→PR comparison.
2. `main`: the `ASV dashboard` workflow restores the orphan branch
   `asv-results`, runs history, and publishes with asv-tachyon to
   [eondocs.org/bench](https://eondocs.org/bench/) and
   [bench.eondocs.org](https://bench.eondocs.org/).

### Local

Same sequence as `benchmarks/README.md`:

```bash
pixi run meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype release
pixi run meson install -C bbdir
pixi run bash -c 'pip install asv asv-tachyon'
pixi run asv machine --yes
pixi run asv run -E "existing:$(which python)" --quick
```

To compare a vesin-NL branch to `main` without waiting on GHA, build and
install both clients (or install sequentially) and run the same ASV class
names. Pair comparison for review uses the PR workflow.

## Fortran

eOn compiles no Fortran. The kernels and the vendored vesin Fortran interface
both live in rgpot. Every kernel reaches vesin the same way:
`rgpot_neighbors` builds a CSR full neighbour list with pair vectors and
distances (Verlet `skin` on the C API), and the kernels read rows out of it.
The per-pot CSR helper modules that used to sit beside the fixed-form kernels
here (`vesin_al`, `vesin_cuh2`, `vesin_fehe`) moved out with those kernels.

Fortran pots take the box diagonal only (orthorhombic). Porting rules and the
neighbour-table cell convention are in
`docs/orgmode/reference/fortran_potentials.org` in rgpot.
