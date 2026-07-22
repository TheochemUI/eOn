---
myst:
  html_meta:
    "description": "Real atomic systems used in eOn documentation tutorials: Morse Pt NEB, LJ cluster minimization, and how they connect to rgpycrumbs plotting."
    "keywords": "eOn systems, Morse Pt, LJ cluster, NEB tutorial, minimization tutorial, documentation examples"
---

# Real systems for documentation

These tutorials ship **self-contained geometries and configs** that use only
**built-in eOn potentials** (no Hugging Face download, no external DFT). They
are the recommended path for:

- learning client jobs end-to-end (`eonclient` + outputs)
- exercising the current
  [rgpycrumbs](https://rgpycrumbs.rgoswami.me) / [chemparseplot](https://chemparseplot.rgoswami.me)
  visualization stack (`plt-neb`, `plt-min`)
- CI-stable MyST-NB execution (docs build does not depend on metatomic models)

For a metatomic / PET-MAD walkthrough (oxadiazole, vinyl alcohol), see the
[atomistic-cookbook eon-pet-neb](https://atomistic-cookbook.org/examples/eon-pet-neb/eon-pet-neb.html)
recipe and the in-tree {doc}`../visualization` tutorial.

## Systems

| System | Potential | Job | Docs page | Data under `systems/data/` |
|--------|-----------|-----|-----------|----------------------------|
| Pt adatom / surface hop (few-atom Morse) | `morse_pt` | NEB | {doc}`morse_pt_neb` | `morse_pt_neb/` |
| Lennard-Jones cluster | `ljcluster` | minimization | {doc}`lj_minimization` | `lj_min/` |
| Pt saddle (dimer) | `morse_pt` | saddle search | {doc}`pt_saddle` | `pt_saddle/` |

Geometries are copied from the in-repo `benchmarks/data/` and `examples/` trees
so the same inputs power unit/ASV benchmarks and the docs.

## Plotting conventions (current stack)

These conventions match the production `rgpycrumbs` / `chemparseplot` tools:

**NEB (`plt-neb`)**

- Prefer full history: `--landscape-path all`, `--show-pts`, `--highlight-last`
- Structure strip: `--plot-structures all` (every band image) or `crit_points`
  (R / SP / P only for a lighter figure)
- Landscape: `--project-path` so the panel is true **1:1 Å** in reaction-valley
  coordinates (`Δs = Δd`); colorbar and strip sit under a content-sized figure
- SP marker: gold star on the 1D profile (and landscape) when the final band
  and/or `sp.con` are present

**Minimization (`plt-min`)**

- Require movies: `[Debug] write_movies = true` (and optionally
  `write_deprecated_outs = true` for legacy `.dat` sidecars)
- **One 2D landscape per endpoint** — do not overlay reactant and product on a
  shared `(s, d)` frame; each trajectory defines its own RMSD basis
- Use `--label reactant` / `--label product` so titles read
  *Reactant minimization* / *Product minimization* with **initial** /
  **minimized** endpoint captions
- Relative energy on the colorbar (readable tick labels)

```{toctree}
:maxdepth: 1
:caption: Systems

morse_pt_neb
lj_minimization
pt_saddle
```
