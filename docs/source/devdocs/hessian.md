---
myst:
  html_meta:
    "description": "When the finite-difference Hessian is colored, and how a potential opts in to a Hessian range."
    "keywords": "eOn Hessian, hessianRange, coloring, rgpot, many-body"
---

# Hessian range

The Hessian job is part of eOn
{cite:t}`dhess-goswamiEfficientExplorationChemical2025`.
The older description of that client is Chill et al.
{cite:t}`dhess-chillEONSoftwareLong2014`.

## When coloring applies

`Hessian::calculate` colors the mobile atoms when all of the following hold.

- `Potential::hessianRange()` is finite and strictly positive.
- `[Main] remove_net_force` is off, or some atom is fixed, or the system has one atom. Net-force removal shifts every atom, so a column is no longer confined to the neighborhood.
- `[Hessian] checkpoint_path` is empty. A checkpoint is stored one coordinate at a time.

The reported value is the distance beyond which every Hessian block between two atoms is exactly zero. Mobile atoms are greedy-colored on the square of that neighborhood. One displacement of a color along x, y, or z fills the columns of that color. Entries outside the neighborhood are written as zero, which is exact only for that distance.

Each row of that square is deduplicated while it is filled, with the same epoch-marked visited array the greedy coloring uses. The temporary storage is the edges, at most one slot per mobile pair. Pushing every toucher pair and sorting afterwards holds `O(sum t_k^2)` integers first, which is the whole cell when a many-body range reaches across it.

When the coloring has as many colors as mobile atoms, `calculateColored` returns false and `calculate` takes the batched path, or the serial path when the potential does not batch. A color per atom does not reduce the force count, and the colored loop would skip batching. Catch2 `Coloring falls back when every mobile atom needs its own color` is that case: a four-atom cluster whose diameter is below the range takes 25 force calls, the serial central count. The eight-atom embedding chain takes 49 force calls total on the central scheme.

`0` and any non-finite value mean the range is unknown. The matrix is then one column per coordinate. A potential that evaluates batches takes the batched path instead, still one sample per coordinate.

## How a potential opts in

Override `Potential::hessianRange()` and return that distance. The default is `0`.

An rgpot kernel wrapped by `RgpotAdapter` does not override the adapter. It specializes `eonc::HessianRangeOf<Kernel>` in `include/eon/potentials/RgpotAdapter/HessianRange.h`. The primary template is `HessianLocality::Unknown` and returns `0`. The adapter calls that trait and passes the result through `acceptedHessianRange`, which turns a non-positive or non-finite value into `0`.

Override `Potential::hessianRange()` directly for a potential that is not an rgpot kernel. `finiteCutoff()` is deprecated and stays as the fallback that default `hessianRange()` calls, so an existing override still colors.

The specialization declares a locality class and returns the matching distance.

| Class | Range |
| --- | --- |
| Pair | the interaction cutoff |
| Atom-centred many-body | twice the interaction cutoff |
| Molecular site model | the intermolecular cutoff plus twice the largest intramolecular extent, or `0` when that extent is not a hard bound |
| Message passing | `n_layers` times the edge cutoff when the layer count is in the kernel, otherwise `0` |

Do not infer the range from `config().cutoff`. That field is the neighbour-list cutoff. For a many-body kernel it is too short, and for water and for UMA it is not a Hessian range at all.

A kernel added in rgpot follows [Porting potentials](project:porting_potentials.md) and then adds this specialization in the same change that wires the factory arm. Leaving the primary template in place keeps coloring off.

## Kernels

The classes below were checked at the current rgpot pin and at rgpot release `v3.5.0`. "Wrapped" means `makePotential` builds an `RgpotAdapter` for that kernel. A kernel that exists in rgpot but is not wrapped still has a specialization when its header is present, so a later factory arm cannot turn coloring on by having a cutoff field.

| Kernel | Pin | 3.5.0 | Class | Reported range |
| --- | --- | --- | --- | --- |
| `LJPot` | wrapped | wrapped | pair | `config().cutoff`, default 15 Å |
| `LJClusterPot` | wrapped | wrapped | pair | `config().cutoff`, default 15 Å |
| `MorsePot` | wrapped | wrapped | pair | `config().cutoff`, default 9.5 Å |
| `ZBLPot` | wrapped | wrapped | pair | `config().cut_global`, default 2.5 Å |
| `EAMAlPot` | wrapped | wrapped | atom-centred | `2 * 7.1` Å |
| `SWPot` | wrapped | wrapped | atom-centred | `2 * 1.8 * 2.0951` Å |
| `TersoffPot` | wrapped | wrapped | atom-centred | `2 * 3.0` Å |
| `EDIPPot` | wrapped | wrapped | atom-centred | `2 * 3.1213820` Å |
| `LenoskyPot` | wrapped | wrapped | atom-centred | `2 * 4.5` Å |
| `FeHePot` | wrapped | wrapped | atom-centred | `8.4` Å |
| `CuH2Pot` | wrapped | wrapped | atom-centred | `2 * 6.1` Å |
| `WaterHPot` (TIP4P-H) | factory arm when water is built | same | unknown | disabled |
| `XTBPot` | wrapped when built | same | unknown | disabled |
| `D3Pot` | wrapped when built | same | unknown | disabled |
| `D4Pot` | wrapped when built | same | unknown | disabled |
| `MOPACPot` | wrapped | wrapped | unknown | disabled |
| `ExprPot` | wrapped when built | same | unknown | disabled |
| `UmaPot` | header, not a factory arm | same, and `config().cutoff` | message passing | disabled |
| `TIP4PPot` | not in this pin | header | molecular site | disabled |
| `SPCEPot` | not in this pin | header | molecular site | disabled |
| `TIP4PPtPot` | not in this pin | header | molecular site | disabled |
| `EAMPot` (cell list) | not in this pin | header | atom-centred | `2 * 5.555` Å |
| `EMTPot` | not in this pin | header | atom-centred | disabled |
| `NWChemPot`, `CPMDPot`, `TBLitePot`, `SkalaPot`, `MetatomicPot` | dlopen engines | same | unknown | disabled |

`NWChemPot`, `CPMDPot`, `TBLitePot`, `SkalaPot` and `MetatomicPot` are reached through `RgpotPot` or an engine plugin, not through `RgpotAdapter`. They do not override `hessianRange`, so they stay at `0`.

Potentials that are not rgpot kernels also stay at `0` unless they override `hessianRange`. That includes the in-tree EMT potential and the in-tree TIP4P, SPC/E and TIP4P-Pt potentials on this pin. Those water potentials are the same flexible models as the 3.5.0 kernels below, and the same reason keeps coloring off.

### Where the numbers come from

Pair kernels. `LJPot`, `LJClusterPot` and `MorsePot` return `config().cutoff`. `ZBLPot` returns `cut_global`. The ZBL switch makes the energy and the force zero at that radius (`ZBLPot.cc`, the switching polynomial at `cut_global`).

Atom-centred Fortran kernels. The C++ face has no cutoff field. The defaults below are the values the neighbour table is built at, and the block range is twice that.

- `EAMAlPot`: `rcut = real(7.1, real32)` in `rgpot_eam_al.f90`, the EON r2549 default-real literal. Pairs and densities are dropped by `truncation_radius`, which is inside `rcut`.
- `SWPot`: `alpha * sigma` with `sigma = 2.0951` and `alpha = 1.8` in `rgpot_sw.f90`. The three-body term couples both ends of a triplet.
- `TersoffPot`: taper end `s = 3.0` in `rgpot_tersoff.f90`. `cutoff_taper` is zero at and beyond `S`.
- `EDIPPot`: `max(a, b)` with `a = b = 3.1213820` in `rgpot_edip.f90`.
- `LenoskyPot`: `default_cutoff = 4.5` in `rgpot_lenosky.f90`, shared by the pair and density splines.
- `FeHePot`: in `rgpot_fehe.f90` the pair cutoffs are Fe-Fe 5.3, Fe-He 3.9028 and He-He 5.4. d-band density vanishes at 4.2 and only Fe embeds it. s-band density vanishes at 4.10 and both species embed it. The range is the larger of the pair cutoffs and twice the largest density cutoff, `max(5.4, 2 * 4.2) = 8.4`.
- `CuH2Pot`: `rcut = 6.1` in `rgpot_cuh2.f90`. `phi` and `rho` are shifted to vanish there.

`WaterHPot`. `rgpot_water_h.f90` says there is no cutoff. Every site interacts with the lone H at its nearest periodic image, however far away. The range stays `0`.

UMA, on the pin and on `v3.5.0`. `UmaConfig::cutoff` (default 6 Å, `UmaConfig.hpp`) is the vesin edge cutoff. `UmaPot::forceImpl` builds that edge list and runs the compiled graph (`UmaPot.cc`, the `vesin_fairchem_edges` call). The package metadata can replace the cutoff. The number of message-passing layers is not a field of the kernel or of `UmaConfig`. `n_layers * cutoff` is not derivable, so the specialization is `MessagePassing` and the range is `0`.

Water on `v3.5.0`. `WaterConfig::cutoff` (default 8.5 Å, `WaterPots.hpp`) is the intermolecular cutoff passed into TIP4P, SPC/E and TIP4P-Pt. The CCL intramolecular term (`ccl.hpp`, `intramolecular`) is a polynomial in the bond stretch and the angle. It has no distance at which it is identically zero. SPC/E also exposes `coulombFull` (`spce_ccl.hpp`), a Coulomb path with no cutoff. A molecule is flexible, so there is no hard intramolecular extent to add to the cutoff. The range stays `0`.

Cell-list EAM on `v3.5.0`. `EAMCellConfig` is empty. `EAM/Parameters.h` lists `r_cut` of 5.555 Å (Al), 5.412 Å (Pd) and 5.5155 Å (Au). The specialization reports twice the largest of those cutoffs, `2 * 5.555` Å.

EMT on `v3.5.0`. `EMTConfig` is a parameter-set flag and has no cutoff. `EMTDefaultParameterProvider::calc_cutoff` builds a cutoff from the elements in the system, and the weight is a Fermi function. There is no single distance that is exact for every system, so the class is recorded and the range stays `0`.

`XTBPot`, `D3Pot`, `D4Pot`, `MOPACPot` and `ExprPot` use the primary template. D3 and D4 do not expose a cutoff on their config, and an expression is only as local as its terms.

## Checked zig-zags

Catch2 `Opted-in kernels match the serial matrix at the declared range` runs each opted-in kernel on its own zig-zag. The colored matrix matches the per-column matrix. A block past the interaction cutoff is nonzero for every many-body kernel, and every block past the declared range is 0. The short-range column is the negative test: SW_SI and EAM_AL are wrapped so the reported range is `r_c` instead of `2 r_c`, and that matrix misses the per-column matrix by more than `1e-2` of max `|H|`.

| Kernel | Atoms | Colors | Serial calls | Colored calls | max `|H|` | Largest nonzero block (Å) | Colored error | Short-range error |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EAM_AL | 14 | 11 | 85 | 67 | 0.730150 | 9.907 | 0 | 0.013011 |
| SW_SI | 10 | 7 | 61 | 43 | 16.16306 | 4.070 | 0 | 0.458243 |
| TERSOFF_SI | 8 | 5 | 49 | 31 | 14.54654 | 4.070 | 0 | not run |
| EDIP | 10 | 7 | 61 | 43 | 12.68665 | 4.070 | 0 | not run |
| LENOSKY_SI | 10 | 7 | 61 | 43 | 13.38487 | 4.613 | 0 | not run |
| FEHE | 12 | 9 | 73 | 55 | 3.670133 | 8.188 | 0 | not run |
| CUH2 | 14 | 11 | 85 | 67 | 2.062919 | 8.833 | 0 | not run |
| ZBL | 8 | 3 | 49 | 19 | 137.3325 | 1.500 | 0 | not run |

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: DHS_
keyprefix: dhess-
---
```
