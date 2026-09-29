# rgpot kernel inventory

This is the current split of kernels between eOn and rgpot, not a
claim that the move is finished.

## Already in rgpot (used via headers / FortranPots)

LJ, LJ cluster, Morse, ZBL, D3, D4, Expr, MOPAC, and the Fortran
kernels (SW, EDIP, Tersoff, Lenosky, EAM-Al, FeHe, CuH2, TIP4P-H)
when `with_fortran_pots` is on.

Highway SIMD for Morse/LJ belongs in **rgpot**, not `client/potentials`.

## Still in-tree (engines / adapters)

| Tree | Why it stays |
|---|---|
| ASE / ASE_NWCHEM / ASE_ORCA | Python calculator adapters |
| AMS / AMS_IO / VASP / LAMMPS / MPIPot / ExtPot / SocketNWChem | subprocess / file engines |
| Metatomic / XTBPot | optional heavy engines (also have rgpot engine .so loaders) |
| EMT / EAM / Water / Water_Pt / GPR / CatLearn | not yet in rgpot |

Adding a new empirical kernel should go to rgpot first.

## Neighbor lists

`eonc::VesinNeighbors` stays in eOn. Metatomic calls `vesin_neighbors`
directly and does not use that wrapper, so Metatomic is not the only
consumer of vesin headers. Wrap builds take vesin from the rgpot
subproject. Builds against an installed rgpot still compile eOn's
vendored translation unit, because an installed `rgpot.pc` does not
export vesin headers. Moving the wrapper into rgpot would be an rgpot
API change, not a deletion of `client/thirdparty/vesin`.

## Adding a kernel

New empirical kernels go into rgpot, then an `RgpotAdapter` arm in eOn.
See [Porting potentials](project:porting_potentials.md).
