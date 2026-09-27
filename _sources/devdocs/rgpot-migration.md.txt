# rgpot kernel inventory

`eOn-6yzk` is the migration epic. This is the current split, not a
claim that the epic is done.

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
