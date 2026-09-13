# eOn 3.2.1

A patch on `v3.2.0`. See monorepo `CHANGELOG.md` section 3.2.1 for the
fragment list.

## IDPP honors frozen atoms

IDPP, collective IDPP, and SIDPP no longer write the full 3N vector.
Frozen sites next to a mover stay put.

## pyeonclient optimizer

`Parameters` exposes `opt_method`, `neb_opt_method`,
`refine_opt_method`, and `refine_threshold`.

## rgpot 3.1

Wrap pin is **v3.1.2**. `potential = dftd3` / `dftd4` and
`potential = expr` with `[ExprPot]` terms (`0.5*lj + d3`). D3/D4
wraps are off on Windows (MSVC rustc). MOPAC through rgpot is not in
this cut.

## ASE / Matter python surface

`from_ase` sets PBC before wrapping. Only `FixAtoms` fully freeze.
`setCell` / `setForces` dirty caches. Positions views are not writeable
in place. `to_ase` attaches SinglePoint energy/forces. NEB keeps the
GIL when the pot is not thread-safe.
